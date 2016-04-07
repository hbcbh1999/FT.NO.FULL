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
*				dprint.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains statistics/printout driver routines.
*/


#include <driver/ddecs.h>
#include <util/cdecs.h>

enum	_TIME_STEP_SET_BY { TIME_STEP_NOT_SET,
			    TIME_STEP_SET_BY_FRONT,
			    TIME_STEP_SET_BY_WAVE,
			    TIME_STEP_SET_BY_PREVIOUS,
			    TIME_STEP_SET_BY_USER_LIMIT,
			    TIME_STEP_SET_BY_PRINTING
};
typedef enum _TIME_STEP_SET_BY TIME_STEP_SET_BY;

#if defined(__GD__)
struct 	_GD_POINT {
	double	xcrds;
	double	left_ycrds;
	double	right_ycrds;
	int	wave_type;
};
typedef	struct	_GD_POINT GD_POINT;

struct	_GD_GRAPH_DATA {
	int	imax;
	double 	xmin,xmax;
	double 	ymin,ymax;
	double	*xcrds;
	double	*ycrds;
	GD_POINT *points;
	int	num_points;
	struct _GD_GRAPH_DATA *prev;
	struct _GD_GRAPH_DATA *next;
	double	time;
};
typedef	struct	_GD_GRAPH_DATA GD_GRAPH_DATA;
#endif /* defined(__GD__) */

	/* LOCAL Function Declarations */
LOCAL	boolean	print_fronts_or_states(Grid*,Printplot*);
LOCAL	const char *TimeStepSetByString(TIME_STEP_SET_BY);
LOCAL	double	check_for_ts_restriction(OUTPUT_DATA*,double,double);
LOCAL	double	comm_time_step(double,double*,int,TIME_STEP_SET_BY*);
LOCAL	void	ensure_printout(OUTPUT_DATA*,Grid*);
LOCAL	void	print_TIME_DATA_stamp(FILE*,CHART*);
LOCAL	void	print_Time_Data_stamp(FILE*,Grid*);
LOCAL	void	print_front_and_wave(FILE*,Grid*,Wave*,Front*,Printplot*);
LOCAL	void	print_next_dt(Grid*,Wave*,Front*,Printplot*,int,FILE*);
LOCAL	void	print_solution(FILE*,Grid*,Wave*,Front*,Printplot*);
LOCAL	void	user_output(Grid*,Wave*,Front*,Printplot*,boolean);
LOCAL	void	wall_time_dump(CHART*,Printplot*,boolean,int);
LOCAL   void    print_repart_solution(CHART*,Printplot*,Grid*,Wave*,Front*);
LOCAL   int     count_digits_in_string(char* str); /* Determine how many digits are in the string */
LOCAL   boolean    file_exists(const char * filename); /* Determine if a file exists */

#if defined(ONED)
LOCAL	void	print_1d_solution_line(FILE*,double*,COMPONENT,Locstate,
				       Front*,Wave*,Printplot*,char*);
#endif /* defined(ONED) */

LOCAL	void	plot_prostar_data(Wave*,Front*,PROSTAR_plot_data*);
#if defined(__GD__)
LOCAL	void	plot_gd_data(Wave*,Front*,GD_plot_data*);
LOCAL	void 	gd_make_movie(GD_GRAPH_DATA*,char*,char*);
#endif /* defined(__GD__) */

#if defined(USE_AMR)
LOCAL   double   find_single_proc_time_step(Grid*,Wave*,Front*,double*,TIME_STEP_SET_BY*);
#endif /* if defined(USE_AMR) */

/* needed for VTK */
LOCAL   boolean    subsampling_is_valid(Wave*,Front*,int);
LOCAL   void    fill_vtk_values2d(Wave*,Front*,VTK_plot_data*);
LOCAL   void    fill_vtk_values2d_sub(Wave*,Front*,VTK_plot_data*);
LOCAL   void    fill_vtk_values3d(Wave*,Front*,VTK_plot_data*);
LOCAL   void	fill_vtk_values3d_sub(Wave*,Front*,VTK_plot_data*);
LOCAL   void    print_vtk_states(Wave*,Front*,VTK_plot_data*);
LOCAL   void    print_intfc_visit_file(int,int,VTK_PRINT_OPTIONS*);
LOCAL   void    print_variable_visit_file(int,char*,VTK_PRINT_OPTIONS*);
LOCAL	void	print_vtk_header_data(char*,int,int,int,double,double,double,double*,double,
					int,char*,int,int,boolean);
/* end needed for VTK  */



/* Jul 11 2002 : Myoung-Nyoun : Add : print only fuels using component*/
#if defined(TWOD)
LOCAL   void print_density_for_only_one_component(
        Grid		*grid,
	Printplot	*prt)
{
        OUTPUT_SOLN	*os = (prt->output_soln)[0];
	OUTPUT_VALUE	*(*solution)(OUTPUT_SOLN*,double*,int*) = os->solution;
	OUTPUT_VALUE	*value;
	RECT_GRID	*gr = grid->rect_grid;
	INTERFACE       *intfc= os->intfc;
	char            nfname[512];
	FILE            *nfile;
	COMPONENT       comp;
	int             xmin=0;
	int             ymin=0;
	int             xmax;
	int             ymax;
	int             ix,iy,icoords[MAXD],i,dim;
	double		coords[MAXD];
	      
	xmax=gr->gmax[0];
	ymax=gr->gmax[1];
	  
	(void) set_output_file_name(pp_mynode(),nfname,prt->outfile,
				grid->step,0);
	(void) sprintf(nfname,"%s-%s",nfname,os->name);
	if ((nfile = fopen(nfname,"w")) == NULL)
	{
	  screen("ERROR in print_density_for_only_one_component(), "
		 "can't open output file %s\n",nfname);
	  clean_up(ERROR);
	}
	
	(void) fprintf(nfile,"\n#ONLY FUEL %s Using COMPONENT DATA\n",os->name);
	for (iy = ymax - 1; iy >= ymin; iy--)
	{
	  icoords[1] = iy;
	  coords[1] = cell_edge(iy,1,gr);
	  for (ix = xmin; ix < xmax; ix++)
	  {
	    icoords[0] = ix;
	    coords[0] = cell_edge(ix,0,gr);
	    comp=component(coords,intfc);
	    if(is_exterior_comp(comp,intfc))
	      comp=nearest_interior_comp(YES,NO_COMP,coords,intfc);
	    if( comp==2 )
	    {
	      value = (*solution)(os,coords,icoords);
	      (void) fprintf(nfile,"%- 15.8g%s",value->uval.fval," ");
	    }
	    else
	      (void) fprintf(nfile,"0 ");
	  }
	  (void) fprintf(nfile,"\n");
	}
	(void) fprintf(nfile,"\n#End ONLY FUEL %s Using COMPONENT DATA\n",os->name);
	(void) Fclose(nfile);
}    /* end print_density_for_only_one_component() */
#endif /* defined(TWOD) */  
/* Jul 11 2002 : Myoung-Nyoun : Add : End */




/*
*			d_print_initial_data():
*/


EXPORT void d_print_initial_data(
	FILE		*file,
	CHART		*chart,
	Printplot	*prt)
{
	Grid	     *grid = chart->grid;
	RECT_GRID    *gr = grid->rect_grid;
	static const char *FORMAT = 
	    "\n      stop_time = %-10g               stop_step = %-10d\n";

	(void) fprintf(file,"\n\n\n");
	(void) foutput(file);
	(void) fprintf(file,"\t\t\tINITIAL DATA:\n\n\n");
	fprint_rectangular_grid(file,gr);
	(void) fprintf(file,FORMAT,stop_time(grid),stop_step(grid));
	(void) fprintf(file,"\n\t\tComputational ");
	switch (gr->Remap.remap)
	{
	case IDENTITY_REMAP:
	    switch (gr->dim)
	    {
	    case 1:
	        (void) fprintf(file,"Length");
	        break;
	    case 2:
	        (void) fprintf(file,"Area");
	        break;
	    case 3:
	        (void) fprintf(file,"Volume");
	        break;
	    }
	    break;
	case CYLINDRICAL_REMAP:
	    switch (gr->dim)
	    {
	    case 1:
	        (void) fprintf(file,"Area");
	        break;
	    case 2:
	        (void) fprintf(file,"Volume");
	        break;
	    case 3:
	    default:
	        screen("ERROR in d_print_initial_data(), "
		       "3D CYLINDRICAL_REMAP not supported\n");
	        clean_up(ERROR);
	        break;
	    }
	    break;
	case SPHERICAL_REMAP:
	    switch (gr->dim)
	    {
	    case 1:
	        (void) fprintf(file,"Volume");
	        break;
	    case 2:
	        screen("ERROR in d_print_initial_data(), "
		       "2D SPHERICAL_REMAP not supported\n");
	        clean_up(ERROR);
	        break;
	    case 3:
	    default:
	        screen("ERROR in d_print_initial_data(), "
		       "3D SPHERICAL_REMAP not supported\n");
	        clean_up(ERROR);
	        break;
	    }
	    break;
	case INVALID_REMAP:
	default:
	    screen("ERROR in d_print_initial_data(), invalid remap %d\n",
		   gr->Remap.remap);
	    clean_up(ERROR);
	    break;
	}
	(void) fprintf(file," = %g\n",gr->Remap.area);
	(void) fprintf(file,"\n\t\tRemap Geometry:  %s\n",
		       gr->Remap.remap_name);

	if (debugging("init_printout"))
	{
	    (void) printf("\t\tBEFORE TIME LOOP:\n");
	    if (prt->printout != NULL)
	    {
	    	int step_sav;

	    	step_sav = grid->step;
	    	if (step_sav > 0)
		    grid->step = step_sav-1;

	    	(*prt->printout)(chart,prt,NO,0);

	    	if (step_sav > 0) grid->step = step_sav;
	    }
	}
}		/*end d_print_initial_data*/


/*
*				d_printout():
*
*	Prints results for the purpose of plotting and data analysis.
*/

/*ARGSUSED*/
EXPORT void d_printout(
	CHART		*chart,
	Printplot 	*prt,
	boolean		about_to_stop,
	int		error)
{
	Grid		*grid = chart->grid;
	Wave		*wave = chart->wave;
	Front		*front = chart->front;
	FILE		*file = stdout;
	char		fname[1024];
	OUTPUT_DATA	**mod = prt->main_output_data;
	boolean		extra_print = NO;
	boolean		print_time = NO;
	boolean		prt_fr_or_sts;
	int		i;

	wall_time_dump(chart,prt,about_to_stop,error);
	interactive_printing_control(prt,grid,&extra_print,&about_to_stop);

	if (extra_print == YES)
	{
	    print_time = YES;
	    prt_fr_or_sts = YES;
	}
	else if (about_to_stop == YES)
	{
	    print_time = YES;
	    prt_fr_or_sts = YES;
	    for (i = 0; i < NUM_OUTPUT_FORMATS; ++i)
	    {
	        if (output_format_on(i,prt))
		    ensure_printout(prt->main_output_data[i],grid);
	    }
	} 
	else
	{
	    prt_fr_or_sts = print_fronts_or_states(grid,prt);
	    for (i = 0; i < NUM_OUTPUT_FORMATS; ++i)
	    {
	        if (is_ts_for_output_format(i,grid,prt))
	    	{
	    	    print_time = YES;
		    break;
	    	}
	    }
	}
        if (!(about_to_stop == YES && grid->repart_at_end_of_run == YES))
       	{
	    if (print_time == NO)
	    {
	        print_Time_Data_stamp(file,grid);
	        print_extreme_values(file,chart,prt);
	    }
	    else
	    { 
		if (prt->outfile != NULL)
	        {
	             /* print time stamp to master output file */
		    print_Time_Data_stamp(file,grid);
                    print_extreme_values(file,chart,prt);

		    if (prt_fr_or_sts == YES)
		    {
	                set_output_file_name(pp_mynode(),fname,
				prt->outfile,grid->step,error);
	                if ((file = fopen(fname,"w")) == NULL)
	                {
		            screen("ERROR in d_printout(), "
	                         "can't open output file %s\n",fname);
	                    clean_up(ERROR);
	                }
         	        if (debugging("nobuf"))
	                    setbuf(file,NULL);
	                record_print_version(file);
	                print_title(file,prt->title);
	                if (error != 0)
		            (void)printf("\n\nCLEAN_UP printout, error = %d\n\n",					error);
	                if (prt->print_initial_data)
		            (*prt->print_initial_data)(file,chart,prt);
		    }
		    else
		        file = NULL;
	        }

	        if (file != NULL)
	        {
	            	/* Print Time Information */

	            print_TIME_DATA_stamp(file,chart);
	            print_extreme_values(file,chart,prt);
	        }
            
	                 /*print interface and states*/
	        print_solution(file,grid,wave,front,prt);
#if defined(USE_AMR)
                if(YES == print_time)
                    /*overture_dump_st_printout(chart, prt); */
#endif /* if defined(USE_AMR) */
	        if (prt->plot_ellip_vars != NULL)
	            (*prt->plot_ellip_vars)(grid);
	   }
		/* Compute Statistics */

	   if ((prt->grid_statistics != NULL) && (error == 0))
	   {
	       start_clock("STATISTICS");
	       (*prt->grid_statistics)(prt->gs_data,grid,wave,front,YES);
	       stop_clock("STATISTICS");
	   }

	   if (debugging("wall_time"))
	       (void) printf("\n\t\tCurrent Wall Time:	%s\n",date_string());

			/* More specific statistics/diagnostics */

	   user_output(grid,wave,front,prt,about_to_stop);
	   (void) printf("\n\n\n");
	   if (print_time == YES)
	   {
	       for (i = 0; i < NUM_OUTPUT_FORMATS; ++i)
	       {
	           if (is_ts_for_output_format(i,grid,prt))
		   {
		       if (real_time_output(Output_mode(mod[i])))
		           Output_next_print_time(mod[i]) +=
			       Output_time_freq(mod[i]);
		       else
		           Output_next_print_step(mod[i]) +=
			       Output_step_freq(mod[i]);
		   }
	       }

	       if (file != NULL)
	       {
	   	   print_next_dt(grid,wave,front,prt,error,file);
	           if (file != stdout)
	           {
		       trace_foutput(file);
		       (void) Fclose(file);
		       if (prt->compress == YES)
		       {
		           char s[256];
		           (void) sprintf(s,"gzip %s &",fname);
			   (void) system(s);
		       }
	           }
	       }
	   }

	   if (about_to_stop == NO)
	   {
	       IMPORT boolean  suppress_prompts;
	       boolean         sav_suppress_prompts;
	       INIT_DATA    *init = grid->init;

	       sav_suppress_prompts = suppress_prompts;
	       suppress_prompts = (interactive_prompting(init) == YES) ? 
	       				NO : YES;

			/*  Call await() at Pause Time */

	       if (real_time_output(pause_mode(grid)) &&
	           (grid->time >= pause_time(grid)))
	       {
	           screen("PAUSE - pause time pause_time(grid) = %g reached\n",
	    	          pause_time(grid));
	           await(grid,front,wave,prt,NO);
	       }
	       else if (mesh_time_output(pause_mode(grid)) &&
	            (grid->step == pause_step(grid)))
	       {
	           screen("PAUSE - pause step pause_step(grid) = %d reached\n",
	    	          pause_step(grid));
	           await(grid,front,wave,prt,NO);
	       }
	       suppress_prompts = sav_suppress_prompts;
	   }
        }
	else
	{
	    print_repart_solution(chart,prt,grid,wave,front);
	}

}   	/*end d_printout*/

LOCAL	void	print_next_dt(
	Grid	  *grid,
	Wave	  *wave,
	Front	  *front,
	Printplot *prt,
	int	  error,
	FILE	  *file)
{
	OUTPUT_DATA	**mod = prt->main_output_data;
	double		next_dt;
	next_dt = (error == 0) ? find_time_step(grid,wave,front,prt,NO) :
				 grid->dt;
	(void) foutput(file);
	(void) fprintf(file,"\tnext_dt:  ");
	if (output_format_on(RECT_STATES,prt) &&
	    (Output_in_binary(mod[RECT_STATES]) == YES))
	{
	    /*
	     * We assume that binary printing of states implies restart
	     * is the main concern, thus a binary next_dt.
	     */

	    (void) fprintf(file,"\f%c",2);
	    (void) fwrite((const void *) &grid->time,sizeof(double),1,file);
	    (void) fwrite((const void *) &next_dt,sizeof(double),1,file);
	}
	else
	    (void) fprintf(file,"%-"FFMT"\n",next_dt);
}		/*end print_next_dt*/

LOCAL	void	ensure_printout(
	OUTPUT_DATA	*data,
	Grid		*grid)
{
	Output_mode(data) = MESH_TIME;
	Output_start_step(data) = 0;
	Output_next_print_step(data) = grid->step;
}		/*end ensure_printout*/

EXPORT	void	set_output_file_name(
	int             my_node,
	char		*fname,
	const char	*bname,
	const int	step,
	const int	error)
{
	char buff[4096],buff2[4096];

	(void) strcpy(buff,bname);
	if (error != 0)
	    (void) strcat(buff,".CLEAN_UP");
	if (step >= 0)
	{
	    (void) sprintf(buff2,"%s.ts%s",buff,
			   right_flush(step,TSTEP_FIELD_WIDTH));
	}
	else
	{
	    (void) sprintf(buff2,"%s.ts%s",buff,
			   right_flush(0,TSTEP_FIELD_WIDTH));
	}
	(void) sprintf(fname,"%s",buff2);
	if (pp_numnodes() > 1)
	    (void) sprintf(fname,"%s-nd%s",buff2,
				 right_flush(my_node,PP_NODE_FIELD_WIDTH));
}		/*end set_output_file_name*/

EXPORT	void	adjoin_node_number(
	char	*s)
{
	int myid = pp_mynode();
	int nn = pp_numnodes();
	int nd;
	for (nd = 0; nn != 0; nn /=10, ++nd);
	(void) sprintf(s,"%s-nd%s",s,right_flush(myid,nd));
}		/*end adjoin_node_number*/

LOCAL	void	wall_time_dump(
	CHART		*chart,
	Printplot 	*prt,
	boolean		about_to_stop,
	int		error)
{
	PRINT_OPTIONS	*pto = wall_time_dump_options(prt);
	OUTPUT_DATA	*prt_save_output_data[NUM_OUTPUT_FORMATS];
	Grid		*grid;
	Wave		*wave;
	Front		*front;
	double		current_time, elapsed_time;
	int		i;
	FILE		*file;
	char		fname[1024];
	boolean		sav_binary;
	int		do_wall_time_print;
	static OUTPUT_DATA **wall_dump_output_data = NULL;

	debug_print("walltime","Entered wall_time_dump()\n");

	if ((error != 0) || (about_to_stop == YES) || (pto == NULL))
	{
	    debug_print("walltime","Left wall_time_dump()\n");
	    return;
	}

	if (print_wall_time_interval(pto) < 0.0)
	{
	    debug_print("walltime","Left wall_time_dump()\n");
	    return;
	}

	current_time = real_time();
	elapsed_time = current_time - last_wall_time_dump(pto);
	do_wall_time_print =
	    (elapsed_time < print_wall_time_interval(pto)) ? 0 : 1;
	pp_global_imax(&do_wall_time_print,1);
	if (do_wall_time_print == 0)
	{
	    debug_print("walltime","Left wall_time_dump()\n");
	    return;
	}

	last_wall_time_dump(pto) = current_time;

	(void) sprintf(fname,"%s%d",
		       print_filename(pto),wall_time_dump_number(pto)++);
	wall_time_dump_number(pto) = wall_time_dump_number(pto)%2;
	if (debugging("walltime"))
	{
	    (void) printf("wall_time_dump_number = %d\n",
			  wall_time_dump_number(pto));
	}
	set_output_file_name(pp_mynode(),fname,fname,-1,0);

	grid = chart->grid;
	wave = chart->wave;
	front = chart->front;
	if (debugging("walltime"))
	    (void) printf("Opening wall time dump file %s\n",fname);

	if ((file = fopen(fname,"w")) == NULL)
	{
	    screen("ERROR in wall_time_dump(), can't open %s\n",fname);
	    clean_up(ERROR);
	}

	if (wall_dump_output_data == NULL)
	{
	    PRINT_OPTIONS Pto;

	    uni_array(&wall_dump_output_data,
		   NUM_OUTPUT_FORMATS,sizeof(OUTPUT_DATA*));

	    set_defaults_for_print_options(&Pto,grid->init);
	    Prt_mode(Pto) = MESH_TIME;
	    Print_step_interval(Pto) = 1;
	    Print_start_step(Pto) = 0;
	    Print_in_binary(Pto) = print_in_binary(pto);

	    scalar(&wall_dump_output_data[PRT_FRONTS],sizeof(OUTPUT_DATA));
	    set_output_data(&Pto,wall_dump_output_data[PRT_FRONTS],grid,
			    prt,NO,NO);
	    scalar(&wall_dump_output_data[RECT_STATES],sizeof(OUTPUT_DATA));
	    set_output_data(&Pto,wall_dump_output_data[RECT_STATES],grid,
			    prt,NO,NO);
	    wall_dump_output_data[TRI_STATES] = NULL;
	    wall_dump_output_data[HDF_STATES] = NULL;
	    wall_dump_output_data[SDS_STATES] = NULL;
	}

	for (i = 0; i < NUM_OUTPUT_FORMATS; ++i)
	{
	    prt_save_output_data[i] = prt->main_output_data[i];
	    prt->main_output_data[i] = wall_dump_output_data[i];
	}
	ensure_printout(prt->main_output_data[PRT_FRONTS],grid);
	ensure_printout(prt->main_output_data[RECT_STATES],grid);

	sav_binary = is_binary_output();
	set_binary_output(print_in_binary(pto));

	record_print_version(file);
	print_title(file,prt->title);
	if (prt->print_initial_data)
	    (*prt->print_initial_data)(file,chart,prt);
	print_TIME_DATA_stamp(file,chart);

	print_solution(file,grid,wave,front,prt);

	print_next_dt(grid,wave,front,prt,error,file);
	trace_foutput(file);
	(void) Fclose(file);

	set_binary_output(sav_binary);
	for (i = 0; i < NUM_OUTPUT_FORMATS; ++i)
	    prt->main_output_data[i] = prt_save_output_data[i];
	debug_print("walltime","Left wall_time_dump()\n");

}		/*end wall_time_dump*/

EXPORT void interactive_printing_control(
	Printplot	*prt,
	Grid		*grid,
	boolean		*extra_print,
	boolean		*about_to_stop)
{
	FILE		*fp;
	OUTPUT_DATA	**mod = prt->main_output_data;
	int		i;
	int		step;
	static	char	filename[NUM_OUTPUT_FORMATS][25];
	static	boolean	first = YES;

	if (first == YES)
	{
	    first = NO;
	    for (i = 0; i < NUM_OUTPUT_FORMATS; ++i)
	        (void) sprintf(filename[i],"%s.ts",output_format_name(i));
	}

	for (i = 0; i < NUM_OUTPUT_FORMATS; ++i)
	{
	    if ((mod[i] != NULL) && (fp = fopen(filename[i],"r")) != NULL)
	    {
		if (fscanf(fp,"%d",&step) != EOF)
		{
		    do
		    {
		        if ((grid->step == step) && 
			(is_print_time(grid->time,Print_control(mod[i])) == NO)
						&&
			(is_print_step(grid->step,Print_control(mod[1])) == NO))
			    *extra_print = YES;
		    }
		    while (fscanf(fp,"%d",&step) != EOF);
	        }
		else
		{
		    *about_to_stop = YES;
		    stop_step(grid) = grid->step;
		    (void) printf("\nRun being stopped from outside\n");
	        }
		(void) fclose(fp);
	    }
	}
}		/*end interactive_printing_control*/


/*
*				user_output():
*
*	Provides a method of recording data at time step intervals.
*	Each output function is loaded in an array in the Printplot.
*	This function steps through the array, executing each function
*	in turn.  Each function has its own data structure of
*	information stored as an array of generic pointers in the
*	Printplot.  The correspondence is one-to-one, so the ith
*	function corresponds to the ith data pointer.  An OUTPUT_DATA
*	structure with suitable values should be the first element of
*	any data structure intended to be used in prt->user_output_data.
*
*	To add a new function to the list, the function
*	add_user_output_function() is provided.
*/

/* ARGSUSED */
LOCAL void user_output(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	boolean		about_to_stop)
{
	int	    i;
	void	    (**f)(Grid*,Wave*,Front*,Printplot*,OUTPUT_DATA*,boolean);
	OUTPUT_DATA *data;

	if (prt->user_output_funcs == NULL)
	    return;

	debug_print("user_output","Entered user_output()\n");
	for (i = 0, f = prt->user_output_funcs; *f != NULL; ++f, ++i)
	{
	    data = prt->user_output_data[i];
	    if ((is_print_time(grid->time,Print_control(data)) == YES) ||
	        (is_print_step(grid->step,Print_control(data)) == YES) ||
		about_to_stop)
	    {
	        (*(*f))(grid,wave,front,prt,data,about_to_stop);
	        if (real_time_output(Output_mode(data)))
		    Output_next_print_time(data) += Output_time_freq(data);
	        else
		    Output_next_print_step(data) += Output_step_freq(data);
	    }
        }
	debug_print("user_output","Left user_output()\n");
}		/*end user_output*/


LOCAL	void print_Time_Data_stamp(
	FILE		*file,
	Grid		*grid)
{
	if (file == NULL)
	    return;
	(void) fprintf(file,"\n\n\n\n");
	(void) foutput(file);
	(void) fprintf(file,"\tTime Data:  t = %-"FFMT" j = %-10d dt = %-"FFMT"\n",
		       grid->time,grid->step,grid->dt);
}		/*end print_Time_Data_stamp*/

LOCAL	void print_TIME_DATA_stamp(
	FILE		*file,
	CHART		*chart)
{
	Grid		*grid = chart->grid;

	if (file == NULL)
	    return;
	(void) fprintf(file,"\n\n\n\n");
	(void) foutput(file);
	(void) fprintf(file,"\tTIME DATA:  t = %-"FFMT" j = %-10d dt = %-"FFMT"\n",
		       grid->time,grid->step,grid->dt);
	(void) foutput(file);
}		/*end print_TIME_DATA_stamp*/

/*
*			print_solution():
*
*	Calls print_front_and_wave() for the front and wave
*	associated with each chart that is below chart_of_front(front)
*	in the tree structure.
*/

LOCAL	void print_solution(
	FILE		*file,
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt)
{
	LEVEL		*lev;
	CHART		*root, *chart;

		/* print information for root */

	print_front_and_wave(file,grid,wave,front,prt);
	root = chart_of_front(front);
	if ((root != NULL) && (root->level != NULL))
	{
	    print_LEVEL(file,root->level);

	    for (lev=root->level->next_finer_level; lev;
	         lev=lev->next_finer_level)
	    {
	        print_LEVEL(file,lev);
	        for (chart = lev->first; chart; chart = chart->next_chart)
	        {
	    	    print_CHART(file,chart);
	    	    print_front_and_wave(file,grid,chart->wave,
		                         chart->front,prt);
	        }
	    }
	}
}		/*end print_solution*/

/*
*       Set up new zoom_grid and print solution with new parallel partition
*/
LOCAL void print_repart_solution(
	CHART           *chart,
	Printplot       *prt,
	Grid		*grid,
	Wave            *wave,
	Front		*front)
{
	FILE            *file;
	Front          	*tempfront;
	RECT_GRID      	dual_grid;
	RECT_GRID      	*zoom_grid = &(grid->new_pp_grid->Zoom_grid);
        char           	fname[1024];
        double          	L[MAXD], U[MAXD], *GL, *GU;
        int            	i,j,k;
	int            	ratio[MAXD];
        int            	dim = front->rect_grid->dim;
        int            	icrds_old[MAXD],icrds_new[MAXD];
	int             gmax[MAXD],lbuf[MAXD], ubuf[MAXD];       
		 
	GL = grid->pp_grid->Global_grid.L;    
	GU = grid->pp_grid->Global_grid.U;

	find_Cartesian_coordinates(pp_mynode(),grid->pp_grid,icrds_old);
	for (i = 0; i < dim; ++i)
            ratio[i] = (grid->new_pp_grid->gmax[i])/(grid->pp_grid->gmax[i]);
	for (i = 0; i < grid->new_pp_grid->nn; ++i)
       	{
	    find_Cartesian_coordinates(i,grid->new_pp_grid,icrds_new);
	    
	    k = 1;
	    for (j = 0; j < dim; ++j)
            {  
                k = k && ((ratio[j]*icrds_old[j] <= icrds_new[j])
			    &&(icrds_new[j] <= ratio[j]*icrds_old[j] + 
		               ratio[j]-1));
            }
            if (k)
            {	
	      
	       /*  set the zoom_grid of new_pp_grid here */	
		for (j = 0; j < dim; ++j)
	        {
	            L[j] = grid->new_pp_grid->dom[j][icrds_new[j]];
	            U[j] = grid->new_pp_grid->dom[j][icrds_new[j] + 1];
	            gmax[j] = irint((U[j] - L[j])/(front->rect_grid->h[j]));

		    switch (dim) /* TODO Unify 2 and 3 D */
	            {
	            case 1:
	            case 2:
	    	        lbuf[j] = (icrds_new[j] > 0) ? 
			             grid->new_pp_grid->buf[j] : 0;
	    	        ubuf[j] = (icrds_new[j] <(grid->new_pp_grid->gmax[j]-1))
			           ? grid->new_pp_grid->buf[j]:0;
	    	        break;
	            case 3:
	    	        lbuf[j] = (icrds_new[j] > 0) 
				? grid->new_pp_grid->buf[j] : 
				front->rect_grid->lbuf[j];

	    	        ubuf[j] = (icrds_new[j] <(grid->new_pp_grid->gmax[j]-1))
				? grid->new_pp_grid->buf[j] :
				front->rect_grid->ubuf[j];
			break;
	            }
	        }
	        set_rect_grid(L,U,GL,GU,lbuf,ubuf,gmax,dim,
		    &front->rect_grid->Remap,&grid->new_pp_grid->Zoom_grid);
                
	           /* set up filename and output front and wave*/
	        set_output_file_name(i,fname,prt->outfile,grid->step,0);
		if ((file = fopen(fname,"w")) == NULL)
		{
	            screen("ERROR: unable to open file '%s'\n",fname);
		    clean_up(ERROR);
		}
	         
		record_print_version(file);
	        print_title(file,prt->title);
		if (prt->print_initial_data)
		    (*prt->print_initial_data)(file,chart,prt);
		print_TIME_DATA_stamp(file,chart);
	        print_extreme_values(file,chart,prt);
		      
		      /*clip the front with the new zoom_grid*/
	          
	        tempfront = copy_front(front);
		set_copy_intfc_states(YES);
	        set_size_of_intfc_state(size_of_state(front->interf));
		tempfront->interf = copy_interface(front->interf);
		for (j = 0; j < dim; ++j)
	        {
		    if (icrds_new[j] > 0)
		    	rect_boundary_type(tempfront->interf,j,0) 
				= SUBDOMAIN_BOUNDARY;
		    if (icrds_new[j] <(grid->new_pp_grid->gmax[j]-1))
		    	rect_boundary_type(tempfront->interf,j,1) 
				= SUBDOMAIN_BOUNDARY;
		}
			
			/*set topological grid here*/
		set_dual_grid(&dual_grid,zoom_grid);	
		for (j = 0; j < dim; ++j)
	            gmax[j] = dual_grid.gmax[j]+dual_grid.lbuf[j]+
		                       dual_grid.ubuf[j];
		set_rect_grid(dual_grid.VL,dual_grid.VU,dual_grid.GL,
		      dual_grid.GU,NOBUF,NOBUF,gmax,dim,&zoom_grid->Remap,
			             &topological_grid(tempfront->interf));
		set_computational_grid(tempfront->interf,zoom_grid); 
			   
			   /*clip the interface to subdomain*/
		clip_front_for_output(tempfront,zoom_grid);
			   /*print repart front and wave*/

			/*set relative coordinates of new nodes */
		for(j = 0; j < dim; ++j)
		    icrds_new[j] = icrds_new[j] - icrds_old[j]*ratio[j];

		for (j = 0; j < prt->n_rect_state_vars; ++j)
		{
		    prt->output_soln[j]->repart_at_end = YES;
		    prt->output_soln[j]->icrds_new = icrds_new;
		}

			/*print repart front and wave*/
		print_front_and_wave(file,grid,wave,tempfront,prt);
		print_next_dt(grid,wave,tempfront,prt,0,file);
		trace_foutput(file);
		(void)Fclose(file);
	        free_front(tempfront);
	   } 
        }
}		/*end print_repart_solution*/

/*
*			print_front_and_wave():
*
*	Provides printout of front->interf and of all state variables
*	associated with wave.
*/

LOCAL void print_front_and_wave(
	FILE		*file,
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt)
{
	int		var;
	int		dim = front->rect_grid->dim;
	POINTER		extra[2];
	OUTPUT_DATA	**mod = prt->main_output_data;
	boolean		sav_binary = is_binary_output();

			/* Initialize */

	if (prt->initialize_for_printout != NULL)
	    (*prt->initialize_for_printout)(grid,wave,front,prt);

	extra[0] = (POINTER) front;
	extra[1] = (POINTER) wave;
	for (var = 0; var < prt->n_rect_state_vars; ++var) 
	{
	    prt->output_soln[var]->intfc = front->interf;
	    prt->output_soln[var]->extra = (POINTER) extra;
	}

	if ((print_fronts_or_states(grid,prt) == YES) && (file != NULL))
	{
	    /* Print Front.  We assume that printing or plotting
	     * states is more useful WITH the fronts. */

	    if      (output_format_on(PRT_FRONTS,prt))
	    	set_binary_output(Output_in_binary(mod[PRT_FRONTS]));
	    else if (output_format_on(RECT_STATES,prt))
	    	set_binary_output(Output_in_binary(mod[RECT_STATES]));
	    else 
	    	set_binary_output(Output_in_binary(mod[TRI_STATES]));

	    (void) fprintf(file,"\n\n\n\n");
	    (void) foutput(file);
	    (void) fprintf(file,"\t\t\tFRONT DATA:\n");
	    fprint_front(front,file);
	    (void) fprintf(file,"\n\n\n\n");
	    (void) foutput(file);
	    (void) fprintf(file,"\t\t\tEND OF FRONT DATA:\n");
	}

	if ((is_ts_for_output_format(RECT_STATES,grid,prt) ||
	     is_ts_for_output_format(TRI_STATES,grid,prt)) && (file != NULL))
	{
	    (void) fprintf(file,"\n\n\n\n");
	    (void) foutput(file);
	    (void) fprintf(file,"\t\t\tSTATE DATA:\n");

	    	/* Printout Components of Locstates in Wave */

	    if (is_ts_for_output_format(RECT_STATES,grid,prt))
	    {
	    	set_binary_output(Output_in_binary(mod[RECT_STATES]));
	    	print_states(file,prt,dim);
	    }

	    	/* Plot Components of Locstates in Wave */

	    if (is_ts_for_output_format(TRI_STATES,grid,prt))
	    {
	    	set_binary_output(Output_in_binary(mod[TRI_STATES]));
	    	plot_states(file,front,wave,prt);
	    }

	    (void) foutput(file);
	    (void) fprintf(file,"\t\t\tWAVE SPEEDS:\n");
	    print_max_wave_speed_info(file,wave);

	    if (debugging("pp_solution"))
	    	(*wave->show_wave_states)(wave);

	    (void) fprintf(file,"\n\n\n\n");
	    (void) foutput(file);
	    (void) fprintf(file,"\t\t\tEND OF STATE DATA\n");
	}
	set_binary_output(sav_binary);

/* needed for VTK */
	if (is_ts_for_output_format(VTK_STATES,grid,prt))
	{
	    for (var = 0; var < prt->n_VTK_vars; ++var)
	        print_vtk_states(wave,front,prt->vtk_data+var);
        }

/* end needed for VTK */
#if defined(USE_HDF)
	if (is_ts_for_output_format(HDF_STATES,grid,prt))
	{
	    for (var = 0; var < prt->n_HDF_vars; ++var)
	    	plot_hdf_data(wave,front,prt->HDF_data+var);
	}
	if (is_ts_for_output_format(SDS_STATES,grid,prt))
	{
	    for (var = 0; var < prt->n_SDS_vars; ++var)
	    	plot_hdf_data(wave,front,prt->SDS_data+var);
	}
#endif /* defined(USE_HDF) */
	if (is_ts_for_output_format(PROSTAR_STATES,grid,prt))
	{
	    for (var = 0; var < prt->n_PROSTAR_vars; var++)
	    	plot_prostar_data(wave,front,prt->PROSTAR_data+var);
	}
#if defined(__GD__)
	if (is_ts_for_output_format(GD_MOVIE,grid,prt))
	{
	    for (var = 0; var < prt->n_GD_vars; var++)
	    	plot_gd_data(wave,front,prt->GD_data+var);
	}
#endif /* defined(__GD__) */
}		/*end print_front_and_wave*/

EXPORT	OUTPUT_DATA	**d_alloc_output_datas(
	int num_datas)
{
	OUTPUT_DATA	**datas;

	bi_array(&datas,num_datas,1,sizeof(OUTPUT_DATA));
	return datas;
}		/*end alloc_output_datas*/


LOCAL	boolean	print_fronts_or_states(
	Grid		*grid,
	Printplot	*prt)
{
	return (is_ts_for_output_format(PRT_FRONTS,grid,prt) ||
	    	is_ts_for_output_format(RECT_STATES,grid,prt) ||
	    	is_ts_for_output_format(TRI_STATES,grid,prt)) ?
		YES : NO;
}		/*end print_fronts_or_states*/


#if defined(ONED)
/*
*			plot_states1d():
*
*	Prints state variables using a graphs format.
*/

/*ARGSUSED*/
EXPORT  void plot_states1d(
	FILE		*file,
	Front		*front,
	Wave		*wave,
	Printplot	*prt)
{
	INTERFACE	*intfc = front->interf;
	CRXING		*crx;
	TRI_GRID	*grid = wave_tri_soln(wave)->tri_grid;
	Table		*T = table_of_interface(grid->grid_intfc);
	RECT_GRID	*gr = &grid->rect_grid;
	COMPONENT	*comps = T->components;
	Locstate	state;
	byte		*storage;
	double		coords[MAXD];
	size_t		sizest = front->sizest;
	size_t		len;
	size_t		max_len;
	int             precision;
	int		var;
	size_t		num_var = prt->n_tri_vars;
	int		ix, xmax;
	int		i, nc, *list;
	char		fmt[80];

	if (num_var == 0)
	    return;

	max_len = 26;
	precision = 20;
	for (var = 0;  var < num_var;  ++var)
	{
	    len = strlen(prt->tri_plot_name[var]);
	    if (len > max_len)
	    	max_len = len;
	}
	(void) sprintf(fmt," %%-%ds",(int)max_len);
	(void) foutput(file);
	(void) fprintf(file,fmt,"POSITION");
	for (var = 0;  var < num_var;  ++var)
	    (void) fprintf(file,fmt,prt->tri_plot_name[var]);
	(void) fprintf(file,"\n");

	xmax = gr->gmax[0];
	storage = grid->rect_state_storage;
	(void) sprintf(fmt," %%-%d.%dg",(int)max_len,precision);
	for (ix = 0; ix < xmax; ++ix)
	{
	    if (!is_excluded_comp(comps[ix],intfc))
	    {
	    	state = (Locstate)(storage + ix*sizest);
	    	coords[0] = cell_edge(ix,0,gr);
	    	print_1d_solution_line(file,coords,comps[ix],state,
				       front,wave,prt,fmt);
	    }
	    nc = T->seg_crx_count[ix];
	    list = T->seg_crx_lists[ix];
	    for (i = 0; i < nc; ++i)
	    {
	    	crx = &(T->crx_store[list[i]]);
	    	coords[0] = Coords(crx->pt)[0];
	    	print_1d_solution_line(file,coords,negative_component(crx->pt),
				       left_state(crx->pt),front,wave,prt,fmt);
		print_1d_solution_line(file,coords,positive_component(crx->pt),
				       right_state(crx->pt),front,wave,prt,fmt);
	    }
	}
}		/*end plot_states1d*/

LOCAL	void	print_1d_solution_line(
	FILE		*file,
	double		*coords,
	COMPONENT	comp,
	Locstate	state,
	Front		*front,
	Wave		*wave,
	Printplot	*prt,
	char		*fmt)
{
	int	var;
	size_t	num_var = prt->n_tri_vars;
	static	double	*fvals = NULL;
	static	size_t	fval_len = 0;

	if (fvals == NULL)
	{
	    fval_len = num_var+1;
	    uni_array(&fvals,fval_len,FLOAT);
	}
	else if (fval_len < (num_var+1))
	{
	    free(fvals);
	    fval_len = num_var+1;
	    uni_array(&fvals,fval_len,FLOAT);
	}
	fvals[0] = coords[0];
	for (var = 0; var < num_var; ++var)
	{
	    fvals[var+1] = (*prt->tri_plot_function[var])(coords,front,
							  wave,comp,state);
	}

	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",(char)(num_var+1));
	    (void) fwrite((const void *) fvals,FLOAT,num_var+1,file);
	}
	else
	{
	    for (var = 0; var <= num_var; ++var)
	    	(void) fprintf(file,fmt,fvals[var]);
	    (void) fprintf(file,"\n");
	}
}		/*end print_1d_solution_line*/
#endif /* defined(ONED) */

#if defined(TWOD)

/*
*			plot_states2d():
*
*	Calls print_tri_soln() to generate input data for the
*	program "tri", which will make discontinuous contour plots.
*/

EXPORT  void plot_states2d(
	FILE		*file,
	Front		*front,
	Wave		*wave,
	Printplot	*prt)
{
	double	   *gl = wave->rect_grid->GL;
	double	   *gu = wave->rect_grid->GU;
	int	   num_sources = wave->num_point_sources;
	int	   i, j, var, dim = wave->rect_grid->dim;
	static const char *CL[3] = {"CXL", "CYL", "CZL"};
	static const char *CU[3] = {"CXU", "CYU", "CZU"};

	/*if (prt->n_tri_vars == 0) */
	  /*  return; */

	(void) foutput(file);
	(void) fprintf(file,"TRI_SOLN:");

	for (var = 0;  var < prt->n_tri_vars;  ++var)
	{
	    (void) fprintf(file,"%s%s"," ",prt->tri_plot_name[var]);
	}
	(void) fprintf(file,"\n");

	(void) fprintf(file,"#Point Source Data\n#%d\n",num_sources);
	for (i = 0;  i < num_sources;  ++i)
	{
	    (void) fprintf(file,"#%d ",wave->source_type[i]);
	    for (j = 0; j < dim; ++j)
	    	(void) fprintf(file,"%g ",wave->srcpt[j][i]);
	    (void) fprintf(file,"%g\n",wave->strength[i]);
	}

	(void) fprintf(file,"#Clipping Boundary\t");
	for (i = 0; i < dim; ++i)
	    (void) fprintf(file,"%s = %g %s = %g ",CL[i],gl[i],CU[i],gu[i]);
	(void) fprintf(file,"\n");

	print_tri_soln(file,front,wave,wave_tri_soln(wave),
		       prt->n_tri_vars,prt->tri_plot_function);

	(void) foutput(file);
	(void) fprintf(file,"END OF TRI_SOLN\n");

}		/*end plot_states2d*/
#endif /* defined(TWOD) */

#if defined(THREED)
/*ARGSUSED*/
EXPORT  void plot_states3d(
	FILE		*file,
	Front		*front,
	Wave		*wave,
	Printplot	*prt)
{
	if (prt->n_tri_vars == 0)
	    return;

	screen("ERROR in plot_states3d(), function not implemented\n");
	clean_up(ERROR);
}		/*end plot_states3d*/
#endif /* defined(THREED) */


LOCAL void plot_prostar_data(
	Wave		  *wave,
	Front		  *front,
	PROSTAR_plot_data *prostar_data)
{

	PROSTAR_PRINT_OPTIONS *opts = &PROSTAR_print_opts(prostar_data);
	char                  *v_file_name = prostar_data->v_file_name;
	char                  *c_file_name = prostar_data->c_file_name;
	char                  *p_file_name = prostar_data->p_file_name;
	char                  *base_file_name = prostar_data->base_file_name;
	FILE                  *vfile, *cfile, *pfile;
	RECT_GRID	      *gr = front->rect_grid;
	int		      nvars = prostar_num_vars(opts);
	int		      i, j, k, dim = prostar_data->dim;
	boolean                  vrt_cel_created = prostar_data->vrt_cel_created;
	int		      var, cell_count;
	int		      *pixels, *vertices;
	int                   width, length, height;
	double		      L[3], U[3];
	double		      coords[3];
	double		      *step=prostar_data->step;
	COMPONENT             comp;
	static Locstate	      state = NULL;
	double                 intfc_jump[1];
	double                 phys_data[6], temp[6];
	int		      stp = front->step;
	const char	      *vert_suffix, *cell_suffix, *phys_suffix;

	vert_suffix = ".vrt";
	cell_suffix = ".cel";
	phys_suffix = ".usr";

	if (debugging("PROSTAR"))
	    printf("Entering plot_prostar_data()\n");

	if (state == NULL)
	    alloc_state(front->interf,&state,front->sizest);

	pixels = prostar_pixels(opts);
	vertices = prostar_vertices(opts);

	for(i = 0; i < 3; i++)
	{
	    vertices[i] = pixels[i] + 1;
	    L[i] = prostar_L0(opts)[i];
	    L[i] = max(L[i],front->rect_grid->GL[i]);
	    U[i] = prostar_U0(opts)[i];	
	    U[i] = min(U[i],front->rect_grid->GU[i]);	
	}	

	for (i = 0; i < dim; i++)
	{
	  if (pixels[i] == 0)
	  {
	      if (debugging("PROSTAR"))
		  (void) printf("pixels[%d] = 0\n",i);
	      break;
	  }
	}

	width = pixels[0];
	length = pixels[1];
	height = pixels[2];

	for(i = 0; i < 3; i++)
	{ 
	    prostar_data->step[i] = prostar_len(opts)[i]/pixels[i];
	    step[i] = prostar_len(opts)[i]/pixels[i];	   	
	}
   	
	(void) sprintf(prostar_data->p_file_name,
		       "%s-ts%s%s",prostar_data->base_file_name,
		       right_flush(stp,TSTEP_FIELD_WIDTH),phys_suffix);

	if (vrt_cel_created == 0)
	{
	    (void) sprintf(prostar_data->v_file_name,
			   "%s%s", prostar_data->base_file_name,
			   vert_suffix);
	    (void) sprintf(prostar_data->c_file_name,
			   "%s%s",prostar_data->base_file_name,
			   cell_suffix);
	    
	    if ((vfile = fopen(v_file_name,"w")) == NULL)
	    {
		(void) screen("WARNING in plot_prostar_data(), "
			      "can't open %s\n",v_file_name);
		return;
	    }
	    if ((cfile = fopen(c_file_name,"w")) == NULL)
	    {
		(void) screen("WARNING in plot_prostar_data(), "
			      "can't open %s\n",c_file_name);
		return;
            }
           
	    
	    /* set up for vertices */
	    for (i = 0; i < 3; i++)
	    {
		uni_array(prostar_data->vert_scale+i+1,
		       prostar_vertices(opts)[i],FLOAT);
	    }
	    for (i = 0; i <= width; i++)
	        prostar_data->vert_scale[1][i] = L[0] + (i)*step[0];
	    for (j = 0; j <= length; j++)
	        prostar_data->vert_scale[2][j] = L[1] + (j)*step[1];
	    for (k = 0; k <= height; k++)
	        prostar_data->vert_scale[3][k] = L[2] + (k)*step[2];
	    
	    cell_count = 1;
	    for (k = 0; k <= height; k++)
	    {
		coords[2] = prostar_data->vert_scale[3][k];
		
		for (j = 0; j <= length; j++)
		{
		    coords[1] = prostar_data->vert_scale[2][j];
		    
		    for (i = 0; i <= width; i++)
	            {
			coords[0] = prostar_data->vert_scale[1][i];
			
			/* print vertex number */
			(void) fprintf(vfile,"%9d      ",
				       k*(width+1)*(length+1)+j*(width+1)+i+1);
			
			/* print coords of vertex */
			if (dim == 3)
		       	    (void) fprintf(vfile,"%16.9lf%16.9lf%16.9lf\n",
					 coords[0], coords[1], 
					 coords[2]);
			/* if dim = 2 need dummy coords for z coordinate */
			if (dim == 2)
			{
			    if (k == 0) 
			        (void) fprintf(vfile,"%16.9lf%16.9lf%16.9lf\n",
					     coords[0], coords[1], 
					     coords[2]);
			    if (k == 1) 
			        /* height will be default=.001 for dim=2*/
			        (void) fprintf(vfile,"%16.9lf%16.9lf%16.9lf\n",
					       coords[0], coords[1], 
					       coords[2]+.001);
			}
			/* print cell data */
			if ((i != width) && (j != length) && ( k != height))
			{
			    /* print cell number */
			    (void) fprintf(cfile,"%9d      ",cell_count++);
			    /* print the numbers of its vertices */
			    (void) fprintf(cfile,
					   "%9d%9d%9d%9d%9d%9d%9d%9d%9d%4d\n",
					   k*(width+1)*(length+1)+j*(width+1)+i+1,
					   k*(width+1)*(length+1)+j*(width+1)+i+2,
					   k*(width+1)*(length+1)+(j+1)*(width+1)
					   +i+2,
					   k*(width+1)*(length+1)+(j+1)*(width+1)
					   +i+1,
					   (k+1)*(width+1)*(length+1)+j*(width+1)
					   +i+1,
					   (k+1)*(width+1)*(length+1)+j*(width+1)
					   +i+2,
					   (k+1)*(width+1)*(length+1)+(j+1)*
					   (width+1)+i+2,
					   (k+1)*(width+1)*(length+1)+(j+1)*
					   (width+1)+i+1,1,1);
			}
		    }
		}
	    }
	    fclose(vfile);
	    fclose(cfile);
	    prostar_data->vrt_cel_created = 1;
	}
	
	if ((pfile = fopen(p_file_name,"w")) == NULL)
	{
	    (void) screen("WARNING in plot_prostar_data(), "
			  "can't open %s\n",p_file_name);
	    return;
	}

	/* set up for  physics, need coords of cell centers */
	for (i = 0; i < width; i++)
	    prostar_data->scale[1][i] = U[0] - (i+.5)*step[0];
	for (j = 0; j < length; j++)
	    prostar_data->scale[2][j] = L[1] + (j+.5)*step[1];
	for (k = 0; k < height; k++)
	    prostar_data->scale[3][k] = L[2] + (k+.5)*step[2];
	
	/* default data */
	for (var = 0; var < 6; var++)
	{
	    phys_data[var] = 0.0; /*default vales */
	    temp[var] = 0.0; 
	}
	cell_count = 1;
	for (k = 0; k < height; k++)
	{
	    coords[2] = prostar_data->scale[3][k];
 
	    for (j = 0; j < length; j++)
            {
		coords[1] = prostar_data->scale[2][j];

		for (i = 0; i < width; i++)
		{
		    coords[0] = prostar_data->scale[1][i];

		    /* find component at coords, needed for hyp_solution */
		    comp = component(coords,front->interf);
		    if (is_exterior_comp(comp,front->interf))
		    {
			comp = nearest_interior_comp(YES,NO_COMP,coords,
						      front->interf);
		    }

		    /* find state data at coords */
		    hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,
				 front,wave,state,NULL);

		    if ((dim == 2)&&(nvars >5))
		        nvars = 5;
		    if ((dim == 3)&&(nvars > 6))
		        nvars = 6;

		    for (var = 0; var < nvars; var++)
		    {
			phys_data[var] =
			  (PROSTAR_frame_plot_function(prostar_data->frame_data[var])
			    (coords,front,wave,comp,state));
		    }
		    if (dim == 2)
		    {
		        for (var = 2; var < nvars; var++)
			{
			    temp[var] = phys_data[var];
			}
			phys_data[2] = 0.0;
		        for (var = 3; var < nvars+1; var++)
			{
			    phys_data[var] = temp[var-1];
			}

		    }
		    (void) fprintf(pfile,"%9d      ",cell_count++);
		    (void) fprintf(pfile,
				   "%16.9lf%16.9lf%16.9lf%16.9lf%16.9lf%16.9lf\n",
				   phys_data[0], phys_data[1], phys_data[2],
				   phys_data[3], phys_data[4], phys_data[5]);
		}
	    }
	}
	fclose(pfile);
	(void) output();

	if (debugging("PROSTAR"))
	    printf("Leaving plot_prostar_data()\n");
}	/* end plot_prostar_data */

LOCAL	void	set_PROSTAR_plotting_range(
	Front		   *front,
	PROSTAR_plot_data  *prostar_data,
	int		   *pixels,
	int                *vertices,
	double		   *L,
	double		   *U)
{
	PROSTAR_PRINT_OPTIONS	*opts = &PROSTAR_print_opts(prostar_data);
	int	i, dim = prostar_data->dim;

	debug_print("PROSTAR","Entered set_PROSTAR_plotting_range().\n");
	for (i = 0; i < dim; i++)
	{
	    pixels[i] = prostar_pixels(opts)[i];
	    vertices[i] = prostar_vertices(opts)[i];
	    vertices[i] = pixels[i] + 1;

	    L[i] = prostar_L0(opts)[i];
	    L[i] = max(L[i],front->rect_grid->GL[i]);
	    U[i] = prostar_U0(opts)[i];	
	    U[i] = min(U[i],front->rect_grid->GU[i]);
	}
	debug_print("PROSTAR","Left set_PROSTAR_plotting_range().\n");
}		/*end set_PROSTAR_plotting_range*/

#if defined(__GD__)
LOCAL void plot_gd_data(
	Wave		  *wave,
	Front		  *front,
	GD_plot_data 	  *gd_data)
{
	GD_PRINT_OPTIONS *opts = &GD_print_opts(gd_data);
	int     	 var,nvars = gd_num_vars(opts);
	RECT_GRID        *gr = wave->rect_grid;
	Locstate	 state;
	COMPONENT	 comp;
	double		 *coords;
	int		 i,xmax,num_pts,icoords[MAXD];
	static	GD_GRAPH_DATA *gd_frame,**first_frame;
	GD_GRAPH_DATA	 *pframe;
	INTERFACE	 *intfc = front->interf;
	POINT		 **pp;

	if (debugging("GD_MOVIE"))
	    printf("Entering plot_gd_data()\n");

	xmax = gr->gmax[0];
	for (pp = intfc->points, num_pts = 0; pp && *pp; ++pp, ++num_pts) ;

	if (first_frame == NULL)
	{
	    uni_array(&first_frame,nvars,sizeof(GD_GRAPH_DATA*));
	    for (var = 0; var < nvars; ++var)
	    	first_frame[var] = NULL;
	}
	uni_array(&gd_frame,nvars,sizeof(GD_GRAPH_DATA));
	for (var = 0; var < nvars; ++var)
	{
	    if (first_frame[var] == NULL)
	    {
	    	first_frame[var] = pframe = gd_frame + var;
	    	first_frame[var]->prev = first_frame[var]->next = NULL;
	    }
	    else
	    {
	    	for (pframe = first_frame[var]; ; pframe = pframe->next)
		    if (pframe->next == NULL) break;
	    	pframe->next = gd_frame + var;
	    	gd_frame[var].prev = pframe;
	    	gd_frame[var].next = NULL;
		pframe = gd_frame + var;
	    }
	    uni_array(&pframe->xcrds,xmax,FLOAT);
	    uni_array(&pframe->ycrds,xmax,FLOAT);
	    uni_array(&pframe->points,num_pts,sizeof(GD_POINT));
	    pframe->imax = xmax;
	    pframe->num_points = num_pts;
	    for (pp = intfc->points, i = 0; pp && *pp; ++pp, ++i)
	    {
		coords = Coords(*pp);
		pframe->points[i].wave_type = wave_type(*pp);
		pframe->points[i].xcrds = coords[0];
		comp = negative_component(*pp);
	    	state = left_state(*pp);
		pframe->points[i].left_ycrds = (GD_frame_plot_function(
				gd_data->frame_opts[var])
	    			(coords,front,wave,comp,state));
		comp = positive_component(*pp);
	    	state = right_state(*pp);
		pframe->points[i].right_ycrds = (GD_frame_plot_function(
				gd_data->frame_opts[var])
	    			(coords,front,wave,comp,state));
	    }
	    for (i = 0; i < xmax; ++i)
	    {
	    	icoords[0] = i;
	    	coords = Rect_coords(icoords,wave);
	    	pframe->xcrds[i] = coords[0];
	    	comp = Rect_comp(icoords,wave);
	    	state = Rect_state(icoords,wave);
	    	pframe->ycrds[i] = (GD_frame_plot_function(
				gd_data->frame_opts[var])
	    			(coords,front,wave,comp,state));
	    }
	    pframe->time = front->time;
	}
	for (var = 0; var < nvars; ++var)
	{
	    GD_FRAME_OPTS *opts = gd_data->frame_opts + var;
	    char *filename = gd_plot_name(opts);
	    char *var_name = gd_var_name(opts);
	    gd_make_movie(first_frame[var],filename,var_name);
	}

	if (debugging("GD_MOVIE"))
	    printf("Leaving plot_gd_data()\n");
}	/* end plot_gd_data */

LOCAL	void gd_make_movie(
	GD_GRAPH_DATA	 *first_frame,
	char *filename,
	char *var_name)
{
	GD_GRAPH_DATA	 *pframe;
	int		 num_frame;
	int		 i,is,ip,iseg,imax,num_pts,num_sts;
	GD_POINT	 *points;
	double		 *x,*y;
	double		 xmin,xmax,ymin,ymax,height;
	char		 title[100];

	if (debugging("GD_MOVIE"))
	    printf("Entering gd_make_movie()\n");

	xmin = ymin =  HUGE_VAL;
	xmax = ymax = -HUGE_VAL;
	imax = first_frame->imax;

	uni_array(&x,imax,FLOAT);
	uni_array(&y,imax,FLOAT);

	for (pframe = first_frame; pframe != NULL; pframe = pframe->next)
	{
	    for (i = 0; i < imax; ++i)
	    {
	    	if (xmin > pframe->xcrds[i])
		    xmin = pframe->xcrds[i];
	    	if (ymin > pframe->ycrds[i])
		    ymin = pframe->ycrds[i];
	    	if (xmax < pframe->xcrds[i])
		    xmax = pframe->xcrds[i];
	    	if (ymax < pframe->ycrds[i])
		    ymax = pframe->ycrds[i];
	    }
	    num_pts = pframe->num_points;
	    for (i = 0; i < num_pts; ++i)
	    {
		if (pframe->points[i].wave_type == SUBDOMAIN_BOUNDARY)
		    continue;
	    	if (xmin > pframe->points[i].xcrds)
		    xmin = pframe->points[i].xcrds;
	    	if (xmax < pframe->points[i].xcrds)
		    xmax = pframe->points[i].xcrds;
	    	if (ymin > pframe->points[i].left_ycrds)
		    ymin = pframe->points[i].left_ycrds;
	    	if (ymax < pframe->points[i].left_ycrds)
		    ymax = pframe->points[i].left_ycrds;
	    	if (ymin > pframe->points[i].right_ycrds)
		    ymin = pframe->points[i].right_ycrds;
	    	if (ymax < pframe->points[i].right_ycrds)
		    ymax = pframe->points[i].right_ycrds;
	    }
	}
	height = ymax - ymin;
	ymin -= 0.15*height;
	ymax += 0.15*height;

	gd_initplot(filename,var_name,xmin,xmax,ymin,ymax,num_pts);

	for (pframe = first_frame; pframe != NULL; pframe = pframe->next)
	{
	    imax = pframe->imax;
	    num_pts = pframe->num_points;
	    points = pframe->points;
	    ip = 0;
	    is = 0;
	    iseg = 0;
	    while (is < imax)
	    {
		num_sts = 0;
	    	if (ip < num_pts)
		{
		    if (points[ip].wave_type != SUBDOMAIN_BOUNDARY &&
		        points[ip].xcrds < pframe->xcrds[is])
		    {
	    		x[num_sts] = points[ip].xcrds;
	    		y[num_sts] = points[ip].right_ycrds;
			++num_sts;
			++ip;
		    }
		    else if (points[ip].xcrds < pframe->xcrds[is])
			++ip;
		}
	    	for (i = is; i < imax; ++i)
		{
		    if (ip < num_pts && pframe->xcrds[i] >= points[ip].xcrds &&
		    	points[ip].wave_type != SUBDOMAIN_BOUNDARY)
		    {
		    	x[num_sts] = points[ip].xcrds;	
		    	y[num_sts] = points[ip].left_ycrds;	
			++num_sts;
			if (pframe->xcrds[i] == points[ip].xcrds) ++i;
			break;
		    }
		    else if (i == imax-1 && pframe->xcrds[i] < points[ip].xcrds
		    	&& points[ip].wave_type != SUBDOMAIN_BOUNDARY)
		    {
		    	x[num_sts] = points[ip].xcrds;	
		    	y[num_sts] = points[ip].left_ycrds;	
			++num_sts;
		    }
		    else
		    {
		    	x[num_sts] = pframe->xcrds[i];
	    		y[num_sts] = pframe->ycrds[i];
			++num_sts;
		    }
		}
		is = i;

	    	gd_plotdata(num_sts,x,y);
		++iseg;
	    }
	    sprintf(title,"Time = %f",pframe->time);
            gd_plotframe(title);
	}
	gd_closeplot();
	free_these(2,x,y);
	if (debugging("GD_MOVIE"))
	    printf("Leaving gd_make_movie()\n");
}	/* end gd_make_movie */
#endif /* defined(__GD__) */

/*
*			find_time_step():
*
*	Computes the new time-step, as large as allowed by the
*	CFL condition and parabolic restrictions. 
*/

EXPORT	double find_time_step(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	boolean		print_set_by)
{
	TIME_STEP_SET_BY TimeStepSetBy;
	double		 max_dt, CFL;
	double		 newdt;
	double	         coords[MAXD], fcrds[MAXD], wcrds[MAXD];
	int	         i, dim = front->rect_grid->dim;

	CFL = Time_step_factor(front);
	if (dim > 1)
	    CFL *= min(Front_spacing(front,GENERAL_WAVE),
		       Front_spacing(front,VECTOR_WAVE));

	newdt = HUGE_VAL;

		/* Maximum Interior Wave Speeds (previous time) */

	TimeStepSetBy = TIME_STEP_NOT_SET;
	if (wave->max_hyp_time_step)
	{
	    max_dt = (*wave->max_hyp_time_step)(wave,wcrds);
#if !defined(_AIX)
	    if (!finite(max_dt) || !finite(-max_dt))
		(void) printf("WARNING in find_time_step(), "
			      "max_hyp_time_step returns infinity\n");
	    else
#endif /* !defined(_AIX) */
	    {
	        max_dt *= CFL;
	        if (max_dt < newdt)
	        {
	    	    newdt = max_dt;
	    	    TimeStepSetBy = TIME_STEP_SET_BY_WAVE;
	    	    for (i = 0; i < dim; ++i)
	    	        coords[i] = wcrds[i];
	        }
	    }
	}

	    /* Maximum Front Speeds (previous time) */

	if (front->max_front_time_step)
	{
	    max_dt = (*front->max_front_time_step)(front,fcrds);
#if !defined(_AIX)
	    if ((!finite(max_dt) || !finite(-max_dt)) &&
	        (front->interf->num_points>0))
		(void) printf("WARNING in find_time_step(), "
			      "max_front_time_step returns infinity\n");
	    else
#endif /* !defined(_AIX) */
	    {
	        max_dt *= CFL;
	        if (max_dt < newdt)
	        {
	    	    newdt = max_dt;
	    	    TimeStepSetBy = TIME_STEP_SET_BY_FRONT;
	    	    for (i = 0; i < dim; ++i)
	    	        coords[i] = fcrds[i];
	        }
	    }
	}
	/* This becomes annoying when the time step was */
        /* forced to be very small to fit the print interval */
	/*
	if ((grid->dt != 0.0) && (grid->step > 0))
	{
	    if (newdt > Max_time_step_modification_factor(front)*grid->dt)
	    {
	        newdt = Max_time_step_modification_factor(front)*grid->dt;
	        TimeStepSetBy = TIME_STEP_SET_BY_PREVIOUS;
	        for (i = 0; i < dim; ++i)
	            coords[i] = HUGE_VAL;
	    }
	}
	*/

        if (debugging("limit_time_step"))
	{
	    if (newdt > dt_lim(grid))
	    {
	    	newdt = dt_lim(grid);
	    	TimeStepSetBy = TIME_STEP_SET_BY_USER_LIMIT;
	    	for (i = 0; i < dim; ++i)
	    	    coords[i] = HUGE_VAL;
	    }
	}

	newdt = comm_time_step(newdt,coords,dim,&TimeStepSetBy);

	if (!stop_run(grid,front,wave))
	{
	    double tmpdt;
	    tmpdt = nonphysics_timestep_reductions(grid,wave,front,prt,newdt);
	    if (tmpdt != newdt)
	    {
	    	newdt = tmpdt;
	    	TimeStepSetBy = TIME_STEP_SET_BY_PRINTING;
	    }
	}

	if (print_set_by == YES)
	{
	    (void) printf("\nTime step set by ");
	    switch (TimeStepSetBy)
	    {
	    case TIME_STEP_SET_BY_WAVE:
		print_general_vector("interior update at ",coords,dim,", ");
		print_general_vector("Maxsp(wave) = ",Maxsp(wave),dim,"\n");
		break;
	    case TIME_STEP_SET_BY_FRONT:
		print_general_vector("front update at ",coords,dim,", ");
		print_general_vector("Spfr(front) = ",Spfr(front),dim+1,"\n");
		break;
	    case TIME_STEP_SET_BY_PREVIOUS:
		(void) printf("%g times previous dt\n",
			Max_time_step_modification_factor(front));
		break;
	    case TIME_STEP_SET_BY_USER_LIMIT:
		(void) printf("user defined limit %g\n",dt_lim(grid));
		break;
	    case TIME_STEP_SET_BY_PRINTING:
		(void) printf("printing restriction\n");
		break;
	    case TIME_STEP_NOT_SET:
	    default:
		screen("ERROR in find_time_step(), time step not set\n");
		clean_up(ERROR);
	    }
	}
	return newdt;
}		/*end find_time_step*/

LOCAL const int NUM_LOOK_AHEAD = 3;/*TOLERANCE*/

LOCAL double check_for_ts_restriction(
	OUTPUT_DATA	*odata,
	double		time,
	double		dt)
{
	int		 i;
	double		newdt;

	if ((odata != NULL) && exact_time_output(Output_mode(odata)))
	{
	    for (i = 1; i <= NUM_LOOK_AHEAD; ++i)
	    {
	        if (is_print_time(time+dt*i,Print_control(odata)) == YES)
		{
		    newdt = Output_next_print_time(odata) - time;
		    if (newdt < dt) dt = newdt;
	        }
	    }
	}
	return dt; /* no restriction found */
}		/*end check_for_ts_restriction*/


/*
*			nonphysics_timestep_reductions():
*
*	Reduces the timestep based on non-physical considerations, allowing
*	the computation to reach (exactly) a particular real time.  The
*	restrictions currently include global print times, pause times, stop
*	times, and user print times.
*
*	This routine looks ahead NUM_LOOK_AHEAD steps to check for target
*	times.  This allows the step sizes up to the target to be adjusted
*	to be approximately equal, avoiding very tiny time steps.
*/

/*ARGSUSED*/
EXPORT double nonphysics_timestep_reductions(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	double		dt)
{
        double		newdt = dt;
	int		i;

	/* If newdt takes us past next printing time, then reduce. */
	for (i = 0; i < NUM_OUTPUT_FORMATS; ++i)
	{
	    newdt = check_for_ts_restriction(prt->main_output_data[i],
					     grid->time,newdt);
	}

	/* If newdt takes us past pause time, then reduce (maybe again). */
	if (exact_time_output(pause_mode(grid)))
	{
	    for (i = 1; i <= NUM_LOOK_AHEAD; ++i)
	    {
		if ((grid->time + newdt*i) >= pause_time(grid))
		{
		    newdt = (pause_time(grid) - grid->time)/i;
		    break;
	        }
	    }
	}

	if (exact_time_output(stop_time_mode(grid)))
	{
	    /* If newdt takes us past stop time, then reduce (maybe again). */
	    for (i = 1; i <= NUM_LOOK_AHEAD; ++i)
	    {
	        if ((grid->time + newdt*i) >= stop_time(grid))
	        {
	    	    newdt = (stop_time(grid) - grid->time)/i;
	    	    break;
	        }
            }
	}

	/* Check user printing routines also. */
	if (prt->user_output_funcs != NULL)
	{
	    int		i;
	    void	(**f)(Grid*,Wave*,Front*,Printplot*,
			      OUTPUT_DATA*,boolean);

	    for (i = 0, f = prt->user_output_funcs; *f != NULL; ++f, ++i)
	    {
	        newdt = check_for_ts_restriction(prt->user_output_data[i],
						 grid->time,newdt);
	    }
        }
	if (grid->user_dt_control)
	    (*grid->user_dt_control)(grid,&newdt);

	return newdt;
}		/*end nonphysics_timestep_reductions*/


/*
*			is_print_time():
*
*	Returns YES if time >= next printing time and NO otherwise.
*/

EXPORT	boolean is_print_time(
	double		time,
	PRINT_CONTROL	*prtctrl)
{
	PRINT_OPTIONS	*pto = &print_options(prtctrl);
	double npt;

	if (!real_time_output(prt_mode(pto)))
	    return NO;
	if (time < (print_start_time(pto) - MACH_EPS)) /*TOLERANCE*/
	    return NO;

	npt = next_print_time(prtctrl);
	npt *= (npt < 0) ?	(1.0 + print_time_tolerance(prtctrl)) :
				(1.0 - print_time_tolerance(prtctrl));
	return (time >= npt) ? YES : NO;
}		/*end is_print_time*/


/*
*			is_print_step():
*
*	Returns YES if step == next printing step and NO otherwise.
*/

EXPORT boolean is_print_step(
	int		step,
	PRINT_CONTROL	*prtctrl)
{
	PRINT_OPTIONS	*pto = &print_options(prtctrl);
	if (mesh_time_output(prt_mode(pto)) &&
	    (step >= print_start_step(pto)) &&
	    (step == next_print_step(prtctrl)))
	    return YES;
	else
	    return NO;
}		/*end is_print_step*/


/*
*				d_await():
*
*	Handles pauses in execution due to exceeding specified pause time.
*
*	Provides indefinite delay before resuming execution; if desired
*	the program may be killed.  It is possible to modify various paramaters
*	governing debugging and printout.
*/

/*ARGSUSED*/
EXPORT void d_await(
	Grid		*grid,
	Front		*front,
	Wave		*wave,
	Printplot	*prt,
	boolean		got_intfc_from_file)
{
	IMPORT boolean  suppress_prompts;
	boolean		sav_suppress_prompts;
	INIT_DATA	*init = grid->init;
	char		s[Gets_BUF_SIZE];

	sav_suppress_prompts = suppress_prompts;
	suppress_prompts = (interactive_prompting(init) == YES) ? NO : YES;

	screen("\n\nRequest continuation or termination of run "
	       "(cont(dflt),term): ");
	(void) Gets(s);
	if ((s[0] == 't') || (s[0] == 'T'))
	{
	    suppress_prompts = sav_suppress_prompts;
	    clean_up(0);
	}

	screen("\nCurrent Time = %g, current step = %d.\n",
	       grid->time,grid->step);
	screen("Specify the new pause time mode, "
	       "exact%s, constant%s, or mesh%s: ",
	       (pause_mode(grid) == EXACT_TIME) ? " (dflt)" : "",
	       (pause_mode(grid) == CONSTANT_TIME) ? " (dflt)" : "",
	       (pause_mode(grid) == MESH_TIME) ? " (dflt)" : "");
	(void) Gets(s);

	if (s[0] == 'm' || s[0] == 'M')
	    pause_mode(grid) = MESH_TIME;
	else if (s[0] == 'c' || s[0] == 'M')
	    pause_mode(grid) = CONSTANT_TIME;
	else if (s[0] == 'e' || s[0] == 'E')
	    pause_mode(grid) = EXACT_TIME;

	if (mesh_time_output(pause_mode(grid)))
	{
	    pause_step(grid) = INT_MAX;
	    screen("Enter the new Pause Time Step (dflt = %d): ",
		   pause_step(grid));
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscanf(s,"%d \n",&pause_step(grid));
	}
	else
	{
	    pause_time(grid) = HUGE_VAL;
	    screen("Enter the new Pause Time (dflt = %g): ",
		   pause_time(grid));
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscan_float(s,&pause_time(grid));
	}

	screen("Type 'm' to modify the debugging: ");
	(void) Gets(s);
	if ((s[0] == 'm') || (s[0] == 'M'))
	{
	    screen("Type 'init', 'debug', or 'nodebug' to initialize, "
	           "continue, or terminate\n\tthe debug mode: ");
	    (void) Gets(s);
	    if      ((s[0] == 'i') || (s[0] == 'I'))
	    	dbparams(prt->init) = init_debug(PROMPT_FOR_DEBUG);
	    else if ((s[0] == 'd') || (s[0] == 'D'))
	    	dbparams(prt->init) = init_debug(SOME);
	    else if ((s[0] == 'n') || (s[0] == 'N'))
	    	dbparams(prt->init) = init_debug(NONE);
	}

	if (debugging("vm_1") || debugging("vmalloc_1"))
	    vmalloc_debug(1);
	if (debugging("vm_2") || debugging("vmalloc_2"))
	    vmalloc_debug(2);
	if (debugging("vm_3") || debugging("vmalloc_3"))
	    vmalloc_debug(3);
	if (debugging("vm_4") || debugging("vmalloc_4"))
	    vmalloc_debug(4);
	if (debugging("vm_5") || debugging("vmalloc_5"))
	    vmalloc_debug(5);

	screen("Type 'm' to modify the printout: ");
	(void) Gets(s);
	if ((s[0] == 'm') || (s[0] == 'M'))
	{
	    init_printplot(prt->init,grid,front,prt);
	    grid->dt = nonphysics_timestep_reductions(grid,wave,
						      front,prt,grid->dt);
	}
	if (Init_redistribution_function(front) != NULL)
	{
	    screen("Type 'm' to modify the redistribution parameters: ");
	    (void) Gets(s);
	    if ((s[0] == 'm') || (s[0] == 'M'))
	    {
	        INTERFACE	*sav_intfc = i_intfc(init);

	        i_intfc(init) = front->interf;
		copy_redistribution_values(init,front);
	        prompt_for_redistribute(init);
	        Init_redistribution(init,front);
	        i_intfc(init) = sav_intfc;
	    }
	}
	suppress_prompts = sav_suppress_prompts;
}		/*end d_await*/

/*
*			give_peep():
*
*	Provides a peep at the solutions.
*/

/* ARGSUSED */
EXPORT void give_peep(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*odata,
	boolean		about_to_stop)
{
	(void) fprintf(stderr,"\n\n\nCurrent Time: %.4g   Time Step: %d\n",
				grid->time,grid->step);
	show_COMP(Output_file(odata),front->interf);
}		/*end give_peep*/


EXPORT	void d_print_stopping_criteria(Grid *grid)
{
	(void) printf("\t\t\tStopping Criteria\n");
	(void) printf("\tstop_time %g\tstop_step %d",
		      stop_time(grid),stop_step(grid));
}		/*end d_print_stopping_criteria*/

EXPORT	void	d_print_D_Front_structure(Front* fr)
{
	(void) printf("\n\n\n\t\tD_Front %p structure\n",(POINTER)fr);
	h_print_H_Front_structure(fr);
	(void) printf("amr info: chart_of_front(fr) = %p\n",
		      (POINTER)chart_of_front(fr));
	(void) printf("\n\n\n\t\tEnd D_Front %p structure\n",(POINTER)fr);
}		/*end d_print_D_Front_structure*/

LOCAL double comm_time_step(
	double			dt,
	double			*coords,
	int			dim,
	TIME_STEP_SET_BY	*TimeStepSetBy)
{
	int		numnodes = pp_numnodes();
	int		i, imin;
	static struct	_COMM_TIME_STEP {
					  double	           dt;
		                          double	           coords[MAXD];
		                          TIME_STEP_SET_BY TimeStepSetBy;
	                                } Cts;
	static	struct   _COMM_TIME_STEP *cts = NULL;

	if (numnodes == 1)
	    return dt;

	debug_print("time_step","Entered comm_time_step()\n");

	if (cts == NULL)
	    uni_array(&cts,numnodes,sizeof(struct _COMM_TIME_STEP));

	Cts.dt = dt;
	Cts.TimeStepSetBy = *TimeStepSetBy;
	for (i = 0; i < dim; ++i)
	    Cts.coords[i] = coords[i];

	pp_all_gather((POINTER)&Cts,sizeof(struct _COMM_TIME_STEP),
		      (POINTER)cts,sizeof(struct _COMM_TIME_STEP));

	if (debugging("time_step"))
	{
	    char buf[20];
	    for (i = 0; i < numnodes; ++i)
	    {
		(void) printf("cts[%d].dt = %g\n",i,cts[i].dt);
		(void) sprintf(buf,"cts[%d].coords = ",i);
		print_general_vector(buf,cts[i].coords,dim,"\n");
		(void) printf("cts[%d].TimeStepSetBy = %s\n",i,
			      TimeStepSetByString(cts[i].TimeStepSetBy));
	    }
	}

	for (imin = 0, i = 1; i < numnodes; ++i)
	    if (cts[imin].dt > cts[i].dt)
		imin = i;

	dt = cts[imin].dt;
	*TimeStepSetBy = cts[imin].TimeStepSetBy;
	for (i = 0; i < dim; ++i)
	    coords[i] = cts[imin].coords[i];

	debug_print("time_step","Left comm_time_step()\n");
	return dt;
}		/*end comm_time_step*/

LOCAL	const char *TimeStepSetByString(
	TIME_STEP_SET_BY TimeStepSetBy)
{
	switch (TimeStepSetBy)
	{
	case TIME_STEP_NOT_SET:
	    return "TIME_STEP_NOT_SET";
	case TIME_STEP_SET_BY_FRONT:
	    return "TIME_STEP_SET_BY_FRONT";
	case TIME_STEP_SET_BY_WAVE:
	    return "TIME_STEP_SET_BY_WAVE";
	case TIME_STEP_SET_BY_PREVIOUS:
	    return "TIME_STEP_SET_BY_PREVIOUS";
	case TIME_STEP_SET_BY_USER_LIMIT:
	    return "TIME_STEP_SET_BY_USER_LIMIT";
	case TIME_STEP_SET_BY_PRINTING:
	    return "TIME_STEP_SET_BY_PRINTING";
	default:
	    return NULL;
	}
}		/*end TimeStepSetByString*/

#if defined(USE_AMR)

EXPORT double find_amr_time_step(
        CHART        **rts,
        int          num_patches,
        Printplot    *prt,
        boolean         print_set_by)
{
        int          i;
        double        dt = HUGE_VAL, tmpdt;
        double        coords[MAXD], crds[MAXD];
        double        newdt;
        int          dim = rts[0]->front->rect_grid->dim;

        TIME_STEP_SET_BY TimeStepSetBy, tmpTimeStepSetBy;

        TimeStepSetBy = TIME_STEP_NOT_SET;

        for(i = 0; i < num_patches; i++)
        {
            if(rts[i]->wave->patch_level ==
                 rts[0]->amrparam->numberOfRefinementLevels-1 or
                 rts[i]->wave->patch_level == 0)
            {
                tmpdt = find_single_proc_time_step(rts[i]->grid,
                      rts[i]->wave,rts[i]->front,crds, &tmpTimeStepSetBy);
                if(tmpdt < dt)
                {
                    dt = tmpdt;
                    ftft_assign(coords,crds,sizeof(double)*MAXD);
                    TimeStepSetBy = tmpTimeStepSetBy;
                }
            }
        }
        newdt = comm_time_step(dt,coords,dim,&TimeStepSetBy);

        /* lack nonphysics_timestep_reductions() here */

        if (print_set_by == YES)
        {
            (void) printf("\nTime step set by ");
            switch (TimeStepSetBy)
            {
            case TIME_STEP_SET_BY_WAVE:
                print_general_vector("interior update at ",coords,dim,", ");
                print_general_vector("Maxsp(wave) = ",Maxsp(rts[0]->wave),dim,"\n");
                break;
            case TIME_STEP_SET_BY_FRONT:
                print_general_vector("front update at ",coords,dim,", ");
                print_general_vector("Spfr(front) = ",Spfr(rts[0]->front),dim+1,"\n");
                break;
            case TIME_STEP_SET_BY_PREVIOUS:
                (void) printf("%g times previous dt\n",
                        Max_time_step_modification_factor(rts[0]->front));
                break;
            case TIME_STEP_SET_BY_USER_LIMIT:
                (void) printf("user defined limit %g\n",dt_lim(rts[0]->grid));
                break;
            case TIME_STEP_SET_BY_PRINTING:
                (void) printf("printing restriction\n");
                break;
            case TIME_STEP_NOT_SET:
            default:
                screen("ERROR in find_time_step(), time step not set\n");
                clean_up(ERROR);
            }
        }
        return newdt;
}

LOCAL double find_single_proc_time_step(
        Grid             *grid,
        Wave             *wave,
        Front            *front,
        double            *coords,
        TIME_STEP_SET_BY *TimeStepSetBy)
{
        double            max_dt, CFL;
        double            newdt;
        double            fcrds[MAXD], wcrds[MAXD];
        int              i, dim = front->rect_grid->dim;

        CFL = Time_step_factor(front);
        if (dim > 1)
            CFL *= min(Front_spacing(front,GENERAL_WAVE),
                       Front_spacing(front,VECTOR_WAVE));

        newdt = HUGE_VAL;
                /* Maximum Interior Wave Speeds (previous time) */

        if (wave->max_hyp_time_step)
        {
            max_dt = (*wave->max_hyp_time_step)(wave,wcrds);
#if !defined(_AIX)
            if (!finite(max_dt) || !finite(-max_dt))
            {
#if defined(USE_AMR)
                if(front->rect_grid->dim == 2)
                {
                    CURVE **c;
                    int   interior = NO;
                    for(c = front->interf->curves; c and *c; c++)
                    {
                        if(wave_type(*c) >= FIRST_PHYSICS_WAVE_TYPE)
                        {
                            interior = YES; break;
                        }
                    }
                    if(YES == interior)
                    {
                        (void) printf("WARNING in find_time_step(), "
                              "max_front_time_step returns infinity\n");
                        (void) printf("From front[%d], level[%d]\n",
                               front->patch_number, front->patch_level);
                    }
                }
#else
                (void) printf("WARNING in find_time_step(), "
                              "max_hyp_time_step returns infinity\n");
#endif /* if defined(USE_AMR) */
            }
            else
#endif /* !defined(_AIX) */
            {
                max_dt *= CFL;
                if (max_dt < newdt)
                {
                    newdt = max_dt;
                    *TimeStepSetBy = TIME_STEP_SET_BY_WAVE;
                    for (i = 0; i < dim; ++i)
                        coords[i] = wcrds[i];
                }
            }
        }

            /* Maximum Front Speeds (previous time) */

        if (front->max_front_time_step)
        {
            max_dt = (*front->max_front_time_step)(front,fcrds);
#if !defined(_AIX)
            if ((!finite(max_dt) || !finite(-max_dt)) &&
                (front->interf->num_points>0))
                (void) printf("WARNING in find_time_step(), "
                              "max_front_time_step returns infinity\n");
            else
#endif /* !defined(_AIX) */
            {
                max_dt *= CFL;
                if (max_dt < newdt)
                {
                    newdt = max_dt;
                    *TimeStepSetBy = TIME_STEP_SET_BY_FRONT;
                    for (i = 0; i < dim; ++i)
                        coords[i] = fcrds[i];
                }
            }
        }
        if ((grid->dt != 0.0) && (grid->step > 0))
        {
            if (newdt > Max_time_step_modification_factor(front)*grid->dt)
            {
                newdt = Max_time_step_modification_factor(front)*grid->dt;
                *TimeStepSetBy = TIME_STEP_SET_BY_PREVIOUS;
                for (i = 0; i < dim; ++i)
                    coords[i] = HUGE_VAL;
            }
        }

        if (debugging("limit_time_step"))
        {
            if (newdt > dt_lim(grid))
            {
                newdt = dt_lim(grid);
                *TimeStepSetBy = TIME_STEP_SET_BY_USER_LIMIT;
                for (i = 0; i < dim; ++i)
                    coords[i] = HUGE_VAL;
            }
        }

        return newdt;
}
#endif /* defined(USE_AMR) */

LOCAL void	print_vtk_states(
	Wave            *wave,
	Front           *front,
	VTK_plot_data   *vtk_data)
{
       	char	dirname[1024];
	char buff[1024];
        VTK_PRINT_OPTIONS *opts = &VTK_print_opts(vtk_data);
	int     sub_factor = opts->subsample_factor;
	boolean    print_in_binary = opts->_print_in_binary;
        static  boolean verified = FALSE;
	int i;
	int dim = front->interf->dim;
	int step = front->step;

	sprintf(buff,"%s",vtk_base(opts));
	sprintf(dirname,"%s-vtk.ts%s",buff,right_flush(step,5));
	if (pp_numnodes() > 1)
	{
	    sprintf(buff,"%s",dirname);
	    sprintf(dirname,"%s-nd%s",buff,right_flush(pp_mynode(),4));
        }

	if(verified || sub_factor == 1 || sub_factor > 1 
	    && subsampling_is_valid(wave,front,sub_factor))
	{
	    if (dim == 1)
	    {
		screen("ERROR in print_vtk_states (dprint.c), code for 1D vtk "
			"output not present.");
		clean_up(ERROR);
	    }
	    if (dim == 2)
	    {
	        if(sub_factor != 1)
                    fill_vtk_values2d_sub(wave,front,vtk_data);
	        else
	            fill_vtk_values2d(wave,front,vtk_data);  

	        vtk_interface_plot(dirname,front->interf,print_in_binary,front->time,step,front->coordinate);
		print_intfc_visit_file(step,dim,opts);
	    }

	    if (dim == 3)
	    {
	        if(sub_factor != 1)
		    fill_vtk_values3d_sub(wave, front, vtk_data);
                else    
		    fill_vtk_values3d(wave,front,vtk_data);

                vtk_interface_plot(dirname,front->interf,print_in_binary,front->time,step,front->coordinate);
		print_intfc_visit_file(step,dim,opts);
	    }
	    verified = TRUE;
	}
	else
	    clean_up(ERROR);
}     	/* end print_vtk_states */

LOCAL   void    print_intfc_visit_file(
	int step,
	int dim,
        VTK_PRINT_OPTIONS *opts)
{
	char intfc_name[1024],buff[1024],buff2[1024];
	int i;
	int numnodes = pp_numnodes();
	int mynode = pp_mynode();
	FILE *intfc_bank;

	/* Open or create intfc.visit */
        sprintf(intfc_name,"%s.intfc.visit",vtk_base(opts));
	if(file_exists(intfc_name) && (numnodes == 1 || mynode == 1))
            intfc_bank = fopen(intfc_name, "a");
        else
        {
            if(numnodes > 1 && mynode == 1)
            {
                intfc_bank = fopen(intfc_name, "w");
                fprintf(intfc_bank, "!NBLOCKS %d\n", numnodes);
            }
            else if(numnodes == 1)
                intfc_bank = fopen(intfc_name, "w");
        }
		
	/* Write to intfc.visit */
        if(numnodes > 1 && mynode == 1)
        {
            for(i = 0; i < numnodes; i++)
            {
                 sprintf(buff,"%s-vtk.ts%s",opts->base_name,
					right_flush(step,5));
	        if (dim == 2)
                    sprintf(buff2,"%s-nd%s/2d-intfc.vtk\n",
					buff,right_flush(i,4));
		if (dim == 3)
                    sprintf(buff2,"%s-nd%s/3d-intfc.vtk\n",buff,
					right_flush(i,4));
                fprintf(intfc_bank, "%s", buff2);
            }
        }
        else if(numnodes == 1)
	    if (dim == 2)
                fprintf(intfc_bank,"%s-vtk.ts%s/2d-intfc.vtk\n",
					opts->base_name,right_flush(step,5));
	    if (dim == 3)
                fprintf(intfc_bank,"%s-vtk.ts%s/3d-intfc.vtk\n",
					opts->base_name,right_flush(step,5));
	
	/* Close intfc.visit */
        if(numnodes == 1 || mynode == 1)
	    fclose(intfc_bank);
}

LOCAL   void    print_variable_visit_file(
	int step,
	char name[512],
        VTK_PRINT_OPTIONS *opts)
{
	char param_name[1024],buff[1024],buff2[1024];
	int i;
	int numnodes = pp_numnodes();
	int mynode = pp_mynode();
	FILE *param_bank;

	/*Create the relevant entry in the relevant .visit file*/
	sprintf(param_name,"%s.%s.visit",vtk_base(opts),name);
	if(file_exists(param_name) && (numnodes == 1 || mynode == 1))
            param_bank = fopen(param_name, "a");
        else
        {
            if(numnodes > 1 && mynode == 1)
            {
                param_bank = fopen(param_name, "w");
                fprintf(param_bank, "!NBLOCKS %d\n", numnodes);
            }
            else if(numnodes == 1)
                param_bank = fopen(param_name, "w");
        }

        if(numnodes > 1 && mynode == 1)
        {
            for(i = 0; i < numnodes; i++)
            {
                sprintf(buff, "%s.%s.scalar.ts%s",opts->base_name,name,right_flush(step,5));
                sprintf(buff2, "%s-nd%s.vtk\n", buff,right_flush(i,4));
                fprintf(param_bank, "%s", buff2);
            }
        }
        else if(numnodes == 1)
            fprintf(param_bank,"%s.%s.scalar.ts%s.vtk\n",opts->base_name,name,right_flush(step,5));

        if(numnodes == 1 || mynode == 1)
	    fclose(param_bank);
}

LOCAL	void	print_vtk_header_data(
	char vtkname[4096],
	int width,
	int length,
	int height,
	double h0,
	double h1,
	double h2,
	double coords[3],
	double time,
	int step,
	char name[512],
	int factor,
	int dim,
	boolean print_in_binary)
{
	FILE *file;
	float value1[1], tmp1;
	int value2[1],tmp2;

        if(print_in_binary)
        {
            /* Print the vtk header for a binary file. */
            char str[128];
	    int length1, length2, length3;

            file = fopen(vtkname, "wb");
            sprintf(str,"# vtk DataFile Version 2.0\n");
            fwrite(str, sizeof(char), 27, file);
            sprintf(str, "comment line\n");
            fwrite(str, sizeof(char), 13, file);
            sprintf(str, "BINARY\n");
            fwrite(str, sizeof(char), 7, file);
            sprintf(str, "DATASET STRUCTURED_POINTS\n");
            fwrite(str, sizeof(char), 26, file);

            sprintf(str, "DIMENSIONS %d %d %d\n", width/factor + 1, length/factor + 1, height/factor + 1);
            length1 = count_digits(width/factor + 1);
            length2 = count_digits(length/factor + 1);
            length3 = count_digits(height/factor + 1);
            fwrite(str, sizeof(char), 14 + length1 + length2 + length3, file);

            sprintf(str, "%g", factor*h0);
            length1 = count_digits_in_string(str);
            sprintf(str, "%g", factor*h1);
            length2 = count_digits_in_string(str);
            sprintf(str, "%g", factor*h2);
            length3 = count_digits_in_string(str);
            sprintf(str,"SPACING %g %g %g\n", factor*h0, factor*h1, factor*h2);
            fwrite(str, sizeof(char), 11 + length1 + length2 + length3, file);
      
            sprintf(str, "%g",coords[0]);
            length1 = count_digits_in_string(str);
            sprintf(str, "%g",coords[1]);
            length2 = count_digits_in_string(str);
            sprintf(str, "%g",coords[2]);
            length3 = count_digits_in_string(str);
            sprintf(str, "ORIGIN %g %g %g\n",coords[0],coords[1],coords[2]);
            fwrite(str, sizeof(char),10 + length1 + length2 + length3 , file);

	    sprintf(str, "FIELD FieldData 2\n");
            fwrite(str, sizeof(char), 18, file);
	    sprintf(str, "TIME 1 1 float\n");
            fwrite(str, sizeof(char), 15, file);

	    value1[1];
	    tmp1 = time;
            if(hardware_is_little_endian())
	        value1[0] = endian_float_swap(tmp1);
            else
		value1[0] = tmp1;
	    fwrite(value1, sizeof(float), 1, file);

	    sprintf(str, "\n");
            fwrite(str, sizeof(char), 1, file);

	    sprintf(str, "CYCLE 1 1 int\n");
            fwrite(str, sizeof(char), 14, file);
	    #if defined(int)
            #undef int
            #define not_int
            #endif
	    tmp2 = step;
            if(hardware_is_little_endian())
	        value2[0] = endian_int_swap(tmp2);
            else
		value2[0] = tmp2;
	    fwrite(value2, sizeof(int), 1, file);
	    #if defined(not_int)
	    #define int double
     	    #undef not_int
	    #endif
	    sprintf(str, "\n");

            fwrite(str, sizeof(char), 1, file);
	    if (dim == 2)
	    {
                length1 = count_digits(width/factor*length/factor);
                sprintf(str, "CELL_DATA %d\n", width/factor*length/factor);
                fwrite(str, sizeof(char), 11 + length1, file);
	    }
	    if (dim == 3)
	    {
                length1 = count_digits(width/factor*length/factor*height/factor);
                sprintf(str, "CELL_DATA %d\n", width/factor*length/factor*height/factor);
                fwrite(str, sizeof(char), 11 + length1, file);
	    }
            length1 = count_digits_in_string(name);
            sprintf(str, "SCALARS %s float 1\n", name);
            fwrite(str, sizeof(char), 17 + length1, file);
            sprintf(str, "LOOKUP_TABLE default\n");
            fwrite(str, sizeof(char), 21, file);
        }
        else
        {
            /* Print the vtk header for an ASCII file. */
            file=fopen(vtkname,"w");
            (void) fprintf(file, "# vtk DataFile Version 2.0\n");
            (void) fprintf(file, "comment line\n");
            (void) fprintf(file, "ASCII\n");
            (void) fprintf(file, "DATASET STRUCTURED_POINTS\n");
            (void) fprintf(file, "DIMENSIONS %d %d %d\n", width/factor + 1, length/factor + 1, height/factor + 1);
            (void) fprintf(file, "SPACING %g %g %g\n", factor*h0, factor*h1, factor*h2);
            (void) fprintf(file, "ORIGIN %g %g %g\n", coords[0],coords[1],coords[2]);
	    (void) fprintf(file, "FIELD FieldData 2\n");
	    (void) fprintf(file, "TIME 1 1 float\n");
	    (void) fprintf(file, "%5.5e\n",time);
	    (void) fprintf(file, "CYCLE 1 1 int\n");
	    (void) fprintf(file, "%d\n",step);
	    if (dim == 2)
                (void) fprintf(file, "CELL_DATA %d\n", width/factor*length/factor);
	    if (dim == 3)
                (void) fprintf(file, "CELL_DATA %d\n", width/factor*length/factor*height/factor);
            (void) fprintf(file, "SCALARS %s float 1\n", name);
            (void) fprintf(file, "LOOKUP_TABLE default\n");
        }
	
	fclose(file);
}

LOCAL	void	fill_vtk_values2d(
	Wave		*wave,
	Front		*front,
	VTK_plot_data	*vtk_data)

{
       	COMPONENT	comp;
        VTK_PRINT_OPTIONS *opts = &VTK_print_opts(vtk_data);
        int             nvars = vtk_num_vars(opts);
	int		width,length,var,i,j,icoords[3];
	double		coords[3];
	double		*coordss;
	char		buff[4096],buff2[4096],vtkname[4096],name[512];
	boolean            print_in_binary = opts->_print_in_binary;
	RECT_GRID       *gr = wave->rect_grid;
	Locstate	state;
	FILE		*file;

	width = gr->gmax[0];
	length =  gr->gmax[1];

	for (var = 0; var < nvars; ++var)
	{
	icoords[0]=0;
	icoords[1]=0;
	icoords[2]=0;	
	coordss = Rect_coords(icoords,wave);
	coords[0]=coordss[0]-.5*gr->h[0];
	coords[1]=coordss[1]-.5*gr->h[1];;
	coords[2]=0;

	sprintf(name,"%s",VTK_frame_plot_name(vtk_data->frame_data[var]));
	sprintf(buff2,"%s.%s.scalar.ts%s",vtk_base(opts),name,right_flush(front->step,5));
	if (pp_numnodes() > 1)
	{
	    sprintf(buff,"%s-nd%s",buff2,right_flush(pp_mynode(),4));
	    sprintf(vtkname,"%s.vtk",buff);
	}
	else
	    sprintf(vtkname,"%s.vtk",buff2);

	print_vtk_header_data(vtkname,width,length,0,gr->h[0],gr->h[1],0.0,
				coords,front->time,front->step,name,1,2,
				print_in_binary);

	if (print_in_binary)
	    file = fopen(vtkname, "ab");
	else
	    file=fopen(vtkname,"a");

   	for (j = 0; j < length; ++j)
	{
	    for (i = 0; i < width; ++i)
	    {
		float tmp;
                icoords[1] = j; icoords[0] = i; 
                coords[0] = Rect_coords(icoords,wave)[0];
                coords[1] = Rect_coords(icoords,wave)[1];
                comp = Rect_comp(icoords,wave); 
                state = Rect_state(icoords,wave);

		tmp= (float)(*VTK_frame_plot_function(
				vtk_data->frame_data[var]))(coords,front,wave,
	      		        comp,state);

		if(print_in_binary)
		{

		    float vals[1];
                    if(hardware_is_little_endian())
		        vals[0] = endian_float_swap(tmp);
                    else
		        vals[0] = tmp;
		    fwrite(vals, sizeof(float), 1, file);
		}
                else
                    fprintf(file,"%f ",tmp);
	    }
	}

	fclose(file);
	print_variable_visit_file(front->step,name,opts);
      	}	
}	/* end fill_vtk_values2d */

LOCAL   void    fill_vtk_values2d_sub(
        Wave            *wave,
        Front           *front,
        VTK_plot_data   *vtk_data)

{
        COMPONENT       comp;
        VTK_PRINT_OPTIONS *opts = &VTK_print_opts(vtk_data);
        int             width,length,var,i,j,icoords[3];
        int             nvars = vtk_num_vars(opts);
        int             factor = opts->subsample_factor;
        double           coords[3];
        double           *coordss;
	char		buff[4096],buff2[4096],name[512],vtkname[4096];
        boolean            print_in_binary = opts->_print_in_binary;
        RECT_GRID       *gr = wave->rect_grid;
        Locstate        state;
        FILE            *file;

        width = gr->gmax[0];
        length =  gr->gmax[1];

        for (var = 0; var < nvars; ++var)
        {
        icoords[0]=0;
        icoords[1]=0;
        icoords[2]=0;
        coordss = Rect_coords(icoords,wave);
        coords[0]=coordss[0]-.5*gr->h[0];
        coords[1]=coordss[1]-.5*gr->h[1];;
        coords[2]=0;

        sprintf(name,"%s",VTK_frame_plot_name(vtk_data->frame_data[var]));
        sprintf(buff,"%s.%s.scalar.ts%s",vtk_base(opts),name,right_flush(front->step,5));
        if (pp_numnodes() > 1)
        {
	    sprintf(buff2,"%s-nd%s",buff,right_flush(pp_mynode(),4));
	    sprintf(vtkname,"%s.vtk",buff2);
	}
	else
            sprintf(vtkname,"%s.vtk",buff);

	print_vtk_header_data(vtkname,width,length,0,gr->h[0],gr->h[1],0.0,
				coords,front->time,front->step,name,factor,
				2,print_in_binary);

	if (print_in_binary)
	    file = fopen(vtkname, "ab");
	else
	    file=fopen(vtkname,"a");

        for (j = 0; j < length; j = j + factor)
        {
            for (i = 0; i < width; i = i + factor)
            {
		float tmp;
                icoords[1] = j; icoords[0] = i;
                coords[0] = Rect_coords(icoords,wave)[0];
                coords[1] = Rect_coords(icoords,wave)[1];
                comp = Rect_comp(icoords,wave);
                state = Rect_state(icoords,wave);

		tmp= (float)(*VTK_frame_plot_function(
				vtk_data->frame_data[var]))(coords,front,wave,
                                comp,state);

		if(print_in_binary)
                {

		    float vals[1];
		    if(hardware_is_little_endian())
                        vals[0] = endian_float_swap(tmp);
                    else
		        vals[0] = tmp;
		    fwrite(vals, sizeof(float), 1, file);
		}
                else
                    fprintf(file,"%f ",tmp);
            }
        }

	fclose(file);
	print_variable_visit_file(front->step,name,opts);
        }
}       /* end fill_vtk_values2d_sub */


LOCAL	void	fill_vtk_values3d(
	Wave		*wave,
	Front		*front,
	VTK_plot_data	*vtk_data)

{
      	COMPONENT	comp;
        VTK_PRINT_OPTIONS *opts = &VTK_print_opts(vtk_data);
	int		width,length,height,var,i,j,k,icoords[3];
        int             nvars = vtk_num_vars(opts);
	double		coords[3];
	double		*coordss;
	char		name[512],vtkname[4096],buff[4096],buff2[4096];
	boolean            print_in_binary = opts->_print_in_binary;
	RECT_GRID       *gr = wave->rect_grid;
	Locstate	state;
	FILE		*file;

	width = gr->gmax[0];
	length =  gr->gmax[1];
	height =  gr->gmax[2];
	for (var = 0; var < nvars; ++var)
	{
	icoords[0]=0;
	icoords[1]=0;
	icoords[2]=0;	
	coordss = Rect_coords(icoords,wave);
	coords[0]=coordss[0]-.5*gr->h[0];
	coords[1]=coordss[1]-.5*gr->h[1];
	coords[2]=coordss[2]-.5*gr->h[2];
	
	sprintf(name,"%s",VTK_frame_plot_name(vtk_data->frame_data[var]));
        sprintf(buff,"%s.%s.scalar.ts%s",vtk_base(opts), name, right_flush(front->step,5));

	if (pp_numnodes() > 1)
	{
	    sprintf(buff2,"%s-nd%s",buff,right_flush(pp_mynode(),4));
	    sprintf(vtkname,"%s.vtk",buff2);
	}
	else
	    sprintf(vtkname,"%s.vtk",buff);

	print_vtk_header_data(vtkname,width,length,height,gr->h[0],gr->h[1],
				gr->h[2],coords,front->time,front->step,name,
				1,3,print_in_binary);

	if (print_in_binary)
	    file = fopen(vtkname, "ab");
	else
	    file=fopen(vtkname,"a");

   	for (k = 0; k < height; ++k)
	{
   	    for (j = 0; j < length; ++j)
	    {
	    	for (i = 0; i < width; ++i)
	    	{
		    float tmp;
                    icoords[2] = k; 
                    icoords[1] = j; 
		    icoords[0] = i; 
                    coords[0] = Rect_coords(icoords,wave)[0];
                    coords[1] = Rect_coords(icoords,wave)[1];
                    coords[2] = Rect_coords(icoords,wave)[2];
                    comp = Rect_comp(icoords,wave); 
                    state = Rect_state(icoords,wave);

		    tmp = (float)(*VTK_frame_plot_function(
				vtk_data->frame_data[var]))(coords,front,
				wave,comp,state);

		    if(print_in_binary)
                    {

			float vals[1];
			if(hardware_is_little_endian())
			    vals[0] = endian_float_swap(tmp);
                        else
			    vals[0] = tmp;
			fwrite(vals, sizeof(float), 1, file);
		    }
                    else
                        fprintf(file,"%f ",tmp);
	    	}
	    }
	}

	print_variable_visit_file(front->step,name,opts);

	fclose(file);
	}	
}	/* end fill_vtk_values3d */

LOCAL   void    fill_vtk_values3d_sub(
        Wave            *wave,
        Front           *front,
        VTK_plot_data   *vtk_data)

{
        COMPONENT       comp;
        VTK_PRINT_OPTIONS *opts = &VTK_print_opts(vtk_data);
        int             width,length,height,var,i,j,k,icoords[3];
        int             nvars = vtk_num_vars(opts);
	int             factor = opts->subsample_factor;
        double           coords[3];
        double           *coordss;
        char            buff[4096],buff2[4096],name[100000],vtkname[4096];
        boolean            print_in_binary = opts->_print_in_binary;
        RECT_GRID       *gr = wave->rect_grid;
        Locstate        state;
        FILE            *file;

        width = gr->gmax[0];
        length =  gr->gmax[1];
        height =  gr->gmax[2];
        for (var = 0; var < nvars; ++var)
        {
        icoords[0]=0;
        icoords[1]=0;
        icoords[2]=0;
        coordss = Rect_coords(icoords,wave);
        coords[0]=(coordss[0]-.5*gr->h[0]);
        coords[1]=(coordss[1]-.5*gr->h[1]);
        coords[2]=(coordss[2]-.5*gr->h[2]);

        sprintf(name,"%s",VTK_frame_plot_name(vtk_data->frame_data[var]));
        sprintf(buff,"%s.%s.scalar.ts%s",vtk_base(opts), name,right_flush(front->step,5));

        if (pp_numnodes() > 1)
	{
            sprintf(buff2,"%s-nd%s",buff,right_flush(pp_mynode(),4));
            sprintf(vtkname,"%s.vtk",buff2);
        }
	else
            sprintf(vtkname,"%s.vtk",buff);

	print_vtk_header_data(vtkname,width,length,height,gr->h[0],gr->h[1],
				gr->h[2],coords,front->time,front->step,name,
				factor,3,print_in_binary);

	if (print_in_binary)
	    file = fopen(vtkname, "ab");
	else
	    file=fopen(vtkname,"a");

        for (k = 0; k < height; k = k + factor)
        {
            for (j = 0; j < length; j = j + factor)
            {
                for (i = 0; i < width; i = i + factor)
                {
		    float tmp;
                    icoords[2] = k;
                    icoords[1] = j;
                    icoords[0] = i;
                    coords[0] = Rect_coords(icoords,wave)[0];
                    coords[1] = Rect_coords(icoords,wave)[1];
                    coords[2] = Rect_coords(icoords,wave)[2];
                    comp = Rect_comp(icoords,wave);
                    state = Rect_state(icoords,wave);

		    tmp = (float)(*VTK_frame_plot_function(
				vtk_data->frame_data[var]))(coords,front,wave,
				comp,state);

		    if(print_in_binary)
                    {
			float vals[1];
			if(hardware_is_little_endian())
                            vals[0] = endian_float_swap(tmp);
			else
			    vals[0] = tmp;
                        fwrite(vals, sizeof(float), 1, file);
		    }
                    else
                        fprintf(file,"%f\n",tmp);
                }
            }
        }

	print_variable_visit_file(front->step,name,opts);

        fclose(file);
        }
}       /* end fill_vtk_values3d_sub */
/* end needed for VTK */

boolean subsampling_is_valid(Wave *wave, Front *front, int factor)
{
        RECT_GRID *gr = wave->rect_grid;
	int dim = front->interf->dim;
	int lengths[3];
	
	if(dim == 2)
	{
	    lengths[0] = gr->gmax[0];
	    lengths[1] = gr->gmax[1];

	    if(lengths[0]%factor == 0 && lengths[1]%factor == 0)
	        return TRUE;
	    return FALSE;
	}
	else
	{  
            lengths[0] = gr->gmax[0];
            lengths[1] = gr->gmax[1];
            lengths[2] = gr->gmax[2];
            
	    if(lengths[0]%factor == 0 && lengths[1]%factor == 0 && lengths[2]%factor == 0)
                return TRUE;
            return FALSE;
	}
}       /* end subsampling_is_valid */

int count_digits_in_string(char* str)
{
        int i = 0;
        while(str[i] != '\0')
            i++;
        return i;
}       /* end count_digits_in_string */

boolean file_exists(const char * filename)
{
       FILE *file;
       if (file = fopen(filename, "r"))
       {
       	   fclose(file);
	   return TRUE;
       }
       return FALSE;
}      /* end file_exists */
