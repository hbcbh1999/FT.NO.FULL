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
*				diprt.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains the driver initializing routines, as called in main,
*	and the principle or default initializing subroutines.
*/


#include <driver/ddecs.h>

	/* LOCAL Function Declarations */
LOCAL	int	count_selected_plot_choices(char*,Plot_choice*);
LOCAL	void	prompt_for_print_options(INIT_DATA*,PRINT_OPTIONS*,
					 boolean,boolean,const char*);
#if defined(THREED)
LOCAL	void	print_geomview_plots(Grid*,Wave*,Front*,Printplot*,
				     OUTPUT_DATA*,boolean);
#endif /* defined(THREED) */

#if defined(USE_HDF)
#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
LOCAL	  double identity_filter(double);
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */
LOCAL	void	prompt_for_HDF_frame_data(HDF_PRINT_OPTIONS*,HDF_FRAME_OPTS*,
					  ScalarPlotItem*);
LOCAL	void	prompt_for_compression_type(HDF_DATA_TYPE,HDF_PRINT_OPTIONS*);
LOCAL	void	set_HDF_frame_data(HDF_plot_data*,int);
LOCAL	void	turn_on_hdf_append(HDF_plot_data*);
#endif /* defined(USE_HDF) */
/* needed for VTK */
LOCAL   void    prompt_for_vtk_frame_data(VTK_PRINT_OPTIONS*, VTK_FRAME_OPTS*, ScalarPlotItem*);
/* end needed for VTK */
LOCAL	void	prompt_for_PROSTAR_frame_data(PROSTAR_PRINT_OPTIONS*,
					      PROSTAR_FRAME_OPTS*,
					      ScalarPlotItem*);

/*
*			d_prompt_for_printing_options():
*
*	Prompts for printing and plotting options.
*/

EXPORT	void	d_prompt_for_printing_options(
	INIT_DATA	*init)
{
	char       *c, s[Gets_BUF_SIZE];
	char       bname[1024];
	static const char *seps = " \t";
	int        i, dim = i_intfc(init)->dim;
	char       *mybname, *mydname;

	/*Set default printing options*/
	set_defaults_for_print_options(&Default_print_options(init),init);
	for (i = 0; i < NUM_OUTPUT_FORMATS; ++i)
	    prt_opts(init)[i] = Default_print_options(init);

	screen("\n\t\tPrinting Control\n\n"
	       "Request main output format(s). Options are\n"
	       "\t\tfront_plots only (F)\n"
	       "\t\tfront_plots plus interior_states (Restart format) (R)\n");
	if (dim == 2)
	{
	   screen("\t\tfront_plots plus tri_plots (T)\n"
	          "\t\tfront_plots, interior_states, and tri_plots (A)\n");
#if defined(USE_HDF)
	   screen("\t\tHDF raster plots (H)\n");
#endif /* defined(USE_HDF) */
	}
#if defined(USE_HDF)
	screen("\t\tSDS files (S)\n");
#endif /* defined(USE_HDF) */
	/* needed for VTK */
	screen("\t\tVTK files (V)\n");
	/* end needed for VTK */
	screen("\t\tPROSTAR plots (P)\n");
#if defined(__GD__)
	if (dim == 1)
	    screen("\t\tGD 1-D Movie (G)\n");
#endif /* defined(__GD__) */
	screen("\t\tsuppress output (dflt)\n"
	       "\tEnter the choices as a space separated list: ");
	(void) Gets(s);

	if (s[0] != '\0')
	{
	    for (c = s; *c != '\0'; ++c)
		*c = tolower(*c);
	    for (c = strtok(s,seps); c != NULL; c = strtok(NULL,seps))
	    {
		switch (c[0])
		{
		case 'f':			/* Frontplots only */
		    Prt_mode(prt_opts(init)[PRT_FRONTS]) = SET_PRT_MODE;
		    break;
		case 'r':			/* Restart output */
		    Prt_mode(prt_opts(init)[PRT_FRONTS]) = SET_PRT_MODE;
		    Prt_mode(prt_opts(init)[RECT_STATES]) = SET_PRT_MODE;
		    break;
		case 't':			/* Tri plot output */
		    Prt_mode(prt_opts(init)[PRT_FRONTS]) = SET_PRT_MODE;
		    Prt_mode(prt_opts(init)[TRI_STATES]) = SET_PRT_MODE;
		    break;
		case 'a':			/* All three above */
		    Prt_mode(prt_opts(init)[PRT_FRONTS]) = SET_PRT_MODE;
		    Prt_mode(prt_opts(init)[RECT_STATES]) = SET_PRT_MODE;
		    if (dim == 2)
		        Prt_mode(prt_opts(init)[TRI_STATES]) = SET_PRT_MODE;
		    break;
		case 'h':			/* HDF raster output */
		    if (dim == 2)
		        Prt_mode(prt_opts(init)[HDF_STATES]) = SET_PRT_MODE;
		    break;
		/* needed for VTK */
	        case 'v':                       /* vtk output */
		        Prt_mode(prt_opts(init)[VTK_STATES]) = SET_PRT_MODE;
		    break;
		/*needed for VTK */

		case 's':			/* HDF raster output */
		    Prt_mode(prt_opts(init)[SDS_STATES]) = SET_PRT_MODE;
		    break;
		case 'p':			/* HDF raster output */
		    Prt_mode(prt_opts(init)[PROSTAR_STATES]) = SET_PRT_MODE;
		    break;
#if defined(__GD__)
		case 'g':			/* HDF raster output */
		    Prt_mode(prt_opts(init)[GD_MOVIE]) = SET_PRT_MODE;
		    break;
#endif /* defined(__GD__) */
		default:
		    break;
		}
	    }
	}

	Prt_mode(Default_print_options(init)) = MESH_TIME;
	for (i = 0; i < NUM_OUTPUT_FORMATS; ++i)
	{
	    if (Prt_mode(prt_opts(init)[i]) == PRINTING_OFF)
		continue;

	    Prt_mode(prt_opts(init)[i]) = MESH_TIME;

	    screen("Prompt for %s printing control.\n",output_format_name(i));
	    prompt_for_print_options(init,prt_opts(init)+i,NO,NO,NULL);
	}
	set_binary_output(NO);	/* default for miscellaneous output */

	geomview_opts(init) = NULL;

	if (i_intfc(init)->dim > 1)
	{
	    screen("To obtain geomview plots type y (default = no): ");
	    (void) Gets(s);

	    if (s[0] == 'y')
	    {
	        screen("Enter a root name for the output directory: ");
	        (void) Gets(s);
	    
	        if (s[0] != '\0')
		    mybname = s;
	        else
	        {
		    base_and_dir_name(output_filename(init), &mydname, &mybname);
                    sprintf(s, "%s/gv%s",dirname(mydname),mybname);
	        }
                printf("geomview output directory is:%s\n",s);
  	        scalar(&geomview_opts(init),sizeof(PRINT_OPTIONS));
	        strcpy(print_directory(geomview_opts(init)),s);
	        strcpy(print_filename(geomview_opts(init)),"NONE");
	        prompt_for_print_options(init,geomview_opts(init),NO,NO,
					 "geomview data");
	    }
	}

	if (output_filename(init)[0] != '\0')
	    (void) sprintf(bname,"%s.lastdump",output_filename(init));
	else
	    (void) sprintf(bname,"lastdump");
	screen("The user can request that restart dumps be "
	       "printed at a specified wall\n\t"
	       "time interval.  These dumps will be named\n\t");
	screen("%s0%s and %s1%s\n\t",
	       bname,(pp_numnodes() > 1) ? "-nd\"node#\"" : "",
	       bname,(pp_numnodes() > 1) ? "-nd\"node#\"" : "");
	screen("and will be alternately overwritten as the run proceeds\n\t");
	screen("The wall time dump frequency can be given in units of "
	       "seconds,\n\t"
	       "minutes (default), or hours.  Indicate the units in the "
	       "obvious way\n\t"
	       "such as 30 minutes,  2 hours, etc.\n");
	screen("To request this option enter the wall time print frequency: ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    PRINT_OPTIONS *pto;
	    double	  freq;
	    char	  *c;

	    scalar(&pto,sizeof(PRINT_OPTIONS));
	    wall_time_dump_opts(init) = pto;
	    c = strtok(s,seps);
	    if (sscan_float(c,&freq) != 1)
	    {
		screen("ERROR in d_prompt_for_printing_options(), "
		       "invalid entry of wall time print frequency\n");
		clean_up(ERROR);
	    }
	    c = strtok(NULL,seps);
	    if (c != NULL)
	    {
		switch (c[0])
		{
		case 's': /*Seconds*/
		case 'S':
		    break;

		case 'h': /*Hours*/
		case 'H':
		    freq *= 3600.0;
		    break;
		    
		case 'd': /*Days*/
		case 'D':
		    freq *= 86400.0;
		    break;

		case 'm': /*Minutes*/
		case 'M':
		default:
		    freq *= 60.0;
		    break;
		}
	    }
	    else
		freq *= 60.0;
	    prt_mode(pto) = WALL_TIME;
	    print_wall_time_interval(pto) = freq;
	    last_wall_time_dump(pto) = real_time();
	    wall_time_dump_number(pto) = 0;
	    print_in_binary(pto) = Print_in_binary(Default_print_options(init));
	    screen("Print wall time dumps in binary (default = %s): ",
		   (print_in_binary(pto)==YES) ? "yes" : "no");
	    (void) Gets(s);
	    switch (s[0])
	    {
	    case 'y':
	    case 'Y':
		print_in_binary(pto) = YES;
		break;
	    case 'n':
	    case 'N':
		print_in_binary(pto) = NO;
		break;
	    default:
		break;
	    }
	    (void) strcpy(print_filename(pto),bname);
	}
	else
	    wall_time_dump_opts(init) = NULL;
}		/*end d_prompt_for_printing_options*/


EXPORT	void	set_defaults_for_print_options(
	PRINT_OPTIONS	*pto,
	INIT_DATA       *init)
{
	prt_mode(pto) = PRINTING_OFF;
	print_time_interval(pto) = 1.0;
	print_start_time(pto) = initial_time(init);
	print_step_interval(pto) = 1;
	print_start_step(pto) = 0;
	print_wall_time_interval(pto) = -HUGE_VAL;
	last_wall_time_dump(pto) = -HUGE_VAL;
	(void) strcpy(print_directory(pto),"");
	(void) strcpy(print_filename(pto),"");
	print_in_binary(pto) = is_binary_output();
}		/*end set_defaults_for_print_options*/

/*
*			d_init_printplot():
*
*	Initializes the non-physics-dependent parts of the printing and
*	plotting structure.  Also called by await().
*/

/*ARGSUSED*/
EXPORT void d_init_printplot(
	INIT_DATA	*init,
	Grid		*grid,
	Front		*front,
	Printplot	*prt)
{
	int			i;
	OUTPUT_DATA		**mod = prt->main_output_data;
	OUTPUT_DATA		**store_mod = prt->store_main_output_data;
	static boolean		first = YES;

	debug_print("Printplot","Entered d_init_printplot()\n");

	if (first == YES)
	{
	    first = NO;
	    prt->user_output_funcs = NULL;
	    prt->user_output_data = NULL;
	    prt->_print_states[0] = d_print_states1d;
	    prt->_print_states[1] = d_print_states;
	    prt->_print_states[2] = d_print_states;
	}

	prt->printout = d_printout;
	prt->print_initial_data = d_print_initial_data;
	prt->initialize_for_printout = NULL;
	prt->plot_ellip_vars = NULL;
	prt->_print_extreme_values = NULL;

	for (i = 0; i < NUM_OUTPUT_FORMATS; ++i)
	{
	    if (store_mod[i] == NULL)
	    {
	    	scalar(&store_mod[i],sizeof(OUTPUT_DATA));
	    }
	    mod[i] = (Prt_mode(prt_opts(init)[i]) != PRINTING_OFF) ?
	        store_mod[i] : NULL;
	}

	for (i = 0; i < NUM_OUTPUT_FORMATS; ++i)
	{
	    if (output_format_on(i,prt))
	        set_output_data(prt_opts(init)+i,mod[i],grid,prt,NO,NO);
	}

#if defined(THREED)
	if (geomview_opts(init) != NULL)
	{
	    OUTPUT_DATA	*data;

	    scalar(&data,sizeof(OUTPUT_DATA));
	    set_output_data(geomview_opts(init),data,grid,prt,NO,NO);
	    add_user_output_function(print_geomview_plots,data,prt);
	}
#endif /* defined(THREED) */

	if (wall_time_dump_opts(init) != NULL)
	    wall_time_dump_options(prt) = wall_time_dump_opts(init);
	
	if (debugging("Printplot"))
	{
	    print_Printplot_structure(prt);
	}
	debug_print("Printplot","Left d_init_printplot()\n");
}		/*end d_init_printplot*/


/*
*			set_next_print_time():
*
*	This function is responsible for computing the next print time/step.
*	See comments before prompt_for_print_options().
*/

EXPORT void set_next_print_time(
	Grid		*grid,
	PRINT_CONTROL	*ptc)
{
	PRINT_OPTIONS	*pto = &print_options(ptc);
	if (real_time_output(prt_mode(pto)))
	{
	    /* Check for potential infinite loop below. */

	    if (print_time_interval(pto) <= 0.0)
	    {
	    	screen("ERROR in set_next_print_time(), "
	    	       "invalid real time print interval: 0.0\n");
	    	clean_up(ERROR);
	    }

	    next_print_time(ptc) = print_start_time(pto);

	    /* This is needed because at initialization, the
	     * solution has not yet been printed, whereas every
	     * other time we get to this point it has. */

	    if (grid->step > 0)
	        while (next_print_time(ptc) <= grid->time)
	    	next_print_time(ptc) += print_time_interval(pto);

	    if (debugging("prt_rst"))
	    {
	    	next_print_time(ptc) = grid->time;
	    	if (print_start_time(pto) > grid->time)
	    	{
	    	    (void) printf("WARNING in set_next_print_time(), "
	    	                  "prt_rst specified with start_time %g "
	    	                  "greater than current time %g.\n"
	    	                  "\tResetting print_start_time.\n",
				  print_start_time(pto),grid->time);
	    	    print_start_time(pto) = grid->time;
	    	}
	    }
	}
	else
	{
	    /* Check for potential infinite loop below. */

	    if (print_step_interval(pto) <= 0)
	    {
	    	screen("ERROR in set_next_print_time(), "
	    	       "invalid mesh time print interval: 0\n");
	    	clean_up(ERROR);
	    }

	    next_print_step(ptc) = print_start_step(pto);

	    /* This is needed because at initialization, the
	     * solution has not yet been printed, whereas every
	     * other time we get to this point it has. */

	    if (grid->step > 0)
	        while (next_print_step(ptc) <= grid->step)
	    	    next_print_step(ptc) += print_step_interval(pto);

	    if (debugging("prt_rst"))
	    {
	    	next_print_step(ptc) = grid->step;
	    	if (print_start_step(pto) == grid->step)
	    	{
	    	    (void) printf("WARNING in set_next_print_time(), "
		                  "prt_rst specified with start_step %d "
		                  "greater than current step %d\n"
		                  "\tResetting print_start_step.\n",
				  print_start_step(pto),grid->step);
		    print_start_step(pto) = grid->step;
		}
	    }
	}
}		/*end set_next_print_time*/


#include <sys/stat.h>

/*
*		prompt_for_print_options():
*
*	Prompts for and reads in the fields in a PRINT_OPTIONS structure.
*	The print modes are as follows:
*
*       Mesh time -- prints out every print_step_interval time steps, beginning
*       at step print_start_step.
*
*	Exact time -- prints out at intervals of print_time_interval in real
*	time, beginning at time print_start_time.  The time step during the
* 	computation is adjusted to precisely land on the intended printing
*	times.
*
*	Constant time -- similar to exact time, except that the time steps are
*	not adjusted to land on the intended printing times; rather, printing
*	occurs when the actual time first exceeds the target printing time.
*/

LOCAL	void	prompt_for_print_options(
	INIT_DATA	*init,
	PRINT_OPTIONS	*pto,
	boolean		prompt_for_file_name,
	boolean		one_io_node,
	const char	*prompt)
{
	char    *c, s[Gets_BUF_SIZE];

	/*Set defaults for printing frequencies to entered data*/
	prt_mode(pto) = Prt_mode(Default_print_options(init));
	print_time_interval(pto) =
	    Print_time_interval(Default_print_options(init));
	print_start_time(pto) =
	    Print_start_time(Default_print_options(init));
	print_step_interval(pto) =
	    Print_step_interval(Default_print_options(init));
	print_start_step(pto) = Print_start_step(Default_print_options(init));
	print_in_binary(pto) = Print_in_binary(Default_print_options(init));

	screen("Specify the interval type for printing");
	if (prompt != NULL && strlen(prompt) != 0)
	    screen(" %s",prompt);
	screen("\n\t");
	screen("[exact%s, constant%s, mesh%s]: ",
	       (exact_time_output(prt_mode(pto)))    ? " (default)" : "",
	       (constant_time_output(prt_mode(pto))) ? " (default)" : "",
	       (mesh_time_output(prt_mode(pto)))     ? " (default)" : "");
	(void) Gets(s);
	if (s[0] == '\0')
	{
	    if (exact_time_output(prt_mode(pto)))
	    	(void) strcpy(s,"e");
	    else if (constant_time_output(prt_mode(pto)))
	    	(void) strcpy(s,"c");
	    else if (mesh_time_output(prt_mode(pto)))
	    	(void) strcpy(s,"m");
	}
	else
	    for (c = s; *c != '\0'; ++c)
		*c = tolower(*c);
	if ((s[0] == 'r') || (s[0] == 'e'))
	    prt_mode(pto) = EXACT_TIME;
	else if ((s[0] == 'c') || (s[0] == 'C'))
	    prt_mode(pto) = CONSTANT_TIME;
	else if ((s[0] == 'm') || (s[0] == 'M'))
	    prt_mode(pto) = MESH_TIME;

	if (real_time_output(prt_mode(pto)))
	{
	    screen("Enter the time interval and first time for printing\n"
	           "\t(default = %g %g): ",print_time_interval(pto),
	           			   print_start_time(pto));
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	    	static const char *fmt = "%lf %lf";
	    	(void) sscanf(s,fmt,&print_time_interval(pto),
	    		            &print_start_time(pto));
	    }
	}
	else
	{
	    screen("Enter the step interval and first step for printing\n"
	           "\t(default = %d %d): ",print_step_interval(pto),
		       			   print_start_step(pto));
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscanf(s,"%d %d",&print_step_interval(pto),
	    		      		&print_start_step(pto));
	}

	screen("Request binary/non-binary output "
	       "%s",(print_in_binary(pto) == YES) ? "[b(dflt),n]: " :
						    "[b,n(dflt)]: ");
	(void) Gets(s);
	if ((s[0] == 'n') || (s[0] == 'N'))
	    print_in_binary(pto) = NO;
	else if ((s[0] == 'b') || (s[0] == 'B'))
	    print_in_binary(pto) = YES;

	if (prompt_for_file_name == YES)
	    prompt_for_output_filenames(pto,one_io_node,prompt);

	(void) printf("\n");

	/*Set defaults for printing frequencies to entered data*/
	Prt_mode(Default_print_options(init)) = prt_mode(pto);
	Print_time_interval(Default_print_options(init)) = 
	    print_time_interval(pto);
	Print_start_time(Default_print_options(init)) = 
	    print_start_time(pto);
	Print_step_interval(Default_print_options(init)) = 
	    print_step_interval(pto);
	Print_start_step(Default_print_options(init)) = print_start_step(pto);
	Print_in_binary(Default_print_options(init)) = print_in_binary(pto);
}		/*end prompt_for_print_options*/

EXPORT	void	prompt_for_output_filenames(
	PRINT_OPTIONS   *pto,
	boolean		one_io_node,
	const char	*prompt)
{
	char    name[Gets_BUF_SIZE];
	char    s[Gets_BUF_SIZE];
	int	nn = pp_numnodes();

	/* Load a file name into name (either use default or input).
	 * Use an empty string if default is "NONE". */

	if (print_filename(pto)[0] == '\0')
	{
	    screen("Enter a file name for the output\n"
	           "\tfor %s (default = stdout): ",prompt);
	    (void) Gets(name);
	}
	else if (strncmp(print_filename(pto),"NO_DEFAULT",10) == 0)
	{
	    screen("Enter a file name for the output\n\tfor %s: ",prompt);
	    (void) Gets(name);
	    if (name[0] == '\0')
	    {
	        screen("ERROR in prompt_for_output_filenames(), "
		       "file name not entered\n");
	        clean_up(ERROR);
	    }
	}
	else if (strncmp(print_filename(pto),"NONE",4) != 0)
	{
	    screen("Enter a file name for the output\n "
	           "\tfor %s (default = %s): ",prompt,print_filename(pto));
	    (void) Gets(name);
	    if (name[0] == '\0')
	        (void) strcpy(name,print_filename(pto));
	}
	else
	{
	    name[0] = '\0';
	}

	/*
	 * Look for a directory component in name by stepping back
	 * from the end of the string looking for a /.
	 */

	if (strstr(name,"/") != NULL)
	{
	    /* File name contains a directory
	     * component which overrides default directory (if given). */
    
	    size_t len = strlen(name);
    
	    if (name[len-1] == '/')
	    {
	        screen("ERROR in prompt_for_output_filenames(), "
    		       "invalid file name %s\n",name);
	        clean_up(ERROR);
	    }
	    (void) strcpy(s,name);
	    for (len--; name[len] != '/'; len--)
	    	name[len] = '\0';
	    (void) strcpy(print_directory(pto),name);
	    if ((nn > 1) && (one_io_node == NO))
	        adjoin_node_number(s);
	    (void) strcpy(print_filename(pto),s);
	}
	else
	{
	    if (print_directory(pto)[0] != '\0')
	    {
	        /* Have a file name with no directory component and a
	         * default directory. */

	        if (strcmp(name,"NONE") != 0)
	        {
	            (void) sprintf(s,"%s/%s",print_directory(pto),name);
	            if ((nn > 1) && (one_io_node == NO))
	                adjoin_node_number(s);
	            (void) strcpy(print_filename(pto),s);
	        }
	    }
	    else
	    {
	        /* No directory component in fname and no default
	         * directory.  Try to use current directory as
	         * default. */

	        if (getcwd(print_directory(pto),64) == NULL)
	        {
		    strcpy(print_directory(pto),".");
	            (void) printf("WARNING in prompt_for_output_filenames(), "
			          "cannot get cwd\n");
	        }
	        screen("Enter a directory name for the output\n"
	               "\tfor %s(default = %s): ",prompt,
		       print_directory(pto));
	        (void) Gets(s);
	        if (s[0] != '\0')
	        {
	            (void) strcpy(print_directory(pto),s);
	        }
	        if (strcmp(name,"NONE") != 0)
	        {
	            (void) sprintf(s,"%s/%s",print_directory(pto),name);
	            if ((nn > 1) && (one_io_node == NO))
	                adjoin_node_number(s);
	            (void) strcpy(print_filename(pto),s);
	        }
	    }
	}
}		/*end prompt_for_output_filenames*/


/*
*			init_output_data():
*
*	This function initializes the OUTPUT_DATA part of a
*	specialized output structure.  Print control is initialized,
*	taking the global value as the default (from prt).
*	For directory and file prompting, see open_data_file().
*/

EXPORT void init_output_data(
	INIT_DATA	*init,
	OUTPUT_DATA	*data,
	Grid		*grid,
	Printplot	*prt,
	const char	*prompt,
	boolean		prompt_for_file_name,
	boolean		open_file,
	boolean		one_io_node)
{
	PRINT_OPTIONS	*pto = &PrtOpts(data);
	prompt_for_print_options(init,pto,prompt_for_file_name,
				 one_io_node,prompt);
	set_output_data(pto,data,grid,prt,open_file,one_io_node);
}		/*end init_output_data*/

EXPORT void set_output_data(
	PRINT_OPTIONS	*pto,
	OUTPUT_DATA	*data,
	Grid		*grid,
	Printplot	*prt,
	boolean		open_file,
	boolean		one_io_node)
{
	PrtOpts(data) = *pto;
	set_next_print_time(grid,Print_control(data));

	Output_file(data) = NULL;
	Output_printplot(data) = prt;

	if ((one_io_node == YES) && !is_io_node(pp_mynode()))
	    open_file = NO;

	if (create_directory(Output_dir(data),NO) == FUNCTION_FAILED)
	{
	    screen("ERROR in set_output_data(), directory "
		   "%s doesn't exist and can't be made\n",Output_dir(data));
	    clean_up(ERROR);
	}

	if (open_file == YES)
	    create_output_data_file(data);

}		/*end set_output_data*/

EXPORT	void	create_output_data_file(
	OUTPUT_DATA	*data)
{
	char s[256];
	boolean io_node = (is_io_node(pp_mynode())) ? YES : NO;
	boolean status;
	if (create_directory(Output_dir(data),NO) == FUNCTION_FAILED)
	{
	    screen("ERROR in create_output_data_file(), directory "
		   "%s doesn't exist and can't be made\n",Output_dir(data));
	    clean_up(ERROR);
	}

	status = YES;
	if ((io_node == YES) && (Output_filename(data)[0] != '\0'))
	{
	    if ((Output_file(data) = fopen(Output_filename(data),"w+")) == NULL)
	    {
		(void) printf("WARNING in create_output_data_file(), "
		              "problem opening %s for writing, errno = %s\n",
			      Output_filename(data),strerror(errno));
		status = NO;
	    }
	    else
	    {
	        /* Buffering doesn't help small files and can cause
	         * useful data to be stored in the buffer for
	         * a long time.  Turn off buffering for these files.
	         */
	        setbuf(Output_file(data),NULL);
	        print_machine_parameters(Output_file(data));
	    }
	}
	if (pp_max_status(status) != YES)
	{
	    (void) sprintf(s,"ERROR in create_output_data_file(), "
		             "unable to open %s\n",Output_filename(data));
	    perror(s);
	    screen(s);
	    clean_up(ERROR);
	}
}		/*end create_output_data_file*/


EXPORT	void	copy_output_data(
	OUTPUT_DATA	*new_data,
	OUTPUT_DATA	*data,
	boolean		open_file,
	char		*filename)
{
	boolean io_node = (is_io_node(pp_mynode())) ? YES : NO;
	boolean status;

	*new_data = *data;
	Output_file(new_data) = NULL;
	if (filename != NULL)
	    (void) strcpy(Output_filename(new_data),filename);
	
	status = YES;
	if ((io_node == YES) && (open_file == YES) &&
	    (Output_filename(new_data) != NULL))
	{
	    Output_file(new_data) = fopen(Output_filename(new_data),"w");
	    if (Output_file(new_data) != NULL)
	    {
	        /* Buffering doesn't help small files and can cause
	         * useful data to be stored in the buffer for
	         * a long time.	 Turn off buffering for these files.
	         */
	         setbuf(Output_file(new_data),NULL);
	         print_machine_parameters(Output_file(new_data));
	    }
	    else
		status = NO;
	}
	if (pp_max_status(status) != YES)
	{
	    screen("ERROR in copy_output_data(), "
		   "unable to open %s\n",Output_filename(new_data));
	    clean_up(ERROR);
	}
}		/*end copy_output_data*/


/*
*				open_data_file():
*
*	This function tries to open a file for output separate from the main
*	output file.  The function allows for specification of default file
*	and directory names.  If either are given, they are used and no
*	prompting is done.  It is possible to include the directory component
*	in the default file name, or if prompted, in the input file name string.
*	If no defaults are provided, and no input is given, stdout is used.  To
*	force use of stdout, pass "NONE" for default_fname.  The final
*	directory and file names are returned in dname and fname.  The
*	final file name will have the form dname/fname (if dname
*	non-NULL), and a parallel suffix if appropriate. If (open_file == NO),
*	no file is opened, but all prompting is done.  The argument one_io_node
*	should be set equal to YES if and only if only one node should write
*	the output.  If one_io_node is YES, and this node is not an I/O node,
*	no file is opened.  The return value is the file pointer (perhaps NULL).
*/

EXPORT FILE *open_data_file(
	Front		*fr,
	const char	*prompt,
	boolean		open_file,
	boolean		one_io_node,
	const char	*default_dname,
	char		**dname,
	const char	*default_fname,
	char		**fname)
{
	FILE    *file = NULL;
	boolean dir_in_fname = NO;
	char    name[Gets_BUF_SIZE];
	char    s[Gets_BUF_SIZE];
	size_t  len;

	*dname = *fname = NULL;
	if ((one_io_node == YES) && !is_io_node(pp_mynode()))
	    open_file = NO;

	/* Load a file name into name (either use default or input).
	 * Use an empty string if default is "NONE". */

	if (default_fname == NULL)
	{
	    screen("Enter a file name for the output\n "
	           "\tfor %s (default = stdout): ",prompt);
	    (void) Gets(name);
	}
	else if (strcmp(default_fname,"NO_DEFAULT") == 0)
	{
	    screen("Enter a file name for the output\n\tfor %s: ",prompt);
	    (void) Gets(name);
	    if (name[0] == '\0')
	    {
	        screen("ERROR in open_data_file(), file name not entered\n");
	        clean_up(ERROR);
	    }
	}
	else if (strcmp(default_fname,"NONE") != 0)
	{
	    screen("Enter a file name for the output\n "
	           "\tfor %s (default = %s): ",prompt,default_fname);
	    (void) Gets(name);
	    if (name[0] == '\0')
	        (void) strcpy(name,default_fname);
	}
	else
	{
	    name[0] = '\0';
	}

	/* Look for a directory component in name by stepping back from the
	 * end of the string looking for a /. */

	if (name[0] != '\0')
	{
	    dir_in_fname = (strstr(name,"/") != NULL) ? YES : NO;
	}

	if (name[0] == '\0')
	{
	    /* Have default fname of "NONE" or empty input. */

	    file = stdout;
	}
	else if (dir_in_fname == YES)
	{
	    /* Default fname or input value contains a directory
	     * component which overrides default directory (if given). */

	    len = strlen(name);

	    if (name[len-1] == '/')
	    {
	        screen("ERROR in open_data_file(), "
    		       "invalid file name %s\n",name);
	        clean_up(ERROR);
	    }
	    (void) strcpy(s,name);
	    for (len--; name[len] != '/'; len--)
	    	name[len] = '\0';
	    uni_array(dname,strlen(name)+1,CHAR);
	    (void) strcpy(*dname,name);
	    if ((fr->pp_grid->nn > 1) && (one_io_node == NO))
	        adjoin_node_number(s);
	    uni_array(fname,(strlen(s)+100),CHAR);
	    (void) strcpy(*fname,s);
	}
	else
	{
	    if (default_dname != NULL)
	    {
	        /* Have a file name with no directory component and a
	         * default directory. */

	        len = strlen(default_dname) + 1;
	        scalar(dname,len*sizeof(char));
	        (void) strcpy(*dname,default_dname);
	        if (strcmp(name,"NONE") != 0)
	        {
	            (void) sprintf(s,"%s/%s",*dname,name);
	            if ((fr->pp_grid->nn > 1) && (one_io_node == NO))
	                adjoin_node_number(s);
	            uni_array(fname,(strlen(s)+100),CHAR);
	            (void) strcpy(*fname,s);
	        }
	    }
	    else
	    {
	        char    cur_dir[64];

	        /* No directory component in fname and no default
	         * directory.  Try to use current directory as
	         * default. */

	        if (getcwd(cur_dir,64) == NULL)
	        {
		    strcpy(cur_dir,".");
	            (void) printf("WARNING in open_data_file(), "
				  "cannot get cwd\n");
	        }
	        screen("Enter a directory name for the output\n"
	               "\tfor %s(default = %s): ",prompt,cur_dir);
	        (void) Gets(s);
	        if (s[0] != '\0')
	        {
	            len = strlen(s) + 1;
	            scalar(dname,len*sizeof(char));
	            (void) strcpy(*dname,s);
	        }
	        else
	        {
	            len = strlen(cur_dir) + 1;
	            scalar(dname,len*sizeof(char));
	            (void) strcpy(*dname,cur_dir);
	        }
	        if (strcmp(name,"NONE") != 0)
	        {
	            (void) sprintf(s,"%s/%s",*dname,name);
	            if ((fr->pp_grid->nn > 1) && (one_io_node == NO))
	                adjoin_node_number(s);
	    	    uni_array(fname,(strlen(s)+100),CHAR);
	            (void) strcpy(*fname,s);
	        }
	    }
	}

	if (create_directory(*dname,NO) == FUNCTION_FAILED)
	{
	    screen("ERROR in open_data_file(), directory %s doesn't exist "
	           "and can't be made\n",*dname);
	    clean_up(ERROR);
	}

	if (open_file == YES)
	{
	    if (*fname != NULL)
	    {
	        file = fopen(*fname,"w+");
	        if (file == NULL)
	        {
	            screen("ERROR in open_data_file(), unable to open %s\n",
			   *fname);
	            clean_up(ERROR);
	        }

	        /* Buffering doesn't help small files and can cause
	         * useful data to be stored in the buffer for
	         * a long time.  Turn off buffering for these files.
	         */
	        setbuf(file,NULL);
	        print_machine_parameters(file);
	    }
	}
	return file;
}		/*end open_data_file*/

/*
*			add_user_output_function():
*
*	This function adds a new function to the array of user statistics
*	functions in the g_User_printplot.  See user_output().
*/

EXPORT	void	add_user_output_function(
	void		(*user_output_func)(Grid*,Wave*,Front*,Printplot*,
					    OUTPUT_DATA*,boolean),
	OUTPUT_DATA	*user_output_data,
	Printplot	*prt)
{
	OUTPUT_DATA	**ptmp;
	int		i, n;
	void		(**ftmp)(Grid*,Wave*,Front*,Printplot*,
				 OUTPUT_DATA*,boolean);

	if (prt->user_output_funcs == NULL)
	{
	    uni_array(&prt->user_output_funcs,2,sizeof(user_output_func));
	    uni_array(&prt->user_output_data,2,sizeof(OUTPUT_DATA*));
	    n = 0;
	}
	else
	{
	    for (n = 0; prt->user_output_funcs[n] != NULL; ++n);
	    ftmp = prt->user_output_funcs;
	    ptmp = prt->user_output_data;
	    uni_array(&prt->user_output_funcs,n+2,sizeof(user_output_func));
	    uni_array(&prt->user_output_data,n+2,sizeof(OUTPUT_DATA*));
	    for (i = 0; i < n; ++i)
	    {
	    	prt->user_output_funcs[i] = ftmp[i];
	    	prt->user_output_data[i] = ptmp[i];
	    }
	    free_these(2,ftmp,ptmp);
	}
	prt->user_output_funcs[n] = user_output_func;
	prt->user_output_funcs[n+1] = NULL;
	prt->user_output_data[n] = user_output_data;
	prt->user_output_data[n+1] = NULL;
}		/*end add_user_output_function*/


EXPORT	void	print_plotting_choices(
	const char	*mesg,
	Plot_choice	*plt_choices)
{
	Plot_choice	*pc;
	size_t		max_len = 60, len = 0;

	if (plt_choices == NULL)
		return;
	screen("%s",mesg);
	screen("The choices are --\n");
	screen("\t\t");
	for (pc = plt_choices; pc != NULL; pc = pc->next)
	{
		if (pc->next == NULL)
		{
			screen("or ");
			len += 3;
		}
		len += strlen(pc->prompt) + strlen(pc->selector) + 5;
		if (len > max_len)
		{
			screen("\n\t\t");
			len = strlen(pc->prompt) + strlen(pc->selector) + 5;
		}
		screen("%s (%s)%s",pc->prompt,pc->selector,
				   (pc->next == NULL) ? "." : ", ");
	}
}		/*end print_plotting_choices*/


LOCAL	int	count_selected_plot_choices(
	char	    *s,
	Plot_choice *plt_choices)
{
	Plot_choice	  *pc;
	char		  *c;
	char		  buf[Gets_BUF_SIZE];
	char		  *usr_requested;
	static const char *separators = " \t";
	int		  before, n_vars;
	size_t		  len;

	if (s == NULL)
	    return 0;
	len = strlen(s)+1;
	if (len > Gets_BUF_SIZE)
	    uni_array(&usr_requested,strlen(s)+1,CHAR);
	else
	    usr_requested = buf;
	(void)strcpy(usr_requested,s);
	for (n_vars = 0, c = strtok(usr_requested,separators); c != NULL;
					c = strtok(NULL,separators))
	{
	    before = n_vars;
	    for (pc = plt_choices; pc != NULL; pc = pc->next)
		n_vars += SelectedScalarItems(c,pc,NULL);
	    if (n_vars == before)
	    {
	        screen("ERROR in count_selected_plot_choices(), "
		       "NO SUCH CHOICE %s\n",c);
	        clean_up(ERROR);
	    }
	}
	if (len > Gets_BUF_SIZE)
	    free(usr_requested);
	return	n_vars;
}		/*end count_selected_plot_choices*/

/*
*			prompt_for_tri_plots():
*/

EXPORT	void	prompt_for_tri_plots(
	INIT_DATA	*init)
{
	print_plotting_choices("\nSpecify variables to be plotted (tri grid). ",
			 plot_choices(init));

	selected_tri_plots(init) = read_plotting_choices(NULL,init);
}		/*end prompt_for_tri_plots*/

#if defined(__GD__)
EXPORT	void	prompt_for_gd_plots(
	INIT_DATA	*init)
{
	GD_PRINT_OPTIONS *opts;
	GD_FRAME_OPTS    *frame_opts;
	PRINT_OPTIONS    *gd_opts;
	Plot_choice      *pc;
	const char       *uc_type_name = "GD_MOVIE";
        const char       *lc_type_name = "gd_movie";
	char	         *c, s[1024];
	int              prt_offset = GD_MOVIE;
	ScalarPlotItem   *spi;
	int              i, var, nvars;
	const char	 *separators = " \t";

	opts = &GD_prt_options(init);

	screen("\n\t\tGD Movie plotting initialization\n");
	print_plotting_choices("\nSpecify variables to be plotted . ",
			 plot_choices(init));
	selected_gd_plots(init) = read_plotting_choices(NULL,init);
	nvars = count_selected_plot_choices(selected_gd_plots(init),
			plot_choices(init));
	gd_num_vars(opts) = nvars;

	if (nvars > 0)
	{
	    gd_opts = prt_opts(init) + prt_offset;
	    if (output_filename(init)[0] != '\0')
	    {
	        char *bname, *dirname;
	        base_and_dir_name((const char*)output_filename(init),
	    		          &dirname,&bname);
	        (void) sprintf(print_filename(gd_opts),"%s%s%s/%s/%s",
	    	                      (strlen(dirname)!=0) ? dirname : "",
	    			      (strlen(dirname)!=0) ? "/" : "",
				      lc_type_name,bname,bname);
	    }
	    else
	        strcpy(print_filename(gd_opts),"NO_DEFAULT");
	    (void) sprintf(s,"%s data",uc_type_name);
	    prompt_for_output_filenames(gd_opts,YES,s);

	    frame_opts = (GD_FRAME_OPTS*)Store(nvars*sizeof(GD_FRAME_OPTS));
	    GD_frame_options(init) = frame_opts;
	    for (var = 0, c = strtok(selected_gd_plots(init),separators); 
	    		c != NULL; c = strtok(NULL,separators), ++var)
	    {
		(void) strcpy(GD_selector(frame_opts[var]),c);
	    }
	    for (var = 0; var < nvars; var++)
	    {
	    	GD_FRAME_OPTS    *fopts = frame_opts+var;
		for (pc = plot_choices(init); pc != NULL; pc = pc->next)
		{
		    if (SelectedScalarItems(GD_selector(frame_opts[var]),
		    	pc,&spi) != 0)
		    {
			for(; spi != NULL; spi = spi->next)
			{
			    gd_plot_function(fopts) = spi->plot_fn;
			    sprintf(gd_var_name(fopts),"%s",spi->name);
			    sprintf(gd_plot_name(fopts),"%s-%s.gif",
					print_filename(gd_opts),spi->name);
			}
			break;
		    }
		    if (pc == NULL)
		    {
		    	(void) printf("WARNING in prompt_for_prostar_plots(), "
					"NO SUCH CHOICE %s\n",c);
		    }
		}
	    }
	}
	screen("\n\t\tEnd GD Movie plotting initialization\n");
}		/*end prompt_for_gd_plots*/
#endif /* defined(__GD__) */

/*
*			prompt_for_prostar_plots():
*/

EXPORT	void	prompt_for_prostar_plots(
	INIT_DATA	*init)
{

	PROSTAR_PRINT_OPTIONS *opts;
	PROSTAR_FRAME_OPTS    *frame_opts;
	PRINT_OPTIONS         *prostar_opts;
	Plot_choice	      *pc;
	ScalarPlotItem	      *spi;
	char	              *c, s[1024];
	char	              *selected_plots;
	const char	      *separators = " \t";
	double	              max_len;
	int                   prt_offset = PROSTAR_STATES;
	const char            *uc_type_name = "PROSTAR";
        const char            *lc_type_name = "prostar";
	int                   nvars;
	int	              i, var, dim = i_intfc(init)->dim;
	const char            *dname[] = {"error", "x", "x and y", 
					  "x, y, and z"};
	const char            *fmt = "%lf %lf %lf";

	if (debugging("PROSTAR"))
	    (void) printf("Entering prompt_for_prostar_plots()\n");

	opts = &PROSTAR_prt_options(init);

	screen("\n\t\tPROSTAR plotting initialization\n");
	zero_scalar(opts,sizeof(PROSTAR_PRINT_OPTIONS));

	screen("\n NOTE- PROSTAR data is written as a vector of 6 variables.\n"
	       "The first three variables must be x-velocity,\n"
	       "y-velocity, and z-velocity(zero will be entered for 2d).\n"
	       "For 2d, you must enter X Y [up to 3 more variables].\n"
	       "For 3d, you must enter X Y Z [up to three more variables]\n");
	print_plotting_choices("\nSpecify variables to be plotted (PROSTAR). ",
			 plot_choices(init));
	selected_plots = read_plotting_choices(NULL,init);
	nvars = count_selected_plot_choices(selected_plots,plot_choices(init));
	prostar_num_vars(opts) = nvars;

	if (nvars > 0)
	{
	    for (i = 0; i < dim; i++)
	    {
	        prostar_L0(opts)[i] = Comp_grid(init).L[i];
	        prostar_U0(opts)[i] = Comp_grid(init).U[i];
	    }

	    screen("\nEnter the coordinates of lower corner\n\t"
	           "of the initial view box (dflt =");
	    for (i = 0; i < dim; i++)
	        screen(" %g",prostar_L0(opts)[i]);
	    screen("): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscanf(s,fmt,prostar_L0(opts),
	    			prostar_L0(opts)+1,
	    			prostar_L0(opts)+2); 
	    screen("Enter the coordinates of upper corner\n\t"
	           "of the initial view box (dflt =");
	    for (i = 0; i < dim; i++)
	        screen(" %g",prostar_U0(opts)[i]);
	    screen("): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscanf(s,fmt,prostar_U0(opts),prostar_U0(opts)+1,
			      prostar_U0(opts)+2);
	
	    max_len = 0;
	    for (i = 0; i < dim; i++)
	    {
	        prostar_len(opts)[i] = prostar_U0(opts)[i] - 
		    prostar_L0(opts)[i];
	        if (prostar_len(opts)[i] > max_len)
	            max_len = prostar_len(opts)[i];
	    }

	    for (i = 0 ; i < 3; i++)
	        prostar_pixels(opts)[i] = 1;

	    screen("Specify the number of pixels in the %s direction%s ",
	    	   dname[dim],(dim>1)?"s":"");
	    screen("(dflt ");
	    for (i = 0; i < dim; i++)
	        screen("%d%s",prostar_pixels(opts)[i],
		       (i==(dim-1)) ? "): " : " ");
	    (void) Gets(s);
	    if (dim == 2)
	    {
		if (s[0] != '\0')
		    (void) sscanf(s,"%d %d %d",prostar_pixels(opts),
				prostar_pixels(opts)+1); 
	    }
	    else if (dim == 3)
	    {
		if (s[0] != '\0')
		    (void) sscanf(s,"%d %d %d",prostar_pixels(opts),
				prostar_pixels(opts)+1,
				prostar_pixels(opts)+2); 
	    }
	    (void) printf("width = %d \t length = %d \t height = %d \n\n", 
		      *prostar_pixels(opts), *(prostar_pixels(opts)+1), 
			  *(prostar_pixels(opts)+2));
	
	    /* Prompt for frame data */

	    frame_opts = (PROSTAR_FRAME_OPTS*)Store(nvars*sizeof
						    (PROSTAR_FRAME_OPTS));
	    PROSTAR_frame_options(init) = frame_opts;

	    for (var = 0, c = strtok(selected_plots,separators); c != NULL; 
		 c = strtok(NULL,separators), var++)
	    {
	        (void) strcpy(PROSTAR_selector(frame_opts[var]),c);
	    }

	    for (var = 0; var < nvars; var++)
	    {
	        PROSTAR_FRAME_OPTS    *fopts = frame_opts+var;
	        for (pc = plot_choices(init); pc != NULL; pc = pc->next)
	        {
	    	    if (SelectedScalarItems(prostar_selector(fopts),pc,&spi) != 0)
	    	    {
	    	        for(; spi != NULL; spi = spi->next)
	    	            prompt_for_PROSTAR_frame_data(opts,fopts,spi);
	    	        break;
	            }
	        }
	        if (pc == NULL)
	        {
	    	    (void) printf("WARNING in prompt_for_prostar_plots(), "
	    	                  "NO SUCH CHOICE %s\n",c);
	        }
	    }

	    (void) sprintf(s,"\nYou will now be prompted for a base "
	                     "file name and optional directory for the "
			     "%s output. Output for each variable is to "
			     "a separate file whose name contains "
	                     "the base name, and the prompt string for "
			     "that variable.\n",uc_type_name);
	    screen_print_long_string(s);

	    prostar_opts = prt_opts(init) + prt_offset;
	    if (output_filename(init)[0] != '\0')
	    {
	        char *bname, *dirname;
	        base_and_dir_name((const char*)output_filename(init),
	    		          &dirname,&bname);
	        (void) sprintf(print_filename(prostar_opts),"%s%s%s/%s/%s",
	    	                      (strlen(dirname)!=0) ? dirname : "",
	    			      (strlen(dirname)!=0) ? "/" : "",
				      lc_type_name,
	    	                      bname,bname);
	    }
	    else
	        strcpy(print_filename(prostar_opts),"NO_DEFAULT");
	    (void) sprintf(s,"%s data",uc_type_name);
	    prompt_for_output_filenames(prostar_opts,YES,s);
	}

	screen("\n\t\tEnd %s plotting initialization\n\n",uc_type_name);	    
	   
	if (debugging("PROSTAR"))
	    (void) printf("Leaving prompt_for_prostar_plots()\n");
}		/*end prompt_for_prostar_plots*/


#if defined(TWOD)
/*
*			init_tri_plots():
*
*	This function parses an input string of letters, loading the name
*	and function fields for tri plots in the Printplot structure.
*
*	TODO: this function will not support multiple calls.  It will delete
*	an existing tri plot initialization, instead of appending.
*/

EXPORT	void	init_tri_plots(
	INIT_DATA	*init,
	Printplot	*prt)
{
	Plot_choice	*plt_choices = plot_choices(init);
	Plot_choice	*pc;
	ScalarPlotItem	*spi;
	int		var;
	char		*s = selected_tri_plots(init);
	static const char	*separators = " \t";
	char		*c;


	/*count number of plotting variables*/
	prt->n_tri_vars = count_selected_plot_choices(s,plt_choices);

	if (prt->n_tri_vars == 0)
	    return;

	uni_array(&prt->tri_plot_name,prt->n_tri_vars,sizeof(char *));
	uni_array(&prt->tri_plot_function,prt->n_tri_vars,
		sizeof(double (*)(double*,Front*,Wave*,COMPONENT,Locstate,
				 double*,double*)));

	for (var = 0, c = strtok(s,separators); c != NULL; 
						c = strtok(NULL,separators))
	{
	    for (pc = plt_choices; pc != NULL; pc = pc->next)
	    {
		if (SelectedScalarItems(c,pc,&spi) > 0)
		{
		    for(; spi != NULL; spi = spi->next)
		    {
		        prt->tri_plot_name[var] = spi->name;
		        prt->tri_plot_function[var++] = spi->plot_fn;
	            }
		    break;
	        }
	    }
	    if (pc == NULL)
	    {
		(void) printf("WARNING in init_tri_plots(), "
		              "NO SUCH CHOICE %s\n",c);
	    }
	}
}		/*end init_tri_plots*/
#endif /* defined(TWOD) */

EXPORT	char	*read_plotting_choices(
	char 		*s,
	INIT_DATA	*init)
{
	struct _line_list { char buf[Gets_BUF_SIZE];
			    struct _line_list *next;} head, *line, *next;
	size_t	len;

	screen("\n\tEnter choices as a space separated list, ");
	screen("using multiple lines if needed.");
	screen("\n\tTerminate all lines EXCEPT THE LAST with a backslash '\\'");
	screen("\n\tEnter choices: ");

	head.next = NULL;
	line = &head;
	(void) Gets(line->buf);
	len = strlen(line->buf);
	if (len != 0)
	{
	    while (line->buf[strlen(line->buf)-1] == '\\')
	    {
	    	scalar(&line->next,sizeof(struct _line_list));
	    	line = line->next;
	    	line->next = NULL;
	    	(void) printf("\t: ");
	    	(void) Gets(line->buf);
	    	if (line->buf[0] == '\0')
	    	    break;
	    	len += strlen(line->buf);
	    }
	}

	if (s == NULL)
	    s = (char*) init_table_Store((len+1)*CHAR,init);
	line = &head;
	next = line->next;
	if (next != NULL)
	    line->buf[strlen(line->buf)-1] = ' ';
	(void) strcpy(s,line->buf);
	for (line = next; line != NULL; line = next)
	{
	    if (line->next != NULL)
	    	line->buf[strlen(line->buf)-1] = ' ';
	    (void) strcat(s,line->buf);
	    next = line->next;
	    free(line);
	}
	return s;
}		/*end read_plotting_choices*/


#if defined(USE_HDF)
/*
*			prompt_for_hdf_plots():
*/

EXPORT	void	prompt_for_hdf_plots(
	HDF_DATA_TYPE     type,
	INIT_DATA	  *init)
{
	HDF_PRINT_OPTIONS *opts;
	HDF_FRAME_OPTS    *frame_opts;
	PRINT_OPTIONS	  *hdf_opts;
	Plot_choice	  *pc;
	ScalarPlotItem	  *spi;
	char	          *c, s[1024];
	const char	  *uc_type_name;
	const char	  *lc_type_name;
	char	          *selected_plots;
	static const char	  *separators = " \t";
	double	          max_len, resolution;
	int               prt_offset;
	int               nvars;
	int	          i, var, dim = i_intfc(init)->dim;
	static const char *dname[] = {"error", "x", "x and y", "x, y, and z"};
	static const char *fmt = "%lf %lf %lf";

	switch (type)
	{
	case HDF_RASTER:
	    if (dim != 2)
		return;
	    opts = &HDF_prt_options(init);
	    uc_type_name = "HDF";
	    lc_type_name = "hdf";
	    prt_offset = HDF_STATES;
	    break;
	case HDF_SDS:
	    opts = &SDS_prt_options(init);
	    uc_type_name = "SDS";
	    lc_type_name = "sds";
	    prt_offset = SDS_STATES;
	    break;
	default:
	    screen("ERROR in prompt_for_hdf_plots(), unknown data type\n");
	    clean_up(ERROR);
	}
	screen("\n\t\t%s plotting initialization\n",uc_type_name);

	zero_scalar(opts,sizeof(HDF_PRINT_OPTIONS));

	(void) sprintf(s,"\nSpecify variables to be plotted (%s). ",
		       uc_type_name);
	print_plotting_choices(s,plot_choices(init));
	selected_plots = read_plotting_choices(NULL,init);

	nvars = count_selected_plot_choices(selected_plots,plot_choices(init));
	hdf_num_vars(opts) = nvars;

	if (nvars > 0)
	{
	    hdf_data_type(opts) = type;
	    for (i = 0; i < dim; ++i)
	    {
	        hdf_L0(opts)[i] = Comp_grid(init).L[i];
	        hdf_U0(opts)[i] = Comp_grid(init).U[i];
	        hdf_V(opts)[i] = 0.0;
	    }

	    screen("\nEnter the coordinates of lower corner\n\t"
	           "of the initial view box (dflt =");
	    for (i = 0; i < dim; ++i)
	        screen(" %g",hdf_L0(opts)[i]);
	    screen("): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscanf(s,fmt,hdf_L0(opts),
	    			hdf_L0(opts)+1,
	    			hdf_L0(opts)+2); 
	    screen("Enter the coordinates of upper corner\n\t"
	           "of the initial view box (dflt =");
	    for (i = 0; i < dim; ++i)
	        screen(" %g",hdf_U0(opts)[i]);
	    screen("): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscanf(s,fmt,hdf_U0(opts),hdf_U0(opts)+1,hdf_U0(opts)+2);
	    screen("Enter the velocity of the view box (dflt =");
	    for (i = 0; i < dim; ++i)
	        screen(" %g",hdf_V(opts)[i]);
	    screen("): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscanf(s,fmt,hdf_V(opts),hdf_V(opts)+1,hdf_V(opts)+2); 

	    max_len = 0;
	    for (i = 0; i < dim; ++i)
	    {
	        hdf_len(opts)[i] = hdf_U0(opts)[i] - hdf_L0(opts)[i];
	        if (hdf_len(opts)[i] > max_len)
	            max_len = hdf_len(opts)[i];
	    }

	    resolution = (dim > 1) ? 600.0 : Comp_grid(init).gmax[0];

	    for (i = 0; i < dim; ++i)
	        hdf_pixels(opts)[i] =
		    irint(resolution*hdf_len(opts)[i]/max_len);
	    for (; i < 3; ++i)
	        hdf_pixels(opts)[i] = 1;

	    screen("Specify the number of pixels in the %s direction%s ",
	    	   dname[dim],(dim>1)?"s":"");
	    screen("(dflt ");
	    for (i = 0; i < dim; ++i)
	        screen("%d%s",hdf_pixels(opts)[i],(i==(dim-1)) ? "): " : " ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscanf(s,"%d %d %d",hdf_pixels(opts),
	    			           hdf_pixels(opts)+1,
	    			           hdf_pixels(opts)+2); 

	    /* Prompt for frame data */

	    frame_opts = (HDF_FRAME_OPTS*)Store(nvars*sizeof(HDF_FRAME_OPTS));

	    switch (type)
	    {
	    case HDF_RASTER:
	        HDF_frame_options(init) = frame_opts;
		break;
	    case HDF_SDS:
	        SDS_frame_options(init) = frame_opts;
		break;
	    }

	    for (var = 0, c = strtok(selected_plots,separators); c != NULL; 
	    	      c = strtok(NULL,separators), ++var)
	    {
	        (void) strcpy(HDF_selector(frame_opts[var]),c);
	    }

	    for (var = 0; var < nvars; ++var)
	    {
	        HDF_FRAME_OPTS    *fopts = frame_opts+var;
	        for (pc = plot_choices(init); pc != NULL; pc = pc->next)
	        {
	    	    if (SelectedScalarItems(hdf_selector(fopts),pc,&spi) != 0)
	    	    {
	    	        for(; spi != NULL; spi = spi->next)
	    	            prompt_for_HDF_frame_data(opts,fopts,spi);
	    	        break;
	            }
	        }
	        if (pc == NULL)
	        {
	    	    (void) printf("WARNING in prompt_for_hdf_plots(), "
	    	                  "NO SUCH CHOICE %s\n",c);
	        }
	    }

	    (void) sprintf(s,"\nYou will now be prompted for a base "
	                     "file name and optional directory for the "
			     "%s output. Output for each variable is to "
			     "a separate file whose name contains "
	                     "the base name, and the prompt string for "
			     "that variable.\n",uc_type_name);
	    screen_print_long_string(s);

	    hdf_opts = prt_opts(init) + prt_offset;
	    if (output_filename(init)[0] != '\0')
	    {
	        char *bname, *dirname;
	        base_and_dir_name((const char*)output_filename(init),
	    		          &dirname,&bname);
	        (void) sprintf(print_filename(hdf_opts),"%s%s%s/%s/%s",
	    	                      (strlen(dirname)!=0) ? dirname : "",
	    			      (strlen(dirname)!=0) ? "/" : "",
				      lc_type_name,
	    	                      bname,bname);
	    }
	    else
	        strcpy(print_filename(hdf_opts),"NO_DEFAULT");
	    (void) sprintf(s,"%s data",uc_type_name);
	    prompt_for_output_filenames(hdf_opts,YES,s);
	    prompt_for_compression_type(type,opts);
	}
	if (prt_offset == HDF_STATES)
	{
	    screen("\nEnter yes to show subdomain partition in the frame: ");
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y')
	    	hdf_subdomain_div(opts) = YES;
	    else
	    	hdf_subdomain_div(opts) = NO;
	}

	screen("\n\t\tEnd %s plotting initialization\n\n",uc_type_name);
}		/*end prompt_for_hdf_plots*/


LOCAL	void	prompt_for_compression_type(
	HDF_DATA_TYPE     type,
	HDF_PRINT_OPTIONS *opts)
{
	char      s[Gets_BUF_SIZE];
	comp_info *c_info = &hdf_compression_info(opts);
	
	sds_compression_type(opts) = COMP_CODE_NONE;
	ras_compression_type(opts) = COMP_RLE;
	screen("Enter the compression type, choices are\n");
	switch (type)
	{
	case HDF_RASTER:
	    screen("\tNone (N%s)\n",
	           (ras_compression_type(opts)==COMP_NONE) ?", default" : "");
	    screen("\tRun length encoding, no data loss (R%s)\n",
	           (ras_compression_type(opts)==COMP_RLE) ? ", default" : "");
	    screen("\tJPEG, some data loss (J%s)\n",
	           (ras_compression_type(opts)==COMP_JPEG) ? ", default" : "");
	    screen("Enter choice: ");
	    (void) Gets(s);
	    switch (s[0])
	    {
	    case 'N':
	    case 'n':
	        ras_compression_type(opts) = COMP_NONE;
	        break;
	    case 'R':
	    case 'r':
	        ras_compression_type(opts) = COMP_RLE;
	        break;
	    case 'J':
	    case 'j':
	        ras_compression_type(opts) = COMP_JPEG;
		c_info->jpeg.quality = 75;
		c_info->jpeg.force_baseline = TRUE;
	        screen("Enter the JPEG quality 0 (worst) <= q <= 100 (best) "
		       "(default = %d): ",c_info->jpeg.quality);
	        (void) Gets(s);
		if (s[0] != '\0')
		{
	            int	q;
		    (void) sscanf(s,"%d",&q);
		    if ((q < 0) || (q > 100))
		    {
			screen("ERROR in prompt_for_compression_type(), "
			       "invalid JPEG quality %d\n",q);
			clean_up(ERROR);
		    }
		    c_info->jpeg.quality = q;
		}
		if (DFR8setcompress(COMP_JPEG,c_info) != SUCCEED)
		{
	            screen("ERROR in prompt_for_compression_type(), "
		           "can't set JPEG compression\n");
		    clean_up(ERROR);
		}
	        break;
	    default:
	        break;
	    }
	    break;
	case HDF_SDS:
	    screen("\tNone (N%s)\n",
	           (sds_compression_type(opts)==COMP_CODE_NONE) ?
		       ", default" : "");
	    screen("\tRun length encoding (R%s)\n",
	           (sds_compression_type(opts)==COMP_CODE_RLE) ?
		       ", default" : "");
	    screen("\tGzip deflation (G%s)\n",
	           (sds_compression_type(opts)==COMP_CODE_DEFLATE) ?
		   ", default" : "");
	    screen("\tAdaptive Huffman algorithm (H%s)\n",
	           (sds_compression_type(opts)==COMP_CODE_SKPHUFF) ?
		   ", default" : "");
	    screen("Enter choice: ");
	    (void) Gets(s);
	    switch (s[0])
	    {
	    case 'N':
	    case 'n':
	        sds_compression_type(opts) = COMP_CODE_NONE;
	        break;
	    case 'R':
	    case 'r':
	        sds_compression_type(opts) = COMP_CODE_RLE;
	        break;
	    case 'G':
	    case 'g':
	        sds_compression_type(opts) = COMP_CODE_DEFLATE;
		c_info->deflate.level = 6;
	        screen("Enter the compression level "
		       "0 (fastest) <= q <= 9 (best) "
		       "(default = %d): ",c_info->deflate.level);
	        (void) Gets(s);
		if (s[0] != '\0')
		{
	            int	l;
		    (void) sscanf(s,"%d",&l);
		    if ((l < 0) || (l > 9))
		    {
			screen("ERROR in prompt_for_compression_type(), "
			       "invalid compression level %d\n",l);
			clean_up(ERROR);
		    }
		    c_info->deflate.level = l;
		}
		break;
	    case 'H':
	    case 'h':
	        sds_compression_type(opts) = COMP_CODE_SKPHUFF;
		break;
	    default:
	        break;
	    }
	    break;
	default:
	    break;
	}
}		/*end prompt_for_compression_type*/

/*
*			init_HDF_plots():
*
*	This function drives the initialization for HDF output.
*
*	TODO: this function will not support multiple calls.  It will delete
*	an existing HDF plot initialization, instead of appending.
*/

EXPORT void init_HDF_plots(
	HDF_DATA_TYPE   type,
	INIT_DATA	*init,
	Printplot	*prt)
{
	HDF_plot_data	  *hdf_data;
	HDF_PRINT_OPTIONS *opts;
	HDF_FRAME_OPTS    *frame_opts;
	OUTPUT_DATA	  *odata;
	const char	  *suffix;
	int     	  prt_offset;
	int		  i, var;
	int		  dim = i_intfc(init)->dim;
	int		  nvars;

	switch (type)
	{
	case HDF_RASTER:
	    opts = &HDF_prt_options(init);
	    nvars = hdf_num_vars(opts);
	    if (nvars == 0)
	    {
	        prt->n_HDF_vars = 0;
	        prt->HDF_data = NULL;
	        return;
	    }
	    prt->n_HDF_vars = 1;
	    scalar(&hdf_data,sizeof(HDF_plot_data));
	    HDF_print_opts(hdf_data) = HDF_prt_options(init);
	    prt->HDF_data = hdf_data;
	    frame_opts = HDF_frame_options(init);
	    prt_offset = HDF_STATES;
	    suffix = ".hdf";
	    break;
	case HDF_SDS:
	    opts = &SDS_prt_options(init);
	    nvars = hdf_num_vars(opts);
	    if (nvars == 0)
	    {
	        prt->n_SDS_vars = 0;
	        prt->SDS_data = NULL;
	        return;
	    }
	    prt->n_SDS_vars = 1;
	    scalar(&hdf_data,sizeof(HDF_plot_data));
	    HDF_print_opts(hdf_data) = SDS_prt_options(init);
	    prt->SDS_data = hdf_data;
	    frame_opts = SDS_frame_options(init);
	    prt_offset = SDS_STATES;
	    suffix = "";
	    break;
	}
	
	hdf_data->dim = dim;

	uni_array(&hdf_data->frame_data,nvars,sizeof(HDF_frame_data));

	for (var = 0; var < nvars; ++var)
	{
	    HDF_frame_opts(hdf_data->frame_data[var]) = frame_opts[var];
	    hdf_data->frame_data[var].first = YES;	
	    hdf_data->frame_data[var].append = NO;	
	    hdf_data->frame_data[var].cumulative_time_min =  HUGE_VAL;
	    hdf_data->frame_data[var].cumulative_time_max = -HUGE_VAL;
	}

	hdf_data->scale[0] = NULL;
	for (i = 0; i < dim; ++i)
	{
	    uni_array(hdf_data->scale+i+1,hdf_pixels(opts)[i],FLOAT);
	    hdf_data->step[i] = hdf_len(opts)[i]/hdf_pixels(opts)[i];
	}
	for (++i; i < 4; ++i)
	    hdf_data->scale[i] = NULL;

	hdf_data->num_values =
	    hdf_pixels(opts)[0]*hdf_pixels(opts)[1]*hdf_pixels(opts)[2];
	for (var = 0; var < nvars; ++var)
	{
	    uni_array(&hdf_data->frame_data[var].values,hdf_data->num_values,
		   FLOAT);
	}
	if (hdf_data_type(opts) == HDF_RASTER)
	{
	    hdf_data->num_raster_data = hdf_pixels(opts)[0]*hdf_pixels(opts)[1];
	    uni_array(&hdf_data->raster_data,hdf_data->num_raster_data,
		   sizeof(uint8));
	}
	else
	    hdf_data->num_raster_data = 0;
	uni_array(&hdf_data->comp,hdf_data->num_values,sizeof(COMPONENT));

	for (var = 0; var < nvars; ++var)
	    set_HDF_frame_data(hdf_data,var);

	odata = prt->main_output_data[prt_offset];
	PrtOpts(odata) = prt_opts(init)[prt_offset];
	create_output_data_file(odata);

	if (is_io_node(pp_mynode()))
	{
	    char *plot_name;
	    for (var = 0; var < nvars; ++var)
	    {
		plot_name = HDF_frame_plot_name(hdf_data->frame_data[var]);
	        uni_array(&hdf_data->frame_data[var].file_name,
	           strlen(plot_name)+strlen(Output_filename(odata))+50,CHAR);
	        (void) sprintf(hdf_data->frame_data[var].file_name,
	                       "%s-%s%s",Output_filename(odata),plot_name,
			       suffix);
	    }
	}
}		/*end init_HDF_plots*/

EXPORT	void	set_hdf_append(
	Printplot	*prt)
{
	turn_on_hdf_append(prt->HDF_data);
	turn_on_hdf_append(prt->SDS_data);
}		/*end set_hdf_append*/


LOCAL	void turn_on_hdf_append(
	HDF_plot_data	*hdf_data)
{
	FILE		*file;
	int		var, nvars;

	if (hdf_data == NULL)
	    return;

	nvars = HDF_num_vars(HDF_print_opts(hdf_data));
	for (var = 0; var < nvars; ++var)
	{
	    if ((file = fopen(hdf_data->frame_data[var].file_name,"r")) != NULL)
	    {
		hdf_data->frame_data[var].append = YES;
		(void) fclose(file);
	    }
	}
}		/*end turn_on_hdf_append*/

struct _HDF_FILTER_PROMPT
{
	const char      *name;
	const char      *selector;
	const char      *title;
	HDF_PLOT_FILTER *filter;
};
typedef struct _HDF_FILTER_PROMPT HDF_FILTER_PROMPT;

LOCAL	void	prompt_for_HDF_frame_data(
	HDF_PRINT_OPTIONS *opts,
	HDF_FRAME_OPTS	  *fopts,
	ScalarPlotItem	  *spi)
{
	char s[Gets_BUF_SIZE];
	int  i;
	static const char *fmt = "%lf %lf";

	static HDF_FILTER_PROMPT Filters[] = {
		        {"Identity", "none",  NULL,    identity_filter},
		        {"Log",      "log",   "LOG",   log},
		        {"Log1p",    "log1p", "LOG1P", log1p},
		        {"Exp",      "exp",   "EXP",   exp},
		        {"Expm1",    "expm1", "EXPM1", expm1},
		        {"Atan",     "atan",  "ATAN",  atan},
		        {"Tan",      "tan",   "TAN",   tan},
		        {NULL,        NULL,   NULL,    NULL}
	               };

	screen("\n");
	(void) strcpy(hdf_plot_name(fopts),spi->name);
	hdf_plot_function(fopts) = spi->plot_fn;
	hdf_plot_filter(fopts) = identity_filter;

	if (hdf_data_type(opts) == HDF_SDS)
	    return;

	screen("Enter an optional plotting filter for %s, choices are\n",
		hdf_plot_name(fopts));
	for (i = 0; Filters[i].name != NULL; ++i)
	{
	    screen("\t%s filter (%s",Filters[i].name,Filters[i].selector);
	    if (Filters[i].filter == hdf_plot_filter(fopts))
		screen(", default");
	    screen(")\n");
	}
	screen("Enter choice: ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    for (i = 0; Filters[i].name != NULL; ++i)
	    {
	        if (strcasecmp(s,Filters[i].selector)==0)
	        {
	    	    hdf_plot_filter(fopts) = Filters[i].filter;
		    if (Filters[i].title != NULL)
		    {
	    	        (void) sprintf(hdf_plot_name(fopts),"%s_%s",
				       Filters[i].title,spi->name);
		    }
		    break;
	        }
	    }
	}

	screen("Enter optional scaling factors for %s: ",hdf_plot_name(fopts));
	(void) Gets(s);
	hdf_dflt_scale(fopts) = YES;	
	if (s[0] != '\0')
	{
	    hdf_dflt_scale(fopts) = NO;
	    if (sscanf(s,fmt,&hdf_scale_min(fopts),&hdf_scale_max(fopts))!=2)
	    {
	    	screen("ERROR in prompt_for_HDF_frame_data(), "
	    	       "invalid input of scale factors\n");
	    	clean_up(ERROR);
	    }
	}

	screen("Color palelette data for "
	       "%s can be entered as either an HDF\n",hdf_plot_name(fopts));
	screen_print_long_string("\tcolor palette file, or as a raw "
	       "color palette consisting of a binary list of unsigned chars "
	       "in the form of 256 red values, followed by 256 green values, "
	       "followed by 256 blue values. A raw format file is indicated "
	       "by appending a blank and the indicator raw after the "
	       "file name.  Otherwise any file entered will be assumes "
	       "to be in HDF palette format.\n");
	screen("Enter an optional color palette file: ");
	(void) Gets(hdf_palette_name(fopts));
}		/*end prompt_for_HDF_frame_data*/

LOCAL	void	set_HDF_frame_data(
	HDF_plot_data	*hdf_data,
	int		var)
{
	HDF_PRINT_OPTIONS	*opts = &HDF_print_opts(hdf_data);
	HDF_frame_data		*frame_data = hdf_data->frame_data + var;
	HDF_FRAME_OPTS		*fopts = &hdf_frame_opts(frame_data);
	int	i;
	static uint8	dflt_palette[3*256] = {
	    0xff,0xff,0xff,0x00,0x00,0x00,0xff,0xff,0x00,0xff,0x00,0xff,
	    0x00,0xff,0xff,0xff,0x00,0x00,0x00,0xff,0x00,0x00,0x00,0xff,
	    0x00,0x00,0xfd,0x01,0x02,0xfb,0x03,0x04,0xf9,0x04,0x06,0xf7,
	    0x06,0x08,0xf5,0x07,0x0a,0xf3,0x09,0x0c,0xf1,0x0a,0x0e,0xef,
	    0x0c,0x10,0xed,0x0e,0x12,0xeb,0x0f,0x14,0xe9,0x11,0x16,0xe7,
	    0x12,0x18,0xe4,0x14,0x1b,0xe2,0x15,0x1d,0xe0,0x17,0x1f,0xde,
	    0x18,0x21,0xdc,0x1a,0x23,0xda,0x1c,0x25,0xd8,0x1d,0x27,0xd6,
	    0x1f,0x29,0xd4,0x20,0x2b,0xd2,0x22,0x2d,0xd0,0x23,0x2f,0xce,
	    0x25,0x31,0xcb,0x27,0x34,0xc9,0x28,0x36,0xc7,0x2a,0x38,0xc5,
	    0x2b,0x3a,0xc3,0x2d,0x3c,0xc1,0x2e,0x3e,0xbf,0x30,0x40,0xbd,
	    0x31,0x42,0xbb,0x33,0x44,0xb9,0x35,0x46,0xb7,0x36,0x48,0xb5,
	    0x38,0x4a,0xb2,0x39,0x4d,0xb0,0x3b,0x4f,0xae,0x3c,0x51,0xac,
	    0x3e,0x53,0xaa,0x3f,0x55,0xa8,0x41,0x57,0xa6,0x43,0x59,0xa4,
	    0x44,0x5b,0xa2,0x46,0x5d,0xa0,0x47,0x5f,0x9e,0x49,0x61,0x9c,
	    0x4a,0x63,0x9a,0x4c,0x65,0x97,0x4e,0x68,0x95,0x4f,0x6a,0x93,
	    0x51,0x6c,0x91,0x52,0x6e,0x8f,0x54,0x70,0x8d,0x55,0x72,0x8b,
	    0x57,0x74,0x89,0x58,0x76,0x87,0x5a,0x78,0x85,0x5c,0x7a,0x83,
	    0x5d,0x7c,0x81,0x5f,0x7e,0x7e,0x60,0x81,0x7c,0x62,0x83,0x7a,
	    0x63,0x85,0x78,0x65,0x87,0x76,0x67,0x89,0x74,0x68,0x8b,0x72,
	    0x6a,0x8d,0x70,0x6b,0x8f,0x6e,0x6d,0x91,0x6c,0x6e,0x93,0x6a,
	    0x70,0x95,0x68,0x71,0x97,0x65,0x73,0x9a,0x63,0x75,0x9c,0x61,
	    0x76,0x9e,0x5f,0x78,0xa0,0x5d,0x79,0xa2,0x5b,0x7b,0xa4,0x59,
	    0x7c,0xa6,0x57,0x7e,0xa8,0x55,0x7f,0xaa,0x53,0x81,0xac,0x51,
	    0x83,0xae,0x4f,0x84,0xb0,0x4d,0x86,0xb2,0x4a,0x87,0xb5,0x48,
	    0x89,0xb7,0x46,0x8a,0xb9,0x44,0x8c,0xbb,0x42,0x8e,0xbd,0x40,
	    0x8f,0xbf,0x3e,0x91,0xc1,0x3c,0x92,0xc3,0x3a,0x94,0xc5,0x38,
	    0x95,0xc7,0x36,0x97,0xc9,0x34,0x98,0xcb,0x31,0x9a,0xce,0x2f,
	    0x9c,0xd0,0x2d,0x9d,0xd2,0x2b,0x9f,0xd4,0x29,0xa0,0xd6,0x27,
	    0xa2,0xd8,0x25,0xa3,0xda,0x23,0xa5,0xdc,0x21,0xa7,0xde,0x1f,
	    0xa8,0xe0,0x1d,0xaa,0xe2,0x1b,0xab,0xe4,0x18,0xad,0xe7,0x16,
	    0xae,0xe9,0x14,0xb0,0xeb,0x12,0xb1,0xed,0x10,0xb3,0xef,0x0e,
	    0xb5,0xf1,0x0c,0xb6,0xf3,0x0a,0xb8,0xf5,0x08,0xb9,0xf7,0x06,
	    0xbb,0xf9,0x04,0xbc,0xfb,0x02,0xbe,0xfd,0x00,0xbe,0xfb,0x00,
	    0xbf,0xf9,0x00,0xbf,0xf7,0x00,0xc0,0xf5,0x00,0xc1,0xf3,0x00,
	    0xc1,0xf1,0x00,0xc2,0xef,0x00,0xc2,0xed,0x00,0xc3,0xeb,0x00,
	    0xc3,0xe9,0x00,0xc4,0xe7,0x00,0xc4,0xe5,0x00,0xc5,0xe3,0x00,
	    0xc5,0xe1,0x00,0xc6,0xde,0x00,0xc6,0xdc,0x00,0xc7,0xda,0x00,
	    0xc7,0xd8,0x00,0xc8,0xd6,0x00,0xc8,0xd4,0x00,0xc9,0xd2,0x00,
	    0xc9,0xd0,0x00,0xca,0xce,0x00,0xca,0xcc,0x00,0xcb,0xca,0x00,
	    0xcb,0xc8,0x00,0xcc,0xc6,0x00,0xcc,0xc4,0x00,0xcd,0xc2,0x00,
	    0xcd,0xbf,0x00,0xce,0xbd,0x00,0xce,0xbb,0x00,0xcf,0xb9,0x00,
	    0xcf,0xb7,0x00,0xd0,0xb5,0x00,0xd1,0xb3,0x00,0xd1,0xb1,0x00,
	    0xd2,0xaf,0x00,0xd2,0xad,0x00,0xd3,0xab,0x00,0xd3,0xa9,0x00,
	    0xd4,0xa7,0x00,0xd4,0xa5,0x00,0xd5,0xa3,0x00,0xd5,0xa1,0x00,
	    0xd6,0x9e,0x00,0xd6,0x9c,0x00,0xd7,0x9a,0x00,0xd7,0x98,0x00,
	    0xd8,0x96,0x00,0xd8,0x94,0x00,0xd9,0x92,0x00,0xd9,0x90,0x00,
	    0xda,0x8e,0x00,0xda,0x8c,0x00,0xdb,0x8a,0x00,0xdb,0x88,0x00,
	    0xdc,0x86,0x00,0xdc,0x84,0x00,0xdd,0x82,0x00,0xdd,0x7f,0x00,
	    0xde,0x7d,0x00,0xde,0x7b,0x00,0xdf,0x79,0x00,0xdf,0x77,0x00,
	    0xe0,0x75,0x00,0xe1,0x73,0x00,0xe1,0x71,0x00,0xe2,0x6f,0x00,
	    0xe2,0x6d,0x00,0xe3,0x6b,0x00,0xe3,0x69,0x00,0xe4,0x67,0x00,
	    0xe4,0x65,0x00,0xe5,0x63,0x00,0xe5,0x61,0x00,0xe6,0x5e,0x00,
	    0xe6,0x5c,0x00,0xe7,0x5a,0x00,0xe7,0x58,0x00,0xe8,0x56,0x00,
	    0xe8,0x54,0x00,0xe9,0x52,0x00,0xe9,0x50,0x00,0xea,0x4e,0x00,
	    0xea,0x4c,0x00,0xeb,0x4a,0x00,0xeb,0x48,0x00,0xec,0x46,0x00,
	    0xec,0x44,0x00,0xed,0x42,0x00,0xed,0x3f,0x00,0xee,0x3d,0x00,
	    0xee,0x3b,0x00,0xef,0x39,0x00,0xef,0x37,0x00,0xf0,0x35,0x00,
	    0xf1,0x33,0x00,0xf1,0x31,0x00,0xf2,0x2f,0x00,0xf2,0x2d,0x00,
	    0xf3,0x2b,0x00,0xf3,0x29,0x00,0xf4,0x27,0x00,0xf4,0x25,0x00,
	    0xf5,0x23,0x00,0xf5,0x21,0x00,0xf6,0x1e,0x00,0xf6,0x1c,0x00,
	    0xf7,0x1a,0x00,0xf7,0x18,0x00,0xf8,0x16,0x00,0xf8,0x14,0x00,
	    0xf9,0x12,0x00,0xf9,0x10,0x00,0xfa,0x0e,0x00,0xfa,0x0c,0x00,
	    0xfb,0x0a,0x00,0xfb,0x08,0x00,0xfc,0x06,0x00,0xfc,0x04,0x00,
	    0xfd,0x02,0x00,0xfd,0x00,0x00,0xff,0xff,0xff,0x00,0x00,0x00
	};

	if (hdf_data_type(opts) == HDF_SDS)
	    return;

	/*Set up color palette */

	frame_data->num_table_colors = 8;
	frame_data->num_colors = 256 - frame_data->num_table_colors - 2;
	frame_data->line_color = 255;

	for (i = 0; i < 768; ++i)
	    frame_data->palette[i] = dflt_palette[i];
	if (hdf_palette_name(fopts)[0] != '\0')
	{
	    char	fname[Gets_BUF_SIZE];
	    char	raw[10];
	    int	n;

	    fname[0] = raw[0] = '\0';
	    n = sscanf(hdf_palette_name(fopts),"%s %s",fname,raw);
	    if (n == 1)
	    {
	    	if (DFPgetpal(fname,frame_data->palette) == FAIL)
	    	{
	    	    screen("ERROR in set_HDF_frame_data(), ");
	    	    screen("can't open palette file %s\n",fname);
	    	    clean_up(ERROR);
	    	}
	    }
	    else if ((n == 2) && (strcmp(raw,"raw") == 0))
	    {
	    	FILE	*fp;
	    	uint8	raw_palette[3*256];
	    
	    	if ((fp = fopen(fname,"r")) == NULL)
	    	{
	    	    screen("ERROR in set_HDF_frame_data(), "
	    	           "can't open %s\n",fname);
	    	    clean_up(ERROR);
	    	}
	    	if (fread(raw_palette,sizeof(uint8),3*256,fp) != 3*256)
	    	{
	    	    screen("ERROR in set_HDF_frame_data(), "
	    	           "invalid color map file\n",fname);
	    	    clean_up(ERROR);
	    	}
	    	for (i = 0; i < 256; ++i)
	    	{
	    	    frame_data->palette[3*i]   = raw_palette[i];
	    	    frame_data->palette[3*i+1] = raw_palette[256+i];
	    	    frame_data->palette[3*i+2] = raw_palette[2*256+i];
	    	}
	    	(void) fclose(fp);
	    }
	    else
	    {
	    	screen("ERROR in set_HDF_frame_data(), "
	    	       "invalid file name %s\n",hdf_palette_name(fopts));
	    	clean_up(ERROR);
	    }
	}
}		/*end set_HDF_frame_data*/

#endif /* defined(USE_HDF) */
#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
LOCAL	double	identity_filter(
	double	x)
{
	return x;
}		/*end identity_filter*/
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */

EXPORT void init_PROSTAR_plots(
	INIT_DATA	*init,
	Printplot	*prt)
{

	PROSTAR_plot_data	  *prostar_data;
	PROSTAR_PRINT_OPTIONS     *opts;
	PROSTAR_FRAME_OPTS        *frame_opts;
	OUTPUT_DATA	          *odata;
	int     	          prt_offset;
	int		          i, var;
	int		          dim = i_intfc(init)->dim;
	int		          nvars;

	opts = &PROSTAR_prt_options(init);
	nvars = prostar_num_vars(opts);

	if (debugging("PROSTAR"))
	    (void) printf("Entering init_PROSTAR_plots()");
	if (nvars == 0)
	{
	    prt->n_PROSTAR_vars = 0;
	    prt->PROSTAR_data = NULL;
	    return;
	}
	prt->n_PROSTAR_vars = 1;
	scalar(&prostar_data,sizeof(PROSTAR_plot_data));
	PROSTAR_print_opts(prostar_data) = PROSTAR_prt_options(init);
	prt->PROSTAR_data = prostar_data;
	frame_opts = PROSTAR_frame_options(init);
	prt_offset = PROSTAR_STATES;
	prostar_data->dim = dim;
	prostar_data->vrt_cel_created = 0;
	uni_array(&prostar_data->frame_data,nvars,sizeof(PROSTAR_frame_data));

	for (var = 0; var < nvars; var++)
	{
	    PROSTAR_frame_opts(prostar_data->frame_data[var]) = 
       	        frame_opts[var];
	    prostar_data->frame_data[var].first = YES;	
	    prostar_data->frame_data[var].append = NO;	
	}

	prostar_data->scale[0] = NULL;
	for (i = 0; i < 3; i++)
	    uni_array(prostar_data->scale+i+1,prostar_pixels(opts)[i],FLOAT);    
	for (i++; i < 4; i++)
	    prostar_data->scale[i] = NULL;

	odata = prt->main_output_data[prt_offset];
	PrtOpts(odata) = prt_opts(init)[prt_offset];
	create_output_data_file(odata);

	if (is_io_node(pp_mynode()))
	{
	    uni_array(&prostar_data->base_file_name,
		   strlen(Output_filename(odata))+50,CHAR);
	    (void) sprintf(prostar_data->base_file_name,
			   "%s",Output_filename(odata));
	    uni_array(&prostar_data->v_file_name,
	           strlen(Output_filename(odata))+50,CHAR);
	    uni_array(&prostar_data->c_file_name,
	           strlen(Output_filename(odata))+50,CHAR);
	    uni_array(&prostar_data->p_file_name,
	           strlen(Output_filename(odata))+50,CHAR);
	}
	if (debugging("PROSTAR"))
	    (void) printf("Leaving init_PROSTAR_plots()");

}	/* end init_PROSTAR_plots */

LOCAL	void	prompt_for_PROSTAR_frame_data(
	PROSTAR_PRINT_OPTIONS     *opts,
	PROSTAR_FRAME_OPTS	  *fopts,
	ScalarPlotItem	          *spi)
{
	char         s[Gets_BUF_SIZE];
	int          i;
	const char *fmt = "%lf %lf";
	
	screen("\n");
	(void) strcpy(prostar_plot_name(fopts),spi->name);
	prostar_plot_function(fopts) = spi->plot_fn;

	prostar_dflt_scale(fopts) = YES;	
	

}		/*end prompt_for_PROSTAR_frame_data*/

#if defined(TWOD)
/*
*			init_cross_sections():
*
*	Prompts for optional cross sectional plots.  A cross section
*	is defined as a line segment joining two points (X0,Y0), (X1,Y1)
*	with a given number of interior points.  State variables can
*	be printed at each point on a cross section in a format
*	to be processed by graphs.  In general this printout
*	should be handled by a physics dependent function called
*	from user_output().
*/

EXPORT void init_cross_sections(
	INIT_DATA	*init,
	Front		*front,
	Grid		*grid,
	Printplot	*prt)
{
	const char          *whiteSpace = " \t";
	CROSS_SECTION	    Cr_sec;
	CROSS_SECTION	    *cr_sec;
	Cross_Sections_data *cr_data;
	char                s[Gets_BUF_SIZE];
	/*
	 * 8 is ONE MORE than the MAXIMUM number of whitespace-separated items
         * which can appear on a single cross section specification line.
         * Update this if changes are made.
         */
	char		    *tokens[8], *crfile;
	int		    dim = front->rect_grid->dim;
	int                 i, j, k, inc_num, num_cr_sec, num_pts, nToks;
	double               L[MAXD], U[MAXD], dir[MAXD], dlen;
	double               len, offset = 0.0;

	if (dim != 2)
	{
	    if (debugging("warning"))
	        (void) printf("WARNING in init_cross_sections(), "
	                  " 3D not implemented.\n");
	    return;
	}
	
	screen("Type 'y' to obtain cross sectional plots: ");
	(void) Gets(s);
	if ((s[0] != 'y') && (s[0] != 'Y'))
	    return;

	screen("Enter the number of cross sectional plots to be initialized: ");
	(void) Scanf("%d\n",&num_cr_sec);
	if (num_cr_sec <= 0)
	    return;

	scalar(&cr_data,sizeof(Cross_Sections_data));
	if (prt->outfile != NULL)
	{
	    char *bname, *dname;
	    base_and_dir_name(prt->outfile,&dname,&bname);
	    crfile = s;
	    (void) sprintf(crfile,"%s%scr-sec/%s",
		                  (strlen(dname)!=0) ? dname : "",
				  (strlen(dname)!=0) ? "/" : "",bname);
	}
	else
	    crfile = NULL;
	init_output_data(init,&cr_data->odata,grid,prt,
			 "cross sectional data",YES,YES,NO);
	add_user_output_function(print_cross_sectional_graphs,
				 &cr_data->odata,prt);

	cr_sec = &Cr_sec;
	screen("Enter the %d coordinates of the start and end points\n"
	       "\tof each cross section, followed by the number of\n"
	       "\tpoints on the cross section (at least two)\n"
	       "\tan optional arclength offset from the start point and an\n"
               "\toptional whitespace-free message (<= 20 char.)\n"
               "\tSetting the start (end) point to coincide with the lower\n"
               "\t(upper) mesh limit in a coordinate direction, the\n",dim);
	screen("\tnumber of points equal to the number of mesh cells in that\n"
	       "\tdirection and omitting the offset will produce output at\n"
	       "\tthe centers of the primal mesh cells.  This is usually what\n"               "\tis desired.\n"

               "\tSpecifying an offset yields num_points+1 evenly spaced\n"
	       "\toutput points from start+offset to end, inclusive.\n"
	       "\tThis allows for finer control of the output point\n"
	       "\tlocations within the mesh cells, see init_cross_sections()\n"
	       "\tfor details.\n"
	       "\tNOTE -- One cross section per input line.\n");
	for (i = 0; i < num_cr_sec; ++i)
	{
	    scalar(&cr_sec->next,sizeof(CROSS_SECTION));
	    cr_sec->next->prev = cr_sec;
	    cr_sec = cr_sec->next;
	    cr_sec->next = NULL;
	    screen("Enter cross section data: ");
	    Gets(s);
	    /*
	     * !!!NOTE!!! this loop produces a vector of string pointers
	     * (tokens) into the storage for string variable 's'.
	     * 's' MUST NOT be modified until we are finished processing the
	     * tokens uni_array.
	     */
	    for (nToks = 0, tokens[nToks] = strtok(s,whiteSpace);
	         tokens[nToks] != NULL;
		 nToks +=1, tokens[nToks] = strtok(NULL,whiteSpace))
                /* NULL loop body */ ;

	    for (j = 0; j < dim; ++j)
	    {
		(void) sscan_float(tokens[j],&L[j]);
		(void) sscan_float(tokens[dim+j],&U[j]);
	    }
	    (void) sscanf(tokens[2*dim],"%d",&num_pts);

	    offset = 0.0;
	    if (nToks == 2*dim+1)
	    {
	        offset = -1.0;
                cr_sec->message = "";
            }
	    else
	    {
	        int     twoDp1 = 2*dim+1;

		k = sscan_float(tokens[twoDp1],&offset);
		if ((k  == 0)  || (k == EOF))
		{
		    /* string is the message field and offset is not present */
		    offset = -1.0;
		    cr_sec->message = strdup(tokens[twoDp1]);
		}
		else
		{
		    /* offset IS present */
		    if (offset < 0.0)
		    {
		        (void) printf("\nWARNING in init_cross_sections(),"
			              "\noffset interpreted as arclength ... "
                                      "changing sign to positive.");
	                offset = -offset;
		    }
		    cr_sec->message = (nToks > (twoDp1+1)) ?
		        strdup(tokens[twoDp1+1]) : "";
		}
	    }
	    /* Finished parsing input */

	    inc_num = (offset < 0.0) ? 2 : 1;

	    num_pts += inc_num;
	    cr_sec->num_pts = num_pts;

	    bi_array(&cr_sec->pcoords,num_pts,dim,FLOAT);
	    for (j = 0; j < dim; ++j)
		dir[j] = (U[j] - L[j]);
	    len = mag_vector(dir,dim);
	    if (len == 0.0)
	    {
		screen("ERROR in init_cross_sections(), "
		       "the cross-section length is zero.");
		clean_up(ERROR);
	    }
	    for (j = 0; j < dim; ++j)
		dir[j] /= len;

	    if (offset < 0.0)
	    {
	        /* Original algorithm, offset ignored */

	        /* Define num_pts to be the number of INTERIOR points of the
	         * cross section.  This has the advantage that if
	         * the user chooses the number of points in the cross
	         * section to be the same as the number of grid cells
	         * in that direction, and the cross section covers the
	         * entire domain in that direction, the points of the
	         * cross section will be at the cell centers, evenly spaced,
	         * with half-sized intervals at the endpoints.
	         */

		 dlen = len/(num_pts - inc_num);

		 for (j = 0; j < dim; ++j)
		 {
		     /* Do endpoint and first point (half-spaced) */
		     cr_sec->pcoords[0][j] = L[j];
		     cr_sec->pcoords[1][j] = L[j] + 0.5*dir[j]*dlen;
	         }
		 len = 0.5*dlen;
		 for (k = 2; k < num_pts - 1; ++k)
		 {
		     len += dlen;
		     for (j = 0; j < dim; ++j)
		         cr_sec->pcoords[k][j] = L[j] + dir[j]*len;
	         }
                 for (j = 0; j < dim; ++j)
		     cr_sec->pcoords[num_pts - 1][j] = U[j];
	    }
	    else
	    {
	        /*              Modified algorithm.
		 *
		 * Define num_pts+1 to be evenly-spaced from (start point
		 * + offset) to (end point) inclusive.  E.g., selecting
		 * start and end to agree with the computational mesh in a
		 * direction, an offset of zero and num_pts = 2x(num mesh
		 * cells) will give output at all mesh cell edges and centers.
		 */

	        for (j = 0; j < dim; ++j)
		    L[j] += dir[j]*offset;
	        len -= offset;
	        dlen = len/(num_pts - inc_num);
		    /* Do start point */
	        for (j = 0; j < dim; ++j)
		    cr_sec->pcoords[0][j] = L[j];

	        len = 0.0;
	        for (k = 1; k < num_pts; ++k)
	        {
		    len += dlen;
		    for (j = 0; j < dim; ++j)
	    	        cr_sec->pcoords[k][j] = L[j] + dir[j]*len;
	        }
	    }
	}
	screen("\n");
	cr_data->cross_sections = Cr_sec.next;
	cr_data->cross_sections->prev = NULL;
}		/*end init_cross_sections*/
#endif /* defined(TWOD) */


#if defined(THREED)

/*ARGSUSED*/
LOCAL	void print_geomview_plots(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	char s[512];

	if (prt->outfile == NULL)
	    return;

	(void) sprintf(s,"%s/%s-gv.ts%s",
		       Output_dir(data),basename(prt->outfile),
		       right_flush(grid->step,TSTEP_FIELD_WIDTH));
	if (front->pp_grid->nn > 1)
	    (void) sprintf(s,"%s-nd%d",s,pp_mynode());
	gview_plot_interface(s,front->interf);
}		/*end print_geomview_plots*/
#endif /* defined(THREED) */


EXPORT void init_peep_time(
	INIT_DATA	*init,
	Grid		*grid,
	Printplot	*prt)
{
	OUTPUT_DATA	*odata;
	char		s[Gets_BUF_SIZE];

	screen("Type 'y' to request a periodic glimpse of the solution\n\t");
	screen("via a plot of the component regions: ");
	(void) Gets(s);
	if ((s[0] == 'y') || (s[0] == 'Y'))
	{
	    scalar(&odata,sizeof(OUTPUT_DATA));
	    init_output_data(init,odata,grid,prt,"the peep data",YES,YES,NO);
	    add_user_output_function(give_peep,odata,prt);
	}
}		/*end init_peep_time*/


/* needed for VTK */
EXPORT	void	prompt_for_vtk_plots(
	INIT_DATA	  *init)
{
	VTK_PRINT_OPTIONS *opts;
	PRINT_OPTIONS	  *vtk_opts;
	VTK_FRAME_OPTS    *frame_opts;
	VTK_FRAME_OPTS    *fopts;	
	Plot_choice	  *pc;
	ScalarPlotItem	  *spi;
	char	          *c, *cv, s[1024];
	const char	  *uc_type_name;
	const char	  *lc_type_name;
	char	          *selected_plots;
	char		  *selected_vec_plots;
	char		  backupselect[10000],backupvecselect[10000];

	static const char	  *separators = " \t";
	double	          max_len, resolution;
	int               prt_offset;
	int 		  factor;
	int               nvars,nvvars,group,vpos;
	int	          i, var, veccount,newcount; 
	int		  dim = i_intfc(init)->dim;
	static const char *dname[] = {"error", "x", "x and y", 
			  	      "x, y, and z"};
	static const char *fmt = "%lf %lf %lf";
	prt_offset = VTK_STATES;
	opts = &VTK_prt_options(init);
	screen("\n\t\tVTK plotting initialization\n");
	zero_scalar(opts,sizeof(VTK_PRINT_OPTIONS));

	(void) sprintf(s,"\nSpecify variables to be plotted (VTK).");

	print_plotting_choices(s,plot_choices(init));
	selected_plots = read_plotting_choices(NULL,init);
	nvars = count_selected_plot_choices(selected_plots,
			plot_choices(init));
	vtk_num_vars(opts) = nvars;
	sprintf(vtk_vars(opts),selected_plots);
	
	screen("\nFrom the above vtk variables choosen, specify in an"
	       " ordered list the \nvaribles which you would like to "
	       "be put in a vector file.  For example,\nif you would "
	       "like a velocity vector your would type \"VX VY\". "
	       "Please \nremember that you MUST INCLUDE THE VARIBLES "
	       "THAT YOU WANT PRINTED AS\nVECTORS IN THE ABOVE LIST "
	       "AS WELL.");
	selected_vec_plots = read_plotting_choices(NULL,init);
	nvvars= count_selected_plot_choices(selected_vec_plots,
			plot_choices(init));
        
	/*Sub-sampling prompt. */
	screen("\nYou may specify a sub-sampling factor (a factor"
	" of 2 would sample half\n        the data points in each dimension uniformly)."
	" Your sub-sampling factor must\n        divide evenly into the number" 
	" of data points contained in each dimension of\n        each processor." 
	" Failing this criteria will cause the program to exit.\n       "
	" If you do not want to sub-sample, enter 1.");
        screen("\n        Specify a sub-sampling factor: ");
	(void) Scanf("%d\n", &factor);
	opts->subsample_factor = factor;

	if (nvars > 0)
	{
	    sprintf(backupselect,"%s",selected_plots);
	    sprintf(backupvecselect,"%s",selected_vec_plots);
	    frame_opts = (VTK_FRAME_OPTS*)Store(nvars*sizeof(VTK_FRAME_OPTS));
	    VTK_frame_options(init) = frame_opts;
	    for (var = 0, c = strtok(selected_plots,separators); c != NULL;
	                          c = strtok(NULL,separators), ++var)
	   {
		     (void) strcpy(VTK_selector(frame_opts[var]),c);
	   }
	


	  for (var = 0; var < nvars; var++) 
	    {
		    VTK_FRAME_OPTS    *fopts = frame_opts+var;
		for (pc = plot_choices(init); pc != NULL; pc = pc->next)
		{
			if (SelectedScalarItems(VTK_selector(frame_opts[var]),pc,&spi) != 0)
			{
				for(; spi != NULL; spi = spi->next)
				{
			           (void) strcpy(vtk_plot_name(fopts),spi->name);
			           vtk_plot_function(fopts) = spi->plot_fn;
			           vtk_plot_filter(fopts) = identity_filter;
				}
				break;
			}
			if (pc == NULL)
			{
			      (void) printf("WARNING in prompt_for_prostar_plots(), "
			                    "NO SUCH CHOICE %s\n",c);
			}

		}
	    }
	 }
	 
	    (void) sprintf(s,"\nYou will now be prompted for a base "
	                     "file name and optional directory for the "
			     "VTK output. Output for each variable is to "
			     "a separate file whose name contains "
	                     "the base name, and the prompt string for "
			     "that variable.\n");
	    screen_print_long_string(s);

	    vtk_opts = prt_opts(init) + prt_offset;
	    
	    /*Set printing in binary on or off. */
	    if(print_in_binary(vtk_opts))
	         opts->_print_in_binary = TRUE;
            else
	         opts->_print_in_binary = FALSE;
 
	    
	    if (output_filename(init)[0] != '\0')
	    {
	        char *bname, *dirname;
	        base_and_dir_name((const char*)output_filename(init),
	    		          &dirname,&bname);
		strcpy(opts->base_name, bname);
	        (void) sprintf(print_filename(vtk_opts),"%s%s%s/%s/%s",
	    	                      (strlen(dirname)!=0) ? dirname : "",
	    			      (strlen(dirname)!=0) ? "/" : "",
				      "vtk",
	    	                      bname,bname);
	    }
	    else
	    {
	        strcpy(print_filename(vtk_opts),"NO_DEFAULT");
		strcpy(opts->base_name, "");
	    }	
	    (void) sprintf(s,"VTK data");
	    prompt_for_output_filenames(vtk_opts,YES,s);
	    sprintf(vtk_base(opts),print_filename(vtk_opts));
	    screen("\n\t\tEnd VTK plotting initialization\n\n");

}		/*end prompt_for_vtk_plots*/



EXPORT void init_vtk_plots(
	INIT_DATA	*init,
	Printplot	*prt)
{
	VTK_plot_data	  *vtk_data;
	VTK_PRINT_OPTIONS *opts;
	VTK_FRAME_OPTS    *frame_opts;
	OUTPUT_DATA	  *odata;
	const char	  *suffix;
	int     	  prt_offset;
	int		  i, var;
	int		  dim = i_intfc(init)->dim;
	int		  nvars;
	char 		  *plot_name;
	    prt_offset= VTK_STATES;
	    opts = &VTK_prt_options(init);
	    nvars = vtk_num_vars(opts);
	    if (nvars == 0)
	    {
	        prt->n_VTK_vars = 0;
	        prt->vtk_data = NULL;
	        return;
	    }
	    prt->n_VTK_vars = 1;
	    scalar(&vtk_data,sizeof(VTK_plot_data));
	    VTK_print_opts(vtk_data) = VTK_prt_options(init);
	    frame_opts = VTK_frame_options(init);
	    uni_array(&vtk_data->frame_data,nvars,sizeof(VTK_frame_data));
	    for (var = 0; var < nvars; var++)
            {
            	VTK_frame_opts(vtk_data->frame_data[var]) =
                	frame_opts[var];
            }


	    
	    prt->vtk_data = vtk_data;
	    prt_offset = VTK_STATES;
	    suffix = "";
	    vtk_data->dim = dim;
	    odata = prt->main_output_data[prt_offset];
	    PrtOpts(odata) = prt_opts(init)[prt_offset];
	    create_output_data_file(odata);

	    for (var = 0; var < nvars; ++var)
	   {
	   /*plot_name = VTK_frame_plot_name(vtk_data->frame_data[var]); */
	  /*      uni_array(&vtk_data->frame_data[var].file_name, */
	   /*     strlen(plot_name)+strlen(Output_filename(odata))+50,CHAR); */
	    /*    (void) sprintf(vtk_data->frame_data[var].file_name, */
	     /*                  "%s",Output_filename(odata)); */
	   }
}		/*end init_vtk_plots*/



LOCAL	void	prompt_for_vtk_frame_data(
	VTK_PRINT_OPTIONS *opts,
	VTK_FRAME_OPTS	  *fopts,
	ScalarPlotItem	  *spi)
{
	char s[Gets_BUF_SIZE];
	int  i;
	static const char *fmt = "%lf %lf";

	screen("from prompt_for_vtk_frame_data plot name is = %s \n ",spi->name);
	/*(void) strcpy(vtk_plot_name(fopts),spi->name); */
	/*vtk_plot_function(fopts) = spi->plot_fn; */
	/*vtk_plot_filter(fopts) = identity_filter; */


}


#if defined(__GD__)
EXPORT void init_gd_movie_plots(
	INIT_DATA	*init,
	Printplot	*prt)
{
	GD_PRINT_OPTIONS 	*opts;
	GD_FRAME_OPTS        	*frame_opts;
	GD_plot_data 		*GD_data;
	int 			nvars,var,dim;
	int     	  	prt_offset;

	opts = &GD_prt_options(init);
	nvars = gd_num_vars(opts);
	if (nvars == 0)
	{
	    prt->n_GD_vars = 0;
	    prt->GD_data = NULL;
	    return;
	}
	prt->n_GD_vars = 1;
	scalar(&GD_data,sizeof(GD_plot_data));
	GD_print_opts(GD_data) = GD_prt_options(init);
	frame_opts = GD_frame_options(init);
	prt->GD_data = GD_data;
	prt_offset = GD_MOVIE;
	GD_data->dim = dim;
	uni_array(&GD_data->frame_opts,nvars,sizeof(GD_FRAME_OPTS));
	for (var = 0; var < nvars; var++)
	{
	    GD_data->frame_opts[var] = frame_opts[var];
	}
}	/* end init_gd_movie_plots */
#endif /* defined(__GD__) */
