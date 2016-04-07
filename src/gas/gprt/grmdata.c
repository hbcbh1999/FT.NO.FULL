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
*				grmdata.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains diagnostic routines for Richtmyer-Meshkov simulations.
*
*/


#define	DEBUG_STRING	"rmdata"
#include <gdecs/gdecs.h>


struct	_POSN_DATA {
	double	position[MAXD];
	double	velocity;
	double	value;
};
typedef struct _POSN_DATA POSN_DATA;

struct	_VEL_DATA {
	char	  *header;	/* Field header			*/
	POSN_DATA *intfc_max;	/* data at interface max	*/
	POSN_DATA *intfc_min;	/* data at interface min	*/
};
typedef	struct _VEL_DATA VEL_DATA;

struct	_RM_DATA {
	GraphUnits	*units;
	double	T;			/* Time Scale			*/
	double	L;			/* Length Scale			*/
	double	V;			/* V = L/T,  velocity scale	*/

	POSN_DATA	Max;		/* Maximum sampled data		*/
	POSN_DATA	Min;		/* Minimum sampled data		*/
	double	t_max, t_min;		/* Time limits			*/
	double	a_max, a_min;		/* Amplitude limits		*/
	double	adot_max, adot_min;	/* Growth rate max/min		*/

	int	alloc_len;		/* Maximum length of arrays	*/
	int	active_len;		/* Array elements used		*/
	int	dim;
	int	nvars;
	int	layer_index;

	double	(*_PositionMetric)(double*,struct _RM_DATA*);
	double	(*_VelocityProjection)(double*,double*,double*,struct _RM_DATA*);

	double	*record_time;		/* Array of times		*/
	VEL_DATA	*global;	/* Data at interface extrema	*/
	boolean    first;
	long int   limits;
	long int   min_values_offset;
	long int   max_values_offset;
	long int   next_rm_data_offset;
	long int   print_next_rm_data_offset_here;
	double	   *values;
	boolean    skip;
};
typedef struct _RM_DATA	RM_DATA;

struct _RM_OUTPUT_DATA {
	OUTPUT_DATA	Data;
	int		layer_index;
	RM_DATA		*rm_data;
};
typedef struct _RM_OUTPUT_DATA RM_OUTPUT_DATA;

#define	Rm_output_data(data)	((RM_OUTPUT_DATA*)(data))

#define	PositionMetric(coords,rm_data)					\
	(*rm_data->_PositionMetric)(coords,rm_data)

#define	VelocityProjection(v,coords,nor,rm_data)			\
	(*rm_data->_VelocityProjection)(v,coords,nor,rm_data)

#define	interface_min(rm_data)						\
	(rm_data)->global->intfc_min

#define	interface_max(rm_data)						\
	(rm_data)->global->intfc_max

#define	time_scale(rm_data)	(rm_data)->units->_time_scale
#define	initial_tmax(rm_data)	(rm_data)->units->_initial_tmax
#define	initial_tmin(rm_data)	(rm_data)->units->_initial_tmin

#define	length_scale(rm_data)	(rm_data)->units->_length_scale
#define	initial_lmax(rm_data)	(rm_data)->units->_initial_lmax
#define	initial_lmin(rm_data)	(rm_data)->units->_initial_lmin

#define	velocity_scale(rm_data)	(rm_data)->units->_velocity_scale
#define	initial_vmax(rm_data)	(rm_data)->units->_initial_vmax
#define	initial_vmin(rm_data)	(rm_data)->units->_initial_vmin

LOCAL const char *SCALES_HEADER = "SCALES FOR UNITS";
LOCAL const char *MINS_HEADER   = "POSITION AND VALUES OF MINIMUM FIELDS";
LOCAL const char *MAXS_HEADER   = "POSITION AND VALUES OF MAXIMUM FIELDS";
LOCAL const char *FMT	        = "%-"FFMT;

	/* LOCAL Function Declarations */
LOCAL	FILE	*open_rm_data_file(OUTPUT_DATA*);
LOCAL	RM_DATA	*alloc_RM_DATA(RECT_GRID*,int,double(*)(double*,RM_DATA*),
			       double(*)(double*,double*,double*,RM_DATA*));
LOCAL	VEL_DATA	*alloc_VEL_DATA(const char*,RM_DATA*);
LOCAL	boolean	set_global_vel_data(RM_DATA*,Front*);
LOCAL	double	VerticalGrowthRate(double*,double*,double*,RM_DATA*);
LOCAL	double	RadialGrowthRate(double*,double*,double*,RM_DATA*);
LOCAL	double	amplitude(int,RM_DATA*);
LOCAL	double	growth_rate(int,RM_DATA*);
LOCAL	double	height(double*,RM_DATA*);
LOCAL	double	radius(double*,RM_DATA*);
LOCAL	void	communicate_extreme_value(POSN_DATA*,double);
LOCAL	void	print_max_data(FILE*,RM_DATA*,boolean);
LOCAL	void	print_min_data(FILE*,RM_DATA*,boolean);
LOCAL	void	print_rm_data_header(FILE*,RM_DATA*,boolean);
LOCAL	void	print_rm_data(OUTPUT_DATA*,RM_DATA*);
LOCAL	void	print_rm_data_limits(OUTPUT_DATA*,RM_DATA*);
LOCAL	void	print_rm_header_item(FILE*,const char*,RM_DATA*);
LOCAL	void	print_rm_value_row(OUTPUT_DATA*,int,double*);
LOCAL	void	print_vel_data_header(FILE*,VEL_DATA*,RM_DATA*);
LOCAL	void	record_rm_data(Grid*,Front*,OUTPUT_DATA*,boolean,
			       double(*)(double*,RM_DATA*),
			       double(*)(double*,double*,double*,RM_DATA*));

/*ARGSUSED*/
EXPORT void record_rm_amp_and_vel_data(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	DEBUG_ENTER(record_rm_amp_and_vel_data)
	record_rm_data(grid,front,data,about_to_stop,height,VerticalGrowthRate);
	DEBUG_LEAVE(record_rm_amp_and_vel_data)
}		/*end record_rm_amp_and_vel_data*/

/*ARGSUSED*/
EXPORT void record_radial_rm_amp_and_vel_data(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	DEBUG_ENTER(record_radial_rm_amp_and_vel_data)
	record_rm_data(grid,front,data,about_to_stop,radius,RadialGrowthRate);
	DEBUG_LEAVE(record_radial_rm_amp_and_vel_data)
}		/*end record_radial_rm_amp_and_vel_data*/

EXPORT	GraphUnits	*RmGraphUnits(
	OUTPUT_DATA     **datas)
{
	IO_TYPE IO_type;
	FILE	*file;
	static	GraphUnits	*units = NULL;

	if (units != NULL)
	    return units;

	scalar(&units,sizeof(GraphUnits));
	units->_time_scale = 1.0;
	units->_initial_tmin =  HUGE_VAL;
	units->_initial_tmax = -HUGE_VAL;
	units->_length_scale = 1.0;
	units->_initial_lmin =  HUGE_VAL;
	units->_initial_lmax = -HUGE_VAL;
	units->_velocity_scale = 1.0;
	units->_initial_vmin =  HUGE_VAL;
	units->_initial_vmax = -HUGE_VAL;

	if (datas == NULL)
	    return units;

	if (Output_file(datas[0]) != NULL)
	    return units;

	if ((file = fopen(Output_filename(datas[0]),"r")) == NULL)
	{
	    /*file doesn't currently exist, nothing more to do */
	    return units;
	}
	determine_io_type(file,&IO_type);

	if (next_output_line_containing_string(file,SCALES_HEADER) != NULL)
	{
	    units->_time_scale = fread_float("time_scale = ",&IO_type);
	    units->_length_scale = fread_float("length_scale = ",&IO_type);
	    units->_velocity_scale = fread_float("velocity_scale = ",&IO_type);
	}
	Fclose(file);
	return units;
}		/*end RmGraphUnits*/

LOCAL void record_rm_data(
	Grid		*grid,
	Front		*front,
	OUTPUT_DATA	*data,
	boolean		about_to_stop,
	double		(*position_metric)(double*,RM_DATA*),
	double		(*velocity_projection)(double*,double*,double*,RM_DATA*))
{
	FILE	*logfile;
	boolean status;
	boolean io_node = (is_io_node(pp_mynode())) ? YES : NO;
	double	time = grid->time;
	RM_DATA	*rm_data = Rm_output_data(data)->rm_data;

	DEBUG_ENTER(record_rm_data)
	if (rm_data == NULL) 
	{
	    int	max_step = stop_step(grid)+1;
	    if ((io_node == YES) && (Output_file(data) != stdout))
	    	max_step = 1;
	    pp_global_imin(&max_step,1L);
	    rm_data = alloc_RM_DATA(front->rect_grid,max_step,position_metric,
	    			    velocity_projection);
	    Rm_output_data(data)->rm_data = rm_data;
	    rm_data->layer_index = Rm_output_data(data)->layer_index;
	}
	if (pp_max_status(rm_data->skip) == YES)
	{
	    if (DEBUG)
		(void) printf("skipping rm data printout\n");
	    DEBUG_LEAVE(record_rm_data)
	    return;
	}
	if (DEBUG)
	    (void) printf("processing RM data\n");

	if ((logfile = open_rm_data_file(data)) == NULL)
	{
	    (void) printf("WARNING in record_rm_data(), "
			  "can't open rm data file %s\n",Output_filename(data));
	    return;
	}

	status = set_global_vel_data(rm_data,front);
	if (pp_global_status(status) != FUNCTION_SUCCEEDED)
	{
	    (void) printf("WARNING in record_rm_data(), "
	                  "set_global_vel_data() failed\n");
	    rm_data->skip = YES;
	    if (logfile != stdout)
	    {
	        if (is_io_node(pp_mynode()))
	        {
	    	    (void) Fclose(logfile);
	    	    Output_file(data) = NULL;
	        }
	        rm_data->active_len = 0; /*Reuse field on next call*/
	    }
	    DEBUG_LEAVE(record_rm_data)
	    return;
	}

	rm_data->record_time[rm_data->active_len++] = time;
	rm_data->t_max = max(time,rm_data->t_max);
	rm_data->t_min = min(time,rm_data->t_min);

	if (logfile != stdout)
	{
	    if (is_io_node(pp_mynode()))
	    {
	    	print_rm_data(data,rm_data);
	    	if (about_to_stop == YES)
	    	{
	    	    (void) Fclose(logfile);
	    	    Output_file(data) = NULL;
	    	}
	    }
	    rm_data->active_len = 0; /*Reuse field on next call*/
	}
	else if (about_to_stop == YES)
	{
	    print_rm_data(data,rm_data);
	    (void) printf("\n");
	}
	DEBUG_LEAVE(record_rm_data)
}		/*end record_rm_data*/

LOCAL	FILE	*open_rm_data_file(
	OUTPUT_DATA	*data)
{
	RM_DATA	   *rm_data = Rm_output_data(data)->rm_data;
	IO_TYPE    IO_type;
	boolean       newfile, appendable;
	boolean       io_node = (is_io_node(pp_mynode())) ? YES : NO;
	boolean       output_file_exists;
	char	   *fname = Output_filename(data);
	const char *line;
	int	   i, dim = rm_data->dim;
	FILE	   *file;

	DEBUG_ENTER(record_rm_data)
	output_file_exists = NO;
	if ((io_node == YES) && (Output_file(data) != NULL))
	    output_file_exists = YES;
	output_file_exists = pp_max_status(output_file_exists);
	if (output_file_exists == YES)
	{
	    DEBUG_LEAVE(record_rm_data)
	    if (DEBUG)
	        (void) printf("Output file exits, done\n");
	    return Output_file(data);
	}

	if (DEBUG)
	    (void) printf("Output file == NULL searching for existing file\n");

	newfile = YES;
	if (io_node == YES)
	{
	    if ((file = fopen(fname,"r")) == NULL)
	    {
	        /*file doesn't exit open new file for writing*/
	        if ((Output_file(data) = file = fopen(fname,"w+")) == NULL)
		{
	            screen("ERROR in open_rm_data_file(), can't open %s\n",
			   fname);
	            clean_up(ERROR);
	            DEBUG_LEAVE(record_rm_data)
	            return NULL;
		}
	        setbuf(file,NULL);
		print_machine_parameters(Output_file(data));
	    }
	    else
	    {
	        newfile = NO;
	        (void) fclose(file);
	    }
	}
	newfile = pp_min_status(newfile);
	if (newfile == YES)
	{
	    if (DEBUG)
	        (void) printf("new file opened done\n");
	    DEBUG_LEAVE(record_rm_data)
	    return Output_file(data);
	}

	appendable = YES;
	file = NULL;
	if (io_node == YES)
	{
	    file = Output_file(data) = fopen(fname,"r+");
	    if (file == NULL)
	    {
	        screen("ERROR in open_rm_data_file(), "
		       "can't open %s for append\n",fname);
	        clean_up(ERROR);
	        DEBUG_LEAVE(record_rm_data)
	        return NULL;
	    }
	    setbuf(file,NULL);

	    /*Read Max/Min data from file to initialize data fields*/
	    if (!append_output(file))
		appendable = NO;
	}
	appendable = pp_min_status(appendable);
	if (appendable == NO)
	{
	    /*Can't append to existing file*/
	    (void) printf("WARNING in open_rm_data_file(), "
		          "can't append to %s\n",fname);
	    if (file != NULL)
	        (void) fclose(file);
	    Output_file(data) = NULL;
	    rm_data->skip = YES;
	    DEBUG_LEAVE(record_rm_data)
	    return NULL;
	}

	if (io_node == NO)
	{
	    file = Output_file(data) = fopen(fname,"r");
	    if (file == NULL)
	    {
	        screen("ERROR in open_rm_data_file(), "
		       "can't open %s for read on non-io_node\n",fname);
	        clean_up(ERROR);
	        DEBUG_LEAVE(record_rm_data)
	        return NULL;
	    }
	}
	rm_data->first = NO;
	determine_io_type(file,&IO_type);

	(void) rewind_read_file(file,NULL);
	if (next_output_line_containing_string(file,MINS_HEADER) != NULL)
	{
	    rm_data->min_values_offset = ftell(file);
	    rm_data->t_min = fread_float("time = ",&IO_type);
	    rm_data->a_min = fread_float("amplitude = ",&IO_type);
	    rm_data->adot_min = fread_float("growth rate = ",&IO_type);
	    rm_data->Min.position[0] = fread_float("position = ",&IO_type);
	    for (i = 1; i < dim; ++i)
	        rm_data->Min.position[i] = fread_float(" ",&IO_type);
	    rm_data->Min.velocity = fread_float("velocity = ",&IO_type);
	}
	else
	    rm_data->first = YES;
	if (next_output_line_containing_string(file,MAXS_HEADER) != NULL)
	{
	    rm_data->max_values_offset = ftell(file);
	    rm_data->t_max = fread_float("time = ",&IO_type);
	    rm_data->a_max = fread_float("amplitude = ",&IO_type);
	    rm_data->adot_max = fread_float("growth rate = ",&IO_type);
	    rm_data->Max.position[0] = fread_float("position = ",&IO_type);
	    for (i = 1; i < dim; ++i)
	        rm_data->Max.position[i] = fread_float(" ",&IO_type);
	    rm_data->Max.velocity = fread_float("velocity = ",&IO_type);
	}
	else
	    rm_data->first = YES;

	line = next_output_line_containing_string(file,
						  "next_rm_data_offset = ");
	if (line != NULL)
	    (void) sscanf(line,"%*s %*s %ld %*s %*s %ld",
			       &rm_data->next_rm_data_offset,
			       &rm_data->print_next_rm_data_offset_here);
	else
	    rm_data->first = YES;

	if (fgetstring(file,"nvars = ") == FUNCTION_SUCCEEDED)
	    (void) fscanf(file,"%d",&rm_data->nvars);
	else
	    rm_data->first = YES;
	
	if (next_output_line_containing_string(file,"time") != NULL)
	    rm_data->limits = ftell(file);
	else
	    rm_data->first = YES;

	rm_data->first = pp_max_status(rm_data->first);

	if (io_node == YES)
	{
	    if (rm_data->next_rm_data_offset > 0)
	        (void) fseek(file,rm_data->next_rm_data_offset,SEEK_SET);
	}
	else
	    (void) fclose(file);

	DEBUG_LEAVE(record_rm_data)
	return file;
}		/*end open_rm_data_file*/

LOCAL	double	height(
	double	*coords,
	RM_DATA	*rm_data)
{
	return coords[rm_data->dim-1];
}		/*end height*/

/*ARGSUSED*/
LOCAL	double	VerticalGrowthRate(
	double	*v,
	double	*coords,
	double	*nor,
	RM_DATA	*rm_data)
{
	return scalar_product(v,nor,rm_data->dim)*nor[rm_data->dim-1];
}		/*end VerticalGrowthRate*/

LOCAL	double	radius(
	double	*coords,
	RM_DATA	*rm_data)
{
	return mag_vector(coords,rm_data->dim);
}		/*end height*/

/*ARGSUSED*/
LOCAL	double	RadialGrowthRate(
	double	*v,
	double	*coords,
	double	*nor,
	RM_DATA	*rm_data)
{
	int dim = rm_data->dim;
	double	r = radius(coords,rm_data);

	if (r == 0.0)
	    return 0.0;
	return scalar_product(v,nor,dim)*scalar_product(nor,coords,dim)/r;
}		/*end VerticalGrowthRate*/

LOCAL	RM_DATA	*alloc_RM_DATA(
	RECT_GRID	*gr,
	int		alloc_len,
	double		(*position_metric)(double*,RM_DATA*),
	double		(*velocity_projection)(double*,double*,double*,RM_DATA*))
{
	RM_DATA	*rm_data;
	int	i, dim = gr->dim;

	DEBUG_ENTER(alloc_RM_DATA)
	scalar(&rm_data,sizeof(RM_DATA));
	rm_data->first = YES;
	rm_data->limits = -1L;
	rm_data->values = NULL;
	rm_data->units = RmGraphUnits(NULL);
	for (i = 0; i < dim; ++i)
	{
	    rm_data->Max.position[i] = initial_lmax(rm_data);
	    rm_data->Min.position[i] = initial_lmin(rm_data);
	}
	rm_data->Max.velocity = initial_vmax(rm_data);
	rm_data->Min.velocity = initial_vmin(rm_data);
	rm_data->Max.value = 0.0;/*Unused*/
	rm_data->Min.value = 0.0;/*Unused*/
	rm_data->t_max = initial_tmax(rm_data);
	rm_data->t_min = initial_tmin(rm_data);
	rm_data->a_max = initial_lmax(rm_data);
	rm_data->a_min = initial_lmin(rm_data);
	rm_data->adot_max = initial_vmax(rm_data);
	rm_data->adot_min = initial_vmin(rm_data);
	rm_data->alloc_len = alloc_len;
	rm_data->active_len = 0;
	rm_data->dim = dim;
	rm_data->_PositionMetric = position_metric;
	rm_data->_VelocityProjection = velocity_projection;
	rm_data->min_values_offset = -1L;
	rm_data->max_values_offset = -1L;
	rm_data->next_rm_data_offset = -1L;
	rm_data->print_next_rm_data_offset_here = -1L;
	uni_array(&rm_data->record_time,rm_data->alloc_len,FLOAT);
	rm_data->global = alloc_VEL_DATA("intfc",rm_data);
	rm_data->skip = NO;
	DEBUG_LEAVE(alloc_RM_DATA)
	return rm_data;
}		/*end alloc_RM_DATA*/

LOCAL	VEL_DATA	*alloc_VEL_DATA(
	const char	*header,
	RM_DATA		*rm_data)
{
	VEL_DATA	*vel_data;

	DEBUG_ENTER(alloc_VEL_DATA)
	scalar(&vel_data,sizeof(VEL_DATA));
	if (header == NULL)
	    header = "";
	scalar(&vel_data->header,strlen(header)+1);
	(void) strcpy(vel_data->header,header);

	uni_array(&vel_data->intfc_max,rm_data->alloc_len,sizeof(POSN_DATA));
	uni_array(&vel_data->intfc_min,rm_data->alloc_len,sizeof(POSN_DATA));

	DEBUG_LEAVE(alloc_VEL_DATA)
	return vel_data;
}		/*end alloc_VEL_DATA*/


LOCAL	void print_rm_data(
	OUTPUT_DATA	*data,
	RM_DATA		*rm_data)
{
	FILE		*file = Output_file(data);
	double		*values;
	int		n, nvars;

	DEBUG_ENTER(print_rm_data)
	if (rm_data->first == YES)
	{
	    rm_data->first = NO;
	    print_rm_data_header(file,rm_data,Output_in_binary(data));
	    print_rm_data_limits(data,rm_data);
	}
	else
	    print_rm_data_limits(data,rm_data);

	if (rm_data->values == NULL)
	    uni_array(&rm_data->values,rm_data->nvars,FLOAT);

	if (rm_data->next_rm_data_offset > 0)
	{
	    if (!erase_last_foutput(file))
	    {
		screen("ERROR in print_rm_data(), can't erase last foutput\n");
		clean_up(ERROR);
	    }
	    (void) fseek(file,rm_data->next_rm_data_offset,SEEK_SET);
	}

	values = rm_data->values;
	for (n = 0; n < rm_data->active_len; ++n)
	{
	    int i, dim = rm_data->dim;

	    nvars = 0;
	    values[nvars++] = rm_data->record_time[n]/time_scale(rm_data);
	    values[nvars++] = amplitude(n,rm_data)/length_scale(rm_data);
	    values[nvars++] = growth_rate(n,rm_data)/
				       velocity_scale(rm_data);
	    for (i = 0; i < dim; ++i)
	        values[nvars++] = interface_max(rm_data)[n].position[i] /
						length_scale(rm_data);
	    values[nvars++] = interface_max(rm_data)[n].velocity /
						velocity_scale(rm_data);
	    for (i = 0; i < dim; ++i)
	        values[nvars++] = interface_min(rm_data)[n].position[i] /
	    				length_scale(rm_data);
	    values[nvars++] = interface_min(rm_data)[n].velocity /
	    				velocity_scale(rm_data);
	    print_rm_value_row(data,nvars,values);
	}
	rm_data->next_rm_data_offset = ftell(file);
	(void) fseek(file,rm_data->print_next_rm_data_offset_here,SEEK_SET);
	(void) fprintf(file,"%9ld",rm_data->next_rm_data_offset);
	(void) fseek(file,rm_data->next_rm_data_offset,SEEK_SET);

	trace_foutput(file);
	DEBUG_LEAVE(print_rm_data)
}		/*end print_rm_data*/

LOCAL	void	print_rm_value_row(
	OUTPUT_DATA	*data,
	int	nvars,
	double	*values)
{
	FILE	*file = Output_file(data);

	if (Output_in_binary(data) == YES)
	{
	    (void) fprintf(file,"\f%c",nvars);
	    (void) fwrite((const void *)values,FLOAT,nvars,file);
	}
	else
	{
	    int i;
	    for (i = 0; i < nvars; ++i)
	    	(void) fprintf(file,"%-"FFMT,values[i]);
	    (void) fprintf(file,"\n");
	}
}		/*end print_rm_value_row*/


LOCAL	void	print_rm_data_limits(
	OUTPUT_DATA	*data,
	RM_DATA	*rm_data)
{
	FILE	*file = Output_file(data);
	long int current;
	int	i, dim = rm_data->dim, nvars;
	double	*values;

	print_min_data(file,rm_data,Output_in_binary(data));
	print_max_data(file,rm_data,Output_in_binary(data));
	if (rm_data->limits < 0)
	{
	    rm_data->limits = ftell(file);
	    current = -1L;
	}
	else
	{
	    current = ftell(file);
	    (void) fseek(file,rm_data->limits,SEEK_SET);
	}

	if (rm_data->values == NULL)
	    uni_array(&rm_data->values,rm_data->nvars,FLOAT);

	values = rm_data->values;
	nvars = 0;
	values[nvars++] = rm_data->t_min/time_scale(rm_data);
	values[nvars++] = rm_data->a_min/length_scale(rm_data);
	values[nvars++] = rm_data->adot_min/velocity_scale(rm_data);
	for (i = 0; i < dim; ++i)
	    values[nvars++] = rm_data->Min.position[i]/length_scale(rm_data);
	values[nvars++] = rm_data->Min.velocity/velocity_scale(rm_data);
	for (i = 0; i < dim; ++i)
	    values[nvars++] = rm_data->Min.position[i]/length_scale(rm_data);
	values[nvars++] = rm_data->Min.velocity/velocity_scale(rm_data);
	print_rm_value_row(data,nvars,values);

	nvars = 0;
	values[nvars++] = rm_data->t_max/time_scale(rm_data);
	values[nvars++] = rm_data->a_max/length_scale(rm_data);
	values[nvars++] = rm_data->adot_max/velocity_scale(rm_data);
	for (i = 0; i < dim; ++i)
	    values[nvars++] = rm_data->Max.position[i]/length_scale(rm_data);
	values[nvars++] = rm_data->Max.velocity/velocity_scale(rm_data);
	for (i = 0; i < dim; ++i)
	    values[nvars++] = rm_data->Max.position[i]/length_scale(rm_data);
	values[nvars++] = rm_data->Max.velocity/velocity_scale(rm_data);
	print_rm_value_row(data,nvars,values);

	if (current >= 0)
	    (void) fseek(file,current,SEEK_SET);
}		/*end print_rm_data_limits*/

LOCAL	void	print_rm_data_header(
	FILE		*file,
	RM_DATA		*rm_data,
	boolean		bio)
{
	long	current, nvars_loc;
	(void) foutput(file);
	(void) fprintf(file,"%s\n",SCALES_HEADER);
	fwrite_float(file,"time_scale = ",time_scale(rm_data),bio,FMT," ");
	fwrite_float(file,"length_scale = ",length_scale(rm_data),bio,FMT," ");
	fwrite_float(file,"velocity_scale = ",
		     velocity_scale(rm_data),bio,FMT,"\n");

	print_min_data(file,rm_data,bio);
	print_max_data(file,rm_data,bio);

	(void) foutput(file);
	(void) fprintf(file,"next_rm_data_offset = ");
	rm_data->print_next_rm_data_offset_here = ftell(file);
	(void) fprintf(file,"%9ld ",rm_data->next_rm_data_offset);
	(void) fprintf(file,"print_next_rm_data_offset_here = %ld\n",
			    rm_data->print_next_rm_data_offset_here);

	(void) fprintf(file,"nvars = ");
	nvars_loc = ftell(file);
	(void) fprintf(file,"%-9d\n",0);

	rm_data->nvars = 0;
	(void) foutput(file);
	print_rm_header_item(file,"time",rm_data);
	print_vel_data_header(file,rm_data->global,rm_data);
	(void) fprintf(file,"\n");

	current = ftell(file);
	(void) fseek(file,nvars_loc,SEEK_SET);
	(void) fprintf(file,"%-9d\n",rm_data->nvars);
	(void) fseek(file,current,SEEK_SET);

}		/*end print_rm_data_header*/

LOCAL	void	print_min_data(
	FILE	*file,
	RM_DATA	*rm_data,
	boolean	bio)
{
	int i, dim = rm_data->dim;
	long current;

	if (rm_data->min_values_offset < 0)
	{
	    (void) foutput(file);
	    (void) fprintf(file,"%s\n",MINS_HEADER);
	    rm_data->min_values_offset = ftell(file);
	    current = -1L;
	}
	else
	{
	    current = ftell(file);
	    (void) fseek(file,rm_data->min_values_offset,SEEK_SET);
	}
	fwrite_float(file,"time = ",rm_data->t_min,bio,FMT," ");
	fwrite_float(file,"amplitude = ",rm_data->a_min,bio,FMT," ");
	fwrite_float(file,"growth rate = ",rm_data->adot_min,bio,FMT,"\n");
	fwrite_float(file,"position = ",rm_data->Min.position[0],bio,FMT,"");
	for (i = 1; i < dim; ++i)
	    fwrite_float(file," ",rm_data->Min.position[i],bio,FMT,"");
	(void) fprintf(file,"\n");
	fwrite_float(file,"velocity = ",rm_data->Min.velocity,bio,FMT,"");
	(void) fprintf(file,"\n");

	if (current >= 0)
	    (void) fseek(file,current,SEEK_SET);
}		/*end print_min_data*/

LOCAL	void	print_max_data(
	FILE	*file,
	RM_DATA	*rm_data,
	boolean	bio)
{
	int i, dim = rm_data->dim;
	long current;

	if (rm_data->max_values_offset < 0)
	{
	    (void) foutput(file);
	    (void) fprintf(file,"%s\n",MAXS_HEADER);
	    rm_data->max_values_offset = ftell(file);
	    current = -1L;
	}
	else
	{
	    current = ftell(file);
	    (void) fseek(file,rm_data->max_values_offset,SEEK_SET);
	}
	fwrite_float(file,"time = ",rm_data->t_max,bio,FMT," ");
	fwrite_float(file,"amplitude = ",rm_data->a_max,bio,FMT," ");
	fwrite_float(file,"growth rate = ",rm_data->adot_max,bio,FMT,"\n");
	fwrite_float(file,"position = ",rm_data->Max.position[0],bio,FMT,"");
	for (i = 1; i < dim; ++i)
	    fwrite_float(file," ",rm_data->Max.position[i],bio,FMT,"");
	(void) fprintf(file,"\n");
	fwrite_float(file,"velocity = ",rm_data->Max.velocity,bio,FMT,"");
	(void) fprintf(file,"\n");
	if (current > 0)
	    (void) fseek(file,current,SEEK_SET);
}		/*end print_max_data*/

LOCAL	void	print_rm_header_item(
	FILE	   *file,
	const char *s,
	RM_DATA	   *rm_data)
{
	(void) fprintf(file,"%-17s ",s);
	++rm_data->nvars;
}		/*end print_rm_header_item*/

LOCAL	void	print_vel_data_header(
	FILE		*file,
	VEL_DATA	*vel_data,
	RM_DATA		*rm_data)
{
	char	   s[256];
	const char *header = vel_data->header;
	int        i, dim = rm_data->dim;
	static const char *varname[] = {"x", "y", "z"};

	(void) sprintf(s,"%s%s",header,"_amp");
	print_rm_header_item(file,s,rm_data);
	(void) sprintf(s,"%s%s",header,"_adot");
	print_rm_header_item(file,s,rm_data);

	for (i = 0; i < dim; ++i)
	{
	    (void) sprintf(s,"%s_%s_%s",header,"max",varname[i]);
	    print_rm_header_item(file,s,rm_data);
	}
	(void) sprintf(s,"%s_%s_%s",header,"max","vel");
	print_rm_header_item(file,s,rm_data);
	for (i = 0; i < dim; ++i)
	{
	    (void) sprintf(s,"%s_%s_%s",header,"min",varname[i]);
	    print_rm_header_item(file,s,rm_data);
	}
	(void) sprintf(s,"%s_%s_%s",header,"min","vel");
	print_rm_header_item(file,s,rm_data);
}		/*end print_vel_data_header*/

LOCAL	boolean	set_global_vel_data(
	RM_DATA		*rm_data,
	Front		*front)
{
	HYPER_SURF_ELEMENT    *hse, *hse_max, *hse_min;
	HYPER_SURF    *hs,  *hs_max,  *hs_min;
	INTERFACE    *intfc = front->interf;
	POSN_DATA    *intfc_max = interface_max(rm_data) + rm_data->active_len;
	POSN_DATA    *intfc_min = interface_min(rm_data) + rm_data->active_len;
	Locstate    sl, sr;
	POINT	    *p,   *p_max,   *p_min;
	boolean	    status;
	double	    vl[MAXD], vr[MAXD], nor[MAXD];
	double	    a, adot;
	int	    i, dim = front->rect_grid->dim;

	DEBUG_ENTER(set_global_vel_data)
	status = (rm_data->active_len >= rm_data->alloc_len) ?
			FUNCTION_FAILED : FUNCTION_SUCCEEDED;
	if (pp_global_status(status) == FUNCTION_FAILED)
	{
	    screen("ERROR in set_global_vel_data(), "
		   "active_len >= alloc_len\n"
	           "active_len = %d, alloc_len = %d\n",
	    	   rm_data->active_len,rm_data->alloc_len);
	    clean_up(ERROR);
	    DEBUG_LEAVE(set_global_vel_data)
	    return FUNCTION_FAILED;
	}

	p_max = p_min = NULL;
	hse_max = hse_min = NULL;
	hs_max = hs_min = NULL;
	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (!is_scalar_wave(wave_type(hs)))
	    	continue;
	    if (point_in_buffer(Coords(p),front->rect_grid) == YES)
	    	continue;
	    if (layer_index(hs) != rm_data->layer_index)
	    	continue;

	    if ((p_max==NULL) || (PositionMetric(Coords(p_max),rm_data) < 
				 PositionMetric(Coords(p),    rm_data))
	    )
	    {
	    	p_max = p;
	    	hse_max = hse;
	    	hs_max = hs;
	    }
	    if ((p_min==NULL) ||
	    		(PositionMetric(Coords(p_min),rm_data) >
	    		 PositionMetric(Coords(p),    rm_data))
	    )
	    {
	    	p_min = p;
	    	hse_min = hse;
	    	hs_min = hs;
	    }
	}
	if (p_max != NULL)
	{
	    intfc_max->value = PositionMetric(Coords(p_max),rm_data);
	    slsr(p_max,hse_max,hs_max,&sl,&sr);
	    normal(p_max,hse_max,hs_max,nor,front);
	    for (i = 0; i < dim; ++i)
	    {
	    	intfc_max->position[i] = Coords(p_max)[i];
	    	vl[i] = vel(i,sl);
	    	vr[i] = vel(i,sr);
	    }
	    intfc_max->velocity = 0.5*(
			VelocityProjection(vl,Coords(p_max),nor,rm_data) +
			VelocityProjection(vr,Coords(p_max),nor,rm_data));
	}
	else
	{
	    intfc_max->value = -HUGE_VAL;
	    for (i = 0; i < dim; ++i)
	    	intfc_max->position[i] = -HUGE_VAL;
	    intfc_max->velocity = -HUGE_VAL;
	}
	if (p_min != NULL)
	{
	    intfc_min->value = PositionMetric(Coords(p_min),rm_data);
	    slsr(p_min,hse_min,hs_min,&sl,&sr);
	    normal(p_min,hse_min,hs_min,nor,front);
	    for (i = 0; i < dim; ++i)
	    {
	    	intfc_min->position[i] = Coords(p_min)[i];
	    	vl[i] = vel(i,sl);
	    	vr[i] = vel(i,sr);
	    }
	    intfc_min->velocity = 0.5*(
			VelocityProjection(vl,Coords(p_min),nor,rm_data) +
			VelocityProjection(vr,Coords(p_min),nor,rm_data));
	}
	else
	{
	    intfc_min->value = HUGE_VAL;
	    for (i = 0; i < dim; ++i)
	    	intfc_min->position[i] = HUGE_VAL;
	    intfc_min->velocity = HUGE_VAL;
	}

	communicate_extreme_value(intfc_max, 1.0);
	communicate_extreme_value(intfc_min,-1.0);

	if ((intfc_max->value == -HUGE_VAL) || (intfc_min->value == HUGE_VAL))
	{
	    (void) printf("WARNING in set_global_vel_data(), failed to "
			  "compute global max or min for layer %d\n",
	                  rm_data->layer_index);
	    return FUNCTION_FAILED;
	}

	rm_data->Max.velocity = max(intfc_max->velocity,rm_data->Max.velocity);
	rm_data->Max.velocity = max(intfc_min->velocity,rm_data->Max.velocity);
	rm_data->Min.velocity = min(intfc_max->velocity,rm_data->Min.velocity);
	rm_data->Min.velocity = min(intfc_min->velocity,rm_data->Min.velocity);
	for (i = 0; i < dim; ++i)
	{
	    rm_data->Max.position[i] = max(intfc_max->position[i],
	    			           rm_data->Max.position[i]);
	    rm_data->Max.position[i] = max(intfc_min->position[i],
					       rm_data->Max.position[i]);
		rm_data->Min.position[i] = min(intfc_max->position[i],
					       rm_data->Min.position[i]);
		rm_data->Min.position[i] = min(intfc_min->position[i],
					       rm_data->Min.position[i]);
	}
	a = amplitude(rm_data->active_len,rm_data);
	rm_data->a_max = max(a,rm_data->a_max);
	rm_data->a_min = min(a,rm_data->a_min);
	adot = growth_rate(rm_data->active_len,rm_data);
	rm_data->adot_max = max(adot,rm_data->adot_max);
	rm_data->adot_min = min(adot,rm_data->adot_min);

	DEBUG_LEAVE(set_global_vel_data)
	return FUNCTION_SUCCEEDED;
}		/* set_global_vel_data */

LOCAL	double	amplitude(
	int	n,
	RM_DATA	*rm_data)
{
	return	0.5*(
		PositionMetric(interface_max(rm_data)[n].position,rm_data) -
		PositionMetric(interface_min(rm_data)[n].position,rm_data));
	      
}		/*end amplitude*/

LOCAL	double	growth_rate(
	int	n,
	RM_DATA	*rm_data)
{
	return	0.5*(interface_max(rm_data)[n].velocity -
		     interface_min(rm_data)[n].velocity);
}		/*end growth_rate*/


LOCAL	void	communicate_extreme_value(
	POSN_DATA	*p_data,
	double		sign)
{
	int		numnodes = pp_numnodes();
	int		i;
	static	POSN_DATA	*Remote = NULL;

	DEBUG_ENTER(communicate_extreme_value)
	if (numnodes == 1)
	{
	    DEBUG_LEAVE(communicate_extreme_value)
	    return;
	}

	if (Remote == NULL)
	    uni_array(&Remote,numnodes,sizeof(POSN_DATA));

	pp_all_gather((POINTER)p_data,sizeof(POSN_DATA),
		      (POINTER)Remote,sizeof(POSN_DATA));

	for (i = 0; i < numnodes; ++i)
	{
	    if (sign*p_data->value < sign*Remote[i].value)
	    	*p_data = Remote[i];
	}
	DEBUG_LEAVE(communicate_extreme_value)
}		/* set_global_vel_data */

EXPORT	OUTPUT_DATA	**rm_alloc_output_datas(
	int num_layers)
{
	OUTPUT_DATA	**datas;

	bi_array(&datas,num_layers,1,sizeof(RM_OUTPUT_DATA));
	return datas;
}		/*end rm_alloc_output_datas*/

EXPORT	void	set_rm_layer_indices(
	OUTPUT_DATA	**datas,
	int		num_layers)
{
	int		i;

	for (i = 0; i < num_layers; ++i)
	    Rm_output_data(datas[i])->layer_index = i+1;
}		/*end set_rm_layer_indices*/
