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

#include <gprt/glayer.h>

LOCAL   void record_moving_body_data(Grid*,Wave*,Front*,Printplot*,
                                OUTPUT_DATA*,boolean);

/*
*			init_moving_body_data():
*
*	Initializer for computation of moving body data.
*/

/*ARGSUSED*/
EXPORT	void init_moving_body_data(
	INIT_DATA	*init,
	Front		*front,
	Grid		*grid,
	Printplot	*prt)
{
	INTERFACE	*intfc = front->interf;
        const RECT_GRID	*rgrid = front->rect_grid;
        const int	dim = rgrid->dim;
	RGB_DATA	*rgb_data;
	char		s[Gets_BUF_SIZE];
	int		i,num_moving_body;
	CURVE		**c;
	char		*dname,*fname;
	char		torque_name[100],omega_name[100];
	char		force_name[100],com_name[100];

	if (dim == 1)
	    return;

	num_moving_body = -1;
	for (c = intfc->curves; c && *c; ++c)
        {
            if (wave_type(*c) == MOVABLE_BODY_BOUNDARY)
                if (body_index(*c) > num_moving_body)
                    num_moving_body = body_index(*c);
        }
        num_moving_body++;
        pp_global_imax(&num_moving_body,1);
	if (num_moving_body == 0) return;

	screen("Type 'y' to request moving body time data: ");
	(void) Gets(s);
	if ((s[0] != 'Y') && (s[0] != 'y'))
	    return;

	scalar(&rgb_data,sizeof(RGB_DATA));
	zero_scalar(rgb_data,sizeof(RGB_DATA));

	Output_mode(&rgb_data->odata) = EXACT_TIME;
	Output_time_freq(&rgb_data->odata) = 1.0;
	Output_start_time(&rgb_data->odata) = 0;
	Output_in_binary(&rgb_data->odata) = NO;

	init_output_data(init,&rgb_data->odata,grid,prt,
			 "moving body data",YES,NO,YES);
	add_user_output_function(record_moving_body_data,&rgb_data->odata,prt);

	if (pp_mynode() != 0) return;

	uni_array(&rgb_data->torque_files,num_moving_body,sizeof(FILE*));
	uni_array(&rgb_data->omega_files,num_moving_body,sizeof(FILE*));
	uni_array(&rgb_data->force_files,num_moving_body,sizeof(FILE*));
	uni_array(&rgb_data->com_files,num_moving_body,sizeof(FILE*));

	set_default_data_file_names("rgb_data","-torque",&dname,&fname,
						init);
	if (create_directory(dname,NO) == FUNCTION_FAILED)
        {
            screen("ERROR in open_data_file(), directory %s doesn't exist "
                   		"and can't be made\n",dname);
            clean_up(ERROR);
        }
	for (i = 0; i < num_moving_body; ++i)
	{
	    sprintf(torque_name,"%s-%d",fname,i);
	    rgb_data->torque_files[i] = fopen(torque_name,"w+");
	}

	set_default_data_file_names("rgb_data","-omega",&dname,&fname,
						init);
	if (create_directory(dname,NO) == FUNCTION_FAILED)
        {
            screen("ERROR in open_data_file(), directory %s doesn't exist "
                   		"and can't be made\n",dname);
            clean_up(ERROR);
        }
	for (i = 0; i < num_moving_body; ++i)
	{
	    sprintf(omega_name,"%s-%d",fname,i);
	    rgb_data->omega_files[i] = fopen(omega_name,"w+");
	}

	set_default_data_file_names("rgb_data","-force",&dname,&fname,
						init);
	if (create_directory(dname,NO) == FUNCTION_FAILED)
        {
            screen("ERROR in open_data_file(), directory %s doesn't exist "
                   		"and can't be made\n",dname);
            clean_up(ERROR);
        }
	for (i = 0; i < num_moving_body; ++i)
	{
	    sprintf(force_name,"%s-%d",fname,i);
	    rgb_data->force_files[i] = fopen(force_name,"w+");
	}

	set_default_data_file_names("rgb_data","-com",&dname,&fname,
						init);
	if (create_directory(dname,NO) == FUNCTION_FAILED)
        {
            screen("ERROR in open_data_file(), directory %s doesn't exist "
                   		"and can't be made\n",dname);
            clean_up(ERROR);
        }
	for (i = 0; i < num_moving_body; ++i)
	{
	    sprintf(com_name,"%s-%d",fname,i);
	    rgb_data->com_files[i] = fopen(com_name,"w+");
	}
}		/*end init_moving_body_data*/

/*
*			record_moving_body_data():
*
*	Main routine for computation of moving body data.
*/

/*ARGSUSED*/
LOCAL	void record_moving_body_data(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*out,
	boolean		about_to_stop)
{
	int i,j,num_moving_body;
	INTERFACE *intfc = front->interf;
	CURVE **c;
	RGB_DATA *rgb_data = RGB_data(out);
	FILE **torque_files = rgb_data->torque_files;
	FILE **omega_files = rgb_data->omega_files;
	FILE **force_files = rgb_data->force_files;
	FILE **com_files = rgb_data->com_files;
	static boolean first = YES;
	static double *torque,*omega;
	static double **force,**com_velo;
	int dim = intfc->dim;
	double t,f[MAXD];

	printf("Entering record_moving_body_data()\n");
	num_moving_body = -1;
	for (c = intfc->curves; c && *c; ++c)
        {
            if (wave_type(*c) == MOVABLE_BODY_BOUNDARY)
                if (body_index(*c) > num_moving_body)
                    num_moving_body = body_index(*c);
        }
        num_moving_body++;
        pp_global_imax(&num_moving_body,1);
	if (num_moving_body == 0) return;

	if (first)
	{
	    uni_array(&torque,num_moving_body,FLOAT);
	    uni_array(&omega,num_moving_body,FLOAT);
	    bi_array(&force,num_moving_body,MAXD,FLOAT);
	    bi_array(&com_velo,num_moving_body,MAXD,FLOAT);
	}
	for (i = 0; i < num_moving_body; ++i)
        {
            torque[i] = 0.0;
            for (j = 0; j < dim; ++j)
                force[i][j] = 0.0;
        }
	for (c = intfc->curves; c && *c; ++c)
        {
            if (wave_type(*c) == MOVABLE_BODY_BOUNDARY)
	    {
		i = body_index(*c);
		FrontForceAndTorqueOnHs(front,Hyper_surf(*c),front->dt,f,&t);
		torque[i] += t;
		omega[i] = angular_velo(*c);
		for (j = 0; j < dim; ++j)
		{
		    com_velo[i][j] = center_of_mass_velo(*c)[j];
		    force[i][j] += f[j];
		}
	    }
        }
	for (i = 0; i < num_moving_body; ++i)
	{
	    pp_global_sum(force[i],dim);
            pp_global_sum(&torque[i],1);
	}
	if (pp_mynode() != 0) return;

	if (first)
	{
	    first = NO;
	    for (i = 0; i < num_moving_body; ++i)
	    {
		fprintf(torque_files[i],"\"Torque of body %d\"\n",i+1);
		fprintf(force_files[i],"\"Total force on body %d\"\n",i+1);
		fprintf(omega_files[i],"\"Angular velocity of body %d\"\n",i+1);
		fprintf(com_files[i],"\"COM velocity of body %d\"\n",i+1);
	    }
	}
	for (i = 0; i < num_moving_body; ++i)
	{
	    fprintf(torque_files[i],"%f  %f\n",front->time,torque[i]);
	    fprintf(omega_files[i],"%f  %f\n",front->time,omega[i]);
	    fprintf(force_files[i],"%f  ",front->time);
	    fprintf(com_files[i],"%f  ",front->time);
	    for (j = 0; j < dim; ++j)
	    {
		fprintf(force_files[i],"%f  ",force[i][j]);
		fprintf(com_files[i],"%f  ",com_velo[i][j]);
	    }
	    fprintf(force_files[i],"\n");
	    fprintf(com_files[i],"\n");

	    fflush(torque_files[i]);
	    fflush(omega_files[i]);
	    fflush(force_files[i]);
	    fflush(com_files[i]);
	}
}		/*end record_moving_body_data*/
