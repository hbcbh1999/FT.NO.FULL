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
*				hinit.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains routines used for initializing a Wave structure.
*/


#include <hyp/hdecs.h>

	/* LOCAL function prototypes */
LOCAL	void	set_min_buff_width(int,PP_GRID*);

EXPORT	void	set_hyp_hooks(
	INIT_DATA	*init)
{
	set_front_hooks(init);
	use_mv_states(init) = NO;
	h_init_data(init)->_init_wave_mv_state_list = NULL;
	h_init_data(init)->_alloc_wave_mv_state_list = NULL;
	h_init_data(init)->_init_mv_list_states = NULL;
	use_blocking_pp_comp(init) = YES;
}		/*end set_hyp_hooks*/

/*
*			prompt_for_hyp_method():
*
*	Prompts for the choice of hyperbolic method and sets *method
*	appropriately.
*/

EXPORT void prompt_for_hyp_method(
	INIT_DATA	*init)
{
	Hyp_method	*methods = available_hyperbolic_methods(init);
	Hyp_method	*default_method = default_hyperbolic_method(init);
	Hyp_method	*meth, *method;
	char		s[Gets_BUF_SIZE];
	DEBUG_ENTER(prompt_for_hyp_method)

	screen("\n\t\tSpecify interior hyperbolic difference method.\n"
	       "\nSelect the difference method for solving the hyperbolic\n"
	       "system in the interior regions away from the fronts.\n"
	       "Choices are\n");
	for (meth = methods; meth->ptype.prompt != NULL; ++meth)
	    screen("\t%s (%s)\n",meth->ptype.prompt,meth->ptype.select);
	screen("Enter choice here");
	if (default_method != NULL)
	    screen("(default = %s)",default_method->ptype.select);
	screen(": ");
	(void) Gets(s);

	method = default_method;
	if (s[0] != '\0')
	{
	    for (meth = methods; meth->ptype.prompt != NULL; ++meth)
	    {
	    	if (strncasecmp(s,meth->ptype.select,meth->ptype.ncmp) == 0)
	    	{
	    	    method = meth;
	    	    break;
	    	}
	    }
	}

	if (method == NULL)
	{
	    screen("ERROR in prompt_for_hyp_method(), no method selected\n");
	    clean_up(ERROR);
	}
	hyperbolic_method(init) = method;
	prompt_for_hyp_method_options(init);
	DEBUG_LEAVE(prompt_for_hyp_method)
}		/*end prompt_for_hyp_method*/

EXPORT	void	prompt_for_unsplit_options(
	INIT_DATA       *init)
{
	UnsplitStencilOptions *usopts = &USopts(init);
	char	s[Gets_BUF_SIZE];
	usopts->use_hyp_solution = NO;
	screen("Use hyp_solution for irregular stencil states (dflt = %s): ",
	       y_or_n(usopts->use_hyp_solution));
	(void) Gets(s);
	if ((s[0] == 'y') || (s[0] == 'Y'))
	    usopts->use_hyp_solution = YES;
	else if ((s[0] == 'n') || (s[0] == 'N'))
	    usopts->use_hyp_solution = NO;
}		/*end prompt_for_unsplit_options*/

/*
*			set_hyperbolic_method():
*
*	Sets the hyperbolic method parameters in the wave structure.
*/

EXPORT void set_hyperbolic_method(
	INIT_DATA	*init,
	Wave		*wave,
	int		(**hyp)(double,double*,Wave*,Front*))
{
	Hyp_method	*method = hyperbolic_method(init);
	DEBUG_ENTER(set_hyperbolic_method)

	if (method == NULL)
	{
	    screen("ERROR in prompt_for_hyp_method(), no method selected\n");
	    clean_up(ERROR);
	}
	wave->method = method->ptype.type.ctype;
	wave->npts_sten = method->npts_sten;
	wave->npts_vsten = method->npts_vsten;
	wave->sten_rad = method->sten_rad;
	*hyp = method->hyp_driver;

	set_min_buff_width(max(wave->sten_rad+1,3),wave->pp_grid);
	DEBUG_LEAVE(set_hyperbolic_method)
	return;
}		/*end set_hyperbolic_method*/


EXPORT	boolean	h_read_print_max_wave_speed_info(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	Wave          *wave)
{
	FILE    *file = io_type->file;
	int	i, j, dim = wave->rect_grid->dim;
	int	c;
	char	s[80];

	if (!fgetstring(file,"Maximum Wave Speed Information"))
		return FUNCTION_FAILED;
	(void) fgetstring(file,"Maxsp = ");
	if ((c = getc(file)) == '\f') /*BINARY*/
	{
	    (void) getc(file);
	    (void) read_binary_real_array(Maxsp(wave),dim,io_type);
	}
	else
	{
	    (void) ungetc(c,file);
	    (void) fscan_float(file,Maxsp(wave));
	    for (i = 1; i < dim; ++i)
	    {
	    	(void) getc(file);/*Grab comma*/
	    	(void) fscan_float(file,Maxsp(wave)+i);
	    }
	}
	for (i = 0; i < dim; ++i)
	{
	    (void) sprintf(s,"MaxWaveSpeedCoords[%d] = ",i);
	    if (fgetstring(file,s) == FUNCTION_SUCCEEDED)
	    {
	        if ((c = getc(file)) == '\f') /*BINARY*/
	        {
	    	    (void) getc(file);
		    (void) read_binary_real_array(MaxWaveSpeedCoords(wave)[i],
		                           dim,io_type);
	        }
	        else
	        {
	    	    (void) ungetc(c,file);
	    	    (void) fscan_float(file,MaxWaveSpeedCoords(wave)[i]);
		    for (j = 1; j < dim; ++j)
		    {
		        (void) getc(file);/*Grab comma*/
		        (void) fscan_float(file,MaxWaveSpeedCoords(wave)[i]+j);
		    }
	        }
	    }
	    else
	    {
		(void) printf("WARNING in h_read_print_max_wave_speed_info(), "
			      "failed to find printout of "
			      "MaxWaveSpeedCoords[%d]\n",i);
	    }
	}
	if (wave->sizest != 0)
	{
	    for (i = 0; i < dim; ++i)
	    {
		(void) sprintf(s,"MaxWaveSpeedState[%d] = ",i);
		(void) fgetstring(file,s);
		(void) read_print_state_data(init,io_type,
		                             MaxWaveSpeedState(wave)[i],
					     wave_tri_soln(wave)->intfc);
	    }
	}
	return FUNCTION_SUCCEEDED;
}		/*end h_read_print_max_wave_speed_info*/




/*
*		init_point_source_composition():
*	
*	Initializes the compositions of the point sources in a wave
*	structure, i.e. a Locstate's associated with inflow and outflow.
*	(*pt_source_initializer)(wave,k,&state) is called and the state
*	is associated with source k.
*/

EXPORT void init_point_source_composition(
	Wave		*wave,
	void		(*pt_source_initializer)(Wave*,int,Locstate*))
{
	int		k;
	Locstate	state;		/* state for point source */
	byte		*compstore;	/* storage address for
					   wave->composition */

	debug_print("init","Entering init_point_source_composition()\n");

	if (wave->num_point_sources == 0 || pt_source_initializer == NULL)
	{
		wave->composition = NULL;
		debug_print("init","Leaving init_point_source_composition()\n");
		return;
	}

	uni_array(&wave->composition,wave->num_point_sources,sizeof(Locstate));
	uni_array(&compstore,wave->num_point_sources,wave->sizest);
	for (k = 0; k < wave->num_point_sources; ++k)
	{
		(*pt_source_initializer)(wave,k,&state);
		wave->composition[k] = (Locstate) (compstore + k*wave->sizest);
		ft_assign(wave->composition[k],state,wave->sizest);
	}

	if (debugging("show_pt_sources"))
		print_pt_sources(wave);

	if (debugging("show_pt_sources") || debugging("init"))
		(void) printf("Leaving init_point_source_composition()\n");
}		/*end init_point_source_composition*/

LOCAL void set_min_buff_width(
	int		minbuf,
	PP_GRID		*pp_grid)
{
	RECT_GRID	*zoom_gr;
	int		i, dim;
	int		lbuf[MAXD],ubuf[MAXD];

	zoom_gr = &pp_grid->Zoom_grid;
	dim = pp_grid->Global_grid.dim;

	for (i = 0; i < dim; ++i)
	{
	    if (pp_grid->buf[i] < minbuf)
	    	pp_grid->buf[i] = minbuf;
	    ubuf[i] = (zoom_gr->ubuf[i] > 0) ? pp_grid->buf[i] : 0;
	    lbuf[i] = (zoom_gr->lbuf[i] > 0) ? pp_grid->buf[i] : 0;
	}
	set_rect_grid(zoom_gr->L,zoom_gr->U,zoom_gr->GL,zoom_gr->GU,
		      lbuf,ubuf,zoom_gr->gmax,dim,&zoom_gr->Remap,zoom_gr);
}		/*end set_min_buff_width*/

EXPORT	void h_set_interface_hooks(
	int		dim,
	INIT_DATA       *init)
{
	t_set_interface_hooks(dim,init);
}		/*end h_set_interface_hooks*/

/* For the parabolic step */
EXPORT	boolean	h_read_print_max_viscosity_info(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	Wave          *wave)
{
	FILE    *file = io_type->file;
	int	i, j, dim = wave->rect_grid->dim;
	int	c;
	char	s[80];


	if (!fgetstring(file,"Maximum Viscosity Information"))
		return FUNCTION_FAILED;
	(void) fgetstring(file,"Maxvisc = ");
	if ((c = getc(file)) == '\f') /*BINARY*/
	{
	    (void) getc(file);
	    /*(void) read_binary_real_array(Maxvisc(wave),dim,io_type); */
	}
	else
	{
	    (void) ungetc(c,file);
	    /*(void) fscan_float(file,Maxvisc(wave)); */
	    
	    	(void) getc(file);/*Grab comma*/
	    	/*(void) fscan_float(file,Maxvisc(wave)); */
	}
	for (i = 0; i < dim; ++i)
	{
	    (void) sprintf(s,"MaxViscosityCoords[%d] = ",i);
	    if (fgetstring(file,s) == FUNCTION_SUCCEEDED)
	    {
	        if ((c = getc(file)) == '\f') /*BINARY*/
	        {
	    	    (void) getc(file);
		    /*(void) read_binary_real_array(MaxViscosityCoords(wave)[i],
		      dim,io_type);*/
	        }
	        else
	        {
	    	    (void) ungetc(c,file);
	    	    /*(void) fscan_float(file,MaxViscosityCoords(wave)[i]); */
		    (void) getc(file);/*Grab comma*/
		    /*(void) fscan_float(file,MaxViscosityCoords(wave)[i]); */
		    
	        }
	    }
	    else
	    {
		(void) printf("WARNING in h_read_print_max_viscosity_info(), "
			      "failed to find printout of "
			      "MaxViscosityCoords[%d]\n",i);
	    }
	}
	if (wave->sizest != 0)
	{
	    for (i = 0; i < dim; ++i)
	    {
		(void) sprintf(s,"MaxViscosityState[%d] = ",i);
		(void) fgetstring(file,s);
		(void) read_print_state_data(init,io_type,
		                             MaxViscosityState(wave)[i],
					     wave_tri_soln(wave)->intfc);
	    }
	}
	return FUNCTION_SUCCEEDED;
}		/*end h_read_print_max_viscosity_info*/
