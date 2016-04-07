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
*				gireadstate.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains restart initialization routines for gas dynamics.
*
*/

#include <ginit/ginit.h> 

	/* LOCAL Function Declarations */
LOCAL	boolean	new_read_print_avisc_field(const char*,const char*,
					   double*,const IO_TYPE*);
LOCAL	uint64_t read_print_Gas_param(INIT_DATA*,const IO_TYPE*,
                                      Gas_param**,int,size_t);
LOCAL	void	get_params_for_comp(INIT_DATA*,const IO_TYPE*,INTERFACE*,
                                    COMPONENT,uint64_t*,uint64_t*,Gas_param**,
				    int,int,size_t);
LOCAL	void	get_params_from_file_or_list(INIT_DATA*,const IO_TYPE*,
                                             Gas_param**,uint64_t,uint64_t*,
					     Gas_param**,int,int,size_t);
LOCAL	void	old_read_print_Dirichlet_bdry_states(INIT_DATA*,const IO_TYPE*,
                                                     INTERFACE*);
LOCAL	void	set_boundary_state_function(INIT_DATA*,const IO_TYPE*,
                                            const char*,INTERFACE*,
					    BOUNDARY_STATE*);
LOCAL	void  	read_print_time_dep_pre_data(const IO_TYPE*,BOUNDARY_STATE*);
LOCAL   void    read_print_time_dependent_boundary_state(INIT_DATA*,
							 const IO_TYPE*,
							 BOUNDARY_STATE*,
							 INTERFACE*);
#if defined(TWOD)
LOCAL	void	read_print_RP_DATA(INIT_DATA*,const IO_TYPE*,RP_DATA*,
                                   INTERFACE*);
#endif /* defined(TWOD) */
#if defined(COMBUSTION_CODE)
LOCAL	void	get_combustion_params_for_comp(INIT_DATA*,const IO_TYPE*,
                                               INTERFACE*,COMPONENT,uint64_t*,
					       uint64_t*,uint64_t*,
					       Gas_param**,int,int,size_t);
LOCAL	void	read_print_combustion_params(INIT_DATA*,const IO_TYPE*,
                                             Gas_param**);
#endif /* defined(COMBUSTION_CODE) */


/*
*			g_read_print_Gas_param_list():
*/

/*	Params arrays indexed by components */
LOCAL	uint64_t *params_as_read = NULL;
LOCAL	Gas_param **gas_param_list = NULL;
#if defined(COMBUSTION_CODE)
LOCAL	uint64_t *other_params_as_read = NULL;
#endif /* defined(COMBUSTION_CODE) */



EXPORT	void g_read_print_Gas_param_list(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	INTERFACE     *intfc)
{
	FILE      *file = io_type->file;
	int	  max_n_comps = max_num_comps();
	int	  dim = intfc->dim;
	size_t	  sizest = g_sizest();
	int	  nprms;
	uint64_t  *old_prms;
	COMPONENT comp;
	Gas_param **prms_lst;

	debug_print("restart","Entered g_read_print_Gas_param_list()\n");

#if defined(COMBUSTION_CODE)
	if (g_composition_type() != PURE_NON_REACTIVE)
	    uni_array(&other_params_as_read,max_n_comps,sizeof(uint64_t));
#endif /* defined(COMBUSTION_CODE) */
	if(params_as_read == NULL)
	{
	    uni_array(&params_as_read,max_n_comps,sizeof(uint64_t));
	    uni_array(&gas_param_list,max_n_comps,sizeof(Gas_param *));
	}
	if (min_component(intfc) < 0 || max_component(intfc) > max_n_comps-1) 
	{
	    screen("ERROR in g_read_print_Gas_param_list(), "
		   "comps out of range\n"
		   "min_component(intfc) = %d, max_component(intfc) = %d\n"
		   "max_n_comps = %d\n",
		   min_component(intfc),max_component(intfc),max_n_comps);
	    clean_up(ERROR);
	}
	nprms = read_Gas_param_list(init,io_type,&old_prms,
	                            &prms_lst,dim,sizest,NO);

	if (next_output_line_containing_string(file,"Components and Gas_param")
		== NULL)
	{
	    g_free_restart_params_list();
	    debug_print("restart","Left g_read_print_Gas_param_list()\n");
	    return;
	}

        if(dim == 3)
	{
	    /*#bjet2 */
            (void) fgetstring(file,"min_component = ");
            (void) fscanf(file,"%d", &(min_component(intfc)));

            (void) fgetstring(file,"max_component = ");
            (void) fscanf(file,"%d", &(max_component(intfc)));
            /*printf("#minmaxin  %d  %d\n", min_component(intfc), max_component(intfc)); */
	}

	for (comp = min_component(intfc); comp <= max_component(intfc); ++comp)
	{
	    comp_type(comp)->type = RESTART;
	    comp_type(comp)->_get_state = NULL;
	    comp_type(comp)->free_comp_type_extra = NULL;
	    (void) fgetstring(file,"Component = ");
#if defined(COMBUSTION_CODE)
	    if (g_composition_type() == PURE_NON_REACTIVE)
#endif /* defined(COMBUSTION_CODE) */
	    {
	    	get_params_for_comp(init,io_type,intfc,comp,params_as_read,
				    old_prms,prms_lst,nprms,dim,sizest);
	    }
#if defined(COMBUSTION_CODE)
	    else
	    {
	    	get_combustion_params_for_comp(init,io_type,intfc,comp,
				               params_as_read,
					       other_params_as_read,
				               old_prms,prms_lst,nprms,
					       dim,sizest);
	    }
#endif /* defined(COMBUSTION_CODE) */
	    if (debugging("restart"))
	    {
		(void) printf("Gas param read for comp %d\n",comp);
	    	print_Gas_param(comp_type(comp)->params);
	    }
	}

	for (comp = min_component(intfc); comp <= max_component(intfc); ++comp)
	    gas_param_list[comp] = comp_type(comp)->params;

	if (debugging("restart")) 
	{
	    (void) printf("Params ft_assigned to components upon restart\n");
	    for (comp=min_component(intfc); comp<=max_component(intfc); ++comp)
	    {
	    	(void) printf("Params for component %d\n",comp);
	    	print_Gas_param(comp_type(comp)->params);
	    }
	}
#if defined(COMBUSTION_CODE)
	if (other_params_as_read != NULL)
	    free(other_params_as_read);
	other_params_as_read = NULL;
#endif /* defined(COMBUSTION_CODE) */
	debug_print("restart","Left g_read_print_Gas_param_list()\n");
}		/*end g_read_print_Gas_param_list*/

LOCAL	void get_params_for_comp(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	INTERFACE     *intfc,
	COMPONENT     comp,
	uint64_t      *par,
	uint64_t      *old_prms,
	Gas_param     **prms_lst,
	int	      nprms,
	int	      dim,
	size_t	      sizest)
{
	FILE	  *file = io_type->file;
	COMPONENT comp1;

	debug_print("restart","Entered get_params_for_comp()\n");
	(void) fscanf(file,"%*d %*s %*s %llu",par+comp);
	for (comp1 = min_component(intfc); comp1 < comp; ++comp1)
	{
	    if (par[comp1] == par[comp]) 
	    {
	    	comp_type(comp)->params = comp_type(comp1)->params;
	    	get_params_from_file_or_list(init,io_type,NULL,par[comp],
	    				     old_prms,prms_lst,
	    				     nprms,dim,sizest);
	    	break;
	    }
	}
	if (comp1 == comp)
	{
	    get_params_from_file_or_list(init,io_type,&comp_type(comp)->params,
	    			         par[comp],old_prms,prms_lst,
	    			         nprms,dim,sizest);
	}
	if (par[comp] == 0L)
	    comp_type(comp)->params = NULL;
	debug_print("restart","Left get_params_for_comp()\n");
}		/*end get_params_for_comp*/

#if defined(COMBUSTION_CODE)
LOCAL	void get_combustion_params_for_comp(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	INTERFACE     *intfc,
	COMPONENT     comp,
	uint64_t      *par,
	uint64_t      *opar,
	uint64_t      *old_prms,
	Gas_param     **prms_lst,
	int	      nprms,
	int	      dim,
	size_t	      sizest)
{
	FILE      *file = io_type->file;
	COMPONENT comp1;

	(void) fscanf(file,"%*d %*s %*s %llu %*s %*s %llu",par+comp,opar+comp);
	for (comp1 = min_component(intfc); comp1 < comp; ++comp1)
	{
	    if (par[comp1] == par[comp]) 
	    {
	    	comp_type(comp)->params = comp_type(comp1)->params;
	    	get_params_from_file_or_list(init,io_type,NULL,par[comp],
		                             old_prms,prms_lst,nprms,
					     dim,sizest);
		break;
	    }
	}
	if (comp1 == comp)
	{
	    get_params_from_file_or_list(init,io_type,&comp_type(comp)->params,
			                 par[comp],old_prms,prms_lst,
					 nprms,dim,sizest);
	}
	for (comp1 = min_component(intfc); comp1 <= comp; ++comp1)
	{
	    if (par[comp1]==opar[comp] && comp_type(comp)->params)
	    {
	    	comp_type(comp)->params->other_params = 
				comp_type(comp1)->params;
		break;
	    }
	    if ((opar[comp1] == opar[comp]) && comp_type(comp)->params &&
		 comp_type(comp1)->params && (comp1 < comp))
	    {
		comp_type(comp)->params->other_params = 
				comp_type(comp1)->params->other_params;
		break;
	    }
	}
	if (comp1 == comp)
	{
	    get_params_from_file_or_list(init,io_type,
					 &comp_type(comp)->params->other_params,
			                 par[comp],old_prms,prms_lst,nprms,
					 dim,sizest);
	}
	else
	{
	    get_params_from_file_or_list(init,io_type,NULL,par[comp],
					 old_prms,prms_lst,nprms,dim,
					 sizest);
	}
	if (par[comp] == 0L)
	    comp_type(comp)->params = NULL;
	if (opar[comp] == 0L && comp_type(comp)->params)
	    comp_type(comp)->params->other_params = NULL;
}		/*end get_combustion_params_for_comp*/
#endif /* defined(COMBUSTION_CODE) */

LOCAL	void get_params_from_file_or_list(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	Gas_param     **params,
	uint64_t      params_as_read,
	uint64_t      *old_prms,
	Gas_param     **prms_lst,
	int	      nprms,
	int	      dim,
	size_t	      sizest)
{
	int		i;

	if (old_prms == NULL)
	{
	    (void) read_print_Gas_param(init,io_type,params,dim,sizest);
	    return;
	}
	if (params == NULL)
	    return;
	*params = NULL;
	if (params_as_read == 0)
	    return;
	for (i = 0; i < nprms; ++i)
	{
	    if (params_as_read == old_prms[i])
	    {
	    	*params = prms_lst[i];
	    	return;
	    }
	}
	screen("ERROR in get_params_from_file_or_list(), "
	       "params not in list\n");
	clean_up(ERROR);	
}		/*end get_params_from_file_or_list*/

EXPORT	int read_Gas_param_list(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	uint64_t      **pold_prms,
	Gas_param     ***pprms_lst,
	int	      dim,
	size_t	      sizest,
	boolean	      free_store)
{
	FILE             *file;
	static uint64_t  *old_prms = NULL;
	static Gas_param **prms_lst = NULL;
	static int	 nprms = 0;
	int		 i;
	static OUTPUT	 *oput = NULL;

	debug_print("restart","Entered read_Gas_param_list()\n");
	if (free_store == YES)
	{
	    if (old_prms != NULL)
	    {
	    	free(old_prms);
	    	old_prms = NULL;
	    }
	    if (pold_prms != NULL)
	    	*pold_prms = NULL;
	    if (prms_lst != NULL)
	    {
	    	free(prms_lst);
	    	prms_lst = NULL;
	    }
	    if (pprms_lst != NULL)
	    	*pprms_lst = NULL;
	    nprms = 0;
	    return nprms;
	}
	else if (prms_lst != NULL)
	{
	    *pold_prms = old_prms;
	    *pprms_lst = prms_lst;
	    return nprms;
	}

	file = io_type->file;
	oput = save_read_file_variables(file,oput);
	if (next_output_line_containing_string(file,
			"Equation of State Params List") == NULL)
	{
	    debug_print("restart","Left read_Gas_param_list()\n");
	    return 0;
	}
	(void) fgetstring(file,"Number of params = ");
	(void) fscanf(file,"%d",&nprms);
	uni_array(&old_prms,nprms,sizeof(uint64_t));
	uni_array(&prms_lst,nprms,sizeof(Gas_param *));

	for (i = 0; i < nprms; ++i)
	{
	    old_prms[i] = read_print_Gas_param(init,io_type,prms_lst+i,
	                                       dim,sizest);
	}

	*pold_prms = old_prms;
	*pprms_lst = prms_lst;
	reset_read_file_variables(oput);
	debug_print("restart","Left read_Gas_param_list()\n");
	return nprms;
}		/*end read_Gas_param_list*/

EXPORT void reset_artificial_viscosity_and_heat_conduction(
	INIT_DATA	*init)
{
	INTERFACE	*intfc = restart_intfc(init);
	char		s[Gets_BUF_SIZE];
	COMPONENT	comp = NO_COMP, *comp_list = NULL;
	Gas_param	*params = NULL, **params_list = NULL;
#if defined(COMBUSTION_CODE)
	Gas_param	**other_params_list = NULL;
#endif /* defined(COMBUSTION_CODE) */
	int		i, j, size;

	screen("Request change of artificial visocity and heat conduction\n");
	screen("\tfor individual gas param structures (y or no (dflt)): ");
	(void) Gets(s);
	if (s[0] != 'y' && s[0] != 'Y')
	    return;

	size = max_component(intfc) - min_component(intfc) + 1;
	uni_array(&params_list,size,sizeof(Gas_param *));
#if defined(COMBUSTION_CODE)
	if (g_composition_type() != PURE_NON_REACTIVE)
	    uni_array(&other_params_list,size,sizeof(Gas_param *));
#endif /* defined(COMBUSTION_CODE) */
	uni_array(&comp_list,size,sizeof(COMPONENT));
	for (i = 0, comp = min_component(intfc); comp <= max_component(intfc); 
		++i, ++comp)
	{
	    comp_list[i] = comp;
	    params_list[i] = gas_param_list[comp];
#if defined(COMBUSTION_CODE)
	    if (other_params_list != NULL)
	    	other_params_list[i] = params_list[i]->other_params;
#endif /* defined(COMBUSTION_CODE) */
	}

	for (i = 0; i < size; ++i)
	{
	    if (params_list[i] == NULL)
		continue;
	    params = params_list[i];

	    screen("Gas_param %d associated with components ",params);
	    for (j = i; j < size; ++j)
	    {
	    	if (params_list[j] != params)
		    continue;
	    	params_list[j] = NULL;
	    	screen("%d ",comp_list[j]);
	    }
	    screen("\n");
#if defined(COMBUSTION_CODE)
	    if (other_params_list != NULL)
	    {
	    	screen("Gas_param %d ",params);
	    	screen("also is other_params for components ");
	    	for (j = 0; j < size; ++j)
	    	{
	    	    if (other_params_list[j] != params)
			continue;
	    	    other_params_list[j] = NULL;
	    	    screen("%d ",comp_list[j]);
	    	}
	    	screen("\n");
	    }
#endif /* defined(COMBUSTION_CODE) */
	    print_Gas_param(params);
	    screen("Do you wish to reset the artificial viscosity and\n");
	    screen("\theat conduction for this params? [y or n (dflt)]: ");
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y')
	    {
	    	prompt_for_artificial_viscosity_and_heat_conduction(
	    		init,"special ","for this EOS model, ",
	    		YES,&params->avisc);
	    }
	}
#if defined(COMBUSTION_CODE)
	if (other_params_list != NULL)
	{
	    for (i = 0; i < size; ++i)
	    {
	        if (other_params_list[i] == NULL)
		    continue;
		params = other_params_list[i];
		screen("Gas_param %d ",params);
		screen("is other params for components ");
		for (j = i; j < size; ++j)
		{
		    if (other_params_list[j] != params)
			continue;
		    other_params_list[j] = NULL;
		    screen("%d ",comp_list[j]);
		}
		screen("\n");
		print_Gas_param(params);
		screen("Do you wish to reset the ");
		screen("artificial viscosity and\n\theat conduction");
		screen("for this params? [y or n (dflt)]: ");
		(void) Gets(s);
		if (s[0] == 'y' || s[0] == 'Y')
		{
		    prompt_for_artificial_viscosity_and_heat_conduction(
				init,"special ","for this EOS model, ",
				YES,&params->avisc);
		}
	    }
	}
#endif /* defined(COMBUSTION_CODE) */
	(void) printf("\n\n");
	free_these(2,params_list,comp_list);
#if defined(COMBUSTION_CODE)
	if (other_params_list != NULL)
	    free(other_params_list);
#endif /* defined(COMBUSTION_CODE) */
}		/*end reset_artificial_viscosity_and_heat_conduction*/

EXPORT	void g_free_restart_params_list(void)
{
#if defined(COMBUSTION_CODE)
	if (other_params_as_read != NULL)
	{
	    free(other_params_as_read);
	    other_params_as_read = NULL;
	}
#endif /* defined(COMBUSTION_CODE) */

	(void) read_Gas_param_list(NULL,NULL,NULL,NULL,0,0,YES);
	if (params_as_read != NULL)
	{
	    free(params_as_read);
	    params_as_read = NULL;
	}
	if (gas_param_list != NULL)
	{
	    free(gas_param_list);
	    gas_param_list = NULL;
	}
}		/*end g_free_restart_params_list*/

#if defined(TWOD)
EXPORT	void g_read_print_RP_DATA_at_nodes(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	INTERFACE     *intfc)
{
	FILE     *file = io_type->file;
	uint64_t RP_as_read;
	NODE	 **n;

	if (intfc->dim != 2)
	    return;

	if (next_output_line_containing_string(file,
		"RP_DATA information at nodes") == NULL)
	{
	    screen("ERROR in g_read_print_RP_DATA_at_nodes(), "
	           "Unable to find RP_DATA printout\n");
	    clean_up(ERROR);
	}

	for (n = intfc->nodes; n && *n; ++n)
	{
	    (void) fgetstring(file,"Node ");
	    (void) fscanf(file,"%*ld %*s %*s %*s %llu",&RP_as_read);
	    if (RP_as_read != 0)
	    	read_print_RP_DATA(init,io_type,Rp_data(*n),intfc);
	}
	(void) next_output_line_containing_string(file,
				"End of RP_DATA information at nodes");
}		/*end g_read_print_RP_DATA_at_nodes*/

LOCAL	void read_print_RP_DATA(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	RP_DATA	      *RP,
	INTERFACE *intfc)
{
	FILE      *file = io_type->file;
	size_t	  sizest = size_of_state(intfc);
	int	  dim = intfc->dim;
	int	  i, c;
	char	  s[80];
	int	  nprms;
	Gas_param **prms_lst;
	uint64_t  *old_prms;

	nprms = read_Gas_param_list(init,io_type,&old_prms,&prms_lst,
	                            dim,sizest,NO);

	if (!fgetstring(file,"Printout of RP_DATA structure "))
	{
	    screen("ERROR in read_print_RP_DATA(), "
		   "RP_DATA printout not found\n");
	    clean_up(ERROR);
	}

	if (!fgetstring(file,"ang_dir = "))
	{
	    screen("ERROR in read_print_RP_DATA(), "
		   "ang_dir printout not found\n");
	    clean_up(ERROR);
	}
	RP->ang_dir = fread_angle_direction(file);

	if (!fgetstring(file,"Angles:"))
	{
	    screen("ERROR in read_print_RP_DATA(), "
		   "angle printout not found\n");
	    clean_up(ERROR);
	}
	if ((c = getc(file)) == '\f') /* BINARY OUTPUT */
	{
	    (void) getc(file);
	    (void) read_binary_real_array(RP->ang,MAX_N_CURVES,io_type);
	}
	else
	{
	    (void) ungetc(c,file);
	    for (i = 0; i < MAX_N_CURVES; ++i)
	    {
	    	(void) sprintf(s,"ang[%d] = ",i);
	    	if (!fgetstring(file,s))
		{
	            screen("ERROR in read_print_RP_DATA(), "
		           "printout of ang[%d] not found\n",i);
	            clean_up(ERROR);
		}
	    	(void) fscan_float(file,&RP->ang[i]);
	    }
	}


	for (i = 0; i < MAX_N_CURVES; ++i)
	{
	    (void) sprintf(s,"state[%d]:",i);
	    if (!fgetstring(file,s))
	    {
	            screen("ERROR in read_print_RP_DATA(), "
		           "printout of state[%d] not found\n",i);
	            clean_up(ERROR);
	    }
	    (void) sprintf(s,"Mach number[%d] = ",i);
	    if (fgetstring(file,s))
	        RP->M[i] = fread_float(NULL,io_type);
	    read_print_gas_data(io_type,&RP->state[i],NO,RP->stype,sizest,
			        nprms,old_prms,prms_lst,dim);
	}
	if (fgetstring(file,"Turning angles:"))
	{
	    if ((c = getc(file)) == '\f') /* BINARY OUTPUT */
	    {
	        (void) getc(file);
	        (void) read_binary_real_array(RP->theta,MAX_N_CURVES,io_type);
	    }
	    else
	    {
	        (void) ungetc(c,file);
	        for (i = 0; i < MAX_N_CURVES; ++i)
	        {
	    	    (void) sprintf(s,"theta[%d] = ",i);
	    	    if (fgetstring(file,s))
		    {
	    	        (void) fscan_float(file,&RP->theta[i]);
		    }
	        }
	    }
	}
}		/*end read_print_RP_DATA*/

EXPORT	void g_read_print_ContactWallNodeParams(
	const IO_TYPE *io_type,
	INTERFACE     *intfc)
{
	FILE	  *file = io_type->file;
	CWNP *cwnp;
	if (intfc->dim != 2)
	    return;

	if (next_output_line_containing_string(file,
		"CONTACT WALL NODE PARAMS for interface") == NULL)
	    return;

	cwnp = contact_wall_node_params(intfc);

	if (fgetstring(file,"wall_bond_len = ") == FUNCTION_FAILED)
	{
	    screen("ERROR in g_read_print_ContactWallNodeParams(), "
		   "wall_bond_len printout not found\n");
	    clean_up(ERROR);
	}
	cwnp->wall_bond_len = fread_float(NULL,io_type);
	if (fgetstring(file,"first_adjust_time = ") == FUNCTION_FAILED)
	{
	    screen("ERROR in g_read_print_ContactWallNodeParams(), "
		   "first_adjust_time printout not found\n");
	    clean_up(ERROR);
	}
	cwnp->first_adjust_time = fread_float(NULL,io_type);
	if (fgetstring(file,"first_adjust_step = ") == FUNCTION_FAILED)
	{
	    screen("ERROR in g_read_print_ContactWallNodeParams(), "
		   "first_adjust_time printout not found\n");
	    clean_up(ERROR);
	}
	(void) fscanf(file,"%d",&cwnp->first_adjust_step);

	if (fgetstring(file,"adjust = ") == FUNCTION_FAILED)
	{
	    screen("ERROR in g_read_print_ContactWallNodeParams(), "
		   "first_adjust_time printout not found\n");
	    clean_up(ERROR);
	}
	cwnp->adjust = fread_boolean(file);
}		/*end g_read_print_ContactWallNodeParams*/
#endif /* defined(TWOD) */

EXPORT  void    g_read_print_boundary_state_data(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	INTERFACE     *intfc,
	int           index)
{
	BOUNDARY_STATE *bstate;
	const char     *s;

	f_read_print_boundary_state_data(init,io_type,intfc,index);
	bstate = bstate_list(intfc)[index];
	s = bstate->_boundary_state_function_name;
	if (s == NULL)
	    return;

	set_boundary_state_function(init,io_type,s,intfc,bstate);
}		/*end g_read_print_boundary_state_data*/

LOCAL	void	set_boundary_state_function(
	INIT_DATA      *init,
	const IO_TYPE  *io_type,
	const char     *s,
	INTERFACE      *intfc,
	BOUNDARY_STATE *bstate)
{
	if (strcmp(s,"flow_through_boundary_state") == 0)
	    bstate->_boundary_state_function = flow_through_boundary_state;
	else if (strcmp(s,"constant_pressure_flow_through_boundary_state") == 0)
	    bstate->_boundary_state_function =
		constant_pressure_flow_through_boundary_state;
	else if (strcmp(s,"time_dep_pressure_flow_through_boundary_state") == 0)
	{
	    /*#bjet2 */
            char c;
            
	    read_print_time_dep_pre_data(io_type,bstate);
            fgetstring(io_type->file,"Reference state:");
            c = getc(io_type->file);
            if ( c == '\n' )
            {
                FD_DATA  *fd_data;
                Locstate ref_st = NULL;
                
		fd_data = (FD_DATA*)bstate->_boundary_state_data;
                if (ref_st == NULL)
                    alloc_state(intfc,&ref_st,size_of_state(intfc));
                ref_st = read_print_state_data(init,io_type,ref_st,intfc);
                fd_data->state = ref_st;
                verbose_print_state("FD_DATA state",fd_data->state);
            }
            bstate->_fprint_boundary_state_data =
                g_fprint_tdp_boundary_state_data;
            bstate->_boundary_state_function =
                time_dep_pressure_flow_through_boundary_state;
	}
	else if (strcmp(s,"fixed_boundary_state") == 0)
	    bstate->_boundary_state_function = fixed_boundary_state;
	else if (strcmp(s,"g_fixed_boundary_state") == 0)
	    bstate->_boundary_state_function = g_fixed_boundary_state;
	else if (strcmp(s,"random_velocity_inlet") == 0)
	    read_print_random_velocity_inlet_data(init,io_type,bstate,intfc);
	else if (strcmp(s,"g_time_dependent_boundary_state") == 0)
	    read_print_time_dependent_boundary_state(init,io_type,bstate,intfc);
	else
	{
	    screen("ERROR in set_boundary_state_function(), "
		   "unknown boundary state function %s\n",s);
	    clean_up(ERROR);
	}
}		/*end set_boundary_state_function*/

LOCAL	void read_print_time_dep_pre_data(
	const IO_TYPE  *io_type,
	BOUNDARY_STATE *bstate)
{
	FD_DATA		*fd_data;
	char s[100];

	stat_scalar(&bstate->_boundary_state_data,sizeof(FD_DATA));
	fd_data = (FD_DATA*)bstate->_boundary_state_data;

	fd_data->tr  = read_print_float("Rise time = ",0.0,io_type);
	fd_data->tp  = read_print_float("Peak time = ",0.0,io_type);
	fd_data->ts  = read_print_float("Shut-off time = ",0.0,io_type);
	fd_data->pr_a = read_print_float("Ambient pressure = ",0.0,io_type);
	fd_data->pr_p = read_print_float("Peak pressure = ",0.0,io_type);
}	/* end read_print_time_dep_pre_data */


EXPORT	void g_read_print_Dirichlet_bdry_states(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	INTERFACE     *intfc)
{
	FILE       *file = io_type->file;
	HYPER_SURF **hs;
	int	   i;

	if (next_output_line_containing_string(file,
		"Hypersurface Dirichlet boundary state information") == NULL)
	{
	    old_read_print_Dirichlet_bdry_states(init,io_type,intfc);
	    return;
	}
	for (i = 0, hs = intfc->hss; hs && *hs; ++i, ++hs)
	{
	    if (wave_type(*hs) != DIRICHLET_BOUNDARY)
		continue;
	    (void) fgetstring(file,"Boundary state index for");
	    (void) fscanf(file,"%*s %*d %*s %d",&bstate_index(*hs));
	}
	(void) next_output_line_containing_string(file,
		"End hypersurface Dirichlet boundary state information");
}		/*end g_read_print_Dirichlet_bdry_states*/


LOCAL	void old_read_print_Dirichlet_bdry_states(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	INTERFACE     *intfc)
{
	FILE           *file = io_type->file;
	BOUNDARY_STATE Bstate, *bstate;
	HYPER_SURF     **hs;
	char	       s[120];
	char	       ss[256];
	int	       i, j, index;
	int	       dim = intfc->dim;
	size_t	       sizest = size_of_state(intfc);
	static const char     *hsname[] = { "point", "curve", "surface"};

	if (dim != 2)	/* Old style restart only supported for 2D */
	    return;

	if (next_output_line_containing_string(file,
		"Dirichlet boundary state information") == NULL)
	{
	    screen("ERROR in old_read_print_Dirichlet_bdry_states(), "
	           "Unable to find boundary state data\n");
	    clean_up(ERROR);
	}
	
	Bstate._fprint_boundary_state_data = f_fprint_boundary_state_data;
	Bstate._boundary_state_data = NULL;
	(void) sprintf(ss,"Boundary state information for %s",hsname[dim]);
	for (index = 6, hs = intfc->hss; hs && *hs; ++hs)
	{
	    if (wave_type(*hs) != DIRICHLET_BOUNDARY)
	    	continue;
	    (void) fgetstring(file,ss);
	    Bstate._boundary_state = read_print_state_data(init,io_type,NULL,intfc);

	    (void) fgetstring(file,"Boundary state function = ");
	    (void) fscanf(file,"%s",s);
	    if (strcmp(s,"NULL") == 0 || strcmp(s,"null") == 0 ||
	    	strcmp(s,"(NULL)") == 0 || strcmp(s,"(null)") == 0)
	    {
	    	Bstate._boundary_state_function = NULL;
	    	Bstate._boundary_state_function_name = NULL;
	    }
	    else
	    {
	    	Bstate._boundary_state_function_name = s;
	    	set_boundary_state_function(init,io_type,s,intfc,&Bstate);
	    }
	    for (i = 0; i < index; ++i)
	    {
	    	bstate = bstate_list(intfc)[i];
	    	if (memcmp((const void*)Bstate._boundary_state,
	    		   (const void*)bstate->_boundary_state,
	    		   sizest) != 0)
	    		continue;
	    	if (Bstate._boundary_state_function !=
	    			bstate->_boundary_state_function)
	    		continue;
	    	bstate_index(*hs) = i;
	    	break;
	    }
	    if (i == index)
	    {
	    	bstate_index(*hs) = index;
	    	(void) add_bstate_to_list(&Bstate,intfc,index++);
	    }
	}
	for (i = 0; i < dim; ++i)
	{
	    for (j = 0; j < 2; ++j)
	    {
		if (rect_boundary_type(intfc,i,j) != DIRICHLET_BOUNDARY)
		    continue;
		if (rect_bstate(intfc,i,j) != NULL)
		    continue;

		index = -1;
	    	for (hs = intfc->hss; hs && *hs; ++hs)
		{
		    RECT_GRID *tgr = &topological_grid(intfc);
		    int ic, jc;

		    if (wave_type(*hs) != DIRICHLET_BOUNDARY)
		    	continue;
		    (void) rect_bdry_side_for_curve(&ic,&jc,
				                    Curve_of_hs(*hs),tgr);
		    if (ic != i || jc != j) continue;
		    if (index == -1)
		    	index = bstate_index(*hs);
		    if (index != bstate_index(*hs))
		    {
		    	index = -1;
		    	break;
		    }
		}
		if (index == -1)
			continue;
		bstate_list(intfc)[2*i+j] = bstate_list(intfc)[index];
	    }
	}
	(void) next_output_line_containing_string(file,
		"End of Dirichlet boundary state information");
}		/*end old_read_print_Dirichlet_bdry_states*/


EXPORT	Locstate g_read_print_state_data(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	Locstate      state,
	INTERFACE     *intfc)
{
	int	  nprms;
	Gas_param **prms_lst;
	uint64_t  *old_prms;
	int	  dim = intfc->dim;
	size_t	  sizest = size_of_state(intfc);

	nprms = read_Gas_param_list(init,io_type,&old_prms,
				    &prms_lst,dim,sizest,NO);
	read_print_gas_data(io_type,&state,(state==NULL)?YES:NO,
			    GAS_STATE,sizest,nprms,old_prms,prms_lst,dim);
	return state;
}		/*end g_read_print_state_data*/

EXPORT	void read_print_gas_data(
	const IO_TYPE *io_type,
	Locstate      *pstate,
	int	      alloc,
	int	      state_type,
	size_t	      sizest,
	int	      num_eos,
	uint64_t      *par,
	Gas_param      **gpl,
	int	      dim)
{
	FILE		*file = io_type->file;
	Locstate	state;
	char		s[120];
	int		i, c;
	uint64_t	params_read;

	(void) fgetstring(file,"State information for the ");
	(void) fscanf(file,"%s",s);
	switch (s[0])
	{
	case 'N':		/* NULL state */
	    if (alloc)
		*pstate = NULL;
	    return;
	case 'O':		/* Obstacle state */
	    if (alloc)
		*pstate = (Locstate) store(sizest);
	    g_obstacle_state(*pstate,sizest);
	    return;
	default:
	    break;
	}
	if (alloc)
	    *pstate = (Locstate) store(sizest);
	state = *pstate;
	set_type_of_state(state,state_type);
	(void) fgetstring(file,"State Data ");
	if ((c = getc(file)) == '\f') /* BINARY */
	{
	    c = getc(file);
	    if (c == 1) /*Old style*/
	    {
	        (void) printf("WARNING in read_print_gas_data(), "
		              "old style printout can only be read\n"
                              "on machines with the same endian as the "
			      "output machine and can only be\n"
                              "restarted using the same floating "
			      "point precision as the output.\n");
	        (void) fread((void *)state,sizest,1,file);
		params_read = ptr2ull(Params(state));
	        set_restart_params(&Params(state),params_read,num_eos,par,gpl);
	    }
	    else
	    {
	        double *x;
		int   stype, failed;
	        x = &Dens(state);
	        (void) read_binary_real_array(x,2+SMAXD,io_type);
		read_binary_uint64_t_array(&params_read,1,io_type);
	        set_restart_params(&Params(state),params_read,num_eos,par,gpl);
		read_binary_int_array(&stype,1,io_type);
		set_type_of_state(state,stype);
		read_binary_int_array(&failed,1,io_type);
		set_material_failure(state,failed);
#if defined(COMBUSTION_CODE)
	        switch (Composition_type(state))
		{
		case ZND:
		case PTFLAME:
	            (void) read_binary_real_array(pdens(state),1,io_type);
		    break;
		case TWO_CONSTITUENT_REACTIVE:
	            (void) read_binary_real_array(pdens(state),2,io_type);
		    break;
		case PURE_NON_REACTIVE:
		default:
		    break;
		}
#endif /* defined(COMBUSTION_CODE) */
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    if(Params(state) != NULL &&
                       Params(state)->n_comps != 1)
                    {
                        int i;
                        for(i = 0; i < Params(state)->n_comps; i++)
                            (void) read_binary_real_array(&(pdens(state)[i]),1,io_type);
                    }
                }
	    }
	    return;
	}
	else
	{
	    double x, v[3];
#if defined(COMBUSTION_CODE)
	    double r;
#endif /* defined(COMBUSTION_CODE) */

	    (void) ungetc(c,file);
	    Dens(state) = fread_float("density = ",io_type);
	    (void) fscanf(file,"%*s%*s");
	    x = fread_float(NULL,io_type);
	    for (i = 0; i < dim; ++i)
	    {
	        (void) fscanf(file,"%*s%*s");
	    	v[i] = fread_float(NULL,io_type);
	    }
	    (void) fgetstring(file,"Gas_param = ");
	    (void) fscanf(file,"%llu",&params_read);
	    set_restart_params(&Params(state),params_read,num_eos,par,gpl);
	    reset_gamma(state);
#if defined(COMBUSTION_CODE)
	    if (Composition_type(state) != PURE_NON_REACTIVE)
	    {
	        (void) fgetstring(file,"burned = ");
	        (void) fscanf(file,"%s",s);
	        if (Composition_type(state) == PTFLAME)
	        {
	    	    if (s[0] == 'B')
	    	        Set_other_params(state,state);
	        }
	    }
#endif /* defined(COMBUSTION_CODE) */

	    switch (state_type)
	    {
	    case GAS_STATE:
		Energy(state) = x;
	        for (i = 0; i < dim; ++i)
		    Mom(state)[i] = v[i];
	    	reset_gamma(state);
#if defined(COMBUSTION_CODE)
		if ((Composition_type(state) == ZND) ||
		    (Composition_type(state) == TWO_CONSTITUENT_REACTIVE))
		{
	            Prod(state) = fread_float("product density = ",io_type);
	            if (Composition_type(state) == TWO_CONSTITUENT_REACTIVE)
		        Dens1(state) = fread_float(" rho1 = ",io_type);
		}
#endif /* defined(COMBUSTION_CODE) */
	        break;
	
	    case TGAS_STATE:
		Press(state) = x;
	        for (i = 0; i < dim; ++i)
		    Vel(state)[i] = v[i];
	    	reset_gamma(state);
#if defined(COMBUSTION_CODE)
		if ((Composition_type(state) == ZND) ||
		    (Composition_type(state) == TWO_CONSTITUENT_REACTIVE))
		{
	            React(state) = fread_float("reaction progress = ",io_type);
	            if (Composition_type(state) == TWO_CONSTITUENT_REACTIVE)
		        Dens1(state) = fread_float(" rho1 = ",io_type);
		}
#endif /* defined(COMBUSTION_CODE) */
	        break;

	case EGAS_STATE:
		Energy(state) = x;
	        for (i = 0; i < dim; ++i)
		    Vel(state)[i] = v[i];
	    	reset_gamma(state);
#if defined(COMBUSTION_CODE)
		if ((Composition_type(state) == ZND) ||
		    (Composition_type(state) == TWO_CONSTITUENT_REACTIVE))
		{
	            React(state) = fread_float("reaction progress = ",io_type);
	            if (Composition_type(state) == TWO_CONSTITUENT_REACTIVE)
		        Dens1(state) = fread_float(" rho1 = ",io_type);
		}
#endif /* defined(COMBUSTION_CODE) */
	        break;

	    default:
	        screen("Unknown state type in read_print_gas_data()\n");
	        clean_up(ERROR);
	    }
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                if(Params(state) != NULL &&
                   Params(state)->n_comps != 1)
                {
                    int i;
                    char sss[256];
                    for(i = 0; i < Params(state)->n_comps; i++)
                    {
                        sprintf(sss,"partial density[%2d] = ", i);
                        pdens(state)[i] = fread_float(sss,io_type);
                    }
                }
            }
	}
}		/*end read_print_gas_data*/



EXPORT	void set_restart_params(
	Gas_param	**params,
	uint64_t	params_read,
	int		num_eos,
	uint64_t	*par,
	Gas_param	**gpl)
{
	int		i;

	for (i = 0; i < num_eos; ++i)
	{
	    if (params_read == par[i])
	    {
	    	*params = gpl[i];
	    	return;
	    }
	}

	screen("ERROR in set_restart_params(), No such params %llu in list\n",
	       params_read);
	clean_up(ERROR);
}		/*end set_restart_params*/


EXPORT	void read_print_avisc_params(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	int           visc_print_style,
	AVISC         *avisc)
{
	FILE    *file = io_type->file;
	AVISC	Avisc;

	default_artificial_viscosity(&Avisc);
	switch (visc_print_style)
	{
	case 4:
	    if (fgetstring(file,"Artificial Viscosities and Heat Conductions"))
	    {
	        use_lapidus_artificial_viscosity(Avisc) =
	            new_read_print_avisc_field(Use_nlav,Coef_nlav,
	                		       &Avisc.lapidus_visc_coef,
					       io_type);

	        use_linear_artificial_viscosity(Avisc) =
	            new_read_print_avisc_field(Use_lav,Coef_lav,
	                		       &Avisc.linear_visc_coef,
					       io_type);

	        use_upwind_artificial_viscosity(Avisc) =
	            new_read_print_avisc_field(Use_uwav,Coef_uwav,
	                		       &Avisc.upwind_visc_coef,
					       io_type);

	        use_muscl_slope_flattening(Avisc) =
	            new_read_print_avisc_field(Use_msf,Coef_msf_ieta,
					       &Avisc.msf_ieta,io_type);

	        if (fgetstring(file,Coef_msf_ms))
	            Avisc.min_shock_jump = fread_float(NULL,io_type);

	        if (fgetstring(file,Coef_msf_msvj))
	            Avisc.min_sp_vol_jump = fread_float(NULL,io_type);

	        if (fgetstring(file,Coef_hc))
	            Avisc.heat_cond = fread_float(NULL,io_type);

	        if (fgetstring(file,Coef_char_speed_cutoff))
	            Avisc.char_speed_cutoff = fread_float(NULL,io_type);

	        if (fgetstring(file,Coef_dst))
		    Avisc.dynamic_st = fread_float(NULL,io_type);

	        if (fgetstring(file,Coef_sp))
		    Avisc.sp_coef = fread_float(NULL,io_type);
		else if (use_lapidus_artificial_viscosity(Avisc))
	            Avisc.sp_coef =
		        lapidus_stability_factor(Avisc.lapidus_visc_coef);

                if (fgetstring(file,Coef_contact_detector))
                    Avisc.contact_detector = fread_float(NULL,io_type);
	    }
	    break;

	case 0:
	case 1:
	case 2:
	case 3:
	default:
	    (void) printf("WARNING in read_print_avisc_params(), unknown or "
	                  "old print style, visc_print_style = %d, "
			  "not supported\n",visc_print_style);
	    prompt_for_artificial_viscosity_and_heat_conduction(
				init,"","for this EOS model, ",
				YES,&Avisc);
	    break;
	}
	if (avisc != NULL)
	{
	    *avisc = Avisc;
	}
}		/*end read_print_avisc_params*/


LOCAL	boolean new_read_print_avisc_field(
	const char    *s1,
	const char    *s2,
	double	      *coef,
	const IO_TYPE *io_type)
{
	FILE	*file = io_type->file;
	boolean	itmp = NO;

	if (fgetstring(file,s1) == FUNCTION_SUCCEEDED)
	{
	    itmp = fread_boolean(file);
	    if (fgetstring(file,s2) == FUNCTION_SUCCEEDED)
	        *coef = fread_float(NULL,io_type);
	}
	return itmp;
}		/*end new_read_print_avisc_field*/

/*
*			read_print_Gas_param():
*/

LOCAL uint64_t read_print_Gas_param(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	Gas_param     **params,
	int	      dim,
	size_t	      sizest)
{
	FILE             *file = io_type->file;
	F_USER_INTERFACE *fuh;
	uint64_t	 params_as_read = 0L;
	int		 c, visc_print_style;

	if (params != NULL)
	    *params = NULL;
	(void) fgetstring(file,"Gas_param = ");
	if (params == NULL)
	    return 0;
	(void) fscanf(file,"%llu",&params_as_read);
	if (debugging("restart_params")) 
	    (void) printf("params_as_read = %llu\n",params_as_read);
	if (!params_as_read)
	    return params_as_read;
	visc_print_style = 0;
	while ((c = getc(file)) == ' ')
	    ++visc_print_style;
	(void) ungetc(c,file);

	if ((*params=read_print_EOS_data(init,io_type,*params)) == NULL)
	    return params_as_read;

	if (params != NULL && *params != NULL)
	{
	    (*params)->dim = dim;
	    fuh = f_user_hook(dim);
	    (*params)->_alloc_state = fuh->_alloc_state;
	    (*params)->_alloc_intfc_state = fuh->_alloc_intfc_state;
	    (*params)->_clear_state = fuh->_clear_state;
	    (*params)->sizest = sizest;
	    read_print_avisc_params(init,io_type,visc_print_style,
	                            &(*params)->avisc);

#if defined(COMBUSTION_CODE)
	    read_print_combustion_params(init,io_type,params);
#endif /* defined(COMBUSTION_CODE) */
	    read_print_thermodynamic_restrictions(*params,io_type);
	}
	return params_as_read;
}		/*end read_print_Gas_param*/

#if defined(COMBUSTION_CODE)
/*ARGSUSED*/
LOCAL	void read_print_combustion_params(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	Gas_param     **params)
{
	FILE    *file = io_type->file;
	int	c;
	int	composition_type;
	double	critical_temperature;
	int	burned;
	double	q;
	double	rate_mult;
	double	activ_en;

	(void) fgetstring(file,"composition_type = ");
	(void) fscanf(file,"%d",&composition_type);
	if (params != NULL && *params != NULL)
	    (*params)->composition_type = composition_type;
	if (composition_type == PURE_NON_REACTIVE) return;
	critical_temperature = fread_float("critical_temperature = ",io_type);
	(void) fgetstring(file,"burned = ");
	(void) fscanf(file,"%d",&burned);
	q = fread_float("q = ",io_type);
	rate_mult = fread_float("rate_mult = ",io_type);
	activ_en = fread_float("activ_en = ",io_type);
	if (params != NULL && *params != NULL)
	{
	    (*params)->critical_temperature = critical_temperature;
	    (*params)->burned = burned;
	    (*params)->q = q;
	    (*params)->rate_mult = rate_mult;
	    (*params)->activ_en = activ_en;
	    read_print_flame_velocity_params(init,io_type,*params);
	}
}		/*end read_print_combustion_params*/

#endif /* defined(COMBUSTION_CODE) */

EXPORT	void	read_print_thermodynamic_restrictions(
	Gas_param     *params,
	const IO_TYPE *io_type)
{
	if (fgetstring(io_type->file,"min_energy = ") == FUNCTION_FAILED)
	    return; /*Old style output file, fields not printed*/
	params->min_energy =
	    read_print_float(NULL,params->min_energy,io_type);
	params->min_pressure =
	    read_print_float("min_pressure = ",params->min_energy,io_type);
	params->vacuum_dens =
	    read_print_float("vacuum_dens = ",params->vacuum_dens,io_type);
	params->raref_press =
	    read_print_float("raref_press = ",params->raref_press,io_type);
#if defined(COMBUSTION_CODE)
	if (params->composition_type == PURE_NON_REACTIVE)
		return;
	params->tol_alpha =
	    read_print_float("tol_alpha = ",params->tol_alpha,io_type);
	params->tol_press =
	    read_print_float("tol_press = ",params->tol_press,io_type);
#endif /* defined(COMBUSTION_CODE) */
}		/*end read_print_thermodynamic_restrictions*/

LOCAL void read_print_time_dependent_boundary_state(
        INIT_DATA      *init,
        const IO_TYPE  *io_type,
        BOUNDARY_STATE *bstate,
        INTERFACE      *intfc)
{
    uint64_t   params_read;
    Gas_param  **prms_lst;
    uint64_t   *old_prms;
    int        nprms;
    int        dim = intfc->dim;
    size_t     sizest = size_of_state(intfc);
    TD_BSTATE  *tds;

    stat_scalar(&bstate->_boundary_state_data,sizeof(TD_BSTATE));
    tds = (TD_BSTATE*)bstate->_boundary_state_data;
    (void) fgetstring(io_type->file,"Time Dependent Boundary State Data");
    (void) fgetstring(io_type->file,"Data file = ");
    (void) fscanf(io_type->file,"%s",tds->fn);
    (void) fgetstring(io_type->file,"Params = ");
    (void) fscanf(io_type->file,"%llu",&params_read);
    nprms = read_Gas_param_list(init,io_type,&old_prms,&prms_lst,
	                            dim,sizest,NO);
    set_restart_params(&tds->params,params_read,nprms,old_prms,prms_lst);
    read_time_dependent_data(tds);

    bstate->_fprint_boundary_state_data = 
	    g_print_time_dependent_boundary_state;
    bstate->_boundary_state_function =
	g_time_dependent_boundary_state;
}/* read_print_time_dependent_boundary_state*/
