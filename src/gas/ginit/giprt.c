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
*				giprt.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include <ginit/ginit.h>
#include <sys/types.h>
#include <time.h>
#include <sys/stat.h>


	/* LOCAL Function Declarations */
LOCAL	PRINTING_LIST	*g_set_optional_printing_variables(PRINTING_LIST*,
							   PRINTING_LIST*,
							   INIT_PHYSICS*,
							   INIT_DATA*,
							   Printplot*);
LOCAL	Plot_choice	*choices_valid_for_dim(Scalar_Plot_choice*,int);
LOCAL	double	**velocity_gradient(double*,Front*,Wave*,COMPONENT);
LOCAL	double	return_Mach_number(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	return_density(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	return_divergence_velocity(double*,Front*,POINTER,
					   COMPONENT,Locstate);
LOCAL	double	return_energy_density(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	return_enthalpy(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	return_entropy(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	return_flow_speed(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	return_kinetic_energy(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	return_momentum(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	return_pressure(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	return_radial_component_of_velocity(double*,Front*,POINTER,
						    COMPONENT,Locstate);
LOCAL	double	return_specific_internal_energy(double*,Front*,POINTER,
						COMPONENT,Locstate);
LOCAL	double	return_temperature(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	return_vorticity(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	return_x_velocity(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	int	g_SelectedScalarItems(const char*,Plot_choice*,
				      ScalarPlotItem**);
LOCAL	void	fgraph_front_states_shell(Grid*,Wave*,Front*,Printplot*,
					  OUTPUT_DATA*,boolean);
LOCAL	void	init_statistics_for_problem(Front*,Grid*,Printplot*,
					    INIT_DATA*);

LOCAL	char	*set_default_output_file_name(const char*,const char*,char*);
LOCAL	void	init_Lp_Diff(INIT_DATA*,Front*,Grid*,Printplot*,int,
			     void(*)(Locstate,double*,Front*,Wave*,int));

LOCAL	double	return_tangential_component_of_xy_velocity(double*,Front*,
				POINTER,COMPONENT,Locstate);
LOCAL	double	return_ss_Mach_number(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	void	set_rm_graph_scale(OUTPUT_DATA**,int);

LOCAL	double	return_y_velocity(double*,Front*,POINTER,COMPONENT,Locstate);

LOCAL	double	return_z_velocity(double*,Front*,POINTER,COMPONENT,Locstate);

#if defined(COMBUSTION_CODE)
LOCAL	double	return_density1(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	return_product_density(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	return_reaction_progress(double*,Front*,POINTER,
					 COMPONENT,Locstate);
#endif /* defined(COMBUSTION_CODE) */


/*
*			g_prompt_for_printing_and_plotting():
*
*	Prompts for physics dependent plotting options.
*/

EXPORT	void	g_prompt_for_printing_and_plotting(
	INIT_DATA	*init)
{
	(void)  set_plot_choices(init);
	if (Prt_mode(prt_opts(init)[TRI_STATES]) != PRINTING_OFF)
	    prompt_for_tri_plots(init);
#if defined(USE_HDF)
	if (Prt_mode(prt_opts(init)[HDF_STATES]) != PRINTING_OFF)
	    prompt_for_hdf_plots(HDF_RASTER,init);
	if (Prt_mode(prt_opts(init)[SDS_STATES]) != PRINTING_OFF)
	    prompt_for_hdf_plots(HDF_SDS,init);
#endif /* defined(USE_HDF) */
	/* needed for VTK */
	if (Prt_mode(prt_opts(init)[VTK_STATES]) != PRINTING_OFF)
	    prompt_for_vtk_plots(init);
#if defined(__GD__)
	if (Prt_mode(prt_opts(init)[GD_MOVIE]) != PRINTING_OFF)
	    prompt_for_gd_plots(init);
#endif /* defined(__GD__) */
	/* end needed for VTK */


	g_prompt_for_optional_printing_variables(init);
}		/*end g_prompt_for_printing_and_plotting*/

/*
*			g_init_printing_and_plotting():
*
*	Initializes the printing (rect grid and	interface) and plotting
*	(tri grid) of state data.
*/

EXPORT void g_init_printing_and_plotting(
	INIT_PHYSICS	*ip,
	INIT_DATA	*init,
	Printplot	*prt,
	int		nfloats)
{
	Front		*fr =ip->root->front;
	RECT_GRID	*gr = fr->rect_grid;
	PRINTING_LIST	*head, *pl;
	int		var;

	if (output_format_on(TRI_STATES,prt))
	    init_tri_plots(init,prt);

#if defined(USE_HDF)
	if (output_format_on(HDF_STATES,prt))
	    init_HDF_plots(HDF_RASTER,init,prt);
	if (output_format_on(SDS_STATES,prt))
	    init_HDF_plots(HDF_SDS,init,prt);
#endif /* defined(USE_HDF) */

	/* needed for VTK */
	if (output_format_on(VTK_STATES,prt))
	    init_vtk_plots(init,prt);
	/* end needed for VTK */

#if defined(__GD__)
	if (output_format_on(GD_MOVIE,prt))
	    init_gd_movie_plots(init,prt);
#endif /* defined(__GD__) */

	head = set_printing_list(ip,init,prt);
	bi_array(&prt->output_soln,prt->n_rect_state_vars,1,sizeof(OUTPUT_SOLN));

	for (pl = head; pl != NULL; pl = pl->next)
	    prt->output_soln[pl->var]->name = pl->name;

	for (var = 0; var < prt->n_rect_state_vars; ++var)
	{
	    prt->output_soln[var]->fit = LINEAR_FIT;
	    prt->output_soln[var]->smoothness = SINGULAR;
	    prt->output_soln[var]->solution = g_iphys(ip)->_phys_solution;
	    prt->output_soln[var]->grid = gr;
	    prt->output_soln[var]->intfc_solution =
	        g_iphys(ip)->_phys_intfc_solution;
	    prt->output_soln[var]->intfc = NULL;
	    prt->output_soln[var]->var = var;
	    prt->output_soln[var]->extra = NULL;
	}

	prt->n_restart_vars = nfloats;
	prt->restart_soln = g_iphys(ip)->_set_input_solution(prt);
	prt->_print_extreme_values = g_print_extreme_values;
}		/*end g_init_printing_and_plotting*/

EXPORT	PRINTING_LIST	*g_set_printing_list(
	INIT_PHYSICS	*ip,
	INIT_DATA	*init,
	Printplot	*prt)
{
	int	   i, var;
	int	   dim = ip->root->front->rect_grid->dim;
	static const char *momentum[] = { "X-MOMENTUM",
					  "Y-MOMENTUM",
					  "Z-MOMENTUM"};
	PRINTING_LIST	*head, *pl;

	var = 0;
	scalar(&head,sizeof(PRINTING_LIST));
	pl = head;

	pl->name = "DENSITY";
	pl->var  = var++;
	pl->next = NULL;

	scalar(&pl->next,sizeof(PRINTING_LIST));
	pl = pl->next;
	pl->name = "ENERGY";
	pl->var  = var++;
	pl->next = NULL;

	for (i = 0; i < dim; ++i)
	{
	    scalar(&pl->next,sizeof(PRINTING_LIST));
	    pl = pl->next;
	    pl->name = momentum[i];
	    pl->var  = var++;
	    pl->next = NULL;
	}

#if defined(COMBUSTION_CODE)
	{
	    int	composition_type = g_composition_type();
	    if (composition_type == PTFLAME)
	    {
	    	scalar(&pl->next,sizeof(PRINTING_LIST));
	    	pl = pl->next;
	    	pl->name = "REACTION_PROGRESS";
	    	pl->var  = var++;
	    	pl->next = NULL;
	    }
	    if ((composition_type == ZND) ||
	        (composition_type == TWO_CONSTITUENT_REACTIVE))
	    {
	    	scalar(&pl->next,sizeof(PRINTING_LIST));
	    	pl = pl->next;
	    	pl->name = "PRODUCT_DENSITY";
	    	pl->var  = var++;
	    	pl->next = NULL;
	    }
	    if (composition_type == TWO_CONSTITUENT_REACTIVE)
	    {
	    	scalar(&pl->next,sizeof(PRINTING_LIST));
	    	pl = pl->next;
	    	pl->name = "DENSITY_1";
	    	pl->var  = var++;
	    	pl->next = NULL;
	    }
	}
#endif /* defined(COMBUSTION_CODE) */
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            int num_comp;
            char  *prtstring;

            /*** According to g_compute_sizest() ***/
            num_comp = NumberFloats(init) - (dim + 2);
            for(i = 0; i < num_comp; i++)
            {
                scalar(&pl->next,sizeof(PRINTING_LIST));
                uni_array(&prtstring,20,CHAR);
                sprintf(prtstring,"P_DENS_%d",i);
                pl = pl->next;
                pl->name = prtstring;
                pl->var  = var++;
                pl->next = NULL;
            }
        }
	prt->n_rect_state_vars = var;

	return g_set_optional_printing_variables(head,pl,ip,init,prt);
}		/*end g_set_printing_list*/

EXPORT	void	g_prompt_for_optional_printing_variables(
	INIT_DATA	*init)
{
	static Prompt_type g_opvs[] =
	{
	    {
	    	"VELOCITY",
	    	"V",
	    	1,
	    	{ERROR}
	    },
	    {
	    	"RADIAL_COMPONENT_OF_VELOCITY",
	    	"RV",
	    	2,
	    	{ERROR}
	    },
	    {
	    	"PRESSURE",
	    	"P",
	    	1,
	    	{ERROR}
	    },
	    {
	    	"SPECIFIC_ENTROPY",
	    	"SPECIFIC_ENTROPY",
	    	16,
	    	{ERROR}
	    },
	    {
	    	"SOUND_SPEED",
	    	"SOUND_SPEED",
	    	11,
	    	{ERROR}
	    },
	    {
	    	"TEMPERATURE",
	    	"TEMPERATURE",
	    	11,
	    	{ERROR}
	    },
	    {
	    	NULL,
	    	NULL,
	    	0,
	    	{ERROR}
	    },
	    {			/*Extra storage for tangential velocity*/
	    	NULL,
	    	NULL,
	    	0,
	    	{ERROR}
	    },
	    {			/*Extra storage for debugging params*/
	    	NULL,
	    	NULL,
	    	0,
	    	{ERROR}
	    }
	};
	int i;

	if (Prt_mode(prt_opts(init)[RECT_STATES]) == PRINTING_OFF)
	{
	    optional_printing_variables(init) = NULL;
	    optional_printing_prompt_types(init) = NULL;
	    return;
	}

	if (i_intfc(init)->dim == 2)
	{
	    static Prompt_type Ptype = {
	    	    			 "TANGENTIAL_COMPONENT_OF_XY_VELOCITY",
	    	    			 "TV",
	    	    			 2,
	    	    			 {ERROR}
	        			};
	    for (i = 0; g_opvs[i].prompt != NULL; ++i);
	    g_opvs[i] = Ptype;
	}

	optional_printing_prompt_types(init) = g_opvs;
	
	screen("\nSpecify additional GAS printing variables.  ");
	screen("The choices are,\n");
	for (i = 0; g_opvs[i].prompt != NULL; ++i)
	    screen("\t\t%s (%s)%s",g_opvs[i].prompt,
				   g_opvs[i].select,
			           (g_opvs[i+1].prompt==NULL) ? "." : ",\n");
	optional_printing_variables(init) = read_plotting_choices(NULL,init);

	if (debugging("prt-params"))
	{
	    char   *s;
	    size_t len;
	    static Prompt_type PrtParams = {
	    	    			     "EOS-PARAMS",
	    	    			     "EOS-PARAMS",
	    	    			     10,
	    	    			     {ERROR}
	        			   };
	    for (i = 0; g_opvs[i].prompt != NULL; ++i);
	    g_opvs[i] = PrtParams;
	    len = strlen(optional_printing_variables(init));
	    len += 2 + strlen(g_opvs[i].select);
	    s = (char*) init_table_Store(len*CHAR,init);
	    (void) sprintf(s,"%s %s",optional_printing_variables(init),
				     g_opvs[i].select);
	    optional_printing_variables(init) = s;
	}
}		/*end g_prompt_for_optional_printing_variables*/

LOCAL	PRINTING_LIST	*g_set_optional_printing_variables(
	PRINTING_LIST	*head,
	PRINTING_LIST	*pl,
	INIT_PHYSICS	*ip,
	INIT_DATA	*init,
	Printplot	*prt)
{
	static const char *velocity[] = { "X-VELOCITY",
					  "Y-VELOCITY",
					  "Z-VELOCITY"};
	Prompt_type *g_opvs = optional_printing_prompt_types(init);
	PRINTING_LIST	dummy_head;
	char	*s;
	int	i, j, k, len, var;
	int	dim = ip->root->front->rect_grid->dim;

	if (!output_format_on(RECT_STATES,prt))
	{
	    return head;
	}

	var = (pl != NULL) ? pl->var+1 : 0;

	/* Init printing of rect grid and interface states */

	s = optional_printing_variables(init);
	if ((s == NULL) || (s[0] == '\0'))
	{
	    return head;
	}
	len = (int) strlen(s);

	if (pl == NULL)
	    pl = &dummy_head;
	for (i = 0; i < len; ++i)
	{
	    if (isspace(s[i]))
	        continue;
	    for (j = 0; g_opvs[j].prompt != NULL; ++j)
	    {
	        if (strncasecmp(s+i,g_opvs[j].select,g_opvs[j].ncmp) == 0)
	        {
		    i += g_opvs[j].ncmp-1;
	    	    if (strcmp(g_opvs[j].prompt,"VELOCITY")==0)
	    	    {
	                for (k = 0; k < dim; ++k)
	    	    	{
	    	            scalar(&pl->next,sizeof(PRINTING_LIST));
			    pl = pl->next;
	    	            pl->name = velocity[k];
	    	            pl->var  = var++;
	    	            pl->next = NULL;
	    	        }
	    	    }
		    else
		    {
	    	        scalar(&pl->next,sizeof(PRINTING_LIST));
			pl = pl->next;
			pl->name = g_opvs[j].prompt;
			pl->var  = var++;
			pl->next = NULL;
		    }
		}
	    }
	}
	if (head == NULL)
	    head = dummy_head.next;
	prt->n_rect_state_vars = var;
	return head;
}		/*end g_set_optional_printing_variables*/

/*
*			g_set_plot_choices():
*
*	Initializes and returns a list of the plot variables currently
*	supported.  To add an additional plotting option simply
*	enter the appropriate structure elements in the array of
*	Plot_choice choices below.  This new entry should never
*	be the first or last entry in this array.
*/

EXPORT	Plot_choice	*g_set_plot_choices(
	INIT_DATA	*init)
{
	RECT_GRID *gr = &Comp_grid(init);
	int	  dim = gr->dim;

	static Scalar_Plot_choice choices[] = {
	    {
	        {
	            "x-velocity",
	            "X",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "X_VELOCITY",
	            return_x_velocity,
	            NULL
	        },
	        { NO, YES, YES, YES }
	    },
	    {
	        {
	            "y-velocity",
	            "Y",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "Y_VELOCITY",
	            return_y_velocity,
	            NULL
	        },
	        { NO, NO, YES, YES }
	    },
	    {
	        {
	            "z-velocity",
	            "Z",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "Z_VELOCITY",
	            return_z_velocity,
	            NULL
	        },
	        { NO, NO, NO, YES }
	    },
	    {
	        {
	            "flow speed",
	            "U",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "FLOW_SPEED",
	            return_flow_speed,
	            NULL
	        },
	        { NO, YES, YES, YES }
	    },
	    {
	        {
	            "kinetic energy density",
	            "K",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "KINETIC_ENERGY",
	            return_kinetic_energy,
	            NULL
	        },
	        { NO, YES, YES, YES }
	    },
	    {
	        {
	            "momentum",
	            "M",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "MOMENTUM",
	            return_momentum,
	            NULL
	        },
	        { NO, YES, YES, YES }
	    },
	    {
	        {
	            "density",
	            "D",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "DENSITY",
	            return_density,
	            NULL
	        },
	        { NO, YES, YES, YES }
	    },
	    {
	        {
	            "energy density",
	            "E",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "ENERGY_DENSITY",
	            return_energy_density,
	            NULL
	        },
	        { NO, YES, YES, YES }
	    },
	    {
	        {
	            "vorticity",
	            "VORTICITY",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "VORTICITY",
	            return_vorticity,
	            NULL
	        },
		{ NO, NO, YES, NO }
	    },
	    {
	        {
	            "divergence of velocity",
	            "W",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "DIVERGENCE_VELOCITY",
	            return_divergence_velocity,
	            NULL
	        },
		{ NO, YES, YES, YES }
	    },
	    {
	        {
	            "pressure",
	            "P",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "PRESSURE",
	            return_pressure,
	            NULL
	        },
	        { NO, YES, YES, YES }
	    },
	    {
	        {
	            "specific internal energy",
	            "I",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "SPECIFIC_INTERNAL_ENERGY",
	            return_specific_internal_energy,
	            NULL
	        },
	        { NO, YES, YES, YES }
	    },
	    {
	        {
	            "temperature",
	            "T",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "TEMPERATURE",
	            return_temperature,
	            NULL
	        },
	        { NO, YES, YES, YES } 
	    },
	    {
	        {
	            "specific enthalpy",
	            "H",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "ENTHALPY",
	            return_enthalpy,
	            NULL
	        },
	        { NO, YES, YES, YES }
	    },
	    {
	        {
	            "entropy",
	            "S",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "ENTROPY",
	            return_entropy,
	            NULL
	        },
	        { NO, YES, YES, YES }
	    },
	    {
	        {
	            "radial component of velocity",
	            "RV",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "RADIAL_VELOCITY",
	            return_radial_component_of_velocity,
	            NULL
	        },
	        { NO, YES, YES, YES }
	    },
            {
                {
                    "tangential component of xy velocity",
                    "TV",
                    g_SelectedScalarItems,
                    NULL
                },
                {
                    "TANGENTIAL_VELOCITY",
                    return_tangential_component_of_xy_velocity,
                    NULL
                },
                NO, NO, YES, NO
            },
	    {
	        {
	            "self similar Mach number",
	            "L",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "SS_MACH_NUM",
	            return_ss_Mach_number,
	            NULL
	        },
	        { NO, NO, YES, NO }
	    },
	    {
	        {
	            "Mach number",
	            "N",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "MACH_NUMBER",
	            return_Mach_number,
	            NULL
	        },
	        { NO, YES, YES, YES }
	    },
#if defined(COMBUSTION_CODE)
	    {
	        {
	            "reaction product density",
	            "Q",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "PRODUCT_DENSITY",
	            return_product_density,
	            NULL
	        },
	        { NO, YES, YES, YES }
	    },
	    {
	        {
	            "reaction progress",
	            "R",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "REACTION_PROGRESS",
	            return_reaction_progress,
	            NULL
	        },
	        { NO, YES, YES, YES }
	    },
	    {
	        {
	            "density of first constituent",
	            "C",
	            g_SelectedScalarItems,
	            NULL
	        },
	        {
	            "DENSITY_1",
	            return_density1,
	            NULL
	        },
	        { NO, YES, YES, YES }
	    },
#endif /* defined(COMBUSTION_CODE) */
	    /* Mandatory last element of choices array */
	    {
	        {
	            NULL,
	            NULL,
	            NULL,
	            NULL
	        },
	        {
	            NULL,
	            NULL,
	            NULL
	        },
	        { NO, NO, NO, NO }
	    }
	};
	if ((gr->Remap.remap == CYLINDRICAL_REMAP) ||
	    (gr->Remap.remap == SPHERICAL_REMAP))
	{
	    choices[0].Choice.prompt = "r-velocity";
	    choices[0].Choice.selector = "VR";
	    choices[0].Item.name = "R_VELOCITY";
	}
	if (gr->Remap.remap == CYLINDRICAL_REMAP)
	{
	    choices[1].Choice.prompt = "z-velocity";
	    choices[1].Choice.selector = "VZ";
	    choices[1].Item.name = "Z_VELOCITY";
	}

	plot_choices(init) = choices_valid_for_dim(choices,dim);
	return plot_choices(init);
}		/*end g_set_plot_choices*/

LOCAL	Plot_choice	*choices_valid_for_dim(
	Scalar_Plot_choice	*choices,
	int			dim)
{
	int var;
	Plot_choice	Head, *pc;

	Head.next = NULL;
	for (var = 0, pc = &Head; choices[var].Choice.prompt != NULL; ++var)
	{
	    choices[var].Choice.next = NULL;
	    if (choices[var].use_for_dim[dim] == YES)
	    {
	    	pc->next = &choices[var].Choice;
	    	pc = pc->next;
	    }
	}
	return Head.next;
}		/*end choices_valid_for_dim*/


EXPORT	INPUT_SOLN	**g_set_input_solution(
	Printplot	*prt)
{
	INPUT_SOLN	**is;
	int		var;

	bi_array(&is,prt->n_restart_vars,1,sizeof(INPUT_SOLN));

	for (var = 0; var < prt->n_restart_vars; ++var)
	{
	    is[var]->name = prt->output_soln[var]->name;
	    is[var]->fit = LINEAR_FIT;
	    is[var]->smoothness = SINGULAR;
	    is[var]->set_intfc_states = g_restart_set_intfc_states;
	}
	return is;
}		/*end g_set_input_solution*/

EXPORT void g_init_statistics(
	Front		*front,
	Grid		*grid,
	Wave		*wave,
	Printplot	*prt,
	INIT_DATA	*init)
{
	char		s[Gets_BUF_SIZE];
	static boolean	first = YES;

	/* Do not allow modification of stat printout at pause times. */
	if (first == YES)
	    first = NO;
	else
	    return;

	screen("\n\t\tStatistics Control\n\n");

	screen("Type 'y' to initiate prompting for statistics/diagnostics: ");
	(void) Gets(s);
	if ((s[0] != 'y') && (s[0] != 'Y')) return;

	g_init_grid_statistics(init,front,grid,wave,prt);

	init_peep_time(init,grid,prt);

	init_cross_sections(init,front,grid,prt);
	init_intfc_extrema(init,front,grid,prt);
	init_layer_stats(init,front,grid,prt);
	init_intfc_stats(init,front,grid,prt);
	init_rect_state_stats(init,front,grid,prt);
	if (front->rect_grid->dim == 2)
	    init_moving_body_data(init,front,grid,prt);

	init_statistics_for_problem(front,grid,prt,init);

	screen("\n");
}		/*end g_init_statistics*/


/*
*			init_statistics_for_problem():
*
*	Initializes statistical output on a problem dependent basis.  Two
*	types of output exist.  The first is a specialized function to print
*	various problem specific quantities in columnar format.  The second
*	is printing of front states for some or all curves.  If a specialized
*	function exists to select certain curves, it is used as the default.
*	If no specialized function exists, or if desired, printing is turned
*	on for all curves.
*/

LOCAL	void	init_statistics_for_problem(
	Front		*front,
	Grid		*grid,
	Printplot	*prt,
	INIT_DATA	*init)
{
	boolean		stat_one_io_node = NO;
	boolean		fr_one_io_node = NO;
	boolean		stat_prompt_for_file_name = YES;
	boolean		stat_open_file = YES;
	char		s[Gets_BUF_SIZE];
	const char	*stat_prompt;
	const char	*fr_prompt;
	char	        *stat_fname = NULL;
	char	        *fr_fname = NULL;
	int		i, num_layers = 1;
	OUTPUT_DATA	**(*alloc_output_datas)(int) = d_alloc_output_datas;
	void		(*stat_func)(Grid*,Wave*,Front*,Printplot*,
				     OUTPUT_DATA*,boolean) = NULL;
	void		(*fr_func)(Grid*,Wave*,Front*,Printplot*,
				   OUTPUT_DATA*,boolean) = NULL;
	void		(*prob_stat_setup)(OUTPUT_DATA**,int) = NULL;
	char		ffn[256], sfn[256];

	switch(prt_problem_type(prt))
	{
	case ASTROPHYSICAL_JET:
	case INJECTION_INLET:
	    stat_func = record_jet_velocity;
	    stat_prompt = "the jet velocity";
	    break;
	case RICHTMYER_MESHKOV:
	case MESHKOV:
	    alloc_output_datas = rm_alloc_output_datas;
	    num_layers = num_layers(front->interf);
	    stat_func = record_rm_amp_and_vel_data;
	    stat_one_io_node = YES;
	    stat_prompt_for_file_name = YES;
	    stat_open_file = (restart_io_type(init) == NULL) ? YES : NO;
	    stat_prompt = "the amplitude and velocity data";
	    stat_fname = set_default_output_file_name(prt->outfile,
						      "rm-avd",sfn);
	    fr_func = show_front_states_for_rm_problem;
	    fr_prompt = "contacts and rarefaction edges";
	    fr_fname = set_default_output_file_name(prt->outfile,
						    "rm-fronts",ffn);
	    prob_stat_setup = set_rm_graph_scale;
	    break;
	case SHOCKED_THERMAL_LAYER:
	    stat_func = NULL;
	    fr_func = show_front_states_along_lower_wall;
	    fr_prompt = "the lower wall";
	    break;
	case SUPERNOVA:
	case IMPLOSION:
	    alloc_output_datas = rm_alloc_output_datas;
	    num_layers = num_layers(front->interf);
	    stat_one_io_node = YES;
	    stat_prompt_for_file_name = YES;
	    stat_open_file = (restart_io_type(init) == NULL) ? YES : NO;
	    stat_func = record_radial_rm_amp_and_vel_data;
	    stat_prompt = "the radial amplitude and velocity data";
	    stat_fname = set_default_output_file_name(prt->outfile,
						      "rad-rm-avd",sfn);
	    fr_func = show_contact_states;
	    fr_prompt = "the contacts";
	    fr_fname = set_default_output_file_name(prt->outfile,
						    "rad-rm-fronts",ffn);
	    prob_stat_setup = set_rm_graph_scale;
	    break;
	case BUBBLES_DROPS:
	    fr_func = show_contact_and_dirichlet_bdry_states;
            fr_prompt = "\n\tthe contacts and dirichlet boundaries";
            break;
	case EXPANDING_SHELLS:
	case SHOCK_JET:
	case PLANE_FRONT:
	    fr_func = show_contact_states;
	    fr_prompt = "the contacts";
	    break;
	case BOWSHOCK:
	    fr_func = show_front_states_for_bowshocks;
	    fr_prompt = "the bow shock, wedge, and exit";
	    break;
	case RAMP_REFLECTION:
	    stat_func = print_ramp_refl_stats;
	    stat_prompt = "the ramp reflection statistics";
	    fr_func = show_front_states_for_ramp_reflections;
	    fr_prompt = "the physical curves and boundaries inside bubble";
	    break;
	case EXPANDING_SHOCK:
	    fr_func = show_front_states_for_expanding_shocks;
	    fr_prompt = "the expanding shock and walls";
	    break;
	case RANDOM_SURFACE:
	case RAYLEIGH_TAYLOR:
	    init_Lp_Diff(init,front,grid,prt,GAS_STATE,lin_pert_state);
	    if (front->rect_grid->dim == 2)
	    {
	        stat_func = multi_bubble_velocities;
	        stat_prompt = "multi-bubble velocities";
	        stat_fname = set_default_output_file_name(prt->outfile,
							  "rt-mbv",sfn);
	        fr_func = show_contact_states;
	        fr_prompt = "the contacts";
	        fr_fname = set_default_output_file_name(prt->outfile,
							"rt-fronts",ffn);
	    }
	    break;
	case RADIAL_RAYLEIGH_TAYLOR:/*TODO*/
	    stat_func = NULL;
	    fr_func = NULL;
	    stat_prompt = NULL;
	    fr_prompt = NULL;
	    break;
	default:
	    stat_func = NULL;
	    fr_func = NULL;
	    stat_prompt = NULL;
	    fr_prompt = NULL;
	    break;
	}

	if (stat_func != NULL)
	{
	    screen("Type 'y' to request printing of %s: ",stat_prompt);
	    (void) Gets(s);
	    if ((s[0] == 'y') || (s[0] == 'Y'))
	    {
	    	OUTPUT_DATA	**datas, *data;
		boolean		open_file = stat_open_file;
		char		sfname[256];

		datas = (*alloc_output_datas)(num_layers);
		if (num_layers > 1)
		    open_file = NO;

		data = datas[0];
		if (stat_fname != NULL)
		    (void) strcpy(Output_filename(data),stat_fname);
	    	init_output_data(init,data,grid,prt,stat_prompt,
	    			 stat_prompt_for_file_name,open_file,
				 stat_one_io_node);
	    	add_user_output_function(stat_func,data,prt);

		if (num_layers > 1)
		{
		    for (i = 1; i < num_layers; ++i)
		    {
			(void) sprintf(sfname,"%s.%d",
				       Output_filename(data),i+1);
		        copy_output_data(datas[i],data,open_file,sfname);
	    	        add_user_output_function(stat_func,datas[i],prt);
		    }
		    (void) sprintf(Output_filename(data),"%s.1",
				   Output_filename(data));
		    if (stat_open_file == YES)
		    {
		        for (i = 0; i < num_layers; ++i)
			    create_output_data_file(datas[i]);
		    }
		}

	        if (prob_stat_setup != NULL)
		    (*prob_stat_setup)(datas,num_layers);
	    }
	}

	if (fr_func != NULL)
	{
	    screen("Type 'y' to request printing of front states along "
	           "%s: ",fr_prompt);
	    (void) Gets(s);
	    if (s[0] != 'y' && s[0] != 'Y')
	    	fr_func = NULL;
	}

	if ((fr_func == NULL) && (front->_fgraph_front_states != NULL))
	{
	    screen("Type 'y' to request printing of states for all fronts: ");
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y')
	    	fr_func = fgraph_front_states_shell;
	    (void) printf("\n");
	}

	if (fr_func != NULL)
	{
	    OUTPUT_DATA	*data;
	    char *fname = NULL;

	    scalar(&data,sizeof(OUTPUT_DATA));

	    /* defaults */

	    Output_mode(data) = CONSTANT_TIME;
	    Output_time_freq(data) = 1.0;
	    Output_start_time(data) = 0;
	    Output_in_binary(data) = NO;

	    if (fr_fname == NULL)
	    {
	        set_default_data_file_names("front_states","-sts",&fr_fname,
					    &fname,init);
	    }
	    if (fr_fname != NULL)
		(void) strcpy(Output_filename(data),fr_fname);
	    init_output_data(init,data,grid,prt,"front states",YES,YES,
			     fr_one_io_node);
	    add_user_output_function(fr_func,data,prt);
	}

}		/*end init_statistics_for_problem*/


LOCAL	void	set_rm_graph_scale(
	OUTPUT_DATA	**datas,
	int		num_layers)
{
	GraphUnits	*rmunits = RmGraphUnits(datas);
	char	s[120];

	screen("Enter an optional length scale factor (default = %g): ",
		rmunits->_length_scale);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&rmunits->_length_scale);
	screen("Enter an optional initial lower range for lengths ");
	screen("in the given scale\n");
	screen("\t(default = %g): ",rmunits->_initial_lmin);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    (void) sscan_float(s,&rmunits->_initial_lmin);
	    rmunits->_initial_lmin *= rmunits->_length_scale;
	}
	screen("Enter an optional initial upper range for lengths ");
	screen("in the given scale\n");
	screen("\t(default = %g): ",rmunits->_initial_lmax);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    (void) sscan_float(s,&rmunits->_initial_lmax);
	    rmunits->_initial_lmax *= rmunits->_length_scale;
	}

	screen("Enter an optional time scale factor (default = %g): ",
		rmunits->_time_scale);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&rmunits->_time_scale);
	screen("Enter an optional initial lower range for times ");
	screen("in the given scale\n");
	screen("\t(default = %g): ",rmunits->_initial_tmin);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    (void) sscan_float(s,&rmunits->_initial_tmin);
	    rmunits->_initial_tmin *= rmunits->_time_scale;
	}
	screen("Enter an optional initial upper range for times ");
	screen("in the given scale\n");
	screen("\t(default = %g): ",rmunits->_initial_tmax);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    (void) sscan_float(s,&rmunits->_initial_tmax);
	    rmunits->_initial_tmax *= rmunits->_time_scale;
	}

	rmunits->_velocity_scale = rmunits->_length_scale/rmunits->_time_scale;
	screen("Enter an optional velocity scale  factor (default = %g): ",
		rmunits->_velocity_scale);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&rmunits->_velocity_scale);
	screen("Enter an optional initial lower range for velocities ");
	screen("in the given scale\n");
	screen("\t(default = %g): ",rmunits->_initial_vmin);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    (void) sscan_float(s,&rmunits->_initial_vmin);
	    rmunits->_initial_vmin *= rmunits->_velocity_scale;
	}
	screen("Enter an optional initial upper range for velocities ");
	screen("in the given scale\n");
	screen("\t(default = %g): ",rmunits->_initial_vmax);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    (void) sscan_float(s,&rmunits->_initial_vmax);
	    rmunits->_initial_vmax *= rmunits->_velocity_scale;
	}

	set_rm_layer_indices(datas,num_layers);
}		/*end set_rm_graph_scale*/

/*
*			fgraph_front_states_shell():
*
*	This function is simply a shell to call fgraph_front_states(), but
*	taking the correct argument list for the user functions in 
*	g_User_printplot.
*/

/*ARGSUSED*/
LOCAL void fgraph_front_states_shell(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	fgraph_front_states(Output_file(data),front);
}		/*end fgraph_front_states_shell*/


/* ARGSUSED */
LOCAL double return_x_velocity(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	if (is_obstacle_state(state)) return 0.0;
	return vel(0,state);
}		/*end return_x_velocity*/

/* ARGSUSED */
LOCAL double return_y_velocity(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	if (is_obstacle_state(state)) return 0.0;
	return vel(1,state);
}		/*end return_y_velocity*/

/* ARGSUSED */
LOCAL double return_z_velocity(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	if (is_obstacle_state(state)) return 0.0;
	return vel(2,state);
}		/*end return_z_velocity*/

/* ARGSUSED */
LOCAL double return_momentum(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	int		i, dim;
	double		tmp = 0.0;

	if (is_obstacle_state(state)) return 0.0;
	dim = Params(state)->dim;
	for (i = 0; i < dim; ++i)
		tmp += sqr(Mom(state)[i]);
	return sqrt(tmp);
}		/*end return_momentum*/

/* ARGSUSED */
LOCAL double return_density(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	if (is_obstacle_state(state)) return 0.0;
	return Dens(state);
}		/*end return_density*/

/* ARGSUSED */
LOCAL double return_energy_density(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	if (is_obstacle_state(state)) return 0.0;
	return Energy(state);
}		/*end return_energy_density*/

/* ARGSUSED */
LOCAL double return_pressure(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	if (is_obstacle_state(state)) return 0.0;
	return pressure(state);
}		/*end return_pressure*/

/* ARGSUSED */
LOCAL double return_specific_internal_energy(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	if (is_obstacle_state(state)) return 0.0;
	return specific_internal_energy(state);
}		/*end return_specific_internal_energy*/

/* ARGSUSED */
LOCAL double return_temperature(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	if (is_obstacle_state(state)) return 0.0;
	return temperature(state);
}		/*end return_temperature*/

/* ARGSUSED */
LOCAL double return_entropy(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	if (is_obstacle_state(state)) return 0.0;
	return entropy(state);
}		/*end return_entropy*/

/* ARGSUSED */
LOCAL double return_radial_component_of_velocity(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	return RadialComponentOfVelocity(state,coords,fr->rect_grid->dim);
}		/*end return_radial_component_of_velocity*/

/* ARGSUSED */
LOCAL double return_tangential_component_of_xy_velocity(
        double           *coords,
        Front           *fr,
	POINTER		wv,
        COMPONENT       comp,
        Locstate        state)
{
        return TangCmptOfXYVelocity(state,coords,fr->rect_grid->dim);
}               /*end return_tangential_component_of_xy_velocity*/

/* ARGSUSED */
LOCAL double return_enthalpy(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	if (is_obstacle_state(state)) return 0.0;
	return specific_enthalpy(state);
}		/*end return_enthalpy*/


/* ARGSUSED */
LOCAL double return_flow_speed(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	int i, dim = fr->rect_grid->dim;
	double v[MAXD];

	if (is_obstacle_state(state)) return 0.0;

	for (i = 0; i < dim; ++i)
		v[i] = vel(i,state);
	return mag_vector(v,dim);
}		/*end return_flow_speed*/


/* ARGSUSED */
LOCAL double return_kinetic_energy(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	return kinetic_energy(state);
}		/*end return_kinetic_energy*/

/*
*			return_ss_Mach_number():
*
*	Computes the self similar Mach number for a given state relative
*	to the corner position in a ramp reflection problem type.
*
*	WARNING: this diagnostic is specific to RAMP_REFLECTION problems.
*	It should work in parallel.
*/

/* ARGSUSED */
LOCAL	double return_ss_Mach_number(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	double		ss_vel[MAXD];
	double		ss_Mach_num;
	double		t;
	double		ss_origin[MAXD];
	int		i, dim = fr->rect_grid->dim;

	if (is_obstacle_state(state)) return 0.0;
	ramp_reflection_corner_posn(ss_origin,NO,dim);

	t = fr->time + fr->dt;
	if (t == 0) return 0.0;
	for (i = 0; i < dim; ++i)
		ss_vel[i] = vel(i,state) -
			(coords[i] - ss_origin[i]) / t;

	ss_Mach_num = mag_vector(ss_vel,dim) / sound_speed(state);

	return ss_Mach_num;
}		/*end return_ss_Mach_number*/

/* ARGSUSED */
LOCAL	double return_Mach_number(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	double		v[MAXD];
	int		i, dim = fr->rect_grid->dim;

	if (is_obstacle_state(state)) return 0.0;

	for (i = 0; i < dim; ++i)
		v[i] = vel(i,state);
	return mag_vector(v,dim) / sound_speed(state);
}		/*end return_Mach_number*/

/* ARGSUSED */
LOCAL	double return_vorticity(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	double		w[MAXD], ans;
	double		**gradv;
	int		dim = fr->interf->dim;
	int		i;

	if (is_obstacle_state(state) || dim < 2)
	    return 0.0;
	gradv = velocity_gradient(coords,fr,(Wave*)wv,comp);
	if (dim == 2)
	{
	    ans = gradv[1][0] - gradv[0][1];
	}
	else if (dim == 3)
	{
	    ans = 0.0;
	    for (i = 0; i < dim; ++i)
	    {
	    	w[i] = gradv[(i+2)%dim][(i+1)%dim] -
	    			gradv[(i+1)%dim][(i+2)%dim];
	    	ans += sqr(w[i]);
	    }
	    ans = sqrt(ans);
	}
	return ans;
}		/*end return_vorticity*/

/* ARGSUSED */
LOCAL	double return_divergence_velocity(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	RECT_GRID *gr = fr->rect_grid;
	double	  ans;
	double	  **gradv;
	int	  dim = gr->dim;
	int	  i;

	if (is_obstacle_state(state))
	    return 0.0;
	gradv = velocity_gradient(coords,fr,(Wave*)wv,comp);
	ans = 0.0;
	for (i = 0; i < dim; ++i)
	    ans += gradv[i][i];
	if ((dim == 2) && (gr->Remap.remap == CYLINDRICAL_REMAP))
	{
	    double rmin = pos_radius(0.0,gr);
	    double r = pos_radius(coords[0],gr);
	    if (r > rmin)
	        ans += vel(0,state)/r;
	}
	return ans;
}		/*end return_divergence_velocity*/


/*
*			velocity_gradient():
*
*	Returns a finite difference approximation to the velocity gradient
*
*		gradv[i][j] = v[i]_x[j]
*/

LOCAL	double **velocity_gradient(
	double		*coords,
	Front		*fr,
	Wave		*wave,
	COMPONENT	comp)
{
	INTERFACE	*intfc = fr->interf;
	double		*h = fr->rect_grid->h;
	double		crds[MAXD];
	double		v[1<<MAXD][MAXD];
	double		fac;
	COMPONENT	cmp;
	int		i, j, k;
	int		dim = fr->interf->dim;
	int		N = 1<<dim;
	static double	sgn[8][3] = {	{-1, -1, -1},
					{ 1, -1, -1},
					{-1,  1, -1},
					{ 1,  1, -1},
					{-1, -1,  1},
					{ 1, -1,  1},
					{-1,  1,  1},
					{ 1,  1,  1}	};
	static Locstate tmpst = NULL;
	static double	**gradv = NULL;

	if (tmpst == NULL)
	{
	    alloc_state(fr->interf,&tmpst,fr->sizest);
	    bi_array(&gradv,MAXD,MAXD,FLOAT);
	}

	for (k = 0; k < N; ++k)
	{
	    for (i = 0; i < dim; ++i)
	    	crds[i] = coords[i] + 0.5*sgn[k][i]*h[i];
	    cmp = component(crds,intfc);
	    if (is_excluded_comp(cmp,intfc))
	        cmp = comp;
	    hyp_solution(crds,cmp,NULL,UNKNOWN_SIDE,fr,wave,tmpst,NULL);
	    for (i = 0; i < dim; ++i)
	    	v[k][i] = vel(i,tmpst);
	}

	for (j = 0; j < dim; ++j)
	{
	    fac = 2.0/(N*h[j]);
	    for (i = 0; i < dim; ++i)
	    {
	    	for (gradv[i][j] = 0.0, k = 0; k < N; ++k)
	    	    gradv[i][j] += sgn[k][j]*v[k][i];
	    	gradv[i][j] *= fac;
	    }
	}
	return gradv;
}		/*end velocity_gradient*/

#if defined(COMBUSTION_CODE)
/* ARGSUSED */
LOCAL	double return_product_density(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	if (is_obstacle_state(state)) return 0.0;
	if (Composition_type(state) == PURE_NON_REACTIVE) return 0.0;
	if (Composition_type(state) == PTFLAME)
		return (Unburned(state)) ? 0.0 : Dens(state);
	return Prod(state);
}		/*end return_product_density*/

/* ARGSUSED */
LOCAL	double return_reaction_progress(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	if (Composition_type(state) == PURE_NON_REACTIVE) return 0.0;
	if (is_obstacle_state(state)) return 0.0;
	if (Composition_type(state) == PTFLAME)
		return (Unburned(state)) ? 0.0 : 1.0;
	return Prod(state) / Dens(state);
}		/*end return_reaction_progress*/

/* ARGSUSED */
LOCAL	double return_density1(
	double		*coords,
	Front		*fr,
	POINTER		wv,
	COMPONENT	comp,
	Locstate	state)
{
	if (Composition_type(state) != TWO_CONSTITUENT_REACTIVE)
		return 0.0;
	return Dens1(state);
}		/*end return_density1*/
#endif /* defined(COMBUSTION_CODE) */


LOCAL	char	*set_default_output_file_name(
	const char	*root,
	const char	*identifier,
	char		*buf)
{
	char	*bname, *dname;

	if (root == NULL)
	    return NULL;
	if (identifier == NULL)
	    identifier = "";

	if (buf == NULL)
	{
	    static char	  *fname_store = NULL;
	    static size_t len_store = 0;
	    if (len_store == 0)
	    {
	       	len_store = 256 + strlen(root) + strlen(identifier);
	        uni_array(&fname_store,len_store,CHAR);
	    }
	    else if ((strlen(root) + strlen(identifier)) >= len_store)
	    {
	       	free(fname_store);
	        len_store = 256 + strlen(root) + strlen(identifier);
	        uni_array(&fname_store,len_store,CHAR);
	    }
	    buf = fname_store;
	}

	base_and_dir_name(root,&dname,&bname);
	(void) sprintf(buf,"%s%s%s%s%s",
		       (strlen(dname)!=0) ? dname : "",
		       (strlen(dname)!=0) ? "/" : "",
		       identifier,
		       (strlen(identifier)!=0 && strlen(bname)!=0) ? "/" : "",
		       bname);
	if (strlen(buf)==0)
	    return NULL;
	return buf;
}		/*end set_default_output_file_name*/


/*
*			init_Lp_Diff():
*
*	Initializes computation and printing of Lp norm between non-linear
*	and linearized solutions for interface instability problems.  This is
*	a problem dependent initialization, and should only be used for
*	the RANDOM_SURFACE problem type.
*/

LOCAL	void	init_Lp_Diff(
	INIT_DATA	*init,
	Front		*front,
	Grid		*grid,
	Printplot	*prt,
	int		stype,
	void		(*alt_state)(Locstate,double*,Front*,Wave*,int))
{
	SINE_PERT	*pert = Sine_pert(comp_type(COMPB));
	Lp_Diff_data	*lpd_data;
	char		s[Gets_BUF_SIZE];
	char		*lp_diff_file;

	screen("Type 'y' if you wish to compute the Lp norm of the \n ");
	screen("\tdifference between the linearized and nonlinear solutions: ");
	(void) Gets(s);
	if (s[0] != 'y' && s[0] != 'Y')
	    return;

	lp_diff_file = set_default_output_file_name(prt->outfile,"lp_diff",s);
	scalar(&lpd_data,sizeof(Lp_Diff_data));
	lpd_data->pert = alloc_and_copy_sine_pert_structure(pert,front);
	if (lp_diff_file != NULL)
	    (void) strcpy(Output_filename(&lpd_data->odata),lp_diff_file);
	init_output_data(init,&lpd_data->odata,grid,prt,"Lp data",YES,YES,NO);

	lpd_data->stype = stype;
	lpd_data->_alternate_state = alt_state;

	lpd_data->p = 2;
	screen("Enter the number p (p = 1 or 2, default = %d): ",lpd_data->p);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%d",&lpd_data->p);

	lpd_data->ratio = 10;
	screen("Enter the ratio of refinement inside the mixing zone ");
	screen("(default = %d): ",lpd_data->ratio);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%d",&lpd_data->ratio);

	add_user_output_function(Lp_diff,&lpd_data->odata,prt);
}		/*end init_Lp_Diff*/

/*
*				init_probe()
*
*	Initializes a data probe.
*
*	Input:
*		center - if non-null locate probe at coords center
*		         otherwise the code will prompt for the probe center.
*		rad    - probe samples over a rectangle with center center
*			 and half side lengths rad.
*		n      - sample n[i] points per distance rad[i] in direction i
*		ptol   - tolerance for detecting wave events
*		init   - INIT_DATA structure needed as data for the output
*			 functions
*		ip     - INIT_PHYSICS structure containing front and grid
*			 data
*
*	Output:
*		An output data structure is initialized to print probe
*		output data as a function of time.
*/

EXPORT	PROBE	*init_probe(
	double        *center,
	double        *rad,
	int          *n,
	const double  *ptol,
	INIT_DATA    *init,
	INIT_PHYSICS *ip)
{
	Grid      *grid = ip->root->grid;
	Front     *front = ip->root->front;
	INTERFACE *intfc = front->interf;
	Printplot *prt = ip->prt;
	PROBE	  *probe;
	RECT_GRID *gr = front->rect_grid;
	char	  s[Gets_BUF_SIZE];
	double     d[3];
	int       i, j, k, dim = gr->dim;
	int       indx, N;
	size_t    sizest = front->sizest;

	scalar(&probe,sizeof(PROBE));
	probe->line = 0;
	probe->time[0] = probe->time[1] = probe->time[2] = grid->time;
	probe->dim = dim;
	probe->L = mag_vector(gr->h,dim)/sqrt((double)dim);

	for (i = 0; i < 3; ++i)
	    probe->c[i] = 0.0;
	if (center != NULL)
	{
	    for (i = 0; i < dim; ++i)
	        probe->c[i] = center[i];
	}
	else
	{
	    screen("Enter the coordinates of the probe center: ");
	    if (Scanf("%f %f %f\n",probe->c,probe->c+1,probe->c+2) != dim)
	    {
	        screen("ERROR in init_probe(), improper input of "
		       "probe center\n");
		clean_up(ERROR);
	    }
	}
	for (i = 0; i < 3; ++i)
	    probe->r[i] = 0.5;
	if (rad != NULL)
	{
	    for (i = 0; i < dim; ++i)
	        probe->r[i] = rad[i];
	}
	else
	{
	    screen("Enter the radii of the probe box: ");
	    if (Scanf("%f %f %f\n",probe->r,probe->r+1,probe->r+2) != dim)
	    {
	        screen("ERROR in init_probe(), improper input of "
		       "probe radius\n");
		clean_up(ERROR);
	    }
	}
	for (i = 0; i < 3; ++i)
	    probe->n[i] = 0;
	if (n != NULL)
	{
	    for (i = 0; i < dim; ++i)
	        probe->n[i] = n[i];
	}
	else
	{
	    screen("Enter the number of sample points per direction (dflt =");
	    for (i = 0; i < dim; ++i)
	        screen(" %d",probe->n[i]+1);
	    screen("): ");
	    if (Scanf("%d %d %d\n",probe->n,probe->n+1,probe->n+2) != dim)
	    {
	        screen("ERROR in init_probe(), improper input of "
		       "the number of sample points per direction\n");
		clean_up(ERROR);
	    }
	    for (i = 0; i < dim; ++i)
	        probe->n[i] = max(0,probe->n[i]-1);
	}
	probe->wave_tag = 0.0;
	if (ptol != NULL)
	    probe->ptol = *ptol;
	else
	{
	    probe->ptol = 0.05;
	    screen("Enter a dimensionless tolerance to detect wave crossings "
	           "(dflt = %g): ",probe->ptol);
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscan_float(s,&probe->ptol);
	}

	n = probe->n;
	probe->N = N = (2*n[0]+1)*(2*n[1]+1)*(2*n[2]+1);
	uni_array(&probe->pcomp,N,sizeof(COMPONENT));
	uni_array(&probe->pstate,N,sizeof(Locstate));
	uni_array(&probe->pcoords,3*N,FLOAT);
	for (i = 0; i < 3; ++i)
	    d[i] = probe->r[i]/(probe->n[i] + 0.5);
	alloc_state(intfc,probe->st,max(sizest,sizeof(VGas)));
	alloc_state(intfc,probe->st+1,max(sizest,sizeof(VGas)));
	alloc_state(intfc,probe->st+2,max(sizest,sizeof(VGas)));
	alloc_state(intfc,&probe->dst,max(sizest,sizeof(VGas)));
	alloc_state(intfc,&probe->Ndst,max(sizest,sizeof(VGas)));
	for (i = -n[0]; i <= n[0]; ++i)
	for (j = -n[1]; j <= n[1]; ++j)
	for (k = -n[2]; k <= n[2]; ++k)
	{
	    indx = probe_index(i,j,k,probe);
	    alloc_state(intfc,probe->pstate+indx,sizest);
	    probe_coords(i,j,k,probe)[0] = probe->c[0] + i*d[0];
	    probe_coords(i,j,k,probe)[1] = probe->c[1] + j*d[1];
	    probe_coords(i,j,k,probe)[2] = probe->c[2] + k*d[2];
	}
	init_output_data(init,&probe->odata,grid,prt,
	                        "the probe data",YES,NO,NO);
	add_user_output_function(evaluate_probe,&probe->odata,prt);


	return probe;
}		/*end init_probe*/

LOCAL	int	g_SelectedScalarItems(
	const char	*c,
	Plot_choice	*pc,
	ScalarPlotItem	**spi)
{
	if (strcasecmp(c,pc->selector) != 0)
	{
	    if (spi != NULL)
	    	*spi = NULL;
	    return 0;
	}
	if (spi != NULL)
	{
	    *spi = scalar_plot_item(pc);
	    (*spi)->next = NULL;
	}
	return 1;
}		/*end g_SelectedScalarItems*/

EXPORT	void fprint_rarefaction_wave_1d(
	FILE                 *file,
	_RAREFACTION_WAVE_1D *rw1d)
{
	char mesg[1024];
	(void) fprintf(file,"_RAREFACTION_WAVE_1D structure 0x%p\n",rw1d);
	(void) print_wave_family(file,"l_or_r = ",rw1d->l_or_r,"\n");
	(void) fprintf(file,"zbar = %g, tbar = %g\n",rw1d->zbar,rw1d->tbar);
	(void) fprintf(file,"zl = %g, zt = %g\n",rw1d->zl,rw1d->zt);
	(void) fprintf(file,"zmin = %g, stmin = 0x%p\n",rw1d->zmin,rw1d->stmin);
	(void) fprintf(file,"zmax = %g, stmax = 0x%p\n",rw1d->zmax,rw1d->stmax);
	(void) sprintf(mesg,"stl 0x%p at position zl = %g",
		       rw1d->stl,rw1d->zl);
	verbose_fprint_state(file,mesg,rw1d->stl);
	(void) sprintf(mesg,"stt 0x%p at position zt = %g",
		       rw1d->stt,rw1d->zt);
	verbose_fprint_state(file,mesg,rw1d->stt);
	(void) printf("End _RAREFACTION_WAVE_1D structure 0x%p\n",rw1d);
}		/*end fprint_rarefaction_wave_1d*/
