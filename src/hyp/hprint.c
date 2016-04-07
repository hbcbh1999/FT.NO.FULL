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
*			hprint.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Storage allocation/free routines, print routines
*	for the trigrid construction
*/


#include <hyp/hdecs.h>

	/* LOCAL Function Declarations */
LOCAL	void	print_el_integral_structure(EL_INTEGRALS*);
LOCAL	void	print_interpolator_structure(INTERPOLATORS*);
LOCAL	void	print_unsplit_structure(UNSPLIT*);
#if defined(TWOD)
LOCAL	void	nearest_interior_state(double*,COMPONENT,Front*,Wave*,
				       INTERFACE*,double*,Locstate);
#endif /* defined(TWOD) */

#if defined(TWOD)

EXPORT	void	print_tri_soln(
	FILE		*file,
	Front		*front,
	Wave		*wave,
	TRI_SOLN	*soln,
	size_t		n_vars,
	double		(**fcn)(double*,Front*,POINTER,COMPONENT,Locstate))
{
	TRI_GRID	*grid = soln->tri_grid;
	Table		*T = table_of_interface(grid->grid_intfc);
	BILINEAR_ELEMENT *bilin = grid->bilin_els;
	COMPONENT	*comp = T->components;
	LINEAR_ELEMENT	*lin = grid->lin_els;
	TG_PT		*node = grid->node_points;
	boolean		binary_output = is_binary_output();
	char		field[15];
	double		coords_on[MAXD];
	double		coords[MAXD];
	int		n_reg_nodes, n, m, n_NODES;
	int		i, j;
	register TG_PT	*p_0 = grid->node_points;
	register Locstate *state = grid->states;
	static const int	FIELD_SIZE = 15;
	static const char	*FIELD = "%-15g";
	static Locstate	tmpst = NULL;
	static size_t	num_alloc = 0;
	static double	*node_binary = NULL;
	static char	*node_text = NULL;

	if (front->rect_grid->dim != 2)
		return;

	if (tmpst == NULL)
		alloc_state(front->interf,&tmpst,front->sizest);
	n_reg_nodes = (grid->rect_grid.gmax[0]+1) * (grid->rect_grid.gmax[1]+1);
	n_NODES = n_reg_nodes + 3*grid->n_lin_els;

	if (num_alloc < 2 + n_vars)
	{
		if (num_alloc != 0)
		{
			free(node_binary);
			free(node_text);
		}
		num_alloc = 2+n_vars;
		uni_array(&node_binary,2+n_vars,              sizeof(double));
		uni_array(&node_text, (2+n_vars)*FIELD_SIZE+1,sizeof(char));
	}

	(void) fprintf(file,"NODES %d 1 %d 1\n",n_NODES,binary_output);

		/* Print out regular NODE data */

	if (binary_output == YES)
	{
	    for (n = 0;  n < n_reg_nodes;  ++n, ++node, ++state, ++comp)
	    {
		if (is_exterior_comp(*comp,grid->grid_intfc))
		{
		    coords[0] = Coords(node)[0];
		    coords[1] = Coords(node)[1];
		    nearest_interior_state(coords,*comp,front,wave,
					   soln->intfc,coords_on,tmpst);
		    node_binary[0] = coords_on[0];
		    node_binary[1] = coords_on[1];
		    for (i = 0; i < n_vars; ++i)
		    {
			node_binary[2+i] = (*fcn[i])(coords_on,front,
						  (POINTER)wave,*comp,tmpst);
		    }
	        }
		else
		{
		    node_binary[0] = coords[0] = Coords(node)[0];
		    node_binary[1] = coords[1] = Coords(node)[1];
		    for (i = 0; i < n_vars; ++i)
		    {
			node_binary[2+i] = (*fcn[i])(coords,front,
						  (POINTER)wave,*comp,*state);
		    }
	        }
		(void) fwrite((const void *) node_binary,sizeof(double),
			      2+n_vars,file);
	    }
	}
	else
	{
	    for (n = 0;  n < n_reg_nodes;  ++n, ++node, ++state, ++comp)
	    {
		if (is_exterior_comp(*comp,grid->grid_intfc))
		{
		    coords[0] = Coords(node)[0];
		    coords[1] = Coords(node)[1];
		    nearest_interior_state(coords,*comp,front,wave,
					   soln->intfc,coords_on,tmpst);
		    (void) sprintf(node_text,FIELD,coords_on[0]);
		    (void) sprintf(field,FIELD,coords_on[1]);
		    node_text = strcat(node_text,field);
		    for (i = 0; i < n_vars; ++i)
		    {
			(void) sprintf(field,FIELD,(*fcn[i])(coords_on,front,
					     (POINTER)wave,*comp,tmpst));
			node_text = strcat(node_text,field);
		    }
	        }
		else
		{
		    coords[0] = Coords(node)[0];
		    coords[1] = Coords(node)[1];
		    (void) sprintf(node_text,FIELD,coords[0]);
		    (void) sprintf(field,FIELD,coords[1]);
		    node_text = strcat(node_text,field);
		    for (i = 0; i < n_vars; ++i)
		    {
			(void) sprintf(field,FIELD,(*fcn[i])(coords,front,
					     (POINTER)wave,*comp,*state));
			node_text = strcat(node_text,field);
		    }
	        }
		(void) fprintf(file,"%s\n",node_text);
	    }
	}

		/* Print out irregular NODE data */

	if (binary_output == YES)
	{
	    for (n = 0;  n < grid->n_lin_els;  ++n, ++lin)
	    {
		for (j = 0; j < 3; ++j)
		{
		    node_binary[0] = Coords(lin->p[j])[0];
		    node_binary[1] = Coords(lin->p[j])[1];
		    for (i = 0; i < n_vars; ++i)
		    {
			node_binary[2+i] = (*fcn[i])(Coords(lin->p[j]),front,
					  (POINTER)wave,lin->comp,lin->s[j]);
		    }
		    (void) fwrite((const void *) node_binary,sizeof(double),
				  2+n_vars,file);
		}
	    }
	}
	else
	{
	    for (n = 0;  n < grid->n_lin_els;  ++n, ++lin)
	    {
		for (j = 0; j < 3; ++j)
		{
		    (void) sprintf(node_text,FIELD,Coords(lin->p[j])[0]);
		    (void) sprintf(field,FIELD,Coords(lin->p[j])[1]);
		    node_text = strcat(node_text,field);
		    for (i = 0; i < n_vars; ++i)
		    {
			(void) sprintf(field,FIELD,
				       (*fcn[i])(Coords(lin->p[j]),front,
					 (POINTER)wave,lin->comp,lin->s[j]));
			node_text = strcat(node_text,field);
		    }
		    (void) fprintf(file,"%s\n",node_text);
		}
	    }
	}

		/* Print out RECTANGLE indices */

	(void) fprintf(file,"RECTANGLES %d 0 %d 1\n",
		       grid->n_bilin_els,binary_output);
	if (binary_output == YES)
	{
	    for (n = 0;  n < grid->n_bilin_els;  ++n, ++bilin)
	    {
	    	int is[4];
    
		is[0] = (int)(bilin->p[0]-p_0);
		is[1] = (int)(bilin->p[1]-p_0);
		is[2] = (int)(bilin->p[3]-p_0);
		is[3] = (int)(bilin->p[2]-p_0);
		(void) fwrite((const void *) is,sizeof(int),4,file);
	    }
	}
	else
	{
	    for (n = 0;  n < grid->n_bilin_els;  ++n, ++bilin)
	    	(void) fprintf(file,"%llu %llu %llu %llu\n",
	    		       ptr2ull(bilin->p[0]-p_0),
	    		       ptr2ull(bilin->p[1]-p_0),
	    		       ptr2ull(bilin->p[3]-p_0),
	    		       ptr2ull(bilin->p[2]-p_0));
	}

		/* Print out TRIANGLE indices */

	lin = grid->lin_els;
	m   = n_reg_nodes;
	(void) fprintf(file,"TRIANGLES %d 0 %d 1\n",
		       grid->n_lin_els,binary_output);
	if (binary_output == YES)
	{
	    for (n = 0;  n < grid->n_lin_els;  ++n, ++lin, m+=3)
	    {
	    	int is[3];

		is[0] = (int)(m);
		is[1] = (int)(m+1);
		is[2] = (int)(m+2);
		(void) fwrite((const void *) is,sizeof(int),3,file);
	    }
	}
	else
	{
	    for (n = 0;  n < grid->n_lin_els;  ++n, ++lin, m+=3)
		(void) fprintf(file,"%d %d %d\n",m,m+1,m+2);
	}

	(void) fprintf(file,"END\n");
}		/*end print_tri_soln*/

LOCAL	void nearest_interior_state(
	double		*coords,
	COMPONENT	ext_comp,
	Front		*front,
	Wave		*wave,
	INTERFACE	*intfc,
	double		*coords_on,
	Locstate	state)
{
	COMPONENT	int_comp;
	HYPER_SURF_ELEMENT	*hse;
	HYPER_SURF	*hs;
	double		t[MAXD];

	if (nearest_interface_point(coords,ext_comp,intfc,INCLUDE_BOUNDARIES,
			            NULL,coords_on,t,&hse,&hs) != YES)
	{
	    screen("ERROR in nearest_interior_state(), "
	           "nearest_interface_point() failed\n");
	    clean_up(ERROR);
	}

	int_comp = (ext_comp == positive_component(hs)) ?
			negative_component(hs) : positive_component(hs);

	switch(wave_type(hs))
	{
	case SUBDOMAIN_BOUNDARY:
	    /*
	    *  You can't get states from parallel
	    *  or subdomain boundaries
	    */
	    hyp_solution(coords,int_comp,hs,
			 (ext_comp == positive_component(hs)) ?
	    	            NEGATIVE_SIDE : POSITIVE_SIDE,
	    	         front,wave,state,NULL);
	    return;
	default:
	    break;
	}

	state_along_hypersurface_element(int_comp,t,hse,hs,state);
}		/*end nearest_interior_state*/
#endif /* defined(TWOD) */

EXPORT	void	h_fprint_max_wave_speed_info(
	FILE	*file,
	Wave	*wave)
{
	int	i, j, dim = wave->rect_grid->dim;

	if (strcmp(wave->method,"ADVANCE_FRONTS_ONLY") == 0)
		return;		/* nothing to print out */

	(void) fprintf(file,"Maximum Wave Speed Information\n");
	(void) fprintf(file,"Maxsp = ");
	if (is_binary_output())
	{
	    (void) fprintf(file,"\f%c",dim);
	    (void) fwrite((const void *) Maxsp(wave),sizeof(double),dim,file);
	}
	else
	{
	    (void) fprintf(file,"%"FFMT,Maxsp(wave)[0]);
	    for (i = 1; i < dim; ++i)
	    	(void) fprintf(file,", %"FFMT,Maxsp(wave)[i]);
	}
	(void) fprintf(file,"\n");
	for (i = 0; i < dim; ++i)
	{
	    (void) fprintf(file,"MaxWaveSpeedCoords[%d] = ",i);
	    if (is_binary_output() == YES)
	    {
		(void) fprintf(file,"\f%c",dim);
	    	(void) fwrite((const void *) MaxWaveSpeedCoords(wave)[i],
			       sizeof(double),dim,file);
	    }
	    else
	    {
	    	(void) fprintf(file,"%"FFMT,MaxWaveSpeedCoords(wave)[i][0]);
		for (j = 1; j < dim; ++j)
		    (void) fprintf(file,", %"FFMT,
				   MaxWaveSpeedCoords(wave)[i][j]);
	    }
	    (void) fprintf(file,"\n");
	}


	for (i = 0; i < dim; ++i)
	{
	    (void) fprintf(file,"MaxWaveSpeedState[%d] = ",i);
	    fprint_state_data(file,MaxWaveSpeedState(wave)[i],
			      wave_tri_soln(wave)->intfc);
	}
}		/*end h_fprint_max_wave_speed_info*/

/* For parabolic step */
EXPORT	void	h_fprint_max_viscosity_info(
	FILE	*file,
	Wave	*wave)
{
	int	i, j, dim = wave->rect_grid->dim;

	if (strcmp(wave->method,"ADVANCE_FRONTS_ONLY") == 0)
		return;		/* nothing to print out */

	(void) fprintf(file,"Maximum Viscosity Information\n");
	(void) fprintf(file,"Maxvisc = ");
	if (is_binary_output() == YES)
	{
	    /* Jun 18 2003: Myoung-Nyoun: fixed */
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void *) & Maxvisc(wave),
			  sizeof(double),1,file);
	}
	else
	    (void) fprintf(file,"%FFMT",Maxvisc(wave));

	(void) fprintf(file,"\n");
	(void) fprintf(file,"MaxViscosityCoords = ");
	if (is_binary_output() == YES)
	{
	    /* Jun 18 2003: Myoung-Nyoun: fixed */
	    (void) fprintf(file,"\f%c",dim);
	    (void) fwrite((const void *) MaxViscosityCoords(wave),
			  sizeof(double),dim,file);
	}
	else
	{
	    (void) fprintf(file,"%g",MaxViscosityCoords(wave)[0]);
	    for (j = 1; j < dim; j++)
	      (void) fprintf(file,", %g",MaxViscosityCoords(wave)[j]);
	}
	(void) fprintf(file,"\n");

	(void) fprintf(file,"MaxViscosityState = ");
	fprint_state_data(file,MaxViscosityState(wave),
			  wave_tri_soln(wave)->intfc);

}		/*end h_fprint_max_viscosity_info*/

/*
*			print_Wave_structure():
*
*	Prints the contents of a Wave structure.
*/

EXPORT void print_Wave_structure(
	Wave		*wave)
{
	(void) printf("\n\n\n\t\tWave %p structure\n",(POINTER)wave);
	if (wave == NULL)
	{
	    (void) printf("\t\tstructure not yet allocated\n"
	                  "\n\t\tEnd Wave %p structure\n\n",(POINTER)wave);
	    return;
	}

	print_RECT_GRID_structure(wave->rect_grid);

	(void) printf("\nsize of locstate: sizest %d (bytes)\n",
		      (int)wave->sizest);
	(void) printf("number floats in locstate: nfloats %d\n",wave->nfloats);
	(void) printf("print_state() %p show_wave_states() %p\n",
	       wave->print_state,wave->show_wave_states);

	(void) printf("\nHyperbolic solution method: %s\n",wave->method);

	(void) printf("\nnpt_solver() %p\n",wave->_npt_solver);
	print_interpolator_structure(&wave->interpolator);
	print_el_integral_structure(&wave->el_integral);
	print_unsplit_structure(&wave->unsplit);
	(void) printf("max_wave_speed %p\n",wave->max_wave_speed);

	print_pt_sources(wave);

	(void) printf("\n");
	print_max_wave_speed_info(stdout,wave);
	(void) printf("\n");
	(void) printf("tri_soln %p\n",wave_tri_soln(wave));
	(void) printf("areas %p  min_comp %d\n",
	       wave_areas(wave),wave_min_comp(wave));


	(void) printf("\n\t\tEnd Wave %p structure\n\n",wave);
}			/*end print_Wave_structure*/



EXPORT	void print_pt_sources(
	Wave		*wave)
{
	int		i, j, dim = wave->rect_grid->dim;

	(void) printf("\n\t\tnum_point_sources %d\n",wave->num_point_sources);
	if (wave->num_point_sources == 0) return;

	(void) printf("type\txpoint\typoint\tstrength  composition\n");
	for (i = 0; i < wave->num_point_sources; ++i)
	{
	    (void) printf("%13s ",(wave->source_type[i] == SOURCE) ? "SOURCE"
	    		: ((wave->source_type[i] == SINK) ? "SINK"
	    			: "NOT_SPECIFIED"));
	    for (j = 0; j < dim; ++j)
	    	(void) printf("%.5f ",wave->srcpt[j][i]);
	    for (j = 0; j < dim; ++j)
	    	(void) printf("%.5f ",wave->pt_src_diam[j][i]);
	    (void) printf("%.5f     ",wave->strength[i]);

	    if (wave->composition != NULL)
	    	(*wave->print_state)(wave->composition[i]);
	    else
	    	(void) printf("NULL\n");
	    (void) printf("\n");
	}
}			/*end print_pt_sources*/

EXPORT	void	h_print_H_Front_structure(
	Front		*fr)
{
	(void) printf("\n\n\n\t\tH_Front %p structure\n",(POINTER)fr);
	f_print_Front_structure(fr);
	(void) printf("wave_of_front(fr) = %p\n",(POINTER)wave_of_front(fr));
	(void) printf("\n\n\n\t\tEnd H_Front %p structure\n",(POINTER)fr);
}		/*end h_print_H_Front_structure*/

EXPORT	void print_Stencil(
	Stencil		*sten)
{
	int		i, imax, dim;
#if defined(TWOD)
	int		j;
#endif /* defined(TWOD) */

	if (sten == NULL)
	{
	    (void) printf("Stencil NULL\n");
	    return;
	}
	(void) printf("reg_stencil %s, prev_reg_stencil %s\n",
		      y_or_n(sten->reg_stencil),
		      y_or_n(sten->prev_reg_stencil));
	(void) printf("npts = %d\n",sten->npts);
	dim = sten->fr->rect_grid->dim;
	imax = (sten->npts % 2) ? sten->npts/2 + 1 : sten->npts/2;
	for (i = -sten->npts/2; i < imax; ++i)
	{
	    (void) printf("sten->hs[%d] = %llu, ",i,
	                  hypersurface_number(sten->hs[i]));/*TODO REMOVE*/
	    (void) printf("sten->crx[%d] = %p, ",i,(POINTER)sten->crx[i]);
	    (void) printf("sten->nc[%d] = %d,\n",i,sten->nc[i]);
#if defined(TWOD)
	    /* TODO Generalize CRXING to multiple dimensions */
	    for (j = 0; j < sten->nc[i]; ++j)
	    {
	    	(void) printf("\t\tsten->crx[%d][%d] = %p,\n",
	    		      i,j,(POINTER)sten->crx[i][j]);
	    	print_crxings(sten->crx[i][j],NO);
	    }
#endif /* defined(TWOD) */
	    (void) printf("\tsten->p[%d] = %llu",i,point_number(sten->p[i]));
	    if (sten->p[i] != NULL)
	    {
	    	(void) printf(", ");
	    	print_general_vector(NULL,Coords(sten->p[i]),dim,"");
	    }
	    (void) printf("\n");
	}


	(void) printf("fr %p newfr %p wave %p newwave %p\n",
		      (POINTER)sten->fr,(POINTER)sten->newfr,
		      (POINTER)sten->wave,(POINTER)sten->newwave);
	(void) printf("\n");
}			/*end print_Stencil*/

LOCAL	void	print_interpolator_structure(
	INTERPOLATORS *interpolator)
{
	(void) printf("INTERPOLATORS\n");
	if (interpolator == NULL)
	{
		(void) printf("\tNONE\n");
		return;
	}

	(void) printf("\tlinear_cell() %p\tbilinear_cell() %p\n",
		interpolator->linear_cell,interpolator->bilinear_cell);
	(void) printf("\tgrad_linear_cell() %p\tgrad_bilinear_cell() %p\n",
		interpolator->grad_linear_cell,
		interpolator->grad_bilinear_cell);
	(void) printf("\ngrad_bond() %p\n",interpolator->grad_bond);
}		/*end print_interpolator_structure*/

LOCAL	void print_el_integral_structure(
	EL_INTEGRALS *el_integral)
{
	(void) printf("ELEMENT INTEGRALS\n");
	if (el_integral == NULL)
	{
		(void) printf("\tNONE\n");
		return;
	}
	(void) printf("\tlinear_cell() %p\tbilinear_cell() %p\n",
		el_integral->linear_cell,el_integral->bilinear_cell);
}		/*end print_el_integral_structure*/

LOCAL	void print_unsplit_structure(
	UNSPLIT *unsplit)
{
	(void) printf("UNSPLIT OPERATORS\n");
	if (unsplit == NULL)
	{
		(void) printf("\tNONE\n");
		return;
	}
	(void) printf("\tflux() %p\tflux_obl() %p",
	       unsplit->flux,unsplit->flux_obl);
	(void) printf("\tsources() %p\n",unsplit->sources);
}		/*end print_unsplit_structure*/
