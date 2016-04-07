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
*				gsesinout.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#if defined(SESAME_CODE) && defined(TWOD)

#define DEBUG_STRING    "ses_hyp"

#include <geos/sesame.h>

	/* LOCAL Function Declarations */
LOCAL	OUTPUT_VALUE	*ses_ps_prt_soln(OUTPUT_SOLN*,double*,int*);
LOCAL	OUTPUT_VALUE	*ses_re_prt_soln(OUTPUT_SOLN*,double*,int*);
LOCAL	OUTPUT_VALUE	*ses_rs_prt_soln(OUTPUT_SOLN*,double*,int*);
LOCAL	OUTPUT_VALUE	*ses_rt_prt_soln(OUTPUT_SOLN*,double*,int*);
LOCAL	void	init_ses_restart_interior_states(SESAME_TABLE_TYPE,Wave*,Front*,
						 void (*)(int*,COMPONENT,
							 INPUT_SOLN**,Locstate),
						 INPUT_SOLN**);
LOCAL	void	init_ses_states(SESAME_TABLE_TYPE,Wave*,Front*,Printplot*,
				void (*)(int*,COMPONENT,INPUT_SOLN**,Locstate),
				const IO_TYPE*,INTERFACE*);
LOCAL	void	read_print_ses_table(INIT_DATA*,SESAME_TABLE_TYPE,Front*,Wave*,
				     Printplot*,const char*,const IO_TYPE*,
				     void (*)(int*,COMPONENT,
					      INPUT_SOLN**,Locstate));
LOCAL	void	ses_ps_prt_intfc_soln(OUTPUT_SOLN*,POINT*,HYPER_SURF_ELEMENT*,
				      HYPER_SURF*,OUTPUT_VALUE*,OUTPUT_VALUE*);
LOCAL	void	ses_ps_restart_initializer(int*,COMPONENT,
					   INPUT_SOLN**,Locstate);
LOCAL	void	ses_ps_set_restart_intfc_states(double*,double*,int,POINT*,
					       HYPER_SURF_ELEMENT*,HYPER_SURF*);
LOCAL	void	ses_re_prt_intfc_soln(OUTPUT_SOLN*,POINT*,HYPER_SURF_ELEMENT*,
				      HYPER_SURF*,OUTPUT_VALUE*,OUTPUT_VALUE*);
LOCAL	void	ses_re_restart_initializer(int*,COMPONENT,
					   INPUT_SOLN**,Locstate);
LOCAL	void	ses_re_set_restart_intfc_states(double*,double*,int,POINT*,
					      HYPER_SURF_ELEMENT*,HYPER_SURF*);
LOCAL	void	ses_rs_prt_intfc_soln(OUTPUT_SOLN*,POINT*,HYPER_SURF_ELEMENT*,
				      HYPER_SURF*,OUTPUT_VALUE*,OUTPUT_VALUE*);
LOCAL	void	ses_rs_restart_initializer(int*,COMPONENT,
					   INPUT_SOLN**,Locstate);
LOCAL	void	ses_rs_set_restart_intfc_states(double*,double*,int,POINT*,
					       HYPER_SURF_ELEMENT*,HYPER_SURF*);
LOCAL	void	ses_rt_prt_intfc_soln(OUTPUT_SOLN*,POINT*,HYPER_SURF_ELEMENT*,
				      HYPER_SURF*,OUTPUT_VALUE*,OUTPUT_VALUE*);
LOCAL	void	ses_rt_restart_initializer(int*,COMPONENT,
					   INPUT_SOLN**,Locstate);
LOCAL	void	ses_rt_set_restart_intfc_states(double*,double*,int,POINT*,
					       HYPER_SURF_ELEMENT*,HYPER_SURF*);

/*
*			init_ses_prt():
*
*/

EXPORT void init_ses_prt(
	Printplot	*prt,
	Front		*fr,
	int		sesame_table)
{
	RECT_GRID	*rect_grid = fr->rect_grid;
	int		nfloats = fr->nfloats;
	int		var;
	void		(*set_restart_intfc_states)(double*,double*,int,POINT*,
						    HYPER_SURF_ELEMENT*,
						    HYPER_SURF*);
	static boolean	first = YES;

	if (first == YES)
	{
	    first = NO;
	    bi_array(&prt->output_soln,nfloats,1,sizeof(OUTPUT_SOLN));
	    bi_array(&prt->restart_soln,nfloats,1,sizeof(INPUT_SOLN));
	}
	prt->printout = ses_printout;
	prt->print_initial_data = NULL;
	prt->initialize_for_printout = NULL;
	prt->n_rect_state_vars = nfloats;

	switch (sesame_table)
	{
	case SESAME_RHO_TEMP:
	    prt->output_soln[0]->name = "COLD_PRESSURE";
	    prt->output_soln[1]->name = "COLD_INT_ENERGY";
	    prt->output_soln[2]->name = "REDUCED_PRESSURE";
	    prt->output_soln[3]->name = "REDUCED_INT_ENERGY";
	    prt->output_soln[4]->name = "ENTROPY";
	    prt->output_soln[5]->name = "ADIABATIC_EXPONENT";
	    prt->output_soln[6]->name = "GRUNEISEN_COEFF";
	    prt->output_soln[7]->name = "INTDPOCR";
	    set_restart_intfc_states = ses_rt_set_restart_intfc_states;
	    break;

	case SESAME_RHO_ENERGY:
	    prt->output_soln[0]->name = "COLD_PRESSURE";
	    prt->output_soln[1]->name = "REDUCED_PRESSURE";
	    prt->output_soln[2]->name = "LOG_TEMPERATURE";
	    prt->output_soln[3]->name = "ENTROPY";
	    prt->output_soln[4]->name = "ADIABATIC_EXPONENT";
	    prt->output_soln[5]->name = "GRUNEISEN_COEFF";
	    prt->output_soln[6]->name = "INTDPOCR";
	    set_restart_intfc_states = ses_re_set_restart_intfc_states;
	    break;

	case SESAME_RHO_ENTROPY:
	    prt->output_soln[0]->name = "COLD_PRESSURE";
	    prt->output_soln[1]->name = "COLD_INT_ENERGY";
	    prt->output_soln[2]->name = "REDUCED_PRESSURE";
	    prt->output_soln[3]->name = "REDUCED_INT_ENERGY";
	    prt->output_soln[4]->name = "LOG_TEMPERATURE";
	    prt->output_soln[5]->name = "ADIABATIC_EXPONENT";
	    prt->output_soln[6]->name = "GRUNEISEN_COEFF";
	    prt->output_soln[7]->name = "INTDPOCR";
	    set_restart_intfc_states = ses_rs_set_restart_intfc_states;
	    break;

	case SESAME_PRESS_ENTROPY:
	    prt->output_soln[0]->name = "COLD_INT_ENERGY";
	    prt->output_soln[1]->name = "REDUCED_INT_ENERGY";
	    prt->output_soln[2]->name = "LOG_TEMPERATURE";
	    prt->output_soln[3]->name = "LOG_DENSITY";
	    prt->output_soln[4]->name = "ADIABATIC_EXPONENT";
	    prt->output_soln[5]->name = "GRUNEISEN_COEFF";
	    prt->output_soln[6]->name = "INTDPOCR";
	    set_restart_intfc_states = ses_ps_set_restart_intfc_states;
	    break;
	}

	for (var = 0; var < prt->n_rect_state_vars; var++)
	{
	    prt->output_soln[var]->fit = LINEAR_FIT;
	    prt->output_soln[var]->smoothness = SINGULAR;
	    prt->output_soln[var]->grid = rect_grid;
	    prt->output_soln[var]->intfc = fr->interf;
	    prt->output_soln[var]->var = var;
	    prt->output_soln[var]->extra = NULL;
	    switch (sesame_table)
	    {
	    case SESAME_RHO_TEMP:
	    	prt->output_soln[var]->solution = ses_rt_prt_soln;
	    	prt->output_soln[var]->intfc_solution = ses_rt_prt_intfc_soln;
		break;
	
	    case SESAME_RHO_ENERGY:
	    	prt->output_soln[var]->solution = ses_re_prt_soln;
	    	prt->output_soln[var]->intfc_solution = ses_re_prt_intfc_soln;
		break;
	
	    case SESAME_RHO_ENTROPY:
		prt->output_soln[var]->solution = ses_rs_prt_soln;
		prt->output_soln[var]->intfc_solution = ses_rs_prt_intfc_soln;
		break;

	    case SESAME_PRESS_ENTROPY:
	    	prt->output_soln[var]->solution = ses_ps_prt_soln;
	    	prt->output_soln[var]->intfc_solution = ses_ps_prt_intfc_soln;
		break;
	    }
	}

	prt->n_restart_vars = nfloats;

	for (var = 0; var < prt->n_restart_vars; var++) 
	{
	    prt->restart_soln[var]->name = prt->output_soln[var]->name;
	    prt->restart_soln[var]->fit = LINEAR_FIT;
 	    prt->restart_soln[var]->smoothness = SINGULAR;
	    prt->restart_soln[var]->set_intfc_states = set_restart_intfc_states;
	}
	prt->n_tri_vars = 0;
	prt->tri_plot_name = NULL;
	prt->tri_plot_function = NULL;
}		/*end init_ses_prt*/

/*
*
*	Routines for printout.  
*/

/*ARGSUSED*/
LOCAL OUTPUT_VALUE *ses_rt_prt_soln(
	OUTPUT_SOLN	*os,
	double		*coords,
	int		*icoords)
{
	static OUTPUT_VALUE Sol;
	Wave		*wave = (Wave *) ((POINTER*) os->extra)[1];
	Locstate	state = Rect_state(icoords,wave);
	
	Sol.utype = Float;
	switch (os->var)
	{
	case 0:
		Sol.uval.fval = ses_rt_coldp(state);
		break;
	case 1:
		Sol.uval.fval = ses_rt_colde(state);
		break;
	case 2:
		Sol.uval.fval = ses_rt_redp(state);
		break;
	case 3:
		Sol.uval.fval = ses_rt_rede(state);
		break;
	case 4:
		Sol.uval.fval = ses_rt_S(state);
		break;
	case 5:
		Sol.uval.fval = ses_rt_adb_gam(state);
		break;
	case 6:
		Sol.uval.fval = ses_rt_gru_gam(state);
		break;
	case 7:
		Sol.uval.fval = ses_rt_riv(state);
		break;
	default:
		screen("ERROR: unknown var in ses_rt_prt_soln()\n");
		Sol.uval.fval = ERROR_FLOAT;
		clean_up(ERROR);
	}
	return &Sol;
}		/*end ses_rt_prt_soln*/


/*ARGSUSED*/
LOCAL OUTPUT_VALUE *ses_re_prt_soln(
	OUTPUT_SOLN 	*os,
	double		*coords,
	int		*icoords)
{
	static OUTPUT_VALUE Sol;
	Wave		*wave = (Wave *) ((POINTER*) os->extra)[1];
	Locstate	state = Rect_state(icoords,wave);
	
	Sol.utype = Float;
	switch (os->var)
	{
	case 0:
		Sol.uval.fval =  ses_re_coldp(state);
		break;
	case 1:
		Sol.uval.fval = ses_re_redp(state);
		break;
	case 2:
		Sol.uval.fval = ses_re_Tvar(state);
		break;
	case 3:
		Sol.uval.fval = ses_re_S(state);
		break;
	case 4:
		Sol.uval.fval = ses_re_adb_gam(state);
		break;
	case 5:
		Sol.uval.fval = ses_re_gru_gam(state);
		break;
	case 6:
		Sol.uval.fval = ses_re_riv(state);
		break;
	default:
		screen("ERROR: unknown var in ses_re_prt_soln()\n");
		Sol.uval.fval = ERROR_FLOAT;
		clean_up(ERROR);
	}
	return &Sol;
}		/*end ses_re_prt_soln*/

/*ARGSUSED*/
LOCAL OUTPUT_VALUE *ses_rs_prt_soln(
	OUTPUT_SOLN	*os,
	double		*coords,
	int		*icoords)
{
	static OUTPUT_VALUE Sol;
	Wave		*wave = (Wave *) ((POINTER*) os->extra)[1];
	Locstate	state = Rect_state(icoords,wave);
	
	Sol.utype = Float;
	switch (os->var)
	{
	case 0:
		Sol.uval.fval = ses_rs_coldp(state);
		break;
	case 1:
		Sol.uval.fval = ses_rs_colde(state);
		break;
	case 2:
		Sol.uval.fval = ses_rs_redp(state);
		break;
	case 3:
		Sol.uval.fval = ses_rs_rede(state);
		break;
	case 4:
		Sol.uval.fval = ses_rs_Tvar(state);
		break;
	case 5:
		Sol.uval.fval = ses_rs_adb_gam(state);
		break;
	case 6:
		Sol.uval.fval = ses_rs_gru_gam(state);
		break;
	case 7:
		Sol.uval.fval = ses_rs_riv(state);
		break;
	default:
		screen("ERROR: unknown var in ses_rs_prt_soln()\n");
		Sol.uval.fval = ERROR_FLOAT;
		clean_up(ERROR);
	}
	return &Sol;
}		/*end ses_rs_prt_soln*/


/*ARGSUSED*/
LOCAL OUTPUT_VALUE *ses_ps_prt_soln(
	OUTPUT_SOLN	*os,
	double		*coords,
	int		*icoords)
{
	static OUTPUT_VALUE Sol;
	Wave		*wave = (Wave *) ((POINTER*) os->extra)[1];
	Locstate	state = Rect_state(icoords,wave);
	
	Sol.utype = Float;
	switch (os->var)
	{
	case 0:
		Sol.uval.fval = ses_ps_colde(state);
		break;
	case 1:
		Sol.uval.fval = ses_ps_rede(state);
		break;
	case 2:
		Sol.uval.fval = ses_ps_Tvar(state);
		break;
	case 3:
		Sol.uval.fval = ses_ps_rho_var(state);
		break;
	case 4:
		Sol.uval.fval = ses_ps_adb_gam(state);
		break;
	case 5:
		Sol.uval.fval = ses_ps_gru_gam(state);
		break;
	case 6:
		Sol.uval.fval = ses_ps_riv(state);
		break;
	default:
		screen("ERROR: unknown var in ses_ps_prt_soln()\n");
		Sol.uval.fval = ERROR_FLOAT;
		clean_up(ERROR);
	}
	return &Sol;
}		/*end ses_ps_prt_soln*/

/* ARGSUSED */
LOCAL void ses_rt_prt_intfc_soln(
	OUTPUT_SOLN		*os,
	POINT			*p,
	HYPER_SURF_ELEMENT 	*hse,
	HYPER_SURF		*hs,
	OUTPUT_VALUE		*left,
	OUTPUT_VALUE		*right)
{
	Locstate		stl, str;
	BOND			*b = Bond_of_hse(hse);
	CURVE			*c = Curve_of_hs(hs);

	stl = left_state_at_point_on_curve(p,b,c);
	str = right_state_at_point_on_curve(p,b,c);

	left->utype = right->utype = Float;
	switch (os->var)
	{
	case 0:
		left->uval.fval = ses_rt_coldp(stl);
		right->uval.fval = ses_rt_coldp(str);
		break;
	case 1:
		left->uval.fval = ses_rt_colde(stl);
		right->uval.fval = ses_rt_colde(str);
		break;
	case 2:
		left->uval.fval = ses_rt_redp(stl);
		right->uval.fval = ses_rt_redp(str);
		break;
	case 3:
		left->uval.fval = ses_rt_rede(stl);
		right->uval.fval = ses_rt_rede(str);
		break;
	case 4:
		left->uval.fval = ses_rt_S(stl);
		right->uval.fval = ses_rt_S(str);
		break;
	case 5:
		left->uval.fval = ses_rt_adb_gam(stl);
		right->uval.fval = ses_rt_adb_gam(str);
		break;
	case 6:
		left->uval.fval = ses_rt_gru_gam(stl);
		right->uval.fval = ses_rt_gru_gam(str);
		break;
	case 7:
		left->uval.fval = ses_rt_riv(stl);
		right->uval.fval = ses_rt_riv(str);
		break;
	default:
		screen("ERROR: unknown var in ses_rt_prt_intfc_soln()\n");
		clean_up(ERROR);
	}
}		/*end ses_rt_prt_intfc_soln*/


/*ARGSUSED*/
LOCAL void ses_re_prt_intfc_soln(
	OUTPUT_SOLN		*os,
	POINT			*p,
	HYPER_SURF_ELEMENT	*hse,
	HYPER_SURF		*hs,
	OUTPUT_VALUE		*left,
	OUTPUT_VALUE		*right)
{
	BOND			*b = Bond_of_hse(hse);
	CURVE			*c = Curve_of_hs(hs);
	Locstate		stl, str;

	stl = left_state_at_point_on_curve(p,b,c);
	str = right_state_at_point_on_curve(p,b,c);

	left->utype = right->utype = Float;
	switch (os->var)
	{
	case 0:
		left->uval.fval = ses_re_coldp(stl);
		right->uval.fval = ses_re_coldp(str);
		break;
	case 1:
		left->uval.fval = ses_re_redp(stl);
		right->uval.fval = ses_re_redp(str);
		break;
	case 2:
		left->uval.fval = ses_re_Tvar(stl);
		right->uval.fval = ses_re_Tvar(str);
		break;
	case 3:
		left->uval.fval = ses_re_S(stl);
		right->uval.fval = ses_re_S(str);
		break;
	case 4:
		left->uval.fval = ses_re_adb_gam(stl);
		right->uval.fval = ses_re_adb_gam(str);
		break;
	case 5:
		left->uval.fval = ses_re_gru_gam(stl);
		right->uval.fval = ses_re_gru_gam(str);
		break;
	case 6:
		left->uval.fval = ses_re_riv(stl);
		right->uval.fval = ses_re_riv(str);
		break;
	default:
		screen("ERROR: unknown var in ses_re_prt_intfc_soln()\n");
		clean_up(ERROR);
	}
}		/*end ses_re_prt_intfc_soln*/

/*ARGSUSED*/
LOCAL void ses_rs_prt_intfc_soln(
	OUTPUT_SOLN		*os,
	POINT			*p,
	HYPER_SURF_ELEMENT	*hse,
	HYPER_SURF		*hs,
	OUTPUT_VALUE		*left,
	OUTPUT_VALUE		*right)
{
	BOND			*b = Bond_of_hse(hse);
	CURVE			*c = Curve_of_hs(hs);
	Locstate		stl, str;

	stl = left_state_at_point_on_curve(p,b,c);
	str = right_state_at_point_on_curve(p,b,c);

	left->utype = right->utype = Float;
	switch (os->var)
	{
	case 0:
		left->uval.fval = ses_rs_coldp(stl);
		right->uval.fval = ses_rs_coldp(str);
		break;
	case 1:
		left->uval.fval = ses_rs_colde(stl);
		right->uval.fval = ses_rs_colde(str);
		break;
	case 2:
		left->uval.fval = ses_rs_redp(stl);
		right->uval.fval = ses_rs_redp(str);
		break;
	case 3:
		left->uval.fval = ses_rs_rede(stl);
		right->uval.fval = ses_rs_rede(str);
		break;
	case 4:
		left->uval.fval = ses_rs_Tvar(stl);
		right->uval.fval = ses_rs_Tvar(str);
		break;
	case 5:
		left->uval.fval = ses_rs_adb_gam(stl);
		right->uval.fval = ses_rs_adb_gam(str);
		break;
	case 6:
		left->uval.fval = ses_rs_gru_gam(stl);
		right->uval.fval = ses_rs_gru_gam(str);
		break;
	case 7:
		left->uval.fval = ses_rs_riv(stl);
		right->uval.fval = ses_rs_riv(str);
		break;
	default:
		screen("ERROR: unknown var in ses_rs_prt_intfc_soln()\n");
		clean_up(ERROR);
	}
}		/*end ses_rs_prt_intfc_soln*/

/*ARGSUSED*/
LOCAL void ses_ps_prt_intfc_soln(
	OUTPUT_SOLN		*os,
	POINT			*p,
	HYPER_SURF_ELEMENT	*hse,
	HYPER_SURF		*hs,
	OUTPUT_VALUE		*left,
	OUTPUT_VALUE		*right)
{
	BOND			*b = Bond_of_hse(hse);
	CURVE			*c = Curve_of_hs(hs);
	Locstate		stl, str;

	stl = left_state_at_point_on_curve(p,b,c);
	str = right_state_at_point_on_curve(p,b,c);

	left->utype = right->utype = Float;
	switch (os->var)
	{
	case 0:
		left->uval.fval = ses_ps_colde(stl);
		right->uval.fval = ses_ps_colde(str);
		break;
	case 1:
		left->uval.fval = ses_ps_rede(stl);
		right->uval.fval = ses_ps_rede(str);
		break;
	case 2:
		left->uval.fval = ses_ps_Tvar(stl);
		right->uval.fval = ses_ps_Tvar(str);
		break;
	case 3:
		left->uval.fval = ses_ps_rho_var(stl);
		right->uval.fval = ses_ps_rho_var(str);
		break;
	case 4:
		left->uval.fval = ses_ps_adb_gam(stl);
		right->uval.fval = ses_ps_adb_gam(str);
		break;
	case 5:
		left->uval.fval = ses_ps_gru_gam(stl);
		right->uval.fval = ses_ps_gru_gam(str);
		break;
	case 6:
		left->uval.fval = ses_ps_riv(stl);
		right->uval.fval = ses_ps_riv(str);
		break;
	default:
		screen("ERROR: unknown var in ses_ps_prt_intfc_soln()\n");
		clean_up(ERROR);
	}
}		/*end ses_ps_prt_intfc_soln*/

/*
*			ses_printout():
*
*	Provides printout of rect_grid, front->interf, and of all state 
*	variables associated with wave.
*/

/*ARGSUSED*/
EXPORT void ses_printout(
	CHART		*root,
	Printplot	*prt,
	boolean		unused1,
	int		unused2)
{
	Front		*front = root->front;
	Wave		*wave = root->wave;
	RECT_GRID	*grid = front->rect_grid;
	FILE		*file = prt->file;
	int		var;
	POINTER		extra[2];

	if (file == NULL) file = stdout;
	fprint_rectangular_grid(file,grid);

			/* Initialize */

	extra[0] = (POINTER) front;
	extra[1] = (POINTER) wave;
	for (var = 0; var < prt->n_rect_state_vars; var++)
	{
	    prt->output_soln[var]->intfc = front->interf;
	    prt->output_soln[var]->extra = (POINTER) extra;
	}

			/* Print Front */

	(void) fprintf(file,"\n\n\n\n");
	(void) foutput(file);
	(void) fprintf(file,"\t\t\tFRONT DATA:\n");
	fprint_interface(file,front->interf);
	
	(void) fprintf(file,"\n\n");
	(void) foutput(file);

	(void) fprintf(file,"\n\n\n\n");
	(void) foutput(file);
	(void) fprintf(file,"\t\t\tSTATE DATA:\n");

			/* Printout Components of Locstates in Wave */

	print_states(file,prt,front->rect_grid->dim);


	(void) fprintf(file,"\n\n\n\n");
	(void) foutput(file);
	(void) fprintf(file,"\t\t\tEND OF STATE DATA\n");
}		/*end ses_printout*/

EXPORT boolean restart_sesame(
	INIT_DATA     *init,
	SESAME_EOS    *seos)
{
	const IO_TYPE     *io_type = seos->restart_io_type;
	FILE		  *file = io_type->file;
	Printplot	  *prt;
	void		  (*restart_initializer)(int*,COMPONENT,
					         INPUT_SOLN**,Locstate);
	int		  sizest, i;
	CHART		  **root = seos->root;
	Front		  **fr = seos->fr;
	Wave		  **wave = seos->wave;
	SESAME_TABLE_TYPE tbl;
	int               ntables = 4;
	const SESAME_TABLE_TYPE *tables = sesame_table_numbers();
	const char **ses_header = sesame_headers();

	if (!next_output_line_containing_string(file,
		                                "SESAME RESTART INFORMATION"))
	    return FUNCTION_FAILED;

	g_preserve_user_hooks(2,SAVE_HOOKS);/*no return before matching call */
	set_user_hooks_for_sesame();
	determine_read_version(file);
	read_print_SESAME_params(seos,io_type);
	(void) fgetstring(file,"Number of tables");
	(void) fscanf(file,"%*s %d",&ntables);
	scalar(&prt,sizeof(Printplot));
	for(i = 0; i < ntables; i++)
	{
	    tbl = tables[i];
	    scalar(&root[tbl],sizeof(CHART));
	    scalar(&fr[tbl],sizeof(Ses_Front));
	    scalar(&wave[tbl],sizeof(Ses_Wave));
	    Ses_front_seos(fr[tbl]) = Ses_wave_seos(wave[tbl]) = seos;
	    root[tbl]->front = fr[tbl];
	    root[tbl]->wave = wave[tbl];
	    scalar(&fr[tbl]->rect_grid,sizeof(RECT_GRID));
	    wave[tbl]->rect_grid = fr[tbl]->rect_grid;
	    switch (tbl)
	    {
	    case SESAME_RHO_TEMP:
	    	sizest = sizeof(SES_RT_STATE);
	    	restart_initializer = ses_rt_restart_initializer;
	    	break;
	    case SESAME_RHO_ENERGY:
	    	sizest = sizeof(SES_RE_STATE);
	    	restart_initializer = ses_re_restart_initializer;
	    	break;
	    case SESAME_RHO_ENTROPY:
	    	sizest = sizeof(SES_RS_STATE);
	    	restart_initializer = ses_rs_restart_initializer;
	    	break;
	    case SESAME_PRESS_ENTROPY:
	    	sizest = sizeof(SES_PS_STATE);
	    	restart_initializer = ses_ps_restart_initializer;
	    	break;
	    }
	    set_default_ses_wave_and_front(init,fr[tbl],wave[tbl],sizest,NO);
	    read_print_ses_table(init,tbl,fr[tbl],wave[tbl],prt,
			         ses_header[tbl],io_type,
				 restart_initializer);
	    if (ntables == 1)
	    {
	    	set_ses_inv_hyp_solns(seos);
	    	break;
	    }
	}
	g_preserve_user_hooks(2,RESTORE_HOOKS);
	return FUNCTION_SUCCEEDED;
}		/*end restart_sesame*/

LOCAL	void read_print_ses_table(
	INIT_DATA         *init,
	SESAME_TABLE_TYPE table,
	Front		  *fr,
	Wave		  *wave,
	Printplot	  *prt,
	const char	  *header,
	const IO_TYPE     *io_type,
	void		  (*restart_initializer)(int*,COMPONENT,
					         INPUT_SOLN**,Locstate))
{
	FILE		 *file = io_type->file;
	F_USER_INTERFACE *fuh;
	INTERFACE	 *restart_intfc;
	REMAP            Remap;
	int		 grids_set;

	if (next_output_line_containing_string(file,header) == NULL)
	{
	    screen("ERROR in read_print_ses_table(), "
	           "unable to find header \"%s\" in ses_rstrt_file\n",
		   header);
	    clean_up(ERROR);
	}
	set_remap(2,IDENTITY_REMAP,&Remap);
	read_rectangular_grid(io_type,fr->rect_grid,YES,&Remap);
	set_size_of_intfc_state(fr->sizest);
	restart_intfc = read_print_interface(init,io_type,NO,&grids_set);
	if (restart_intfc == NULL)
	{
	    screen("ERROR in read_print_ses_table(), "
	           "Can't read interface for table %d, header \"%s\"\n",
		   table,header);
	    clean_up(ERROR);
	}
	if (grids_set == NO)
	{
	    set_topological_grid(restart_intfc,fr->rect_grid);
	    (void) adjust_top_grid_for_square(&topological_grid(restart_intfc),
				              fr->rect_grid);
	    set_computational_grid(restart_intfc,fr->rect_grid);
	}
	fr->interf = restart_intfc;
	fuh = &f_user_interface(fr->interf);
	fuh->_bi_interpolate_intfc_states = fr->_state_interpolator;
	fuh->_tri_interpolate_intfc_states = fr->_tri_state_interpolator;
	init_ses_prt(prt,fr,table);
	init_ses_states(table,wave,fr,prt,restart_initializer,
			io_type,restart_intfc);
}		/*end read_print_ses_table*/


/*
*			init_ses_states():
*
*	Initializes state pointers and values in the tracked case.
*/

LOCAL void init_ses_states(
	SESAME_TABLE_TYPE table,
	Wave		  *wave,
	Front		  *front,
	Printplot	  *prt,
	void		  (*restart_initializer)(int*,COMPONENT,
					         INPUT_SOLN**,Locstate),
	const IO_TYPE     *io_type,
	INTERFACE	  *restart_intfc)
{
	int  var;

	DEBUG_ENTER(init_ses_states)

	if (!read_state_variables(io_type,prt->n_restart_vars,
			          restart_intfc,prt->restart_soln,2)) 
	{
	    screen("ERROR in init_ses_states(), "
	           "read_state_variables failed\n");
	    clean_up(ERROR);
	}

	init_ses_restart_interior_states(table,wave,front,restart_initializer,
				         prt->restart_soln);


	for (var = 0; var < prt->n_restart_vars; var++)
	    free_input_soln(prt->restart_soln[var]);


	DEBUG_LEAVE(init_ses_states)
}		/*end init_ses_states*/


/*
*			init_ses_restart_interior_states():
*
*	Initializes the states in a wave structure by calling
*
*		(*restart_initializer)(icoords,comp,restart_soln,state)
*
*	at the centers of the grid blocks of wave->rect_grid.
*/


LOCAL void init_ses_restart_interior_states(
	SESAME_TABLE_TYPE table,
	Wave		  *wave,
	Front		  *front,
	void		  (*restart_initializer)(int*,COMPONENT,
					         INPUT_SOLN**,Locstate),
	INPUT_SOLN	  **restart_soln)
{
	INTERFACE	*intfc;
	int		ix, iy;
	int		icoords[MAXD];
	int		xmax = wave->rect_grid->gmax[0];
	int		ymax = wave->rect_grid->gmax[1];
	size_t		sizest = front->sizest;
	Locstate	state;
	COMPONENT	comp;
	int		status;

	DEBUG_ENTER(init_ses_restart_interior_states)
	if ((wave->sizest == 0) || (restart_initializer == NULL))
	{
		DEBUG_LEAVE(init_ses_restart_interior_states)
		return;
	}

	set_ses_intrp_flag_all(table);
	status = init_ses_hyp_soln_func(wave,front);
	if (status != GOOD_STEP) 
	{
		screen("ERROR: init_ses_hyp_soln_func() failed\n");
		clean_up(ERROR);
	}
	intfc = front->interf;

	for (iy = 0; iy < ymax; iy++)
	{
		icoords[1] = iy;
		for (ix = 0; ix < xmax; ix++) 
		{
			icoords[0] = ix;
			comp = Rect_comp(icoords,wave);
			state = Rect_state(icoords,wave);
			if (is_exterior_comp(comp,intfc))
			{
				clear_state(front->interf,state,sizest);
			}
			else
				(*restart_initializer)(icoords,comp,
						       restart_soln,state);
		}
	}
	DEBUG_LEAVE(init_ses_restart_interior_states)
}		/*end init_ses_restart_interior_states*/

/*ARGSUSED*/
LOCAL void ses_rt_restart_initializer(
	int		*icoords,
	COMPONENT	comp,
	INPUT_SOLN	**restart_soln,
	Locstate	state)
{
	int		ix = icoords[0];
	int		iy = icoords[1];

	ses_rt_coldp(state)  = restart_soln[0]->states2d[iy][ix];
	ses_rt_colde(state)  = restart_soln[1]->states2d[iy][ix];
	ses_rt_redp(state)   = restart_soln[2]->states2d[iy][ix];
	ses_rt_rede(state)   = restart_soln[3]->states2d[iy][ix];
	ses_rt_S(state)      = restart_soln[4]->states2d[iy][ix];
	ses_rt_adb_gam(state)= restart_soln[5]->states2d[iy][ix];
	ses_rt_gru_gam(state)= restart_soln[6]->states2d[iy][ix];
	ses_rt_riv(state)    = restart_soln[7]->states2d[iy][ix];
}		/*end ses_rt_restart_initializer*/


/*ARGSUSED*/
LOCAL void ses_re_restart_initializer(
	int		*icoords,
	COMPONENT	comp,
	INPUT_SOLN	**restart_soln,
	Locstate	state)
{
	int		ix = icoords[0];
	int		iy = icoords[1];

	ses_re_coldp(state)  = restart_soln[0]->states2d[iy][ix];
	ses_re_redp(state)   = restart_soln[1]->states2d[iy][ix];
	ses_re_Tvar(state)   = restart_soln[2]->states2d[iy][ix];
	ses_re_S(state)      = restart_soln[3]->states2d[iy][ix];
	ses_re_adb_gam(state)= restart_soln[4]->states2d[iy][ix];
	ses_re_gru_gam(state)= restart_soln[5]->states2d[iy][ix];
	ses_re_riv(state)    = restart_soln[6]->states2d[iy][ix];
}		/*end ses_re_restart_initializer*/

/*ARGSUSED*/
LOCAL void ses_rs_restart_initializer(
	int		*icoords,
	COMPONENT	comp,
	INPUT_SOLN	**restart_soln,
	Locstate	state)
{
	int		ix = icoords[0];
	int		iy = icoords[1];

	ses_rs_coldp(state)   = restart_soln[0]->states2d[iy][ix];
	ses_rs_colde(state)   = restart_soln[1]->states2d[iy][ix];
	ses_rs_redp(state)    = restart_soln[2]->states2d[iy][ix];
	ses_rs_rede(state)    = restart_soln[3]->states2d[iy][ix];
	ses_rs_Tvar(state)    = restart_soln[4]->states2d[iy][ix];
	ses_rs_adb_gam(state) = restart_soln[5]->states2d[iy][ix];
	ses_rs_gru_gam(state) = restart_soln[6]->states2d[iy][ix];
	ses_rs_riv(state)     = restart_soln[7]->states2d[iy][ix];
}		/*end ses_rs_restart_initializer*/

/*ARGSUSED*/
LOCAL void ses_ps_restart_initializer(
	int		*icoords,
	COMPONENT	comp,
	INPUT_SOLN	**restart_soln,
	Locstate	state)
{
	int		ix = icoords[0];
	int		iy = icoords[1];

	ses_ps_colde(state)   = restart_soln[0]->states2d[iy][ix];
	ses_ps_rede(state)    = restart_soln[1]->states2d[iy][ix];
	ses_ps_Tvar(state)    = restart_soln[2]->states2d[iy][ix];
	ses_ps_rho_var(state) = restart_soln[3]->states2d[iy][ix];
	ses_ps_adb_gam(state) = restart_soln[4]->states2d[iy][ix];
	ses_ps_gru_gam(state) = restart_soln[5]->states2d[iy][ix];
	ses_ps_riv(state)     = restart_soln[6]->states2d[iy][ix];
}		/*end ses_ps_restart_initializer*/

/*ARGSUSED*/
LOCAL void ses_rt_set_restart_intfc_states(
	double			*sl,
	double			*sr,
	int			var,
	POINT			*p,
	HYPER_SURF_ELEMENT	*hse,
	HYPER_SURF		*hs)
{
	BOND		*b = Bond_of_hse(hse);
	CURVE		*c = Curve_of_hs(hs);
	Locstate	lstate, rstate;

	lstate = left_state_at_point_on_curve(p,b,c);
	rstate = right_state_at_point_on_curve(p,b,c);

	switch (var)
	{
	case 0:
		ses_rt_coldp(lstate) = *sl;	ses_rt_coldp(rstate) = *sr;
		break;
	case 1:
		ses_rt_colde(lstate) = *sl;	ses_rt_colde(rstate) = *sr;
		break;
	case 2:
		ses_rt_redp(lstate) = *sl;	ses_rt_redp(rstate) = *sr;
		break;
	case 3:
		ses_rt_rede(lstate) = *sl;	ses_rt_rede(rstate) = *sr;
		break;
	case 4:
		ses_rt_S(lstate) = *sl;		ses_rt_S(rstate) = *sr;
		break;
	case 5:
		ses_rt_adb_gam(lstate) = *sl;	ses_rt_adb_gam(rstate) = *sr;
		break;
	case 6:
		ses_rt_gru_gam(lstate) = *sl;	ses_rt_gru_gam(rstate) = *sr;
		break;
	case 7:
		ses_rt_riv(lstate) = *sl;	ses_rt_riv(rstate) = *sr;
	}

}		/*end ses_rt_set_restart_intfc_states*/

/*ARGSUSED*/
LOCAL void ses_re_set_restart_intfc_states(
	double			*sl,
	double			*sr,
	int			var,
	POINT			*p,
	HYPER_SURF_ELEMENT	*hse,
	HYPER_SURF		*hs)
{
	BOND		*b = Bond_of_hse(hse);
	CURVE		*c = Curve_of_hs(hs);
	Locstate	lstate, rstate;

	lstate = left_state_at_point_on_curve(p,b,c);
	rstate = right_state_at_point_on_curve(p,b,c);

	switch (var)
	{
	case 0:
		ses_re_coldp(lstate) = *sl;	ses_re_coldp(rstate) = *sr;
		break;
	case 1:
		ses_re_redp(lstate) = *sl;	ses_re_redp(rstate) = *sr;
		break;
	case 2:
		ses_re_Tvar(lstate) = *sl;	ses_re_Tvar(rstate) = *sr;
		break;
	case 3:
		ses_re_S(lstate) = *sl;		ses_re_S(rstate) = *sr;
		break;
	case 4:
		ses_re_adb_gam(lstate) = *sl;	ses_re_adb_gam(rstate) = *sr;
		break;
	case 5:
		ses_re_gru_gam(lstate) = *sl;	ses_re_gru_gam(rstate) = *sr;
		break;
	case 6:
		ses_re_riv(lstate) = *sl;	ses_re_riv(rstate) = *sr;
	}
}		/*end ses_re_set_restart_intfc_states*/

/*ARGSUSED*/
LOCAL void ses_rs_set_restart_intfc_states(
	double			*sl,
	double			*sr,
	int			var,
	POINT			*p,
	HYPER_SURF_ELEMENT	*hse,
	HYPER_SURF		*hs)
{
	BOND		*b = Bond_of_hse(hse);
	CURVE		*c = Curve_of_hs(hs);
	Locstate	lstate, rstate;

	lstate = left_state_at_point_on_curve(p,b,c);
	rstate = right_state_at_point_on_curve(p,b,c);

	switch (var)
	{
	case 0:
		ses_rs_coldp(lstate) = *sl;	ses_rs_coldp(rstate) = *sr;
		break;
	case 1:
		ses_rs_colde(lstate) = *sl;	ses_rs_colde(rstate) = *sr;
		break;
	case 2:
		ses_rs_redp(lstate) = *sl;	ses_rs_redp(rstate) = *sr;
		break;
	case 3:
		ses_rs_rede(lstate) = *sl;	ses_rs_rede(rstate) = *sr;
		break;
	case 4:
		ses_rs_Tvar(lstate) = *sl;	ses_rs_Tvar(rstate) = *sr;
		break;
	case 5:
		ses_rs_adb_gam(lstate) = *sl;	ses_rs_adb_gam(rstate) = *sr;
		break;
	case 6:
		ses_rs_gru_gam(lstate) = *sl;	ses_rs_gru_gam(rstate) = *sr;
		break;
	case 7:
		ses_rs_riv(lstate) = *sl;	ses_rs_riv(rstate) = *sr;
	}

}		/*end ses_rs_set_restart_intfc_states*/


/*ARGSUSED*/
LOCAL void ses_ps_set_restart_intfc_states(
	double			*sl,
	double			*sr,
	int			var,
	POINT			*p,
	HYPER_SURF_ELEMENT	*hse,
	HYPER_SURF		*hs)
{
	BOND		*b = Bond_of_hse(hse);
	CURVE		*c = Curve_of_hs(hs);
	Locstate	lstate, rstate;

	lstate = left_state_at_point_on_curve(p,b,c);
	rstate = right_state_at_point_on_curve(p,b,c);

	switch (var)
	{
	case 0:
		ses_ps_colde(lstate) = *sl;	ses_ps_colde(rstate) = *sr;
		break;
	case 1:
		ses_ps_rede(lstate) = *sl;	ses_ps_rede(rstate) = *sr;
		break;
	case 2:
		ses_ps_Tvar(lstate) = *sl;	ses_ps_Tvar(rstate) = *sr;
		break;
	case 3:
		ses_ps_rho_var(lstate) = *sl;	ses_ps_rho_var(rstate) = *sr;
		break;
	case 4:
		ses_ps_adb_gam(lstate) = *sl;	ses_ps_adb_gam(rstate) = *sr;
		break;
	case 5:
		ses_ps_gru_gam(lstate) = *sl;	ses_ps_gru_gam(rstate) = *sr;
		break;
	case 6:
		ses_ps_riv(lstate) = *sl;	ses_ps_riv(rstate) = *sr;
	}
}		/*end ses_ps_set_restart_intfc_states*/

#endif /* defined(SESAME_CODE) && defined(TWOD) */
