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
*				gseshyp.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Uses the SESAME tabular equation of state to initialize
*	hyp solution functions for the evaluation of the equation
*	of state for gas dynamics.
*/

#if defined(SESAME_CODE) && defined(TWOD)

#define DEBUG_STRING    "ses_hyp"

#include <geos/sesame.h>

struct _SPOLY_FIT {
	double GAM;
	double Pinf;
	double Einf;
	double cv;
	double et;
};
typedef struct _SPOLY_FIT SPOLY_FIT;

	/* LOCAL Function Declarations */
LOCAL	boolean	PS_bisect(double,double*,POINTER);
LOCAL	boolean	RE_bisect(double,double*,POINTER);
LOCAL	boolean	RS_bisect(double,double*,POINTER);
LOCAL	boolean	VP_bisect(double,double*,POINTER);
LOCAL	boolean	remap_ses_intfc(Front*,Front*,
				void (*)(POINT*,BOND*,CURVE*,POINT*,
					 BOND*,CURVE*,boolean,RECT_GRID*,
					 POINTER),
				void (*)(INTERFACE*,INTERFACE*,
				 	 void (*)(POINT*,BOND*,CURVE*,
					  	  POINT*,BOND*,CURVE*,
					  	  boolean,RECT_GRID*,
						  POINTER),
				 	 POINTER),
				SESAME_EOS*,
				boolean,SESAME_TABLE_TYPE);
LOCAL	boolean	ses_replace_unphys_loop(NNLIST*,NNLIST**,CURVE**,Front*,
					int,double,int);
LOCAL	const char *eos_wave_type_as_string(int);
LOCAL	void	least_squares_spoly_fit(FILE*,const char*,double*,
					double*,double**,double**,int,int);
LOCAL	void	map_RT_to_PS(POINT*,BOND*,CURVE*,POINT*,BOND*,CURVE*,
			     boolean,RECT_GRID*,POINTER);
LOCAL	void	map_RT_to_RE(POINT*,BOND*,CURVE*,POINT*,BOND*,CURVE*,
			     boolean,RECT_GRID*,POINTER);
LOCAL	void	map_RT_to_RS(POINT*,BOND*,CURVE*,POINT*,BOND*,CURVE*,
			     boolean,RECT_GRID*,POINTER);
LOCAL	void	map_RT_to_VP(POINT*,BOND*,CURVE*,POINT*,BOND*,CURVE*,
			     boolean,RECT_GRID*,POINTER);
LOCAL	void	map_rt_intfc_states_to_tri_grid_intfc_states(Wave*,Front*);
LOCAL	int	read_eos_wave_type_from_string(const char*);
LOCAL	int	read_eos_node_type_from_string(const char*,INTERFACE*);
LOCAL	void	point_not_in_RE_window(double,double,COMPONENT,double,double,
				       double,double,SESAME_EOS*);
LOCAL	void	point_not_in_RS_window(double,double,COMPONENT,double,double,
				       double,double,SESAME_EOS*);
LOCAL	void	point_not_in_VP_window(double,double,COMPONENT,double,double,
				       double,double,SESAME_EOS*);
LOCAL	void	set_PS_rect_grid(INTERFACE*,INTERFACE*,
				 void (*)(POINT*,BOND*,CURVE*,POINT*,BOND*,
					  CURVE*,boolean,RECT_GRID*,POINTER),
				 POINTER);
LOCAL	void	set_RE_rect_grid(INTERFACE*,INTERFACE*,
				 void (*)(POINT*,BOND*,CURVE*,POINT*,BOND*,
					  CURVE*,boolean,RECT_GRID*,POINTER),
				 POINTER);
LOCAL	void	set_RS_rect_grid(INTERFACE*,INTERFACE*,
				 void (*)(POINT*,BOND*,CURVE*,POINT*,BOND*,
					  CURVE*,boolean,RECT_GRID*,POINTER),
				 POINTER);
LOCAL	void	set_VP_rect_grid(INTERFACE*,INTERFACE*,
				 void (*)(POINT*,BOND*,CURVE*,POINT*,BOND*,
					  CURVE*,boolean,RECT_GRID*,POINTER),
				 POINTER);
LOCAL	void	fprint_eos_node_type(FILE*,const char*,int,const char*,
				     INTERFACE*);
LOCAL	void	get_phase_ent_state(double,double,COMPONENT,Front*,Wave*,
				    double*,double*);
LOCAL	void	get_rho_on_phase_bdry_at_entropy(double,double,double,
						 double*,CURVE*);
LOCAL	void	get_riv_int(Front*,Wave*,double*,double,COMPONENT*,
			    double*,double*,double*,double*);
LOCAL	void	get_state_on_ref_curve(Front*,Wave*,double,double*,double*,
				       COMPONENT*);
LOCAL	void	init_PS_interior_states(SESAME_EOS*);
LOCAL	void	init_RS_interior_states(SESAME_EOS*);
LOCAL	void	init_rT_hyp_soln(Front*,Wave*,SESAME_EOS*,POINTER,POINTER);
LOCAL	void	init_ses_interior_states(Wave*,Front*,
					 void (*)(double*,COMPONENT,Locstate,
						  SESAME_EOS*),
					 SESAME_EOS*,SESAME_TABLE_TYPE);
LOCAL	void	linear_spoly_fit(SESAME_EOS*,double,double,SPOLY_FIT*);
LOCAL	void	print_sesame_interface(SESAME_EOS*,INTERFACE*,
				       SESAME_TABLE_TYPE);
LOCAL	void	print_spoly_fit(FILE*,SESAME_EOS*);
LOCAL	void	reset_rgrid_after_untan(Front*,Front*,
					void (*)(POINT*,BOND*,CURVE*,POINT*,
						 BOND*,CURVE*,boolean,
						 RECT_GRID*,POINTER),
					void (*)(INTERFACE*,INTERFACE*,
				     	         void (*)(POINT*,BOND*,CURVE*,
					                  POINT*,BOND*,CURVE*,
					                  boolean,RECT_GRID*,
							  POINTER),
				                 POINTER),
					SESAME_EOS*,SESAME_TABLE_TYPE);
LOCAL	void	rho_T_initializer(double*,COMPONENT,Locstate,SESAME_EOS*);
LOCAL	void	s_ses_solution(double*,COMPONENT,HYPER_SURF*,SIDE,
			       Front*,POINTER,Locstate,Locstate);
LOCAL	void	ses_PS_initializer(double*,COMPONENT,Locstate,SESAME_EOS*);
LOCAL	void	ses_RE_initializer(double*,COMPONENT,Locstate,SESAME_EOS*);
LOCAL	void	ses_RS_initializer(double*,COMPONENT,Locstate,SESAME_EOS*);
LOCAL	void	ses_VP_initializer(double*,COMPONENT,Locstate,SESAME_EOS*);
LOCAL	void	set_RT_riv(Front*,Wave*,SESAME_EOS*);
LOCAL	void	set_ph_RS_riv_state(SESAME_EOS*);
LOCAL	void	set_rho_T_state(Locstate,double,double,SESAME_EOS*);
LOCAL	void	set_up_multiphase_rt_table(Front*,Wave*,SESAME_EOS*,
					   POINTER,POINTER);
LOCAL	void	set_up_rt_table(Front*,Wave*,SESAME_EOS*);
LOCAL	void	set_smax(SESAME_EOS*);

EXPORT	const SESAME_TABLE_TYPE *sesame_table_numbers(void)
{
	static const SESAME_TABLE_TYPE tables[] = {
						    SESAME_RHO_TEMP,
						    SESAME_RHO_ENERGY,
						    SESAME_RHO_ENTROPY,
						    SESAME_PRESS_ENTROPY,
						    SESAME_VOLUME_PRESSURE
						   };
	return tables;
}		/*end sesame_table_numbers*/

EXPORT	const size_t *sesame_size_of_states(void)
{
	static const size_t sizests[] = {
					  sizeof(SES_RT_STATE),
					  sizeof(SES_RE_STATE),
					  sizeof(SES_RS_STATE),
					  sizeof(SES_PS_STATE),
					  sizeof(SES_VP_STATE)
				         };
	return sizests;
}		/*end sesame_size_of_states*/

EXPORT	const char **sesame_table_names(void)
{
	static const char *ses_table_names[] = {
						 "DENSITY-TEMPERATURE",
	    		                         "DENSITY-ENERGY",
	    		                         "DENSITY-ENTROPY",
	    		                         "PRESSURE-ENTROPY",
	    		                         "SPECIFIC_VOLUME-PRESSURE",
						 NULL
	    		                       };
	return ses_table_names;
}		/*end sesame_table_names*/

EXPORT	const char **sesame_headers(void)
{
	static char **ses_headers = NULL;
	const char **ses_table_names;
	size_t len;
	int   i, nnames;

	if (ses_headers == NULL)
	{
	    ses_table_names = sesame_table_names();
	    for (nnames = 0; ses_table_names[nnames] != NULL; nnames++);
	    uni_array(&ses_headers,nnames+1,sizeof(char*));
	    for (i = 0; i < nnames; i++)
	    {
		len = strlen(ses_table_names[i]) + 7;
		uni_array(ses_headers+i,len,sizeof(char));
		(void) sprintf(ses_headers[i],"%s TABLE",ses_table_names[i]);
	    }
	    ses_headers[nnames] = NULL;
	}
	return (const char **)ses_headers;
}		/*end sesame_headers*/

EXPORT void init_sesame_hyp_tri_solns(
	SESAME_EOS 	*seos,
	POINTER		pbdry,
	POINTER 	ccur)
{
	CHART 		**root = seos->root;
	Front 		**fr = seos->fr;
	Wave 		**wave = seos->wave;
	RECT_GRID	*gr;
	double		p_gr_max, p_gr_min;
	int 		i;
	char 		s[Gets_BUF_SIZE];
	static const char	*fmt = "%lf %lf";

	DEBUG_ENTER(init_sesame_hyp_tri_solns)

	/* Generate the density temperature table */

	init_rT_hyp_soln(fr[SESAME_RHO_TEMP],wave[SESAME_RHO_TEMP],seos,
			 pbdry,ccur);
	gr = fr[SESAME_RHO_TEMP]->rect_grid;
	RS_entropy_scale(seos) =
		(gr->U[0] - gr->L[0])/(Entropy_max(seos)-Entropy_min(seos));
	RS_entropy_shift(seos) = 0.5*(gr->U[0] + gr->L[0] -
		RS_entropy_scale(seos)*(Entropy_max(seos)+Entropy_min(seos)));
	p_gr_min = ses_ps_grid_from_press(Pressure_min(seos),seos);
	p_gr_max = ses_ps_grid_from_press(Pressure_max(seos),seos);
	PS_entropy_scale(seos) =
		(p_gr_max - p_gr_min)/(Entropy_max(seos)-Entropy_min(seos));
	PS_entropy_shift(seos) = 0.5*(p_gr_max + p_gr_min -
		PS_entropy_scale(seos)*(Entropy_max(seos)+Entropy_min(seos)));

	screen("The grid coordinate representing entropy ");
	screen("on the density-entropy table\n");
	screen("is a linear function of entropy S_grid = a_rs + b_rs * S\n");
	screen("Enter a_rs, b_rs (default = %g %g): ",
		RS_entropy_shift(seos),RS_entropy_scale(seos));
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if(sscanf(s,fmt,&RS_entropy_shift(seos),&RS_entropy_scale(seos))!=2)
	    {
		screen("ERROR in init_sesame_hyp_tri_solns(), ");
		screen("Both a_rs AND b_rs must be specified\n");
		clean_up(ERROR);
	    }
	}
	screen("The grid coordinate representing entropy ");
	screen("on the pressure-entropy table\n");
	screen("is a linear function of entropy S_grid = a_ps + b_ps * S\n");
	screen("Enter a_ps, b_ps (default = %g %g): ",
		PS_entropy_shift(seos),PS_entropy_scale(seos));
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if(sscanf(s,fmt,&PS_entropy_shift(seos),&PS_entropy_scale(seos))!=2)
	    {
		screen("ERROR in init_sesame_hyp_tri_solns(), ");
		screen("Both a_ps AND b_ps must be specified\n");
		clean_up(ERROR);
	    }
	}

	if (DEBUG)
	{
	    (void) printf("Initialized %s interface\n",
			  ses_table_name(seos,SESAME_RHO_TEMP));
	    print_sesame_interface(seos,fr[SESAME_RHO_TEMP]->interf,
				   SESAME_RHO_TEMP);
	}

	/* Remap the density-temperature interface onto the 	*/
	/* 		other interfaces			*/

	set_ses_inv_hyp_solns(seos); 

#if defined(DEBUG_SES_HYP)
	if (DEBUG)
	{
	    for (i = 1; i < NUMBER_SESAME_TABLES; i++)
	    {
	        print_sesame_interface(seos,fr[i]->interf,i);
	        (void) printf("Remapped %s rect grid\n",s);
	        print_RECT_GRID_structure(fr[i]->rect_grid);
	    }
	}
#endif /* defined(DEBUG_SES_HYP) */

 	screen("Do you want to save the inversions (default = no)?: ");
	(void) Gets(s);
	if (s[0] == 'Y' || s[0] == 'y')
	{
	    FILE	*file = stdout;
	    boolean	bio = is_binary_output();
	    Printplot	Prt, *prt = &Prt;
	    int		ntables = 4;
	    const SESAME_TABLE_TYPE *tables;

	    tables = sesame_table_numbers();
	    screen("Print density temperature table only(dflt = NO): ");
	    (void) Gets(s);
	    if (s[0] == 'Y' || s[0] == 'y')
		ntables = 1;
	    screen("Enter a filename for the inverted tables (dflt=stdout): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
		if ((file = fopen(s,"w")) == NULL)
		{
		     screen("ERROR in init_sesame_hyp_tri_solns(), "
		            "can't open %s\n",s);
		     clean_up(ERROR);
		}
		record_print_version(file);
		print_title_for_sesame(file,seos);
	        (void) fprintf(file,"\n\n");
	    }
	    screen("Print tables in binary (dflt = %s): ",
		   (bio==YES)?"yes":"no");
	    (void) Gets(s);
	    if (s[0] == 'N' || s[0] == 'n')
		set_binary_output(NO);
	    else if (s[0] == 'Y' || s[0] == 'y')
		set_binary_output(YES);

	    zero_scalar(prt,sizeof(Printplot));
	    prt->file = file;
	    prt->_print_states[1] = d_print_states;

	    (void) foutput(file);
	    (void) fprintf(file,"\t\t\tSESAME RESTART INFORMATION\n\n");
	    record_print_version(file);
	    (void) fprintf(file,"\tSESAME_EOS variables\n");
	    fprint_SESAME_params(file,seos);

	    (void) fprintf(file,"\tNumber of tables    = %12d\n",ntables);

	    if (ntables > 1 && debugging("graph_phase"))
		ntables = 5;

	    for (i = 0; i < ntables; i++)
	    {
	    	(void) foutput(file);
	    	(void) fprintf(file,"\t\t\t%s TABLE\n",
			ses_table_name(seos,tables[i]));
	    	init_ses_prt(prt,fr[tables[i]],tables[i]);
	    	ses_printout(root[tables[i]],prt,NO,0);
	    }

	    (void) foutput(file);
	    (void) fprintf(file,"\t\t\tEnd SESAME RESTART INFORMATION\n\n");
	    set_binary_output(bio);
	}
	print_spoly_fit(stdout,seos);
	DEBUG_LEAVE(init_sesame_hyp_tri_solns)
}		/*end init_sesame_hyp_tri_solns*/

EXPORT	void set_ses_inv_hyp_solns(
	SESAME_EOS    *seos)
{
	Front 		**fr = seos->fr;
	Wave 		**wave = seos->wave;

	DEBUG_ENTER(set_ses_inv_hyp_solns)

	if (remap_ses_intfc(fr[SESAME_RHO_ENERGY],fr[SESAME_RHO_TEMP],
		map_RT_to_RE,set_RE_rect_grid,seos,YES,SESAME_RHO_ENERGY) == 
								FUNCTION_FAILED)
	{
	    screen("ERROR in set_ses_inv_hyp_solns(), "
	           "remap_ses_intfc() failed for %s table\n",
		   ses_table_name(seos,SESAME_RHO_TEMP));
	    print_sesame_interface(seos,fr[SESAME_RHO_ENERGY]->interf,
				   SESAME_RHO_ENERGY);
	    clean_up(ERROR);
	}
	init_ses_interior_states(wave[SESAME_RHO_ENERGY],fr[SESAME_RHO_ENERGY],
				 ses_RE_initializer,seos,SESAME_RHO_ENERGY);

	if (remap_ses_intfc(fr[SESAME_RHO_ENTROPY],fr[SESAME_RHO_TEMP],
		map_RT_to_RS,set_RS_rect_grid,seos,YES,SESAME_RHO_ENTROPY) ==
								FUNCTION_FAILED)
	{
	    screen("ERROR in set_ses_inv_hyp_solns(), "
	           "remap_ses_intfc() failed for %s table\n",
		   ses_table_name(seos,SESAME_RHO_TEMP));
	    print_sesame_interface(seos,fr[SESAME_RHO_ENTROPY]->interf,
				   SESAME_RHO_ENTROPY);
	    clean_up(ERROR);
	}
	init_RS_interior_states(seos);
#if defined(PHASE_CODE)
	if (multiphase_eos(seos) == YES)
		set_ph_RS_riv_state(seos);
#endif /* defined(PHASE_CODE) */

	if (remap_ses_intfc(fr[SESAME_PRESS_ENTROPY],fr[SESAME_RHO_TEMP],
		map_RT_to_PS,set_PS_rect_grid,seos,NO,SESAME_PRESS_ENTROPY) ==
								FUNCTION_FAILED)
	{
	    screen("ERROR in set_ses_inv_hyp_solns(), "
	           "remap_ses_intfc() failed for %s table\n",
		   ses_table_name(seos,SESAME_PRESS_ENTROPY));
	    print_sesame_interface(seos,fr[SESAME_PRESS_ENTROPY]->interf,
				   SESAME_PRESS_ENTROPY);
	    clean_up(ERROR);
	}
	init_PS_interior_states(seos);

	if (debugging("graph_phase"))
	{
	    if (remap_ses_intfc(fr[SESAME_VOLUME_PRESSURE],fr[SESAME_RHO_TEMP],
		map_RT_to_VP,set_VP_rect_grid,seos,NO,SESAME_VOLUME_PRESSURE) ==
			FUNCTION_FAILED)
	    {
		screen("ERROR in set_ses_inv_hyp_solns(), "
		       "remap_ses_intfc() failed for %s table\n",
		       ses_table_name(seos,SESAME_VOLUME_PRESSURE));
	        print_sesame_interface(seos,fr[SESAME_VOLUME_PRESSURE]->interf,
				       SESAME_VOLUME_PRESSURE);
		clean_up(ERROR);
	    }
	    init_ses_interior_states(wave[SESAME_VOLUME_PRESSURE],
		fr[SESAME_VOLUME_PRESSURE],ses_VP_initializer,seos,
		SESAME_VOLUME_PRESSURE);
	}
	DEBUG_LEAVE(set_ses_inv_hyp_solns)
}		/*end set_ses_inv_hyp_solns*/

LOCAL	boolean	remap_ses_intfc(
	Front		  *newfr,
	Front		  *fr,
	void		  (*remap)(POINT*,BOND*,CURVE*,POINT*,BOND*,
				   CURVE*,boolean,RECT_GRID*,POINTER),
	void		  (*set_rgrid)(INTERFACE*,INTERFACE*,
				       void (*)(POINT*,BOND*,CURVE*,
					        POINT*,BOND*,CURVE*,
					        boolean,RECT_GRID*,POINTER),
				       POINTER),
	SESAME_EOS	  *seos,
	boolean		  density_table,
	SESAME_TABLE_TYPE table)
{
	CROSS	         *cross;
	int	         flag = LAST_ATTEMPT_TO_UNTANGLE;
	CURVE	         **c;
	NODE	         **n;
	boolean	         sav_intrp;

	DEBUG_ENTER(remap_ses_intfc)

	set_size_of_intfc_state(newfr->sizest);
	newfr->interf = remap_interface(fr->interf,remap,set_rgrid,
					(POINTER)seos);
	if (newfr->interf == NULL)
	{
	    (void) printf("WARNING in remap_ses_intfc(), "
	                  "remap_interface() failed\n");
	    DEBUG_LEAVE(remap_ses_intfc)
	    return FUNCTION_FAILED;
	}
	set_ses_hooks(newfr,table);
	copy_rect_grid(newfr->rect_grid,computational_grid(newfr->interf));
	set_computational_grid(newfr->interf,newfr->rect_grid);
	if (DEBUG)
	{
	    (void) printf("Interface returned from remap_interface()\n");
	    print_sesame_interface(seos,newfr->interf,table);
	    if (debugging("ses_states"))
		show_intfc_states(newfr->interf);
	}
	sav_intrp = interpolate_intfc_states(newfr->interf);
	interpolate_intfc_states(newfr->interf) = YES;
	for (c = newfr->interf->curves; c && *c; c++)
	{
	    if ((density_table == YES) && (wave_type(*c)==ISOCHOR_BOUNDARY))
	    {
		if ((*c)->num_points > 2)
	        {
		    BOND *b;
		    double ys = Coords((*c)->start->posn)[1];
		    double ye = Coords((*c)->end->posn)[1];
		    for (b = (*c)->first->next; b != NULL; b = b->next)
		    {
		        if (!Between(Coords(b->start)[1],ys,ye))
			    (void) delete_start_of_bond(b,*c);
		    }
	        }
	    }
	    (void) equi_curve_redistribute(newfr,*c,YES);
	}

	if (intfc_delete_fold_back_bonds(newfr) == NO)
	{
	    (void) printf("WARNING in remap_ses_intfc(), "
	                  "intfc_delete_fold_back_bonds() failed\n");
	    print_sesame_interface(seos,newfr->interf,table);
	    if (debugging("ses_states"))
		show_intfc_states(newfr->interf);
	    interpolate_intfc_states(newfr->interf) = sav_intrp;
	    DEBUG_LEAVE(remap_ses_intfc)
	    return FUNCTION_FAILED;
	}
	for (n = newfr->interf->nodes; n && *n; n++)
	{
#if defined(PHASE_CODE)
	    if (node_type(*n) == PHASE_BDRY_NODE)
	    	continue;
#endif /* defined(PHASE_CODE) */
	    node_type(*n) = EOS_BOUNDARY_NODE;
	}
	if (intersections(newfr->interf,&cross,YES) == FUNCTION_FAILED)
	{
	    (void) printf("WARNING in remap_ses_intfc(), "
	                  "intersections() failed\n"
	                  "Interface after remap_interface()\n");
	    print_sesame_interface(seos,newfr->interf,table);
	    if (debugging("ses_states"))
		show_intfc_states(newfr->interf);
	    interpolate_intfc_states(newfr->interf) = sav_intrp;
	    DEBUG_LEAVE(remap_ses_intfc)
	    return FUNCTION_FAILED;
	}
	if (cross != NULL)
	{
	    int status;
	    char	name[80];

	    (void) sprintf(name,"%s inside remap_ses_intfc() ",
			   ses_table_name(seos,table));
	    (void) print_number_of_tangles(name,newfr->interf,cross);
	    status = scalar_unravel(newfr,&cross,flag);
	    if (status != CURVES_UNTANGLED)
	    {
	        screen("ERROR in remap_ses_intfc(), interface is tangled\n");
	        print_sesame_interface(seos,newfr->interf,table);
		if (debugging("ses_states"))
			show_intfc_states(newfr->interf);
	        clean_up(ERROR);
	        interpolate_intfc_states(newfr->interf) = sav_intrp;
	        DEBUG_LEAVE(remap_ses_intfc)
	        return FUNCTION_FAILED;
	    }
	    reset_rgrid_after_untan(fr,newfr,remap,set_rgrid,seos,table);
	}
	interpolate_intfc_states(newfr->interf) = sav_intrp;
	if (DEBUG)
	{
	    (void) printf("Interface after remap_interface()\n");
	    print_sesame_interface(seos,newfr->interf,table);
	    if (debugging("ses_states"))
	    	show_intfc_states(newfr->interf);
	}
	DEBUG_LEAVE(remap_ses_intfc)
	return FUNCTION_SUCCEEDED;
}		/*end remap_ses_intfc*/

EXPORT	void set_user_hooks_for_sesame(void)
{
	F_USER_INTERFACE *fuh = f_user_hook(2);
	G_USER_INTERFACE *guh = g_user_hook(2);

	zero_scalar(guh,sizeof(G_USER_INTERFACE));
	guh->_intfc_type = EOS_INTERFACE;

	fuh->_wave_type_as_string = eos_wave_type_as_string;
	fuh->_fprint_hsbdry_type = fprint_eos_node_type;
	fuh->_read_wave_type_from_string = read_eos_wave_type_from_string;
	fuh->_read_hsbdry_type_from_string = read_eos_node_type_from_string;
	fuh->_alloc_state = f_alloc_state;
	fuh->_alloc_intfc_state = f_alloc_intfc_state;
	fuh->_clear_state = f_clear_state;
}		/*end set_user_hooks_for_sesame*/

EXPORT	void	set_ses_hooks(
	Front             *front,
	SESAME_TABLE_TYPE table)
{
	F_USER_INTERFACE *fuh;
	static void (*prst[])(FILE*,Locstate,INTERFACE*) = {
						ses_rt_fprint_state_data,
						ses_re_fprint_state_data,
						ses_rs_fprint_state_data,
						ses_ps_fprint_state_data,
						ses_vp_fprint_state_data
					     };

	fuh = &f_user_interface(front->interf);
	fuh->_bi_interpolate_intfc_states = front->_state_interpolator;
	fuh->_tri_interpolate_intfc_states = front->_tri_state_interpolator;
	fuh->_fprint_state_data = prst[table];
	fuh->_fprint_intfc_state = prst[table];
}		/*end set_ses_hooks*/

LOCAL	void	reset_rgrid_after_untan(
	Front		  *fr,
	Front		  *newfr,
	void		  (*remap)(POINT*,BOND*,CURVE*,POINT*,BOND*,
				   CURVE*,boolean,RECT_GRID*,POINTER),
	void		  (*set_rgrid)(INTERFACE*,INTERFACE*,
			  	       void (*)(POINT*,BOND*,CURVE*,
			  		        POINT*,BOND*,CURVE*,
			  		        boolean,RECT_GRID*,POINTER),
			  	     POINTER),
	SESAME_EOS	  *seos,
	SESAME_TABLE_TYPE table)
{
	POINT			*p;
	HYPER_SURF_ELEMENT	*hse;
	HYPER_SURF		*hs;
	INTERFACE		*intfc = newfr->interf;
	RECT_GRID		*cgr = computational_grid(intfc);
	RECT_GRID		*tgr = &topological_grid(intfc);
	double			L[3], U[3];
	int			dim = intfc->dim;
	int			i;

	for (i = 0; i < dim; i++)
	{
	    L[i] =  HUGE_VAL;
	    U[i] = -HUGE_VAL;
	}
	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs) != NO)
	{
	    for (i = 0; i < dim; i++)
	    {
	    	L[i] = min(L[i],Coords(p)[i]);
	    	U[i] = max(U[i],Coords(p)[i]);
	    }
	}
	for (i = 0; i < dim; i++)
	{
	    if (L[i] > newfr->rect_grid->L[i])
	    	break;
	    if (U[i] < newfr->rect_grid->U[i])
	    	break;
	}
	if (i == dim)
		return;
	for (i = 0; i < dim; i++)
	{
	    cgr->L[i] = tgr->L[i] = L[i];
	    cgr->U[i] = tgr->U[i] = U[i];
	}
	
	screen("The %s interface was tangled,  after untangling the ",
		ses_table_name(seos,table));
	screen("\tthe new grid limits need to be reset\n");
	(*set_rgrid)(fr->interf,intfc,remap,(POINTER)seos);
	copy_rect_grid(newfr->rect_grid,computational_grid(newfr->interf));
	set_computational_grid(newfr->interf,newfr->rect_grid);
}		/*end reset_rgrid_after_untan*/

LOCAL void init_rT_hyp_soln(
	Front 		*fr,
	Wave 		*wave,
	SESAME_EOS 	*seos,
	POINTER 	pbdry,
	POINTER 	ccur)
{
	RECT_GRID	*gr;
	REMAP           Remap;
	size_t 		sizest = fr->sizest;
	double 		L[MAXD], U[MAXD];
	int 		nhyp[2];
	int 		nrho_hyp = Nrho_hyp(seos), nT_hyp = Ntemp_hyp(seos);

	DEBUG_ENTER(init_rT_hyp_soln)

	set_remap(2,IDENTITY_REMAP,&Remap);
	L[0] = ses_rt_grid_from_rho(Rho_min(seos),seos);
	U[0] = ses_rt_grid_from_rho(Rho_max(seos),seos);
	L[1] = ses_rt_grid_from_temp(Temp_min(seos),seos);
	U[1] = ses_rt_grid_from_temp(Temp_max(seos),seos);
	nhyp[0] = nrho_hyp;
	nhyp[1] = nT_hyp;

	if (DEBUG)
	{
		(void) printf("sizest = %d, nrho_hyp = %d, nT_hyp = %d\n",
			      (int)sizest,nrho_hyp,nT_hyp);
		(void) printf("rho_min = %g, rhomin_grid = %g\n",
			      Rho_min(seos),L[0]);
		(void) printf("rho_max = %g, rhomax_grid = %g\n",
			      Rho_max(seos),U[0]);
		(void) printf("Tmin = %g, Tmin_grid = %g\n",
			      Temp_min(seos),L[1]);
		(void) printf("Tmax = %g, Tmax_grid = %g\n",
			      Temp_max(seos),U[1]);
	}


	/*   initialize interface   */

	set_size_of_intfc_state(sizest);
	fr->interf = make_interface(2);
	set_ses_hooks(fr,SESAME_RHO_TEMP);

	/*   initialize the rect grid   */

	gr = fr->rect_grid;
	set_rect_grid(L,U,L,U,NOBUF,NOBUF,nhyp,2,&Remap,gr);
	set_topological_grid(fr->interf,gr);
	(void) adjust_top_grid_for_square(&topological_grid(fr->interf),gr);
	set_computational_grid(fr->interf,gr);

	seos->terr0 = 0.001*hypot(gr->h[0],gr->h[1]);/*TOLERANCE*/
	Gru_gam_abs_max(seos) = -HUGE_VAL;

#if defined(PHASE_CODE)
	if (multiphase_eos(seos) == YES)
	    set_up_multiphase_rt_table(fr,wave,seos,pbdry,ccur);

	else
#endif /* defined(PHASE_CODE) */
	    set_up_rt_table(fr,wave,seos);

	DEBUG_LEAVE(init_rT_hyp_soln)
}		/*end init_rT_hyp_soln*/

LOCAL	void	set_up_rt_table(
	Front		*fr,
	Wave		*wave,
	SESAME_EOS	*seos)
{
	NODE 	*node[4];
	CURVE 	*curve[4];
	POINT 	*pt;
	double 	*L, *U;
	double 	coords[MAXD];
	double 	drho, dT;
	double 	rho_grid, T_grid;
	int 	i;
	int 	nrho_hyp = Nrho_hyp(seos), nT_hyp = Ntemp_hyp(seos);
	size_t	sizest = fr->sizest;

	DEBUG_ENTER(set_up_rt_table)
	drho = fr->rect_grid->h[0];
	dT = fr->rect_grid->h[1];
	L = fr->rect_grid->L;
	U = fr->rect_grid->U;

	coords[0] = L[0];	coords[1] = L[1];
	node[0] = make_node(Point(coords));
	coords[0] = L[0];	coords[1] = U[1];
	node[1] = make_node(Point(coords));
	coords[0] = U[0];	coords[1] = U[1];
	node[2] = make_node(Point(coords));
	coords[0] = U[0];	coords[1] = L[1];
	node[3] = make_node(Point(coords));

	for (i = 0; i < 4; i++)
	{
	    node_type(node[i]) = FIXED_NODE;
	    curve[i] = make_curve(EXTERIOR_COMP,COMP_PURE_PHASE,
				  node[i], node[(i+1)%4]);
	    wave_type(curve[i]) = (i%2) ? ISOTHERM_BOUNDARY :
	    			          ISOCHOR_BOUNDARY;
	    start_status(curve[i]) = INCIDENT;
	    end_status(curve[i]) = INCIDENT;
	    set_is_bdry(curve[i]);

	    set_rho_T_state(left_start_state(curve[i]),
	    		    Coords(curve[i]->start->posn)[0],
	    		    Coords(curve[i]->start->posn)[1],seos);
	    ft_assign(right_start_state(curve[i]),
		   left_start_state(curve[i]),sizest);
	    set_rho_T_state(left_end_state(curve[i]),
	    		    Coords(curve[i]->end->posn)[0],
	    		    Coords(curve[i]->end->posn)[1],seos);
	    ft_assign(right_end_state(curve[i]),left_end_state(curve[i]),sizest);
	}
	interpolate_intfc_states(fr->interf) = NO;
	for (i = 1, T_grid = L[1] + dT; i < nT_hyp; i++, T_grid += dT)
	{
	    coords[0] = L[0];	coords[1] = T_grid;
	    pt = Point(coords);
	    set_rho_T_state(left_state(pt),Coords(pt)[0],
	    		    Coords(pt)[1],seos);
	    ft_assign(right_state(pt),left_state(pt),sizest);
	    if (insert_point_in_bond(pt,curve[0]->last,curve[0]) !=
		FUNCTION_SUCCEEDED)
	    {
	        screen("ERROR in set_up_rt_table(), "
		       "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }
	    coords[0] = U[0];
	    pt = Point(coords);
	    set_rho_T_state(left_state(pt),Coords(pt)[0],
	    		    Coords(pt)[1],seos);
	    ft_assign(right_state(pt),left_state(pt),sizest);
	    if (insert_point_in_bond(pt,curve[2]->first,curve[2]) !=
		FUNCTION_SUCCEEDED)
	    {
	        screen("ERROR in set_up_rt_table(), "
		       "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }
	}
	for (i=1, rho_grid=L[0]+drho; i < nrho_hyp; i++, rho_grid += drho)
	{
	    coords[0] = rho_grid;	coords[1] = U[1];
	    pt = Point(coords);
	    set_rho_T_state(left_state(pt),Coords(pt)[0],
	    		Coords(pt)[1],seos);
	    ft_assign(right_state(pt),left_state(pt),sizest);
	    if (insert_point_in_bond(pt,curve[1]->last,curve[1]) !=
		FUNCTION_SUCCEEDED)
	    {
	        screen("ERROR in set_up_rt_table(), "
		       "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }
	    coords[1] = L[1];
	    pt = Point(coords);
	    set_rho_T_state(left_state(pt),Coords(pt)[0],Coords(pt)[1],seos);
	    ft_assign(right_state(pt),left_state(pt),sizest);
	    if (insert_point_in_bond(pt,curve[3]->first,curve[3]) !=
		FUNCTION_SUCCEEDED)
	    {
	        screen("ERROR in set_up_rt_table(), "
		       "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }
	}

	/*  init_ses_interior_states calls init_ses_hyp_soln_func */

	init_ses_interior_states(wave,fr,rho_T_initializer,seos,
				 SESAME_RHO_TEMP);
	set_RT_entropy(fr,wave,seos);
	set_RT_riv(fr,wave,seos);
	map_rt_intfc_states_to_tri_grid_intfc_states(wave,fr);
	DEBUG_LEAVE(set_up_rt_table)
}		/*end set_up_rt_table*/

LOCAL	void	map_rt_intfc_states_to_tri_grid_intfc_states(
	Wave	*wave,
	Front	*fr)
{
	POINT				*p;
	HYPER_SURF			*hs, *hs_on;
	HYPER_SURF_ELEMENT		*hse;
	INTERFACE	*intfc = fr->interf;
	INTERFACE	*tg_intfc = wave_tri_soln(wave)->tri_grid->grid_intfc;
	Locstate	sl, sr, state;
	double		coords_on[3];

	scalar(&state,sizeof(SES_RT_STATE));
	(void) next_point(tg_intfc,NULL,NULL,NULL);
	while (next_point(tg_intfc,&p,&hse,&hs) != NO)
	{
	    if (nearest_intfc_state(Coords(p),positive_component(hs),
				intfc,state,coords_on,&hs_on) != YES)
	    {
		screen("ERROR in "
		       "map_rt_intfc_states_to_tri_grid_intfc_states(), "
		       "nearest_intfc_state failed\n");
		clean_up(ERROR);
	    }
	    slsr(p,hse,hs,&sl,&sr);
	    ses_rt_S(sl) = ses_rt_S(sr) = ses_rt_S(state);
	    ses_rt_riv(sl) = ses_rt_riv(sr) = ses_rt_riv(state);
	}
	free(state);
}		/*end map_rt_intfc_states_to_tri_grid_intfc_states*/

/*ARGSUSED*/
LOCAL void rho_T_initializer(
	double 		*coords,
	COMPONENT 	comp,
	Locstate 	state,
	SESAME_EOS 	*seos)
{
	double 		rho_grid = coords[0], T_grid = coords[1];
	size_t 		sizest = sizeof(SES_RT_STATE);

	if (comp == EXTERIOR_COMP)
		clear_state(seos->fr[SESAME_RHO_TEMP]->interf,state,sizest);

	set_rho_T_state(state,rho_grid,T_grid,seos);
}		/*end rho_T_initializer*/


typedef struct {
	Front		*fr;
	Wave		*wv;
	COMPONENT	comp;
	double		rho_grid;
	SESAME_EOS	*seos;
} SESAME_INVERT_PARAMS;

/*
*			ses_RE_initializer():
*	Called from init_interior_states at each mesh point in the
*	density-specific internal energy interface.  Uses the (rho, T)
*	hyp function to invert e = e(rho,T) to get T.
*/

/*ARGSUSED*/
LOCAL void ses_RE_initializer(
	double 		*coords,
	COMPONENT 	comp,
	Locstate 	state,
	SESAME_EOS 	*seos)
{
	double 		rho_grid = coords[0], e_grid = coords[1];
	double		BIS_EPS = seos->BIS_EPS;
	double		ABS_SES_EPS = seos->ABS_SES_EPS;
	SESAME_INVERT_PARAMS re_inv_params;
	Front 		*rt_fr = re_inv_params.fr = seos->fr[SESAME_RHO_TEMP];
	Wave 		*rt_wv = re_inv_params.wv = seos->wave[SESAME_RHO_TEMP];
	double 		Tmin_grid, Tmax_grid;
	double 		T_grid, var[NUM_SES_VAR];
	double 		epsilon, delta;
	double 		T_min, T_max, emin, emax, en;
	double 		lcoords[MAXD];

	re_inv_params.seos = seos;
	Tmin_grid = rt_fr->rect_grid->L[1];
	Tmax_grid = rt_fr->rect_grid->U[1];
	epsilon = seos->reler0*fabs(e_grid);

#if defined(PHASE_CODE)
	if (multiphase_eos(seos))
	{
		double T, E, S;

		if (get_phase_hyp_state(rho_grid,seos,&T,&E,&S))
		{
			E  = ses_re_grid_from_engy(E,seos);
			if (comp == COMP_MIXED_PHASE)	Tmax_grid = T;
			if (comp == COMP_PURE_PHASE)	Tmin_grid = T;
		}
	}
#endif /* defined(PHASE_CODE) */
	delta = seos->reler0*0.5*(Tmax_grid + Tmin_grid);

	re_inv_params.rho_grid = lcoords[0] = coords[0];
	re_inv_params.comp = comp;
	en = ses_re_engy_from_grid(e_grid,seos);
	set_ses_intrp_flag(EVALUATE_ENERGY,SESAME_RHO_TEMP);
	T_min =  ses_rt_temp_from_grid(Tmin_grid,seos);
	lcoords[1] = Tmin_grid;
	ses_solution(lcoords,comp,NULL,POSITIVE_SIDE,rt_fr,rt_wv,var);
	emin = var[RT_CE] + T_min*var[RT_RE];
	T_max = ses_rt_temp_from_grid(Tmax_grid,seos);
	lcoords[1] = Tmax_grid;
	ses_solution(lcoords,comp,NULL,POSITIVE_SIDE,rt_fr,rt_wv,var);
	emax = var[RT_CE] + T_max*var[RT_RE];
	if (en <= emin)
	{
	    if ((emin - en) <= (fabs(en)*BIS_EPS + ABS_SES_EPS))
	    	T_grid = Tmin_grid;
	    else
	    {
	        screen("ERROR in ses_RE_initializer(), en < emin\n");
		point_not_in_RE_window(ses_re_rho_from_grid(rho_grid,seos),en,
				       comp,T_min,T_max,emin,emax,seos);
	        clean_up (ERROR);
	    }
	}
	else if (en >= emax)
	{
	    if(en - emax <= fabs(en)*BIS_EPS)
	        T_grid = Tmax_grid;
	    else
	    {
	        screen("ERROR in ses_RE_initializer(), en > emax\n");
		point_not_in_RE_window(ses_re_rho_from_grid(rho_grid,seos),en,
				       comp,T_min,T_max,emin,emax,seos);
	        clean_up (ERROR);
	    }
	}
	else
	{
	    if (find_root(RE_bisect,(POINTER)&re_inv_params,e_grid,&T_grid,
			  Tmin_grid,Tmax_grid,epsilon,delta) == FUNCTION_FAILED) 
	    {
		if (DEBUG)
		{
	    	     (void) printf("WARNING in ses_RE_initializer(), "
	    	                   "Unable to invert function.\n");
		}
	    }
	}
	set_ses_intrp_flag_all(SESAME_RHO_TEMP);
	lcoords[1] = T_grid;
	ses_solution(lcoords,comp,NULL,POSITIVE_SIDE,rt_fr,rt_wv,var);
	ses_re_coldp(state) = var[RT_CP];
	ses_re_redp(state) = var[RT_RP];
	ses_re_Tvar(state) = T_grid;
	ses_re_S(state) = var[RT_S];
	ses_re_adb_gam(state) = var[RT_AG];
	ses_re_gru_gam(state) = var[RT_GG];
	ses_re_riv(state) = var[RT_RF];
}		/*end ses_RE_initializer*/


LOCAL boolean RE_bisect(
	double 		T_grid,
	double 		*e_grid,
	POINTER 	params)
{
	double 		T, e, var[NUM_SES_VAR];
	Front 		*rt_fr = ((SESAME_INVERT_PARAMS *) params)->fr;
	Wave 		*rt_wv = ((SESAME_INVERT_PARAMS *) params)->wv;
	COMPONENT 	comp = ((SESAME_INVERT_PARAMS *) params)->comp;
	SESAME_EOS	*seos = ((SESAME_INVERT_PARAMS *) params)->seos;
	double 		coords[MAXD];

	coords[0] = ((SESAME_INVERT_PARAMS *) params)->rho_grid;
	coords[1] = T_grid;

	ses_solution(coords,comp,NULL,POSITIVE_SIDE,rt_fr,rt_wv,var);
	T = ses_rt_temp_from_grid(T_grid,seos);
	e = var[RT_CE] + T*var[RT_RE];
	*e_grid = ses_re_grid_from_engy(e,seos);
	return FUNCTION_SUCCEEDED;
}		/*end RE_bisect*/

typedef struct {
	Front		*fr;
	Wave		*wv;
	COMPONENT	comp;
	double		vol_grid;
	SESAME_EOS	*seos;
} SESAME_VP_INVERT_PARAMS;

/*
*			ses_VP_initializer():
*	Called from init_interior_states at each mesh point in the
*	volume - pressure interface.  Uses the (rho, T)
*	hyp function to invert P = P(vol,T) to get T.
*/

/*ARGSUSED*/
LOCAL void ses_VP_initializer(
	double 		*coords,
	COMPONENT 	comp,
	Locstate 	state,
	SESAME_EOS 	*seos)
{
	double		BIS_EPS = seos->BIS_EPS;
	double		ABS_SES_EPS = seos->ABS_SES_EPS;
	double 		vol_grid = coords[0], p_grid = coords[1];
	SESAME_VP_INVERT_PARAMS vp_inv_params;
	Front 		*rt_fr = vp_inv_params.fr = seos->fr[SESAME_RHO_TEMP];
	Wave 		*rt_wv = vp_inv_params.wv = seos->wave[SESAME_RHO_TEMP];
	double 		Tmin_grid, Tmax_grid;
	double 		T_grid, var[NUM_SES_VAR];
	double 		epsilon, delta;
	double 		T_min, T_max, pmin, pmax, pr;
	double 		lcoords[MAXD];
	double 		rho_grid;
		
	vp_inv_params.seos = seos;
	rho_grid = 1.0/ses_vp_vol_from_grid(vol_grid,seos);
	rho_grid = ses_rt_grid_from_rho(rho_grid,seos);
	(void) printf("%g\n",rho_grid);

	Tmin_grid = rt_fr->rect_grid->L[1];
	Tmax_grid = rt_fr->rect_grid->U[1];
	epsilon = seos->reler0*fabs(p_grid);
#if defined(PHASE_CODE)
	if (multiphase_eos(seos))
	{
		double T, E, S;

		if (get_phase_hyp_state(rho_grid,seos,&T,&E,&S))
		{
			E  = ses_re_grid_from_engy(E,seos);
			if (comp == COMP_MIXED_PHASE)	Tmax_grid = T;
			if (comp == COMP_PURE_PHASE)	Tmin_grid = T;
		}
	}
#endif /* defined(PHASE_CODE) */
	delta = seos->reler0*0.5*(Tmax_grid + Tmin_grid);

	vp_inv_params.vol_grid = coords[0];
	vp_inv_params.comp = comp;
	pr = ses_vp_press_from_grid(p_grid,seos);
	set_ses_intrp_flag(EVALUATE_PRESSURE,SESAME_RHO_TEMP);
	T_min =  ses_rt_temp_from_grid(Tmin_grid,seos);
	lcoords[0] = rho_grid;
	lcoords[1] = Tmin_grid;
	ses_solution(lcoords,comp,NULL,POSITIVE_SIDE,rt_fr,rt_wv,var);
	pmin = var[RT_CP] + ses_rt_rho_from_grid(rho_grid,seos)*
		T_min*var[RT_RP];
	T_max = ses_rt_temp_from_grid(Tmax_grid,seos);
	lcoords[1] = Tmax_grid;
	ses_solution(lcoords,comp,NULL,POSITIVE_SIDE,rt_fr,rt_wv,var);
	pmax = var[RT_CP] + ses_rt_rho_from_grid(rho_grid,seos)*
		T_max*var[RT_RP];
	if (pr <= pmin)
	{
	    if ((pmin - pr) <= (fabs(pr)*BIS_EPS + ABS_SES_EPS))
	        T_grid = Tmin_grid;
	    else
	    {
	        screen("ERROR in ses_VP_initializer(), ");
	        screen ("(vol,p) point not in window.\n");
		point_not_in_VP_window(ses_vp_vol_from_grid(vol_grid,seos),pr,
				       comp,T_min,T_max,pmin,pmax,seos);
	        clean_up (ERROR);
	    }
	}
	else if (pr >= pmax)
	{
	    if (pr - pmax <= fabs(pr)*BIS_EPS)
	        T_grid = Tmax_grid;
	    else
	    {
	        screen("ERROR in ses_VP_initializer(), ");
	        screen ("(vol,p) point not in window.\n");
		point_not_in_VP_window(ses_vp_vol_from_grid(vol_grid,seos),pr,
				       comp,T_min,T_max,pmin,pmax,seos);
	        clean_up (ERROR);
	    }
	}
	else
	{
	    if (find_root(VP_bisect,(POINTER)&vp_inv_params,p_grid,&T_grid,
			  Tmin_grid,Tmax_grid,epsilon,delta) == FUNCTION_FAILED) 
	    {
		if (DEBUG)
		{
	            (void) printf("WARNING in ses_VP_initializer(), "
	                          "Unable to invert function.\n");
		}
	    }
	}
	set_ses_intrp_flag_all(SESAME_RHO_TEMP);
	lcoords[0] = rho_grid;
	lcoords[1] = T_grid;
	ses_solution(lcoords,comp,NULL,POSITIVE_SIDE,rt_fr,rt_wv,var);
	ses_vp_colde(state) = var[RT_CE];
	ses_vp_rede(state) = var[RT_RE];
	ses_vp_Tvar(state) = T_grid;
	ses_vp_S(state) = var[RT_S];
	ses_vp_adb_gam(state) = var[RT_AG];
	ses_vp_gru_gam(state) = var[RT_GG];
	ses_vp_riv(state) = var[RT_RF];
}		/*end ses_VP_initializer*/


LOCAL boolean VP_bisect(
	double 		T_grid,
	double 		*p_grid,
	POINTER 	params)
{
	double 		T, rho, p, var[NUM_SES_VAR];
	Front 		*rt_fr = ((SESAME_VP_INVERT_PARAMS *) params)->fr;
	Wave 		*rt_wv = ((SESAME_VP_INVERT_PARAMS *) params)->wv;
	COMPONENT 	comp = ((SESAME_VP_INVERT_PARAMS *) params)->comp;
	SESAME_EOS	*seos = ((SESAME_VP_INVERT_PARAMS *) params)->seos;
	double 		coords[MAXD];

	coords[0] = ((SESAME_VP_INVERT_PARAMS *) params)->vol_grid;
	coords[0] = ses_vp_vol_from_grid(coords[0],seos);
	coords[0] = 1.0/coords[0];
	coords[0] = ses_rt_grid_from_rho(coords[0],seos);
	coords[1] = T_grid;

	ses_solution(coords,comp,NULL,NEGATIVE_SIDE,rt_fr,rt_wv,var);
	T = ses_rt_temp_from_grid(T_grid,seos);
	rho = ses_rt_rho_from_grid(coords[0],seos);
	p = var[RT_CP] + rho*T*var[RT_RP];
	*p_grid = ses_vp_grid_from_press(p,seos);
	return FUNCTION_SUCCEEDED;
}		/*end VP_bisect*/


/*ARGSUSED*/
LOCAL void ses_RS_initializer(
	double 		*coords,
	COMPONENT 	comp,
	Locstate 	state,
	SESAME_EOS 	*seos)
{
	double		ABS_SES_EPS = seos->ABS_SES_EPS;
	double		BIS_EPS = seos->BIS_EPS;
	SESAME_INVERT_PARAMS rs_inv_params;
	Front 		*fr = rs_inv_params.fr = seos->fr[SESAME_RHO_TEMP];
	Wave 		*wv = rs_inv_params.wv = seos->wave[SESAME_RHO_TEMP];
	double 		rho_grid = coords[0], S_grid = coords[1];
	double 		T_grid, var[NUM_SES_VAR];
	double 		T_min, T_max, Smin, Smax, S;
	double 		lcoords[MAXD];
	double 		Tmin_grid, Tmax_grid;
	double 		epsilon, delta;

	rs_inv_params.seos = seos;
	Tmin_grid = fr->rect_grid->L[1];
	Tmax_grid = fr->rect_grid->U[1];
	epsilon = seos->reler0*fabs(S_grid);

#if defined(PHASE_CODE)
	if (multiphase_eos(seos))
	{
		double T, E;
		if (get_phase_hyp_state(rho_grid,seos,&T,&E,&S))
		{
			if(S_grid < S)	Tmax_grid = T;
			if(S_grid >= S)	Tmin_grid = T;
		}
	}	
#endif /* defined(PHASE_CODE) */
	delta = seos->reler0*0.5*(Tmax_grid + Tmin_grid);
	rs_inv_params.rho_grid = lcoords[0] = rho_grid;
	rs_inv_params.comp = comp;
	S = ses_rs_entpy_from_grid(S_grid,seos);
	set_ses_intrp_flag(EVALUATE_ENTROPY,SESAME_RHO_TEMP);
	T_min = ses_rt_temp_from_grid(Tmin_grid,seos);
	lcoords[1] = Tmin_grid;
	ses_solution(lcoords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	Smin = var[RT_S];
	T_max = ses_rt_temp_from_grid(Tmax_grid,seos);
	lcoords[1] = Tmax_grid;
	ses_solution(lcoords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	Smax = var[RT_S];
	if (S <= Smin)
	{
	    if ((Smin - S) <= (fabs(S)*BIS_EPS + ABS_SES_EPS))
	        T_grid = Tmin_grid;
	    else
	    {
		screen("ERROR in ses_RS_initializer(), S < Smin.\n");
		point_not_in_RS_window(ses_rs_rho_from_grid(rho_grid,seos),
				       S,comp,T_min,T_max,Smin,Smax,seos);
	        clean_up (ERROR);
	    }
	}
	else if (S >= Smax)
	{
	    if ((S - Smax) <= fabs(S)*BIS_EPS)
	        T_grid = Tmax_grid;
	    else
	    {
	        screen("ERROR in ses_RS_initializer(), S > Smax.\n");
		point_not_in_RS_window(ses_rs_rho_from_grid(rho_grid,seos),
				       S,comp,T_min,T_max,Smin,Smax,seos);
	        clean_up (ERROR);
	    }
	}
	else
	{
	    if (find_root(RS_bisect,(POINTER)&rs_inv_params,S_grid,&T_grid,
			  Tmin_grid,Tmax_grid,epsilon,delta) == FUNCTION_FAILED) 
	    {
		if (DEBUG)
		{
	            (void) printf("WARNING in ses_RS_initializer(), "
	                          "Unable to invert function.\n");
		}
	    }
	}
	set_ses_intrp_flag_all(SESAME_RHO_TEMP);
	lcoords[1] = T_grid;
	ses_solution(lcoords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	ses_rs_coldp(state) = var[RT_CP];
	ses_rs_colde(state) = var[RT_CE];
	ses_rs_redp(state) = var[RT_RP];
	ses_rs_rede(state) = var[RT_RE];
	ses_rs_Tvar(state) = T_grid;
	ses_rs_adb_gam(state) = var[RT_AG];
	ses_rs_gru_gam(state) = var[RT_GG];
	ses_rs_riv(state) = var[RT_RF];
}		/*end ses_RS_initializer*/

LOCAL	void	point_not_in_RE_window(
	double		rho,
	double		en,
	COMPONENT	comp,
	double		T_min,
	double		T_max,
	double		emin,
	double		emax,
	SESAME_EOS	*seos)
{
	FILE	*file;
	char	s[80];;
	(void) printf("rho = %g, en = %g\n",rho,en);
	(void) printf("T_min = %g, e(rho,T_min) = %g\n",T_min,emin);
	(void) printf("T_max = %g, e(rho,T_max) = %g\n",T_max,emax);
	(void) printf("comp = %d\n",comp);
	print_sesame_interface(seos,seos->fr[SESAME_RHO_ENERGY]->interf,
			       SESAME_RHO_ENERGY);
	(void) sprintf(s,"rt%d.ts",seos->ids2);
	if ((file = fopen(s,"w")) != NULL)
	{
	    print_rt_tri_soln(file,seos);
	    trace_foutput(file);
	    (void) fclose(file);
	}
}		/*end point_not_in_RE_window*/

LOCAL	void	point_not_in_RS_window(
	double		rho,
	double		S,
	COMPONENT	comp,
	double		T_min,
	double		T_max,
	double		S_min,
	double		S_max,
	SESAME_EOS	*seos)
{
	FILE	*file;
	char	s[80];;
	(void) printf("rho=%g, S=%g\n",rho,S);
	(void) printf("T_min = %g, S(rho,T_min) = %g\n",T_min,S_min);
	(void) printf("T_max = %g, S(rho,T_max) = %g\n",T_max,S_max);
	(void) printf("comp = %d\n",comp);
	print_sesame_interface(seos,seos->fr[SESAME_RHO_ENTROPY]->interf,
			       SESAME_RHO_ENTROPY);
	(void) sprintf(s,"rt%d.ts",seos->ids2);
	if ((file = fopen(s,"w")) != NULL)
	{
	    print_rt_tri_soln(file,seos);
	    trace_foutput(file);
	    (void) fclose(file);
	}
}		/*end point_not_in_RS_window*/

LOCAL	void	point_not_in_VP_window(
	double		vol,
	double		pr,
	COMPONENT	comp,
	double		T_min,
	double		T_max,
	double		pmin,
	double		pmax,
	SESAME_EOS	*seos)
{
	FILE	*file;
	char	s[80];;

	(void) printf("vol = %g, pr = %g\n",vol,pr);
	(void) printf("T_min = %g, p(vol,T_min) = %g\n",T_min,pmin);
	(void) printf("T_max = %g, p(vol,T_max) = %g\n",T_max,pmax);
	(void) printf("comp = %d\n",comp);

	print_sesame_interface(seos,seos->fr[SESAME_VOLUME_PRESSURE]->interf,
			       SESAME_VOLUME_PRESSURE);
	(void) sprintf(s,"rt%d.ts",seos->ids2);
	if ((file = fopen(s,"w")) != NULL)
	{
	    print_rt_tri_soln(file,seos);
	    trace_foutput(file);
	    (void) fclose(file);
	}
}		/*end point_not_in_VP_window*/

LOCAL boolean RS_bisect(
	double 		T_grid,
	double 		*S_grid,
	POINTER 	params)
{
	double 		var[NUM_SES_VAR];
	Front 		*fr = ((SESAME_INVERT_PARAMS *) params)->fr;
	Wave 		*wv = ((SESAME_INVERT_PARAMS *) params)->wv;
	COMPONENT 	comp = ((SESAME_INVERT_PARAMS *) params)->comp;
	SESAME_EOS	*seos = ((SESAME_INVERT_PARAMS *) params)->seos;
	double 		coords[MAXD];

	coords[0] = ((SESAME_INVERT_PARAMS *) params)->rho_grid;
	coords[1] = T_grid;
	ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	*S_grid = ses_rs_grid_from_entpy(var[RT_S],seos);
	return FUNCTION_SUCCEEDED;
}		/*end RS_bisect*/

typedef struct {
	Front *fr;
	Wave *wv;
	COMPONENT comp;
	double S_grid;
} SES_PS_INVERT_PARAMS;

/*ARGSUSED*/
LOCAL void ses_PS_initializer(
	double 		*coords,
	COMPONENT 	comp,
	Locstate 	state,
	SESAME_EOS 	*seos)
{
	SES_PS_INVERT_PARAMS ps_inv_params;
	Front 		*fr = ps_inv_params.fr = seos->fr[SESAME_RHO_ENTROPY];
	Wave 		*wv = ps_inv_params.wv = seos->wave[SESAME_RHO_ENTROPY];
	double 		p_grid = coords[0];
	double		S = ses_ps_entpy_from_grid(coords[1],seos);
	double 		rho_grid, var[NUM_SES_VAR];
	double 		rhomin, rhomax, pmin, pmax, press;
	double 		T_min, T_max;
	double 		lcoords[MAXD];
	double 		epsilon, delta;
	double 		rmin_grid, rmax_grid;

	rmin_grid = fr->rect_grid->L[0];
	rmax_grid = fr->rect_grid->U[0];
	epsilon = seos->reler0*fabs(coords[0]);
	delta = seos->reler0*0.5*(fabs(rmin_grid) + fabs(rmax_grid));
#if defined(PHASE_CODE)
	if (multiphase_eos(seos))
	{
		double rho0_grid, rho1_grid;

		rho0_grid = rmin_grid;
		rho1_grid = rmax_grid;
		get_phase_ent_state(S,p_grid,comp,
			seos->fr[SESAME_RHO_ENTROPY],
			seos->wave[SESAME_RHO_ENTROPY],
			&rho0_grid,&rho1_grid);
		rmin_grid = rho0_grid;
		rmax_grid = rho1_grid;
		delta = seos->reler0*fabs(rmin_grid);
	}
#endif /* defined(PHASE_CODE) */
	lcoords[0] = rmin_grid;
	ps_inv_params.S_grid = lcoords[1] = ses_rs_grid_from_entpy(S,seos);
	ps_inv_params.comp = comp;
	set_ses_intrp_flag(EVALUATE_PRESSURE,SESAME_RHO_ENTROPY);
	ses_solution(lcoords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	press = ses_ps_press_from_grid(p_grid,seos);
	rhomin = ses_rs_rho_from_grid(rmin_grid,seos);
	T_min = ses_rs_temp_from_var(var[RS_T],seos);
	pmin = rhomin*T_min*var[RS_RP] + var[RS_CP];
	lcoords[0] = rmax_grid;
	ses_solution(lcoords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	rhomax = ses_rs_rho_from_grid(rmax_grid,seos);
	T_max = ses_rs_temp_from_var(var[RS_T],seos);
	pmax = rhomax*T_max*var[RS_RP] + var[RS_CP];
	if(press <= pmin)
	{
	    rho_grid = rmin_grid;
	}
	else if(press >= pmax)
	{
	    rho_grid = rmax_grid;
	}
	else
	{
	    if (find_root(PS_bisect,(POINTER)&ps_inv_params,p_grid,&rho_grid,
			  rmin_grid,rmax_grid,epsilon,delta) == FUNCTION_FAILED) 
	    {
	        double rho_last, T_last, p_last;

	        lcoords[0] = rho_grid;
	        ses_solution(lcoords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	        rho_last = ses_rs_rho_from_grid(rho_grid,seos);
	        T_last = ses_rs_temp_from_var(var[RS_T],seos);
	        p_last = rho_last*T_last*var[RS_RP] +var[RS_CP];
		if (DEBUG)
		{
	            (void) printf("WARNING in ses_PS_initializer(), "
	                          "Unable to invert function.\n");
	            (void) printf("P = %g,S = %g\n",press,
	               	          ses_ps_entpy_from_grid(S,seos));
	            (void) printf("epsilon = %g, delta = %g\n",epsilon,delta);
	            (void) printf("rmax = %g, P(rmax,S) = %g, Tmax = %g\n",
				  rhomax,pmax,T_max);
	            (void) printf("rmin = %g, P(rmin,S) = %g, Tmin = %g\n",
				  rhomin,pmin,T_min);
	            (void) printf("last iter: rho = %g, P = %g, T = %g\n",
	                          rho_last,p_last,T_last);
		 }
	    }
	}
	set_ses_intrp_flag_all(SESAME_RHO_ENTROPY);
	lcoords[0] = rho_grid;
	ses_solution(lcoords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	ses_ps_colde(state) = var[RS_CE];
	ses_ps_rede(state) = var[RS_RE];
	ses_ps_Tvar(state) = var[RS_T];
	ses_ps_rho_var(state) = rho_grid;
	ses_ps_adb_gam(state) = var[RS_AG];
	ses_ps_gru_gam(state) = var[RS_GG];
	ses_ps_riv(state) = var[RS_RF];
}		/*end ses_PS_initializer*/

#if defined(PHASE_CODE)

LOCAL	void	set_up_multiphase_rt_table(
	Front 		*fr,
	Wave 		*wave,
	SESAME_EOS 	*seos,
	POINTER 	pbdry,
	POINTER 	ccur)
{
	BDRY_SIDE  side;
	PHASE_BDRY *phase_bound = (PHASE_BDRY *) pbdry;
	COLD_CURVE *cold_curve = (COLD_CURVE *) ccur;
	CURVE	   **c;
	NODE	   **n;
	RECT_GRID  *gr = fr->rect_grid;
	boolean	   sav_intrp = interpolate_intfc_states(fr->interf);
	double 	   *tbls = seos->sestab.tbls;
	int	   nr = (int)(tbls[2]), nt = (int)(tbls[3]);

	DEBUG_ENTER(set_up_multiphase_rt_table)
	interpolate_intfc_states(fr->interf) = NO;

	seos->de_params.len = max(nt,nr)+4;
	uni_array(&seos->de_params.R,seos->de_params.len,FLOAT);
	uni_array(&seos->de_params.E,seos->de_params.len,FLOAT);
	uni_array(&seos->de_params.P,seos->de_params.len,FLOAT);
	uni_array(&seos->de_params.slopeE,seos->de_params.len,FLOAT);
	uni_array(&seos->de_params.slopeP,seos->de_params.len,FLOAT);
	uni_array(&seos->de_params.T,seos->de_params.len,FLOAT);
		/* Initialize the phase boundary */

	init_new_phase_bound(wave,fr,seos,cold_curve,phase_bound);
	if(debugging("ses_print_hyp"))
	    verbose_ses_show_intfc_states(fr->interf,SESAME_RHO_TEMP,seos);

	ses_phase_states(fr,seos,phase_bound);

	/* Set the boundary */

	if (DEBUG)
	    (void) printf("calling set_boundary\n");
	if (set_boundary(fr->interf,gr,COMP_PURE_PHASE,grid_tolerance(gr)) != FUNCTION_SUCCEEDED)
	{
	    screen("ERROR in set_up_multiphase_rt_table(), "
		   "set_boundary failed\n");
	    clean_up(ERROR);
	}
	for (c = fr->interf->curves; c && *c; c++)
	{
	    if (is_bdry(*c) && (wave_type(*c) == UNKNOWN_WAVE_TYPE))
	    {
		side = rect_bdry_side_for_curve(NULL,NULL,*c,gr);
		switch (side)
		{
		case LEFT_BDRY:
		case RIGHT_BDRY:
		    wave_type(*c) = ISOCHOR_BOUNDARY;
		    break;
	        case LOWER_BDRY:
		case UPPER_BDRY:
		    wave_type(*c) = ISOTHERM_BOUNDARY;
		    break;
		default:
		    screen("ERROR in set_up_multiphase_rt_table(), "
		           "unsupported boundary side %d\n",side);
		    clean_up(ERROR);
		}
	    }
	}
	for (n = fr->interf->nodes; n && *n; n++)
	    if (is_bdry(*n) && (node_type(*n) == UNKNOWN_NODE_TYPE))
	        node_type(*n) = FIXED_NODE;
	set_boundary_states(fr,seos,phase_bound,cold_curve);

	if(debugging("ses_print_hyp"))
	    verbose_ses_show_intfc_states(fr->interf,SESAME_RHO_TEMP,seos);

	/* Set the interior states */

	init_RT_interior_states(wave,fr,seos,cold_curve,phase_bound);

	free_these(6,seos->de_params.R,seos->de_params.E,seos->de_params.P,
		     seos->de_params.T,seos->de_params.slopeE,
		     seos->de_params.slopeP);
	seos->de_params.R = NULL;
	seos->de_params.E = NULL;
	seos->de_params.P = NULL;
	seos->de_params.T = NULL;
	seos->de_params.slopeE = NULL;
	seos->de_params.slopeP = NULL;;
	seos->de_params.len = 0;
	interpolate_intfc_states(fr->interf) = sav_intrp;
	DEBUG_LEAVE(set_up_multiphase_rt_table)
}		/*end set_up_multiphase_rt_table*/


LOCAL	void get_phase_ent_state(
	double 		S_grid,
	double 		p_grid,
	COMPONENT 	comp,
	Front 		*fr,
	Wave 		*wv,
	double 		*rmin_grid,
	double 		*rmax_grid)
{
	CURVE 		**cur;
	BOND 		*b;
	INTERFACE 	*intfc = fr->interf;
	SESAME_EOS	*seos = Ses_front_seos(fr);
	double 		tmp;
	double 		p1, p1_grid, T1, rho1;
	double 		var[NUM_SES_VAR];
	double 		coords[MAXD];
	double 		r, rend;
	int 		flag[2];


	r  = *rmin_grid;	rend = *rmax_grid;
	flag[1] = NO; flag[0] = NO;
	coords[1] = S_grid;
	for (cur = intfc->curves; cur && *cur; cur++)
	{
	    if(wave_type(*cur) == PHASE_BOUNDARY)
	    {
	        for (b = (*cur)->first; b != NULL; b = b->next) 
	        {
	            if (!Between(S_grid,Coords(b->start)[1],
	                           Coords(b->end)[1]))
	                continue;

	            tmp = (S_grid - Coords(b->start)[1]) /
	                 (Coords(b->start)[1] - Coords(b->end)[1]);
	            coords[0] = (Coords(b->start)[0] - 
	                     Coords(b->end)[0])*tmp + 
	                         Coords(b->start)[0];
	            rho1 = ses_rs_rho_from_grid(coords[0],seos);
	            set_ses_intrp_flag(EVALUATE_PRESSURE,
	                SESAME_RHO_ENTROPY);
	            ses_solution(coords,comp,Hyper_surf(*cur),
				 POSITIVE_SIDE,fr,wv,var);
	            T1 = ses_rs_temp_from_var(var[RS_T],seos);
	            p1 = var[RS_CP] + rho1*T1*var[RS_RP];
	            p1_grid = ses_ps_grid_from_press(p1,seos);
	            if (p_grid < p1_grid)
	            {
	                if (flag[1] == NO || 
	                (flag[1] == YES && coords[0] < rend))
	                {
	                    rend = coords[0];
	                    flag[1] = YES;
	                }
	            }
	            else if (p_grid >= p1_grid)
	            {
	                if (flag[0] == NO || 
	                (flag[0] == YES && coords[0] > r))
	                {
	                    r = coords[0];
	                    flag[0] = YES;
	                }
	            }
	        }
	    }
	}
	if (flag[0] == YES)
	{
	    if (comp == COMP_MIXED_PHASE)	*rmin_grid = r;
	    if (comp == COMP_PURE_PHASE)	*rmin_grid = r;
	}
	if (flag[1] == YES)
	{
	    if (comp == COMP_MIXED_PHASE)	*rmax_grid = rend;
	    if (comp == COMP_PURE_PHASE)	*rmax_grid = rend;
	}
}		/*end get_phase_ent_states*/
#endif /* defined(PHASE_CODE) */

LOCAL boolean PS_bisect(
	double 		rho_grid,
	double 		*p_grid,
	POINTER 	params)
{
	double 		p, rho, T, var[NUM_SES_VAR];
	Front 		*rs_fr = ((SES_PS_INVERT_PARAMS *) params)->fr;
	Wave 		*rs_wv = ((SES_PS_INVERT_PARAMS *) params)->wv;
	SESAME_EOS	*seos = Ses_front_seos(rs_fr);
	COMPONENT 	comp = ((SES_PS_INVERT_PARAMS *) params)->comp;
	double 		coords[MAXD];

	coords[0] = rho_grid;
	coords[1] = ((SES_PS_INVERT_PARAMS *) params)->S_grid;
	ses_solution(coords,comp,NULL,POSITIVE_SIDE,rs_fr,rs_wv,var);
	rho = ses_rs_rho_from_grid(rho_grid,seos);
	T = ses_rs_temp_from_var(var[RS_T],seos);
	p = rho*T*var[RS_RP] + var[RS_CP];
	*p_grid = ses_ps_grid_from_press(p,seos);
	return FUNCTION_SUCCEEDED;
}		/*end PS_bisect*/

LOCAL void set_rho_T_state(
	Locstate 	state,
	double 		rho_grid,
	double 		T_grid,
	SESAME_EOS 	*seos)
{
	double 		rho, T;

	rho = ses_rt_rho_from_grid(rho_grid,seos);
	T = ses_rt_temp_from_grid(T_grid,seos);
	setrt(rho,T,state,seos);
}		/*end set_rho_T_state*/

/*ARGSUSED*/
LOCAL void map_RT_to_RE(
	POINT 		*op,
	BOND 		*ob,
	CURVE 		*oldc,
	POINT 		*np,
	BOND 		*nb,
	CURVE 		*newc,
	boolean		first_call,
	RECT_GRID 	*rect_grid,
	POINTER 	params)
{
	SESAME_EOS	*seos = (SESAME_EOS*)params;
	Locstate 	osr, osl, nsl, nsr;
	double 		temp, el, er;

	osr = right_state_at_point_on_curve(op,ob,oldc);
	osl = left_state_at_point_on_curve(op,ob,oldc);
	nsr = right_state_at_point_on_curve(np,nb,newc);
	nsl = left_state_at_point_on_curve(np,nb,newc);

	temp = ses_rt_temp_from_grid(Coords(op)[1],seos);
	Coords(np)[0] = Coords(op)[0];
	er = ses_rt_colde(osr) + temp*ses_rt_rede(osr);
	el = ses_rt_colde(osl) + temp*ses_rt_rede(osl);
	Coords(np)[1] = 0.5*(ses_re_grid_from_engy(el,seos) +
			     ses_re_grid_from_engy(er,seos));

	/* Left States */
	ses_re_coldp(nsl) = ses_rt_coldp(osl);
	ses_re_redp(nsl) = ses_rt_redp(osl);
	ses_re_Tvar(nsl) = Coords(op)[1];
	ses_re_S(nsl) = ses_rt_S(osl);
	ses_re_adb_gam(nsl) = ses_rt_adb_gam(osl);
	ses_re_gru_gam(nsl) = ses_rt_gru_gam(osl);
	ses_re_riv(nsl) = ses_rt_riv(osl);

	/* Right States */
	ses_re_coldp(nsr) = ses_rt_coldp(osr);
	ses_re_redp(nsr) = ses_rt_redp(osr);
	ses_re_Tvar(nsr) = Coords(op)[1];
	ses_re_S(nsr) = ses_rt_S(osr);
	ses_re_adb_gam(nsr) = ses_rt_adb_gam(osr);
	ses_re_gru_gam(nsr) = ses_rt_gru_gam(osr);
	ses_re_riv(nsr) = ses_rt_riv(osr);

	if (first_call == YES)
	{
	    rect_grid->U[0] = rect_grid->L[0] = Coords(np)[0];
	    rect_grid->U[1] = rect_grid->L[1] = Coords(np)[1];
	}
	else
	{
	    rect_grid->U[0] = max(rect_grid->U[0],Coords(np)[0]);
	    rect_grid->L[0] = min(rect_grid->L[0],Coords(np)[0]);
	    rect_grid->U[1] = max(rect_grid->U[1],Coords(np)[1]);
	    rect_grid->L[1] = min(rect_grid->L[1],Coords(np)[1]);
	}
}		/*end map_RT_to_RE*/

/*ARGSUSED*/
LOCAL void map_RT_to_RS(
	POINT 		*op,
	BOND 		*ob,
	CURVE 		*oldc,
	POINT 		*np,
	BOND 		*nb,
	CURVE 		*newc,
	boolean		first_call,
	RECT_GRID 	*rect_grid,
	POINTER 	params)
{
	SESAME_EOS	*seos = (SESAME_EOS*)params;
	Locstate	osr, osl, nsl, nsr;

	osr = right_state_at_point_on_curve(op,ob,oldc);
	osl = left_state_at_point_on_curve(op,ob,oldc);
	nsr = right_state_at_point_on_curve(np,nb,newc);
	nsl = left_state_at_point_on_curve(np,nb,newc);

	Coords(np)[0] = Coords(op)[0];
	Coords(np)[1] = 0.5*(ses_rs_grid_from_entpy(ses_rt_S(osl),seos) +
			     ses_rs_grid_from_entpy(ses_rt_S(osr),seos));

	/* Set left states */
	ses_rs_coldp(nsl) = ses_rt_coldp(osl);
	ses_rs_colde(nsl) = ses_rt_colde(osl);
	ses_rs_redp(nsl) = ses_rt_redp(osl);
	ses_rs_rede(nsl) = ses_rt_rede(osl);
	ses_rs_Tvar(nsl) = Coords(op)[1];
	ses_rs_adb_gam(nsl) = ses_rt_adb_gam(osl);
	ses_rs_gru_gam(nsl) = ses_rt_gru_gam(osl);
	ses_rs_riv(nsl) = ses_rt_riv(osl);

	/* Set right states */
	ses_rs_coldp(nsr) = ses_rt_coldp(osr);
	ses_rs_colde(nsr) = ses_rt_colde(osr);
	ses_rs_redp(nsr) = ses_rt_redp(osr);
	ses_rs_rede(nsr) = ses_rt_rede(osr);
	ses_rs_Tvar(nsr) = Coords(op)[1];
	ses_rs_adb_gam(nsr) = ses_rt_adb_gam(osr);
	ses_rs_gru_gam(nsr) = ses_rt_gru_gam(osr);
	ses_rs_riv(nsr) = ses_rt_riv(osr);

	if (first_call == YES)
	{
	    rect_grid->U[0] = rect_grid->L[0] = Coords(np)[0];
	    rect_grid->U[1] = rect_grid->L[1] = Coords(np)[1];
	}
	else
	{
	    rect_grid->U[0] = max(rect_grid->U[0],Coords(np)[0]);
	    rect_grid->L[0] = min(rect_grid->L[0],Coords(np)[0]);
	    rect_grid->U[1] = max(rect_grid->U[1],Coords(np)[1]);
	    rect_grid->L[1] = min(rect_grid->L[1],Coords(np)[1]);
	}
}		/*end map_RT_to_RS*/


/*ARGSUSED*/
LOCAL void map_RT_to_VP(
	POINT 		*op,
	BOND 		*ob,
	CURVE 		*oldc,
	POINT 		*np,
	BOND 		*nb,
	CURVE 		*newc,
	boolean 	first_call,
	RECT_GRID 	*rect_grid,
	POINTER 	params)
{
	SESAME_EOS	*seos = (SESAME_EOS*)params;
	Locstate 	osr, osl, nsl, nsr;
	double 		temp, pl, pr;
	double 		V, rho;

	osr = right_state_at_point_on_curve(op,ob,oldc);
	osl = left_state_at_point_on_curve(op,ob,oldc);
	nsr = right_state_at_point_on_curve(np,nb,newc);
	nsl = left_state_at_point_on_curve(np,nb,newc);
	negative_component(newc) = positive_component(oldc);
	positive_component(newc) = negative_component(oldc);

	temp = ses_rt_temp_from_grid(Coords(op)[1],seos);
	rho = ses_rt_rho_from_grid(Coords(op)[0],seos);
	V = 1.0/ses_rt_rho_from_grid(Coords(op)[0],seos);
	Coords(np)[0] = ses_vp_grid_from_vol(V,seos);
	pr = ses_rt_coldp(osr) + rho*temp*ses_rt_redp(osr);
	pl = ses_rt_coldp(osl) + rho*temp*ses_rt_redp(osl);
	Coords(np)[1] = 0.5*(ses_vp_grid_from_press(pl,seos) + 
			     ses_vp_grid_from_press(pr,seos));

	/* Left States */
	ses_vp_colde(nsl) = ses_rt_colde(osr);
	ses_vp_rede(nsl) = ses_rt_rede(osr);
	ses_vp_Tvar(nsl) = Coords(op)[1];
	ses_vp_S(nsl) = ses_rt_S(osr);
	ses_vp_adb_gam(nsl) = ses_rt_adb_gam(osr);
	ses_vp_gru_gam(nsl) = ses_rt_gru_gam(osr);
	ses_vp_riv(nsl) = ses_rt_riv(osr);

	/* Right States */
	ses_vp_colde(nsr) = ses_rt_colde(osl);
	ses_vp_rede(nsr) = ses_rt_rede(osl);
	ses_vp_Tvar(nsr) = Coords(op)[1];
	ses_vp_S(nsr) = ses_rt_S(osl);
	ses_vp_adb_gam(nsr) = ses_rt_adb_gam(osl);
	ses_vp_gru_gam(nsr) = ses_rt_gru_gam(osl);
	ses_vp_riv(nsr) = ses_rt_riv(osl);


	if (first_call == YES)
	{
	    rect_grid->U[0] = rect_grid->L[0] = Coords(np)[0];
	    rect_grid->U[1] = rect_grid->L[1] = Coords(np)[1];
	}
	else
	{
	    rect_grid->U[0] = max(rect_grid->U[0],Coords(np)[0]);
	    rect_grid->L[0] = min(rect_grid->L[0],Coords(np)[0]);
	    rect_grid->U[1] = max(rect_grid->U[1],Coords(np)[1]);
	    rect_grid->L[1] = min(rect_grid->L[1],Coords(np)[1]);
	}
}		/*end map_RT_to_VP*/

/*ARGSUSED*/
LOCAL void map_RT_to_PS(
	POINT 		*op,
	BOND 		*ob,
	CURVE 		*oldc,
	POINT 		*np,
	BOND 		*nb,
	CURVE 		*newc,
	boolean 	first_call,
	RECT_GRID 	*rect_grid,
	POINTER 	params)
{
	SESAME_EOS	*seos = (SESAME_EOS*)params;
	Locstate 	osr, osl, nsl, nsr;
	double 		rho, temp, p;

	osr = right_state_at_point_on_curve(op,ob,oldc);
	osl = left_state_at_point_on_curve(op,ob,oldc);
	nsr = right_state_at_point_on_curve(np,nb,newc);
	nsl = left_state_at_point_on_curve(np,nb,newc);

	rho = ses_rt_rho_from_grid(Coords(op)[0],seos);
	temp = ses_rt_temp_from_grid(Coords(op)[1],seos);

	p = 0.5*(ses_rt_coldp(osr) + rho*temp*ses_rt_redp(osr) +
			ses_rt_coldp(osl) + rho*temp*ses_rt_redp(osl));
	Coords(np)[0] = ses_ps_grid_from_press(p,seos);
	Coords(np)[1] = 0.5*(ses_ps_grid_from_entpy(ses_rt_S(osl),seos) +
			     ses_ps_grid_from_entpy(ses_rt_S(osr),seos));

	/* Set left state */
	ses_ps_colde(nsl) = ses_rt_colde(osl);
	ses_ps_rede(nsl) = ses_rt_rede(osl);
	ses_ps_Tvar(nsl) = Coords(op)[1];
	ses_ps_rho_var(nsl) = Coords(op)[0];
	ses_ps_adb_gam(nsl) = ses_rt_adb_gam(osl);
	ses_ps_gru_gam(nsl) = ses_rt_gru_gam(osl);
	ses_ps_riv(nsl) = ses_rt_riv(osl);

	/* Set right state */
	ses_ps_colde(nsr) = ses_rt_colde(osr);
	ses_ps_rede(nsr) = ses_rt_rede(osr);
	ses_ps_Tvar(nsr) = Coords(op)[1];
	ses_ps_rho_var(nsr) = Coords(op)[0];
	ses_ps_adb_gam(nsr) = ses_rt_adb_gam(osr);
	ses_ps_gru_gam(nsr) = ses_rt_gru_gam(osr);
	ses_ps_riv(nsr) = ses_rt_riv(osr);

	if (first_call == YES)
	{
	    rect_grid->U[0] = rect_grid->L[0] = Coords(np)[0];
	    rect_grid->U[1] = rect_grid->L[1] = Coords(np)[1];
	}
	else
	{
	    rect_grid->U[0] = max(rect_grid->U[0],Coords(np)[0]);
	    rect_grid->L[0] = min(rect_grid->L[0],Coords(np)[0]);
	    rect_grid->U[1] = max(rect_grid->U[1],Coords(np)[1]);
	    rect_grid->L[1] = min(rect_grid->L[1],Coords(np)[1]);
	}
}		/*end map_RT_to_PS*/

EXPORT	void set_default_ses_wave_and_front(
	INIT_DATA *init,
	Front 	  *fr,
	Wave 	  *wave,
	size_t 	  sizest,
	boolean   set_sesame_hooks)
{
	DEBUG_ENTER(set_default_ses_wave_and_front)

	if (set_sesame_hooks == YES)
	{
	    g_preserve_user_hooks(2,SAVE_HOOKS);
	    set_user_hooks_for_sesame();
	}

	/* These are needed by f_set_default_front_parameters(). */
	fr->rect_grid->dim = 2; /* The sesame tables are inherently 2d. */
	fr->sizest = sizest;

	f_set_default_front_parameters(NULL,fr);
	fr->_state_interpolator = ses_lin_comb_states;
	fr->_tri_state_interpolator = ses_tri_lin_comb_states;
	fr->_replace_unphys_loop = ses_replace_unphys_loop;
	fr->_hyp_solution = s_ses_solution;
	Clear_redistribution_parameters(fr);
	set_dflt_cur_redist_params(fr);
	Frequency_of_redistribution(fr,GENERAL_NODE) = 1;
	Frequency_of_redistribution(fr,GENERAL_WAVE) = 1;
	Frequency_of_redistribution(fr,VECTOR_WAVE) = -1;
	Front_spacing(fr,GENERAL_WAVE) = Front_spacing(fr,VECTOR_WAVE) = 0.5;
	Cosine_big_angle(fr,GENERAL_WAVE) = -1.0;
	Cosine_big_angle(fr,VECTOR_WAVE) = -1.0;
	fr->nfloats = (int) (sizest/FLOAT);

	/* initialize Wave */

	wave->rect_grid = fr->rect_grid;
	wave->sizest = sizest;
	wave->nfloats = fr->nfloats;

	h_set_default_wave_parameters(init,wave);
	wave->print_state = NULL;
	wave->show_wave_states = NULL;
	wave->interpolator.linear_cell = ses_tri_interpolator;
	wave->interpolator.bilinear_cell = ses_quad_interpolator;
	wave->interpolator.grad_linear_cell = NULL;
	wave->interpolator.grad_bond = NULL;
	wave->interpolator.grad_bilinear_cell = NULL;
	wave->el_integral.bilinear_cell = NULL;
	wave->el_integral.linear_cell = NULL;
	wave->unsplit.flux = NULL;
	wave->unsplit.flux_obl = NULL;
	wave->unsplit.sources = NULL;
	wave->max_wave_speed = NULL;
	wave->Tri_grid_hooks._construct_tri_grid = t_construct_tri_grid;
	wave->Tri_grid_hooks._set_components = set_components2d;
	wave->Tri_grid_hooks._triangulate_mesh = triangulate_mesh2d;
	wave->Tri_grid_hooks._method = COMPLETE_TRIANGULATION;
	wave->Tri_grid_hooks._pt_in_lin_el = point_in_triangle;
	wave->Tri_grid_hooks._blk_triangulate = exact_blk_triangulate;

	if (set_sesame_hooks == YES)
	    g_preserve_user_hooks(2,RESTORE_HOOKS);

	DEBUG_LEAVE(set_default_ses_wave_and_front)
}		/*end set_default_ses_wave_and_front*/

LOCAL	boolean	ses_replace_unphys_loop(
	NNLIST	*nl,
	NNLIST	**new_node_list,
	CURVE	**newc,
	Front	*fr,
	int	i,
	double	min_area,
	int	flag)
{
	NODE	*m;
	boolean status;

	m = nl->m;
	status = f_replace_unphys_loop(nl,new_node_list,
				       newc,fr,i,min_area,flag);
	if (status == FUNCTION_FAILED)
	    return status;

	if (
	    (m->in_curves != NULL) &&
	    (m->in_curves[0] != NULL) &&
	    (m->in_curves[1] == NULL) &&
	    (m->out_curves != NULL) &&
	    (m->out_curves[0] != NULL) &&
	    (m->out_curves[1] == NULL) &&
	    ((node_type(m) == ERROR) || (node_type(m) == UNKNOWN_NODE_TYPE)) &&
	    (wave_type(m->in_curves[0]) != wave_type(m->out_curves[0]))
	)
	{
		node_type(m) = EOS_BOUNDARY_NODE;
	}
	return status;
}		/*end ses_replace_unphys_loop*/

/*ARGSUSED*/
LOCAL	void set_RE_rect_grid(
	INTERFACE 	*rt_intfc,
	INTERFACE 	*re_intfc,
	void		(*remap)(POINT*,BOND*,CURVE*,POINT*,BOND*,CURVE*,
				 boolean,RECT_GRID*,POINTER),
	POINTER 	params)
{
	SESAME_EOS *seos = (SESAME_EOS*)params;
	RECT_GRID  *tgr = &topological_grid(re_intfc);
	RECT_GRID  *cgr = computational_grid(re_intfc);
	RECT_GRID  *rtgr = seos->fr[SESAME_RHO_TEMP]->rect_grid;
	double      *U, *L, *rtU, *rtL;
	int        *gmax, *rtgmax;

	DEBUG_ENTER(set_RE_rect_grid)

	copy_rect_grid(cgr,tgr);
	U = cgr->U;
	L = cgr->L;
	gmax = cgr->gmax;
	rtU = rtgr->U;
	rtL = rtgr->L;
	rtgmax = rtgr->gmax;
	if (seos->restart_io_type == NULL)
	{
	    char s[Gets_BUF_SIZE];
	    screen("Limits for the density-energy table\n");
	    screen("\tminimum density = %g g/cc, "
	           "minimum density grid = %g\n",
	           ses_re_rho_from_grid(L[0],seos),L[0]);
	    screen("\tmaximum density = %g g/cc, "
	           "maximum density grid = %g\n",
	           ses_re_rho_from_grid(U[0],seos),U[0]);
	    screen("\tminimum energy = %g kJ/g, "
	           "minimum energy grid = %g\n",
	           ses_re_engy_from_grid(L[1],seos),L[1]);
	    screen("\tmaximum energy = %g kJ/g, "
	           "maximum energy grid = %g\n",
	           ses_re_engy_from_grid(U[1],seos),U[1]);
	    gmax[0] = rtgmax[0];
	    gmax[1] = (int)(rtgmax[1]*((U[1] - L[1])/(rtU[1] - rtL[1])));
	    set_rect_grid(L,U,L,U,NOBUF,NOBUF,gmax,cgr->dim,&rtgr->Remap,cgr);
	    (void) adjust_top_grid_for_square(cgr,cgr);
	    screen("Enter the number of grid zones\n\t");
	    screen("for density and energy (dflt = %d %d): ",gmax[0],gmax[1]);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscanf(s,"%d %d",gmax,gmax+1);
	}
	else
	{
	    FILE *restart_file = seos->restart_io_type->file;
	    (void) fgetstring(restart_file,"density mesh = ");
	    (void) fscanf(restart_file,"%d%*s%*s%d",gmax,gmax+1);
	}
	cgr->h[0] = (U[0] - L[0])/gmax[0];
	cgr->h[1] = (U[1] - L[1])/gmax[1];
	copy_rect_grid(tgr,cgr);
	(void) adjust_top_grid_for_square(tgr,cgr);
	re_intfc->table->fixed_grid = YES;
	DEBUG_LEAVE(set_RE_rect_grid)
}		/*end set_RE_rect_grid*/

/*ARGSUSED*/
LOCAL	void set_RS_rect_grid(
	INTERFACE 	*rt_intfc,
	INTERFACE 	*rs_intfc,
	void		(*remap)(POINT*,BOND*,CURVE*,POINT*,BOND*,CURVE*,
				 boolean,RECT_GRID*,POINTER),
	POINTER		params)
{
	SESAME_EOS *seos = (SESAME_EOS*)params;
	RECT_GRID  *tgr = &topological_grid(rs_intfc);
	RECT_GRID  *cgr = computational_grid(rs_intfc);
	RECT_GRID  *rtgr = seos->fr[SESAME_RHO_TEMP]->rect_grid;
	double      *U, *L, *rtU, *rtL;
	int        *gmax, *rtgmax;

	DEBUG_ENTER(set_RS_rect_grid)

	copy_rect_grid(cgr,tgr);
	U = cgr->U;
	L = cgr->L;
	gmax = cgr->gmax;
	rtU = rtgr->U;
	rtL = rtgr->L;
	rtgmax = rtgr->gmax;
	if (seos->restart_io_type == NULL)
	{
	    char s[Gets_BUF_SIZE];
	    screen("Limits for the density-entropy table\n");
	    screen("\tminimum density = %g g/cc, "
	           "minimum density grid = %g\n",
	           ses_rs_rho_from_grid(L[0],seos),L[0]);
	    screen("\tmaximum density = %g g/cc, "
	           "maximum density grid = %g\n",
	           ses_rs_rho_from_grid(U[0],seos),U[0]);
	    screen("\tminimum entropy = %g, minimum entropy grid = %g\n",
	    	ses_rs_entpy_from_grid(L[1],seos),L[1]);
	    screen("\tmaximum entropy = %g, "
	           "maximum entropy grid = %g\n",
	           ses_rs_entpy_from_grid(U[1],seos),U[1]);
	    gmax[0] = rtgmax[0];
	    gmax[1] = (int)(rtgmax[1]*((U[1] - L[1])/(rtU[1] - rtL[1])));
	    set_rect_grid(L,U,L,U,NOBUF,NOBUF,gmax,cgr->dim,&rtgr->Remap,cgr);
	    (void) adjust_top_grid_for_square(cgr,cgr);
	    screen("Enter the number of grid zones\n\t");
	    screen("for density and entropy (dflt = %d %d): ",gmax[0],gmax[1]);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscanf(s,"%d %d",gmax,gmax+1);
	}
	else
	{
	    FILE *restart_file = seos->restart_io_type->file;
	    (void) fgetstring(restart_file,"density mesh = ");
	    (void) fscanf(restart_file,"%d%*s%*s%d",gmax,gmax+1);
	}
	cgr->h[0] = (U[0] - L[0])/gmax[0];
	cgr->h[1] = (U[1] - L[1])/gmax[1];
	copy_rect_grid(tgr,cgr);
	(void) adjust_top_grid_for_square(tgr,cgr);
	rs_intfc->table->fixed_grid = YES;
	DEBUG_LEAVE(set_RS_rect_grid)
}		/*end set_RS_rect_grid*/

/*ARGSUSED*/
LOCAL	void set_PS_rect_grid(
	INTERFACE 	*rt_intfc,
	INTERFACE 	*ps_intfc,
	void		(*remap)(POINT*,BOND*,CURVE*,POINT*,BOND*,CURVE*,
				 boolean,RECT_GRID*,POINTER),
	POINTER		params)
{
	SESAME_EOS *seos = (SESAME_EOS*)params;
	RECT_GRID  *tgr = &topological_grid(ps_intfc);
	RECT_GRID  *cgr = computational_grid(ps_intfc);
	RECT_GRID  *rtgr = seos->fr[SESAME_RHO_TEMP]->rect_grid;
	double      *U, *L, *rtU, *rtL;
	int        *gmax, *rtgmax;

	DEBUG_ENTER(set_PS_rect_grid)
	copy_rect_grid(cgr,tgr);
	U = cgr->U;
	L = cgr->L;
	gmax = cgr->gmax;
	rtU = rtgr->U;
	rtL = rtgr->L;
	rtgmax = rtgr->gmax;
	if (seos->restart_io_type == NULL)
	{
	    char s[Gets_BUF_SIZE];
	    screen("Limits for the pressure-entropy table\n");
	    screen("\tminimum pressure = %g, minimum pressure grid = %g\n",
	    	ses_ps_press_from_grid(L[0],seos),L[0]);
	    screen("\tmaximum pressure = %g, maximum pressure grid = %g\n",
	    	ses_ps_press_from_grid(U[0],seos),U[0]);
	    screen("\tminimum entropy = %g, minimum entropy grid = %g\n",
	    	ses_ps_entpy_from_grid(L[1],seos),L[1]);
	    screen("\tmaximum entropy = %g, maximum entropy grid = %g\n",
	    	ses_ps_entpy_from_grid(U[1],seos),U[1]);
	    gmax[0] = (int)(rtgmax[0]*((U[0]-L[0])/(rtU[0]-rtL[0])));
	    gmax[1] = (int)(rtgmax[1]*((U[1]-L[1])/(rtU[1]-rtL[1])));
	    set_rect_grid(L,U,L,U,NOBUF,NOBUF,gmax,cgr->dim,&rtgr->Remap,cgr);
	    (void) adjust_top_grid_for_square(cgr,cgr);
	    screen("Enter the number of grid zones\n\t");
	    screen("for pressure and entropy (dflt = %d %d): ",gmax[0],gmax[1]);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscanf(s,"%d %d",gmax,gmax+1);
	}
	else
	{
	    FILE *restart_file = seos->restart_io_type->file;
	    (void) fgetstring(restart_file,"pressure mesh = ");
	    (void) fscanf(restart_file,"%d%*s%*s%d",gmax,gmax+1);
	}
	cgr->h[0] = (U[0] - L[0])/gmax[0];
	cgr->h[1] = (U[1] - L[1])/gmax[1];
	copy_rect_grid(tgr,cgr);
	(void) adjust_top_grid_for_square(tgr,cgr);
	ps_intfc->table->fixed_grid = YES;
	DEBUG_LEAVE(set_PS_rect_grid)
}		/*end set_PS_rect_grid*/

/*ARGSUSED*/
LOCAL	void set_VP_rect_grid(
	INTERFACE 	*rt_intfc,
	INTERFACE 	*vp_intfc,
	void 		(*remap)(POINT*,BOND*,CURVE*,POINT*,BOND*,CURVE*,
				 boolean,RECT_GRID*,POINTER),
	POINTER 	params)
{
	SESAME_EOS	*seos = (SESAME_EOS*)params;
	RECT_GRID	*tgr = &topological_grid(vp_intfc);
	RECT_GRID 	*cgr = computational_grid(vp_intfc);
	RECT_GRID	*rtgr = seos->fr[SESAME_RHO_TEMP]->rect_grid;
	double      *U, *L, *rtU, *rtL;
	int        *gmax, *rtgmax;

	DEBUG_ENTER(set_VP_rect_grid)
	copy_rect_grid(cgr,tgr);
	U = cgr->U;
	L = cgr->L;
	gmax = cgr->gmax;
	rtU = rtgr->U;
	rtL = rtgr->L;
	rtgmax = rtgr->gmax;
	if (seos->restart_io_type == NULL)
	{
	    char s[Gets_BUF_SIZE];
	    screen("Limits for the volume-pressure table\n");
	    screen("\tminimum volume = %g, minimum volume grid = %g\n",
	    	ses_vp_vol_from_grid(L[0],seos),L[0]);
	    screen("\tmaximum volume = %g, maximum volume grid = %g\n",
	    	ses_vp_vol_from_grid(U[0],seos),U[0]);
	    screen("\tminimum pressure = %g, minimum pressure grid = %g\n",
	    	ses_vp_press_from_grid(L[1],seos),L[1]);
	    screen("\tmaximum pressure = %g, maximum pressure grid = %g\n",
	    	ses_vp_press_from_grid(U[1],seos),U[1]);
	    gmax[0] = rtgmax[0];
	    gmax[1] = (int)(rtgmax[1]*((U[1]-L[1])/(rtU[1]-rtL[1])));
	    set_rect_grid(L,U,L,U,NOBUF,NOBUF,gmax,cgr->dim,&rtgr->Remap,cgr);
	    (void) adjust_top_grid_for_square(cgr,cgr);
	    screen("Enter the number of grid zones\n\t");
	    screen("for specific volume and pressure (dflt = %d %d): ",
	    	   gmax[0],gmax[1]);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscanf(s,"%d %d",gmax,gmax+1);
	}
	else
	{
	    FILE *restart_file = seos->restart_io_type->file;
	    (void) fgetstring(restart_file,"density mesh = ");
	    (void) fscanf(restart_file,"%d%*s%*s%d",gmax,gmax+1);
	}
	cgr->h[0] = (U[0] - L[0])/gmax[0];
	cgr->h[1] = (U[1] - L[1])/gmax[1];
	copy_rect_grid(tgr,cgr);
	(void) adjust_top_grid_for_square(tgr,cgr);
	vp_intfc->table->fixed_grid = YES;
	DEBUG_LEAVE(set_VP_rect_grid)
}		/*end set_VP_rect_grid*/



/*
*			s_ses_solution():
*
*	Calculates the state at an arbitrary point x,y, as considered
*	to belong to a given component comp.  The state is obtained by
*	interpolation between interior states and states on the front
*	and is loaded in the storage pointed to by state.  If comp is
*	exterior_component(front->interf) then boundary conditions are applied.
*/

/* ARGSUSED */
LOCAL void s_ses_solution(
	double 		*coords,
	COMPONENT 	comp,
	HYPER_SURF 	*hs,
	SIDE 		side,
	Front 		*front,
	POINTER		wv,
	Locstate 	state,
	Locstate	dflt_state)
{
	double 		coords_on[MAXD],t;
	HYPER_SURF_ELEMENT *hse;
	INTERFACE 	*intfc = front->interf;
	Wave		*wave = (Wave*)wv;

	if ((is_exterior_comp(comp,intfc)) || (comp != component(coords,intfc)))
	{
	    if (nearest_interface_point(coords,comp,intfc,INCLUDE_BOUNDARIES,
					hs,coords_on,&t,&hse,&hs) != YES)
	    {
	        screen("ERROR in s_ses_solution(), "
	               "nearest_interface_point failed\n");
	        clean_up(ERROR);
	    }
	    right_state_along_bond(t,Bond_of_hse(hse),Curve_of_hs(hs),state);
		return;
	}

	(void) tri_solution(coords,comp,wave_tri_soln(wave),state,dflt_state);
}		/*end s_ses_solution*/



/*
*			init_ses_interior_states():
*			init_RS_interior_states():
*			init_PS_interior_states():
*
*	Initializes the states in a wave structure by calling
*
*		(*initializer)(coords,y,comp,state,params)
*
*	at the centers of the grid blocks of wave->rect_grid.
*/

LOCAL void init_ses_interior_states(
	Wave 		  *wave,
	Front 		  *front,
	void 		  (*initializer)(double*,COMPONENT,Locstate,SESAME_EOS*),
	SESAME_EOS	  *seos,
	SESAME_TABLE_TYPE table)
{
	INTERFACE 	*intfc;
	int 		ix,iy;
	int 		icoords[MAXD];
	int 		xmax = wave->rect_grid->gmax[0];
	int 		ymax = wave->rect_grid->gmax[1];
	size_t 		sizest = front->sizest;
	double 		*coords;
	Locstate 	state;
	COMPONENT 	comp;
	int 		status;

	DEBUG_ENTER(init_ses_interior_states)

	if (wave->sizest == 0 || initializer == NULL)
	{
		DEBUG_LEAVE(init_ses_interior_states)
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
			coords = Rect_coords(icoords,wave);
			comp = Rect_comp(icoords,wave);
			state = Rect_state(icoords,wave);
			if (is_exterior_comp(comp,intfc))
			{
				clear_state(front->interf,state,sizest);
			}
			else
				(*initializer)(coords,comp,state,seos);
		}
	}

	DEBUG_LEAVE(init_ses_interior_states)
}		/*end init_ses_interior_states*/

LOCAL void init_RS_interior_states(
	SESAME_EOS 	*seos)
{
	Wave 		*wave = seos->wave[SESAME_RHO_ENTROPY];
	Front 		*front = seos->fr[SESAME_RHO_ENTROPY];
	COMPONENT 	comp;
	INTERFACE 	*intfc;
	int 		ix,iy;
	int 		icoords[MAXD];
	int 		xmax = wave->rect_grid->gmax[0];
	int 		ymax = wave->rect_grid->gmax[1];
	double 		*coords;
	Locstate 	state;
	int 		status;
	size_t 		sizest = front->sizest;

	DEBUG_ENTER(init_RS_interior_states)
	set_ses_intrp_flag_all(SESAME_RHO_ENTROPY);
	status = init_ses_hyp_soln_func(wave,front);
	if (status != GOOD_STEP) 
	{
	    screen("ERROR in init_RS_interior_states(), ");
	    screen("init_ses_hyp_soln_func(), failed");
	    print_sesame_interface(seos,front->interf,SESAME_RHO_ENTROPY);
	    clean_up(ERROR);
	}
	intfc = front->interf;

	for (iy = 0; iy < ymax; iy++)
	{
		icoords[1] = iy;
		for (ix = 0; ix < xmax; ix++) 
		{
			icoords[0] = ix;
			coords = Rect_coords(icoords,wave);
			comp = Rect_comp(icoords,wave);
			state = Rect_state(icoords,wave);
			if (is_exterior_comp(comp,intfc))
			{
				clear_state(front->interf,state,sizest);
			}
			else
			{
				ses_RS_initializer(coords,comp,state,seos);
			}
		}
	}

	DEBUG_LEAVE(init_RS_interior_states)
}		/*end init_RS_interior_states*/

LOCAL void init_PS_interior_states(
	SESAME_EOS 	*seos)
{
	Wave 		*wave = seos->wave[SESAME_PRESS_ENTROPY];
	Front 		*front = seos->fr[SESAME_PRESS_ENTROPY];
	COMPONENT 	comp;
	INTERFACE 	*intfc;
	int 		ix,iy;
	int 		icoords[MAXD];
	int 		xmax = wave->rect_grid->gmax[0];
	int 		ymax = wave->rect_grid->gmax[1];
	double 		*coords;
	Locstate 	state;
	int 		status;
	size_t 		sizest = front->sizest;

	DEBUG_ENTER(init_PS_interior_states)

	if (sizest == 0 ) 
	{
	    DEBUG_LEAVE(init_PS_interior_states)
	    return;
	}

	set_ses_intrp_flag_all(SESAME_PRESS_ENTROPY);
	status = init_ses_hyp_soln_func(wave,front);
	if (status != GOOD_STEP) 
	{
	    screen("ERROR in init_PS_interior_states(), "
	           "init_ses_hyp_soln_func() failed\n");
	    clean_up(ERROR);
	}
	intfc = front->interf;

	for (iy = 0; iy < ymax; iy++)
	{
	    icoords[1] = iy;
	    for (ix = 0; ix < xmax; ix++) 
	    {
	    	icoords[0] = ix;
	    	coords = Rect_coords(icoords,wave);
	    	comp = Rect_comp(icoords,wave);
	    	state = Rect_state(icoords,wave);
	    	if (is_exterior_comp(comp,intfc))
	    	{
	    	    clear_state(front->interf,state,sizest);
	    	}
	    	else
	    	{
	    	    ses_PS_initializer(coords,comp,state,seos);
	    	}
	    }
	}

	DEBUG_LEAVE(init_PS_interior_states)
}		/*end init_PS_interior_states*/


/*
*			init_ses_hyp_soln_func():
*
*	Performs the necessary initialization to enable use of
*	s_ses_solution():  constructs the triangulated
*	grid and allocates the storage for the interior states.
*	Also initializes wave_areas(wave) for use in Rect_area().
*
*	Note that the interior states need to be loaded following
*	the call to this routine in order for s_ses_solution() to work.
*/

EXPORT int init_ses_hyp_soln_func(
	Wave 		*wave,
	Front 		*front)
{
	RECT_GRID 	Dual_gr, *gr = wave->rect_grid;
	int 		status;

	DEBUG_ENTER(init_ses_hyp_soln_func)
	intfc_delete_very_short_bonds(front);
	if (intfc_delete_fold_back_bonds(front) == FUNCTION_FAILED)
	{
		/*
		* The occurence of fold back bonds is
		* possible with the current EOS code,  at least
		* in the P-S interface.  This can cause errors
		* in the trigrid generation.  The cause of the
		* fold back bonds is currently unknown and
		* is probably a bug.  You should attempt to
		* discover the root cause.
		*/
		screen("ERROR in init_ses_hyp_soln_func(), ");
		screen("intfc_delete_fold_back_bonds() failed\n");
		clean_up(ERROR);
	}

	clear_wave_pointers(wave);
	if (wave->sizest == 0)
	{
		DEBUG_LEAVE(init_ses_hyp_soln_func)
		return GOOD_STEP;
	}

	/*
	*	trisoln.c locates interior states at crossings
	*	of grid lines; thus it should not be given
	*	wave->rect_grid but rather its dual.
	*/

	scalar(&wave_tri_soln(wave),sizeof(TRI_SOLN));
	if (wave_tri_soln(wave) == NULL)
	{
	    (void) printf("WARNING in init_ses_hyp_soln_func(), "
	                  "can't allocate tri solution\n");
	    DEBUG_LEAVE(init_ses_hyp_soln_func)
	    return ERROR_IN_STEP;
	}
	wave_tri_soln(wave)->Tri_grid_hooks = wave->Tri_grid_hooks;

	set_dual_grid(&Dual_gr,gr);

	status = hyp_tri_grid_driver(front,wave,&Dual_gr);

	if (status != GOOD_STEP)
	{
	    (void) printf("WARNING in init_ses_hyp_soln_func(), "
	                  "hyp_tri_grid_driver() failed\n");
	    DEBUG_LEAVE(init_ses_hyp_soln_func)
	    return status;
	}

	front->_hyp_solution = s_ses_solution;
	front->_hyp_grad_solution = NULL;

	DEBUG_LEAVE(init_ses_hyp_soln_func)
	return GOOD_STEP;
}		/*end init_ses_hyp_soln_func*/

LOCAL	const char *eos_wave_type_as_string(
	int	w_type)
{
	switch (w_type)
	{
	case ISOTHERM_BOUNDARY:
	    return "ISOTHERM_BOUNDARY";
	case ISOCHOR_BOUNDARY:
	    return "ISOCHOR_BOUNDARY";
#if defined(PHASE_CODE)
	case PHASE_BOUNDARY:
	    return "PHASE_BOUNDARY";
#endif /* defined(PHASE_CODE) */
	default:
	    return g_wave_type_as_string(w_type);
	}
}		/*end eos_wave_type_as_string*/

LOCAL	int read_eos_wave_type_from_string(
	const char 	*type)
{
	int 	w_type = UNKNOWN_WAVE_TYPE;

	switch(type[0]) 
	{
	case 'I':
	case 'i':
	    if ((strcmp(type,"ISOTHERM_BOUNDARY") == 0) ||
	        (strcmp(type,"isotherm_boundary") == 0))
	        w_type = ISOTHERM_BOUNDARY;
	    else if ((strcmp(type,"ISOCHOR_BOUNDARY") == 0) ||
	             (strcmp(type,"isochor_boundary") == 0))
	        w_type = ISOCHOR_BOUNDARY;
	    break;
#if defined(PHASE_CODE)
	case 'P':
	case 'p':
	    w_type = PHASE_BOUNDARY;
	    break;
#endif /* defined(PHASE_CODE) */

	default:
	    break;
	}
	if (w_type == UNKNOWN_WAVE_TYPE)
	    w_type = g_read_wave_type_from_string(type);
	return w_type;
}	/*end read_eos_wave_type_from_string*/

/*ARGSUSED*/
LOCAL	void fprint_eos_node_type(
	FILE 		*file,
	const char	*mesg1,
	int 		n_type,
	const char 	*mesg2,
	INTERFACE 	*intfc)
{
	switch (n_type)
	{
	case EOS_BOUNDARY_NODE:
	    if (mesg1 != NULL)
	    	(void) fprintf(file,"%s",mesg1);
	    (void) fprintf(file,"EOS_BOUNDARY_NODE");
	    if (mesg2 != NULL)
	    	(void) fprintf(file,"%s",mesg2);
	    return;
#if defined(PHASE_CODE)
	case PHASE_BDRY_NODE:
	    if (mesg1 != NULL)
	    	(void) fprintf(file,"%s",mesg1);
	    (void) fprintf(file,"PHASE_BDRY_NODE");
	    if (mesg2 != NULL)
	    	(void) fprintf(file,"%s",mesg2);
	    return;
#endif /* defined(PHASE_CODE) */
	default:
	    g_fprint_hsbdry_type(file,mesg1,n_type,mesg2,intfc);
	    break;
	}
}		/*end fprint_eos_node_type*/

/*ARGSUSED*/
LOCAL	int read_eos_node_type_from_string(
	const char *type,
	INTERFACE  *intfc)
{
	int 		n_type = UNKNOWN_NODE_TYPE;

	switch(type[0]) 
	{
	case 'E':
	case 'e':
	    n_type = EOS_BOUNDARY_NODE;
	    break;
#if defined(PHASE_CODE)
	case 'P':
	case 'p':
	    n_type = PHASE_BDRY_NODE;
	    break;
#endif /* defined(PHASE_CODE) */

	default:
	    n_type = g_read_hsbdry_type_from_string(type,intfc);
	    break;
	}
	return n_type;
}		/*end read_eos_node_type_from_string*/


/*
*			set_RT_riv():
*
*	Computes the Riemann function integral c/rho drho at constant
*	entropy for the density temperature table.  This function is
*	computed by solving the scalar partial differential equation
*
*		f_rho + (GAMMA*T/rho)*f_T = c/rho
*
*	here GAMMA is the Gruneisen exponent.
*
*	This equation can be written as
*
*		f_t + A*f_x = B
*
*	where
*
*	rho = ses_rt_rho_from_grid(t,seos),
*	T = ses_rt_temp_from_grid(x,seos),
*	A = GAMMA*dlog_rho_drho_grid(rho,seos)*dT_grid_dlogT(T,seos)
*	B = dlog_rho_drho_grid(rho,seos)*c
*
*	This equation is solved
*	using a second order upwind method:
*
*
*	f[n][i] = f[n-1][i] + B[n-1/2][i]*drho
*	            - max(A[n-1/2][i],0)*h[n-1/2][i-1/2]
*		    - min(A[n-1/2][i],0)*h[n-1/2][i+1/2]
*
*	where
*
*	A[n-1/2][i] = 0.5*(A[n][i] + A[n-1][i]);
*	B[n-1/2][i] = 0.5*(B[n][i] + B[n-1][i]);
*
*	and
*
*	h[n-1/2][i+1/2] = (dt/dx)*(df[n-1][i+1/2] + dB[n-1][i+1/2]*dt/2)/
*			          (1.0 + dA[n-1][i+1/2]*(dt/dx)/2)
*
*
*	where
*
*	df[n-1][i+1/2] = f[n-1][i+1] - f[n-1][i]
*	dB[n-1][i+1/2] = B[n-1][i+1] - B[n-1][i]
*	dA[n-1][i+1/2] = A[n-1][i+1] - A[n-1][i]
*
*	Initial conditions are f[0][i] = 0,  and the boundary conditions
*	are given by dropping the term from the update for f that corresponds
*	to signals from outside the domain.
*/

/*ARGSUSED*/
LOCAL	void	set_RT_riv(
	Front		*fr,
	Wave		*wave,
	SESAME_EOS	*seos)
{
	RECT_GRID		*gr = fr->rect_grid;
	POINT			*pt;
	HYPER_SURF		*hs;
	HYPER_SURF_ELEMENT	*hse;
	INTERFACE		*intfc = fr->interf;
	Locstate		s, sl, sr;
	double			GAMMA, GAMMA_max;
	double			***fstore;
	double			**f, **A, **B, **c;
	double			Anm1;
	double			dx, dt;
	double			tmin, tmax, xmin;
	double			rho, T, t, x;
	double			p[3], e[3];
	double			csq;
	double			l;
	double			a, b;
	double			*crds;
	double			cref;
	int			nt, nx;
	int			n, i, icrds[3];

	DEBUG_ENTER(set_RT_riv)
	tmin = gr->L[0];
	tmax = gr->U[0];
	xmin = gr->L[1];
	nx = gr->gmax[1];
	dx = gr->h[1];
	GAMMA_max = Gru_gam_abs_max(seos);
	dt = (2.0/3.0)*dx/GAMMA_max;
	dt = min(gr->h[0],dt);
	nt = irint((tmax-tmin)/dt);
	dt = (tmax-tmin)/nt;
	tri_array(&fstore,4,nt+1,nx+1,FLOAT);
	f = fstore[0];
	A = fstore[1];
	B = fstore[2];
	c = fstore[3];

	if (DEBUG)
	{
	    (void) printf("dt = %g, nt = %d, dx = %g, nx = %d\n",
			  dt,nt,dx,nx);
	    (void) printf("tmin = %g, tmax = %g, xmin = %g, xmax = %g\n",
			  tmin,tmax,xmin,gr->U[1]);
	}

	l = dt/dx;
	for (n = 0; n <= nt; n++)
	{
	    rho = ses_rt_rho_from_grid(tmin + n*dt,seos);
	    for (i = 0; i <= nx; i++)
	    {
		T = ses_rt_temp_from_grid(xmin + i*dx,seos);
		s2eos_lookup(rho,T,p,e,seos);
		GAMMA = p[2]/(e[2]*rho);
		A[n][i] = GAMMA*dlog_rho_drho_grid(rho,seos)*
				dT_grid_dlogT(T,seos);
		csq = sesame_rt_sound_speed_squared(rho,T,p,e);
		c[n][i] = sqrt(csq);
		B[n][i] = c[n][i]*dlog_rho_drho_grid(rho,seos);
	    }
	}
	for (i = 0; i <= nx; i++)
	    f[0][i] = 0.0;
	for (n = 1; n <= nt; n++)
	{
	    f[n][0] = f[n-1][0] + 0.5*(B[n-1][0]+B[n][0])*dt;
	    Anm1 = 0.5*(A[n-1][0]+A[n][0]);
	    if (Anm1 < 0.0)
	    {
		f[n][0] -= Anm1*l*
		    ((f[n-1][1]-f[n-1][0]) + 0.5*(B[n-1][1]-B[n-1][0])*dt)/
		    (1.0 + 0.5*l*(A[n-1][1]-A[n-1][0]));
	    }
	    for (i = 1; i < nx; i++)
	    {
	        f[n][i] = f[n-1][i] + 0.5*(B[n-1][i]+B[n][i])*dt;
	        Anm1 = 0.5*(A[n-1][i]+A[n][i]);
		if (Anm1 < 0.0)
		{
		    f[n][i] -= Anm1*l*((f[n-1][i+1]-f[n-1][i]) +
					0.5*(B[n-1][i+1]-B[n-1][i])*dt)/
		               (1.0 + 0.5*l*(A[n-1][i+1]-A[n-1][i]));
		}
		else if (Anm1 > 0.0)
		{
	            f[n][i] -= Anm1*l*((f[n-1][i]-f[n-1][i-1]) + 
				        0.5*(B[n-1][i]-B[n-1][i-1])*dt)/
		                    (1.0 + 0.5*l*(A[n-1][i]-A[n-1][i-1]));
		}
	    }
	    f[n][nx] = f[n-1][nx] + 0.5*(B[n-1][nx]+B[n][nx])*dt;
	    Anm1 = 0.5*(A[n-1][nx]+A[n][nx]);
	    if (Anm1 > 0.0)
	    {
	        f[n][nx] -= Anm1*l*((f[n-1][nx]-f[n-1][nx-1]) + 
				       0.5*(B[n-1][nx]-B[n-1][nx-1])*dt)/
		                    (1.0 + 0.5*l*(A[n-1][nx]-A[n-1][nx-1]));
	    }
	}

	cref = Reference_sound_speed(seos);
	for (n = 0; n <= nt; n++)
	for (i = 0; i <= nx; i++)
	{
	    f[n][i] /= 0.5*(c[n][i]+cref);
	}

	/*Set front states*/
	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&pt,&hse,&hs) != NO)
	{
	    slsr(pt,hse,hs,&sl,&sr);
	    t = Coords(pt)[0];
	    x = Coords(pt)[1];
	    n = (int)((t - tmin)/dt);
	    n = max(n,0);
	    n = min(n,nt-1);
	    i = (int)((x - xmin)/dx);
	    i = max(i,0);
	    i = min(i,nx-1);
	    a = (t - (tmin+n*dt))/dt;
	    a = max(a,0.0);
	    a = min(a,1.0);
	    b = (x - (xmin+i*dx))/dx;
	    b = max(b,0.0);
	    b = min(b,1.0);
	    ses_rt_riv(sl) = ses_rt_riv(sr) = 
		(1.0-b)*((1.0-a)*f[n][i]   + a*f[n+1][i]) +
		      b*((1.0-a)*f[n][i+1] + a*f[n+1][i+1]);
	}
	/*Set interior states*/
	icrds[2] = 0;
	for (icrds[0] = 0; icrds[0] < wave->rect_grid->gmax[0]; icrds[0]++)
	for (icrds[1] = 0; icrds[1] < wave->rect_grid->gmax[1]; icrds[1]++)
	{
	    crds = Rect_coords(icrds,wave);
	    s = Rect_state(icrds,wave);
	    t = crds[0];
	    x = crds[1];
	    n = (int)((t - tmin)/dt);
	    n = max(n,0);
	    n = min(n,nt-1);
	    i = (int)((x - xmin)/dx);
	    i = max(i,0);
	    i = min(i,nx-1);
	    a = (t - (tmin+i*dt))/dt;
	    a = max(a,0.0);
	    a = min(a,1.0);
	    b = (x - (xmin+i*dx))/dx;
	    b = max(b,0.0);
	    b = min(b,1.0);
	    ses_rt_riv(s) = (1.0-b)*((1.0-a)*f[n][i]   + a*f[n+1][i]) +
		                  b*((1.0-a)*f[n][i+1] + a*f[n+1][i+1]);
	}
	free(fstore);
	DEBUG_LEAVE(set_RT_riv)
}		/*end set_RT_riv*/

LOCAL	void	print_sesame_interface(
	SESAME_EOS	  *seos,
	INTERFACE	  *intfc,
	SESAME_TABLE_TYPE table_no)
{
	(void) output();
	(void) printf("%s INTERFACE\n",ses_table_name(seos,table_no));
	print_interface(intfc);
}		/*end print_sesame_interface*/

/*
*			print_spoly_fit():
*
*	Computes and prints EOS fits to the Stiffened Polytropic EOS,
*	In terms of the sesame representation p = P(rho,T), e = e(rho,T),
*	we have a local fit of the form:
*
*	      de  |
*	C  = ---- |
*	 V    dT  |rho
*
*	         dp  |                      de  |
*	        ---- |         P - rho*rho*---- |
*	         dT  |rho                  drho |T
*	GAM = ------------ = ----------------------
*	        rho*C               rho*T*C
*	             V                     V
*
*	GAM = Gruneisen exponent
*
*	                            dp  |
*	                     GAM*T*---- |
*	        dp  |               dT  |rho
*	c2 =   ---- |    + ------------------
*	       drho | T           rho
*
*
*	c2 = sound speed squared
*
*	GAM = GAM(rho0,T0)
*
*	         1      dp  |                de  |
*	Einf = ----- * ---- |   - e - rho * ---- |
*	        GAM    drho |T              drho |T
*
*
*	        GAM*rho*(e + Einf) - p
*	Pinf = ------------------------ 
*	                GAM + 1
*
*
*	       de  |
*	 C  = ---- |
*	  V    dT  |T
*
*
*	       C  * rho * T - rho*(e + Einf) + Pinf
*	        V
*	et =  --------------------------------------
*                                  GAM+1
*                           GAM*rho
*/

LOCAL	void	print_spoly_fit(
	FILE		*file,
	SESAME_EOS	*seos)
{
	Front     *fr = seos->fr[SESAME_RHO_TEMP];
	Wave 	  *wv = seos->wave[SESAME_RHO_TEMP];
	RECT_GRID *gr = fr->rect_grid;
	SPOLY_FIT Spfl, Spfm, Spfs;
	Locstate  state;
	double	  N;
	int	  gnr, gnt, nr, nt, NR, NT;
	double     *rho, *T, rho_grid, T_grid;
	double	  **p, **e;
	double	  rmin, rmax, tmin, tmax;
	double 	  *tbls = seos->sestab.tbls;
	double	  *rho_a, *T_a, *e_a, *p_a;
	int	  irmin, irmax, jtmin, jtmax;
	int	  i, j, icoords[3];


	zero_scalar(&Spfm,sizeof(SPOLY_FIT));
	zero_scalar(&Spfs,sizeof(SPOLY_FIT));

	gnr = gr->gmax[0], gnt = gr->gmax[1];
	nr = irint(tbls[2]);
	nt = irint(tbls[3]);
	NR = max(nr,gnr);
	NT = max(nt,gnt);
	uni_array(&rho,NR,FLOAT);
	uni_array(&T,NT,FLOAT);
	bi_array(&p,NR,NT,FLOAT);
	bi_array(&e,NR,NT,FLOAT);

	for (i = 0; i < gnr; i++)
	{
	    icoords[0] = i;
	    rho_grid = cell_center(i,0,gr);
	    rho[i] = ses_rt_rho_from_grid(rho_grid,seos);
	    for (j = 0; j < gnt; j++)
	    {
		icoords[1] = j;
		T_grid = cell_center(j,1,gr);
		T[j] = ses_rt_temp_from_grid(T_grid,seos);
		state = Rect_state(icoords,wv);
		p[i][j] = rho[i]*T[j]*ses_rt_redp(state) + ses_rt_coldp(state);
		e[i][j] = T[j]*ses_rt_rede(state) + ses_rt_colde(state);

		linear_spoly_fit(seos,rho[i],T[j],&Spfl);
		Spfm.GAM += Spfl.GAM;   Spfs.GAM += sqr(Spfl.GAM);
		Spfm.Pinf += Spfl.Pinf; Spfs.Pinf += sqr(Spfl.Pinf);
		Spfm.Einf += Spfl.Einf; Spfs.Einf += sqr(Spfl.Einf);
		Spfm.cv += Spfl.cv;     Spfs.cv += sqr(Spfl.cv);
		Spfm.et += Spfl.et;     Spfs.et += sqr(Spfl.et);
	    }
	}
	N = gnr*gnt;

	Spfm.GAM /= N;
	Spfs.GAM = sqrt(Spfs.GAM/(N-1) - (N/(N-1))*Spfm.GAM*Spfm.GAM);

	Spfm.Pinf /= N;
	Spfs.Pinf = sqrt(Spfs.Pinf/(N-1) - (N/(N-1))*Spfm.Pinf*Spfm.Pinf);

	Spfm.Einf /= N;
	Spfs.Einf = sqrt(Spfs.Einf/(N-1) - (N/(N-1))*Spfm.Einf*Spfm.Einf);

	Spfm.cv /= N;
	Spfs.cv = sqrt(Spfs.cv/(N-1) - (N/(N-1))*Spfm.cv*Spfm.cv);

	Spfm.et /= N;
	Spfs.et = sqrt(Spfs.et/(N-1) - (N/(N-1))*Spfm.et*Spfm.et);

	(void) fprintf(file,"\nStiffened polytropic average fit,\n"
	    "    Gruneisen coefficient = %- 14.8g,       variance = %- 14.8g\n"
	    "    P_infinity            = %- 14.8g GPa,   variance = %- 14.8g\n"
	    "    E_infinity            = %- 14.8g kJ/g,  variance = %- 14.8g\n"
	    "    C_V                   = %- 14.8g kJ/gK, variance = %- 14.8g\n"
	    "    et                    = %- 14.8g kJ/g,  variance = %- 14.8g\n"
	    "    R =                   = %- 14.8g kJ/gK\n",
	    Spfm.GAM,Spfs.GAM,
	    Spfm.Pinf,Spfs.Pinf,
	    Spfm.Einf,Spfs.Einf,
	    Spfm.cv,Spfs.cv,
	    Spfm.et,Spfs.et,Spfm.GAM*Spfm.cv);

	least_squares_spoly_fit(file,
				"\nStiffened polytropic least squares fit,\n",
				rho,T,p,e,gnr,gnt);

	/* Least squares fit (raw table windowed to density/temperature domain*/

	rmin = Rho_min(seos);
	rmax = Rho_max(seos);
	tmin = Temp_min(seos);
	tmax = Temp_max(seos);
	rho_a = tbls + 4;
	T_a   = rho_a + nr;
	p_a   = T_a + nt;
	e_a   = p_a + nr*nt;

	/* Find density window */
	for (irmin = 0; irmin < nr; irmin++)
	    if (rho_a[irmin] >= rmin)
		break;
	for (irmax = nr-1; irmax > -1; irmax--)
	    if (rho_a[irmax] <= rmax)
		break;
	/* Find temperature window */
	for (jtmin = 0; jtmin < nt; jtmin++)
	    if (T_a[jtmin] >= tmin)
		break;
	for (jtmax = nt-1; jtmax > -1; jtmax--)
	    if (T_a[jtmax] <= tmax)
		break;

	if ((irmin == nr) || (irmax == -1) || (irmin > irmax))
	{
	    (void) printf("WARNING in SESAME_stiffen_polytropic_fit(), "
			  "no density points in window\n");
	}
	else if ((jtmin == nt) || (jtmax == -1) || (jtmin > jtmax))
	{
	    (void) printf("WARNING in SESAME_stiffen_polytropic_fit(), "
			  "no temperature points in window\n");
	}
	else
	{
	    NR = irmax - irmin +1;
	    NT = jtmax - jtmin + 1;

	    for (i = irmin; i <= irmax; i++)
	    {
		rho[i-irmin] = rho_a[i];
	        for (j = jtmin; j <= jtmax; j++)
	        {
		    T[j-jtmin] = T_a[j];
		    p[i-irmin][j-jtmin] = p_a[nr*j+i];
		    e[i-irmin][j-jtmin] = e_a[nr*j+i];
		}
	    }
	    least_squares_spoly_fit(file,"\nStiffened polytropic raw windowed "
				    "(rho,T) table fit,\n",
				    rho,T,p,e,NR,NT);
	}

	/*compute spoly fit from whole raw table*/
	nr = irint(tbls[2]);
	nt = irint(tbls[3]);
	for (i = 0; i < nr; i++)
	{
	    rho[i] = rho_a[i];
	    for (j = 0; j < nt; j++)
	    {
		T[j] = T_a[j];
		p[i][j] = p_a[nr*j+i];
		e[i][j] = e_a[nr*j+i];
	    }
	}
	least_squares_spoly_fit(file,"\nStiffened polytropic raw whole "
				"(rho,T) table fit,\n",rho,T,p,e,nr,nt);
	(void) fprintf(file,"\n");
	free_these(4,rho,T,p,e);
}		/*end print_spoly_fit*/

LOCAL	void	linear_spoly_fit(
	SESAME_EOS *seos,
	double	   rho,
	double	   T,
	SPOLY_FIT  *spf)
{
	double p, e, pv[3], ev[3];
	double GAM, gam, Einf, Pinf, cv;

	s2eos_lookup(rho,T,pv,ev,seos);

	p = pv[0];
	e = ev[0];
	spf->GAM = GAM = pv[2]/(rho*ev[2]);
	gam = GAM + 1.0;
	spf->Einf = Einf = pv[1]/GAM - e - rho*ev[1];
	spf->Pinf = Pinf = (GAM*rho*(e + Einf) - p)/gam;
	spf->cv = cv = ev[2];
	spf->et = (cv*rho*T - rho*(e + Einf) + Pinf)/pow(rho,gam);
	if (spf->et < 0.0)
	    spf->et = 0.0;
}		/*end linear_spoly_fit*/

LOCAL	void	least_squares_spoly_fit(
	FILE        *file,
	const char  *mesg,
	double       *rho,
	double       *T,
	double       **p,
	double       **e,
	int         nr,
	int         nt)
{
	double RHO, TEMP, P, Pfit, E, Efit, RG;
	double DRHO, DTEMP;
	double GAM, Einf, Pinf, cv, et;
	double pe, ee, npe, nee, rpe, ree;
	double icpe, ricpe, nicpe;
	double err, **perr, **eerr, **rperr, **reerr;
	double **icperr, **ricperr;
	double pe_max, rpe_max, ee_max, ree_max;
	double RHO_at_pe_max, T_at_pe_max, P_at_pe_max, E_at_pe_max;
	double RHO_at_ee_max, T_at_ee_max, P_at_ee_max, E_at_ee_max;
	double RHO_at_rpe_max, T_at_rpe_max, P_at_rpe_max, E_at_rpe_max;
	double RHO_at_ree_max, T_at_ree_max, P_at_ree_max, E_at_ree_max;
	double icpe_max, ricpe_max;
	double RHO_at_icpe_max, T_at_icpe_max, P_at_icpe_max, E_at_icpe_max;
	double RHO_at_ricpe_max, T_at_ricpe_max, P_at_ricpe_max, E_at_ricpe_max;
	double pcell2, pbar, ecell2, ebar;
	double N, D;
	double a00, a01, a02,
	      a10, a11, a12,
	      a20, a21, a22,
	       x0,  x1,  x2,
	       b0,  b1,  b2;
	int   i, j;

	bi_array(&perr,nr,nt,FLOAT);
	bi_array(&rperr,nr,nt,FLOAT);

	bi_array(&eerr,nr,nt,FLOAT);
	bi_array(&reerr,nr,nt,FLOAT);

	bi_array(&icperr,nr,nt,FLOAT);
	bi_array(&ricperr,nr,nt,FLOAT);

	if (rho[0] == 0.0)
	{
	    rho++;
	    nr--;
	}
	if (T[0] == 0.0)
	{
	    T++;
	    nt--;
	}

	a00 = 0.0; a01 = 0.0; a02 = 0.0;
	a10 = 0.0; a11 = 0.0; a12 = 0.0;
	a20 = 0.0; a21 = 0.0; a22 = 0.0;
	 b0 = 0.0;  b1 = 0.0;  b2 = 0.0;;
	for (i = 0; i < nr; i++)
	{
	    RHO = rho[i];
	    for (j = 0; j < nt; j++)
	    {
		E = e[i][j];
		P = p[i][j];
		a00 += sqr(RHO*E); a01 += RHO*RHO*E; a02 += RHO*E;
				   a11 += RHO*RHO;   a12 += RHO;
		b0 += RHO*E*P;     b1 += RHO*P;      b2 += P;
	    }
	}
	N = nr*nt;
	a00  /= N;  a01 /= N; a02  /= N;
	a10 = a01;  a11 /= N; a12  /= N;
	a20 = a02; a21 = a12; a22 = 1.0;
	 b0  /= N;  b1  /= N;  b2  /= N;

	/*Solve Ax = b using Cramer's rule*/
	D = a00*a11*a22 + a01*a12*a20 + a02*a10*a21 -
	    a02*a11*a20 - a00*a12*a21 - a01*a10*a22;
	x0 = (b0*a11*a22 + a01*a12*b2 + a02*b1*a21 -
	      a02*a11*b2 - b0*a12*a21 - a01*b1*a22)/D;
	x1 = (a00*b1*a22 + b0*a12*a20 + a02*a10*b2 -
	      a02*b1*a20 - a00*a12*b2 - b0*a10*a22)/D;
	x2 = (a00*a11*b2 + a01*b1*a20 + b0*a10*a21 -
	      b0*a11*a20 - a00*b1*a21 - a01*a10*b2)/D;
	GAM = x0;
	Einf = x1/GAM;
	Pinf = -x2/(GAM+1.0);

	a00 = 0.0; a01 = 0.0;
	a10 = 0.0; a11 = 0.0;
	a20 = 0.0; a21 = 0.0;
	 b0 = 0.0;  b1 = 0.0;
	for (i = 0; i < nr; i++)
	{
	    RHO = rho[i];
	    RG = pow(RHO,GAM+1.0);
	    for (j = 0; j < nt; j++)
	    {
		TEMP = T[j];
		E = e[i][j];
		P = p[i][j];
		a00 += sqr(RHO*TEMP);        a01 -= RHO*RG*TEMP;
				             a11 += RG*RG;
		b0 += RHO*TEMP*(RHO*(E+Einf) - Pinf);
		b1 += RG*(Pinf - RHO*(E+Einf));
	    }
	}
	a00  /= N; a01 /= N;
	a10 = a01; a11 /= N;
	 b0  /= N;  b1 /= N;

	D = a00*a11 - a10*a01;
	cv = (b0*a11 - b1*a01)/D;
	et = (a00*b1 - a10*b0)/D;

	pe_max = rpe_max = ee_max = ree_max = icpe_max = ricpe_max = 0.0;
	RHO_at_pe_max = T_at_pe_max = P_at_pe_max = E_at_pe_max = 0.0;
	RHO_at_rpe_max = T_at_rpe_max = P_at_rpe_max = E_at_rpe_max = 0.0;
	RHO_at_ee_max = T_at_ee_max = P_at_ee_max = E_at_ee_max = 0.0;
	RHO_at_ree_max = T_at_ree_max = P_at_ree_max = E_at_ree_max = 0.0;
	RHO_at_icpe_max = T_at_icpe_max = P_at_icpe_max = E_at_icpe_max = 0.0;
	RHO_at_ricpe_max = T_at_ricpe_max =
	    P_at_ricpe_max = E_at_ricpe_max = 0.0;
	for (i = 0; i < nr; i++)
	{
	    RHO = rho[i];
	    RG = pow(RHO,GAM);
	    for (j = 0; j < nt; j++)
	    {
		TEMP = T[j];
		E = e[i][j];
		P = p[i][j];
		Pfit = GAM*RHO*(cv*TEMP - et*RG) - Pinf;
		perr[i][j] = fabs(P - Pfit);
		if (perr[i][j] > pe_max)
		{
		    pe_max = perr[i][j];
		    RHO_at_pe_max = RHO;
		    T_at_pe_max = TEMP;
		    P_at_pe_max = P;
		    E_at_pe_max = E;
		}
		rperr[i][j] = (P > 0.0) ? fabs(Pfit/P - 1.0) : 0.0;
		if (rperr[i][j] > rpe_max)
		{
		    rpe_max = rperr[i][j];
		    RHO_at_rpe_max = RHO;
		    T_at_rpe_max = TEMP;
		    P_at_rpe_max = P;
		    E_at_rpe_max = E;
		}

		Efit = (cv*TEMP+Pinf/RHO - et*RG) - Einf;
		eerr[i][j] = fabs(E - Efit);
		if (eerr[i][j] > ee_max)
		{
		    ee_max = eerr[i][j];
		    RHO_at_ee_max = RHO;
		    T_at_ee_max = TEMP;
		    P_at_ee_max = P;
		    E_at_ee_max = E;
		}
		reerr[i][j] = (E > 0.0) ? fabs(Efit/E - 1.0) : 0.0;
		if (reerr[i][j] > ree_max)
		{
		    ree_max = reerr[i][j];
		    RHO_at_ree_max = RHO;
		    T_at_ree_max = TEMP;
		    P_at_ree_max = P;
		    E_at_ree_max = E;
		}

	        Pfit = GAM*RHO*(E + Einf) - (GAM+1.0)*Pinf;
		icperr[i][j] = fabs(P - Pfit);
		if (icperr[i][j] > icpe_max)
		{
		    icpe_max = icperr[i][j];
		    RHO_at_icpe_max = RHO;
		    T_at_icpe_max = TEMP;
		    P_at_icpe_max = P;
		    E_at_icpe_max = E;
		}
		ricperr[i][j] = (P > 0.0) ? fabs(Pfit/P - 1.0) : 0.0;
		if (ricperr[i][j] > ricpe_max)
		{
		    ricpe_max = ricperr[i][j];
		    RHO_at_ricpe_max = RHO;
		    T_at_ricpe_max = TEMP;
		    P_at_ricpe_max = P;
		    E_at_ricpe_max = E;
		}

	    }
	}
	rpe_max *= 100.0;
	ree_max *= 100.0;
	ricpe_max *= 100.0;

	pe = ee = rpe = ree = pbar = ebar = icpe = ricpe = 0.0;
	for (i = 0; i < nr-1; i++)
	{
	    DRHO = rho[i+1] - rho[i];
	    for (j = 0; j < nt-1; j++)
	    {
		DTEMP = T[j+1] - T[j];
		pcell2 = 0.25*(sqr(p[i][j  ]) + sqr(p[i+1][j  ]) +
		               sqr(p[i][j+1]) + sqr(p[i+1][j+1]));
		pbar += DRHO*DTEMP*pcell2;

		err = 0.25*(sqr(perr[i][j  ]) + sqr(perr[i+1][j  ]) +
		            sqr(perr[i][j+1]) + sqr(perr[i+1][j+1]));
		pe += DRHO*DTEMP*err;

		err = 0.25*(sqr(rperr[i][j  ]) + sqr(rperr[i+1][j  ]) +
		            sqr(rperr[i][j+1]) + sqr(rperr[i+1][j+1]));
		rpe += DRHO*DTEMP*err;

		ecell2 = 0.25*(sqr(e[i][j  ]) + sqr(e[i+1][j  ]) +
		               sqr(e[i][j+1]) + sqr(e[i+1][j+1]));

		ebar += DRHO*DTEMP*ecell2;
		err = 0.25*(sqr(eerr[i][j  ]) + sqr(eerr[i+1][j  ]) +
		            sqr(eerr[i][j+1]) + sqr(eerr[i+1][j+1]));
		ee += DRHO*DTEMP*err;

		err = 0.25*(sqr(eerr[i][j  ]) + sqr(eerr[i+1][j  ]) +
		            sqr(eerr[i][j+1]) + sqr(eerr[i+1][j+1]));
		ee += DRHO*DTEMP*err;

		err = 0.25*(sqr(reerr[i][j  ]) + sqr(reerr[i+1][j  ]) +
		            sqr(reerr[i][j+1]) + sqr(reerr[i+1][j+1]));
		ree += DRHO*DTEMP*err;

		err = 0.25*(sqr(icperr[i][j  ]) + sqr(icperr[i+1][j  ]) +
		            sqr(icperr[i][j+1]) + sqr(icperr[i+1][j+1]));
		icpe += DRHO*DTEMP*err;

		err = 0.25*(sqr(ricperr[i][j  ]) + sqr(ricperr[i+1][j  ]) +
		            sqr(ricperr[i][j+1]) + sqr(ricperr[i+1][j+1]));
		ricpe += DRHO*DTEMP*err;
	    }
	}

	DRHO = rho[nr-1]-rho[0];
	DTEMP = T[nt-1] - T[0];
	D = DRHO*DTEMP;

	pe = sqrt(pe/D);
	rpe = 100.0*sqrt(rpe/D);
	pbar = sqrt(pbar/D);
	npe = 100.0*pe/pbar;

	ee = sqrt(ee/D);
	ree = 100.0*sqrt(ree/D);
	ebar = sqrt(ebar/D);
	nee = 100.0*ee/ebar;

	icpe = sqrt(icpe/D);
	ricpe = 100.0*sqrt(ricpe/D);
	nicpe = 100.0*icpe/pbar;

	(void) fprintf(file,"%s\n"
	    "  nr = %d, nt = %d\n"
	    "  Gruneisen coefficient                    = %- 14.8g\n"
	    "  P_infinity                               = %- 14.8g GPa\n"
	    "  E_infinity                               = %- 14.8g kJ/g\n"
	    "  C_V                                      = %- 14.8g kJ/gK\n"
	    "  et                                       = %- 14.8g kJ/g\n"
	    "  R                                        = %- 14.8g kJ/gK\n",
	    mesg,nr,nt,GAM,Pinf,Einf,cv,et,GAM*cv);

	(void) fprintf(file,"\nINCOMPLETE EOS PRESSURE ERROR\n"
	    "  L2 norm            = %- 14.8g GPa\n"
	    "  normalized L2 norm = %- 14.8g %%\n"
	    "  maximum            = %- 14.8g GPa\n"
	    "  at\n"
	    "    rho              = %- 14.8g gram/cc\n"
	    "    T                = %- 14.8g K\n"
	    "    P                = %- 14.8g GPa\n"
	    "    E                = %- 14.8g kJ/g\n",
	    icpe,nicpe,icpe_max,
	    RHO_at_icpe_max,T_at_icpe_max,P_at_icpe_max,E_at_icpe_max);

	(void) fprintf(file,"\nINCOMPLETE EOS RELATIVE PRESSURE ERROR\n"
	    "  L2 norm            = %- 14.8g %%\n"
	    "  maximum            = %- 14.8g %%\n"
	    "  at\n"
	    "    rho              = %- 14.8g gram/cc\n"
	    "    T                = %- 14.8g K\n"
	    "    P                = %- 14.8g GPa\n"
	    "    E                = %- 14.8g kJ/g\n",
	    ricpe,ricpe_max,
	    RHO_at_ricpe_max,T_at_ricpe_max,P_at_ricpe_max,E_at_ricpe_max);

	(void) fprintf(file,"\nPRESSURE ERROR\n"
	    "  L2 norm            = %- 14.8g GPa\n"
	    "  normalized L2 norm = %- 14.8g %%\n"
	    "  maximum            = %- 14.8g GPa\n"
	    "  at\n"
	    "    rho              = %- 14.8g gram/cc\n"
	    "    T                = %- 14.8g K\n"
	    "    P                = %- 14.8g GPa\n"
	    "    E                = %- 14.8g kJ/g\n",
	    pe,npe,pe_max,RHO_at_pe_max,T_at_pe_max,P_at_pe_max,E_at_pe_max);

	(void) fprintf(file,"\nRELATIVE PRESSURE ERROR\n"
	    "  L2 norm            = %- 14.8g %%\n"
	    "  maximum            = %- 14.8g %%\n"
	    "  at\n"
	    "    rho              = %- 14.8g gram/cc\n"
	    "    T                = %- 14.8g K\n"
	    "    P                = %- 14.8g GPa\n"
	    "    E                = %- 14.8g kJ/g\n",
	    rpe,rpe_max,RHO_at_rpe_max,T_at_rpe_max,P_at_rpe_max,E_at_rpe_max);

	(void) fprintf(file,"\nENERGY ERROR\n"
	    "  L2 norm            = %- 14.8g kJ/g\n"
	    "  normalized L2 norm = %- 14.8g %%\n"
	    "  maximum            = %- 14.8g kJ/g\n"
	    "  at\n"
	    "    rho              = %- 14.8g gram/cc\n"
	    "    T                = %- 14.8g K\n"
	    "    P                = %- 14.8g GPa\n"
	    "    E                = %- 14.8g kJ/g\n",
	    ee,nee,ee_max,RHO_at_ee_max,T_at_ee_max,P_at_ee_max,E_at_ee_max);

	(void) fprintf(file,"\nRELATIVE ENERGY ERROR\n"
	    "  L2 norm            = %- 14.8g %%\n"
	    "  maximum            = %- 14.8g %%\n"
	    "  at\n"
	    "    rho              = %- 14.8g gram/cc\n"
	    "    T                = %- 14.8g K\n"
	    "    P                = %- 14.8g GPa\n"
	    "    E                = %- 14.8g kJ/g\n",
	    ree,ree_max,RHO_at_ree_max,T_at_ree_max,P_at_ree_max,E_at_ree_max);

	free_these(6,perr,rperr,eerr,reerr,icperr,ricperr);
}		/*end least_squares_spoly_fit*/

#if defined(PHASE_CODE)

enum { NINT  = 5 }; /* number of regions for integration */

LOCAL	void set_ph_RS_riv_state(
	SESAME_EOS	*seos)
{
	Front 		*fr = seos->fr[SESAME_RHO_ENTROPY];
	Wave 		*wv = seos->wave[SESAME_RHO_ENTROPY];
	INTERFACE 	*intfc = fr->interf;
	CURVE 		**cur;
	BOND 		*b;
	COMPONENT 	comp, compold;
	Locstate 	state, rstate, lstate;
	double		ABS_SES_EPS = seos->ABS_SES_EPS;
	double 		drho;
	double 		S,p,rho,c,riv;
	double 		rho0,cor0;
	double 		rhold,corold;
	double 		sign;
	double 		coords[MAXD];
	double 		dr, cl, cr;
	double		cref = Reference_sound_speed(seos);
	int 		j;
	int 		icoords[MAXD];
	int 		nx = fr->rect_grid->gmax[0];
	int 		ny = fr->rect_grid->gmax[1];
	int 		ix, iy;
	static const	double	RIVEPS = 1.0e-05; /* TOLERANCE */

	DEBUG_ENTER(set_ph_RS_riv_state)

	set_smax(seos);
	drho = (fr->rect_grid->U[0] - fr->rect_grid->L[0])/((nx-1)*NINT);

	/*
	* TODO: the code below may cause problems if the diagonal 
	* (Smax,rmin)->(Smin,rmax) does not remain in the computational domain.
	*/
	

		/*Compute riv on boundary curves*/

	for (cur = intfc->curves; cur && *cur; cur++)
	{
	    for (b = (*cur)->first; b != NULL; b = b->next) 
	    {

	        S = Coords(b->start)[1];    rho = Coords(b->start)[0];

	        lstate = left_state_at_point_on_curve(b->start,b,*cur);
	        rstate= right_state_at_point_on_curve(b->start,b,*cur);

	        get_state_on_ref_curve(fr,wv,S,&rho0,&cor0,&comp);

	    
	        if(fabs(rho - rho0) < RIVEPS*drho)
	        {
	            ses_rs_riv(lstate) = 0.0;
	            ses_rs_riv(rstate) = 0.0;
		    if (DEBUG)
		    {
		        (void) printf("Left/Right state Riv is zero in %s\n\t",
				      "set_ph_RS_riv_state()");
		        (void) printf("rho = %g, rho0 = %g, S = %g\n",
				      rho,rho0,S);
		    }
	            continue;
	        }


	        sign = (rho0 > rho) ? -1.0 : 1.0;
	        rhold = rho0;
	        corold = cor0;
	        compold = comp;
		riv = 0.0;
	        for (j = 0; j < nx*NINT; j++)
	        {
	            rho0 = rho0 + sign*j*drho;
	            if(sign*(rho - rho0) < 0)
	            {
	                /* At end of region of intergration */
	                rho0 = rho;
			p = ses_rs_redp(rstate)*ses_rs_rho_from_grid(rho0,seos)*
				ses_rs_temp_from_var(ses_rs_Tvar(rstate),seos) +
				ses_rs_coldp(rstate);
			if (wave_type(*cur) != PHASE_BOUNDARY)
			{
			    c = ses_rs_adb_gam(rstate)*p/
				ses_rs_rho_from_grid(rho0,seos);
			    c = pow(c,0.5);
			    c = c/ses_rs_rho_from_grid(rho0,seos);
			}
	                else 
	                {
			    if (compold == COMP_MIXED_PHASE)
				c = ses_rs_adb_gam(rstate)*p/
					ses_rs_rho_from_grid(rho0,seos);
			    else
				c = ses_rs_adb_gam(lstate)*p/
					ses_rs_rho_from_grid(rho0,seos);
			    c = pow(c,0.5);
			    c = c/ses_rs_rho_from_grid(rho0,seos);
			}
			dr = ses_rs_rho_from_grid(rho0,seos) - 
				ses_rs_rho_from_grid(rhold,seos);
	                riv = riv + 0.5*(c + corold)*dr;
			c = c*ses_rs_rho_from_grid(rho0,seos);
			if (wave_type(*cur) == PHASE_BOUNDARY)
			{
			    p = ses_rs_redp(rstate)*ses_rs_rho_from_grid(rho0,seos)*
			    	ses_rs_temp_from_var(ses_rs_Tvar(rstate),seos) +
				ses_rs_coldp(rstate);
			    cr = ses_rs_adb_gam(rstate)*p/
					ses_rs_rho_from_grid(rho0,seos);
			    cl = ses_rs_adb_gam(lstate)*p/
					ses_rs_rho_from_grid(rho0,seos);

	                    ses_rs_riv(lstate) = 2.0*riv/(cl+cref);
	                    ses_rs_riv(rstate) = 2.0*riv/(cr+cref);
			}
			else
			{
	                    ses_rs_riv(lstate) = 2.0*riv/(c+cref);
	                    ses_rs_riv(rstate) = ses_rs_riv(lstate);
			}
			if (DEBUG)
			{
			    (void) printf("In set_ph_RS_riv_state()\n\t");
			    (void) printf("right state Riv = %g%s, ",
			                  ses_rs_riv(rstate),
			    	      (fabs(ses_rs_riv(rstate)) < ABS_SES_EPS) ?
			    	          " (zero)" : "");
			    (void) printf("left state Riv = %g%s\n\t",
			                  ses_rs_riv(lstate),
			    	      (fabs(ses_rs_riv(lstate)) < ABS_SES_EPS) ?
			    	         " (zero)" : "");
			    (void) printf("rho = %g, S = %g\n",rho,S);
			}
	                break;
	            }
	            else
	            {
	                coords[0] = rho0;    coords[1] = S;
	                get_riv_int(fr,wv,coords,sign,&compold,
	                                &rho0,&rhold,&corold,&riv);
	            }
	        }
	    }
	    b = (*cur)->last;
	    if(DEBUG)
	    	(void) printf("Setting state at end of curve\n");

	    lstate = left_state_at_point_on_curve(b->end,b,*cur);
	    rstate= right_state_at_point_on_curve(b->end,b,*cur);
	    S = Coords(b->end)[1];    rho = Coords(b->end)[0];

	    get_state_on_ref_curve(fr,wv,S,&rho0,&cor0,&comp);

	    
	    if(fabs(rho - rho0) < RIVEPS*drho)
	    {
		ses_rs_riv(lstate) = 0.0;
		ses_rs_riv(rstate) = 0.0;
		if (DEBUG)
		{
		    (void) printf("Left/Right state Riv is zero in %s\n\t",
		                  "set_ph_RS_riv_state()");
		    (void) printf("rho = %g, rho0 = %g, S = %g\n",rho,rho0,S);
		}
		continue;
	    }


	    sign = (rho0 > rho) ? -1.0 : 1.0;
	    rhold = rho0;
	    corold = cor0;
	    compold = comp;
	    riv = 0.0;
	    for (j = 0; j < nx*NINT; j++)
	    {
	        rho0 = rho0 + sign*j*drho;
	        if(sign*(rho - rho0) < 0)
	        {
	        /* At end of region of intergration */
	            rho0 = rho;
		    p = ses_rs_redp(rstate)*ses_rs_rho_from_grid(rho0,seos)*
				ses_rs_temp_from_var(ses_rs_Tvar(rstate),seos) +
				ses_rs_coldp(rstate);
		    if (wave_type(*cur) != PHASE_BOUNDARY)
		    {
			c = ses_rs_adb_gam(rstate)*p/
				ses_rs_rho_from_grid(rho0,seos);
			c = pow(c,0.5);
			c = c/ses_rs_rho_from_grid(rho0,seos);
		    }
	            else 
	            {
			if (compold == COMP_MIXED_PHASE)
			    c = ses_rs_adb_gam(rstate)*p/
				ses_rs_rho_from_grid(rho0,seos);
			else
			    c = ses_rs_adb_gam(lstate)*p/
				ses_rs_rho_from_grid(rho0,seos);
	                c = pow(c,0.5);
	                c = c/ses_rs_rho_from_grid(rho0,seos);
		    }
		    dr = ses_rs_rho_from_grid(rho0,seos) - 
			ses_rs_rho_from_grid(rhold,seos);
	            riv = riv + 0.5*(c + corold)*dr;
		    c = c*ses_rs_rho_from_grid(rho0,seos);
		    ses_rs_riv(rstate) = 2.0*riv/(c+cref);
		    ses_rs_riv(lstate) = ses_rs_riv(rstate);
		    if (DEBUG)
		    {
		        (void) printf("In set_ph_RS_riv_state()\n\t");
		        (void) printf("right state Riv = %g%s, ",
		                      ses_rs_riv(rstate),
		                      (fabs(ses_rs_riv(rstate)) < ABS_SES_EPS) ?
		        	      " (zero)" : "");
		        (void) printf("left state Riv = %g%s\n\t",
		                      ses_rs_riv(lstate),
		                      (fabs(ses_rs_riv(lstate)) < ABS_SES_EPS) ?
		        	      " (zero)" : "");
		        (void) printf("rho = %g, S = %g\n",rho,S);
		    }
	            break;
		}
	        else
	        {
		    coords[0] = rho0;    coords[1] = S;
	            get_riv_int(fr,wv,coords,sign,&compold,
	                         &rho0,&rhold,&corold,&riv);
	        }
	    }
	}

/* Initialize riv at interior states */
	if (DEBUG) (void) printf("Initializing interior states\n");

	for(ix = 0; ix < nx; ix++)
	{
	    for(iy = 0; iy < ny; iy++)
	    {
	        icoords[0] = ix;    icoords[1] = iy;
	        comp = Rect_comp(icoords,wv);
		rho = Rect_coords(icoords,wv)[0];
		S   = Rect_coords(icoords,wv)[1];
	        state = Rect_state(icoords,wv);
	        if (is_exterior_comp(comp,intfc))
	        {
	            ses_rs_riv(state) = ERROR_FLOAT;
	        }
	        else
	        {
	       	    get_state_on_ref_curve(fr,wv,S,&rho0,&cor0,&comp);

	            set_ses_intrp_flag_all(SESAME_RHO_ENTROPY);
	        
	        
	            if(fabs(rho - rho0) < RIVEPS*drho)
	            {
	                ses_rs_riv(state) = 0.0;
		        if (DEBUG)
		        {
		            (void) printf("Grid state (%d, %d) ",ix,iy);
			    (void) printf("Riv is zero in %s\n\t",
				          "set_ph_RS_riv_state()");
		            (void) printf("rho = %g, rho0 = %g, S = %g\n",
					  rho,rho0,S);
		        }
	                break;
	            }


	            sign = (rho0 > rho) ? -1.0 : 1.0;
	            rhold = rho0;
	            corold = cor0;
	       	    compold = comp;
		    riv = 0.0;
	            for (j = 0; j < nx*NINT; j++)
	            {
	                rho0 = rho0 + sign*j*drho;
	                if(sign*(rho - rho0) < 0)
	                {
	                    /* At end of region of intergration */
	                    rho0 = rho;
			    p = ses_rs_redp(state)*ses_rs_rho_from_grid(rho0,seos)*
			        ses_rs_temp_from_var(ses_rs_Tvar(state),seos) +
				ses_rs_coldp(state);
			    c = ses_rs_adb_gam(state)*p/
				ses_rs_rho_from_grid(rho0,seos);
	                    c = pow(c,0.5);
	                    c = c/ses_rs_rho_from_grid(rho0,seos);
			    dr = ses_rs_rho_from_grid(rho0,seos) - 
				ses_rs_rho_from_grid(rhold,seos);
	                    riv = riv + 0.5*(c + corold)*dr;
			    c = c*ses_rs_rho_from_grid(rho0,seos);
	                    ses_rs_riv(state) = 2.0*riv/(c+cref);
			    if (DEBUG)
			    {
			        (void) printf("In %s, Riv(%d, %d) = %g%s, ",
				              "set_ph_RS_riv_state()",
			                       ix,iy,ses_rs_riv(state),
			    	       (fabs(ses_rs_riv(state)) < ABS_SES_EPS) ?
			    	               " (zero)" : "");
			        (void) printf("rho = %g, S = %g\n",rho,S);
			    }
			    break;
	                }
	                else
	                {
	                    coords[0] = rho0;
	                    coords[1] = S;
	                    get_riv_int(fr,wv,coords,sign,&compold,&rho0,
				&rhold,&corold,&riv);
	                }
	            }
	        }
	    }
	}
	DEBUG_LEAVE(set_ph_RS_riv_state)
}		/*end set_ph_RS_riv_state*/

LOCAL	void get_riv_int(
	Front 		*fr,
	Wave 		*wv,
	double 		*coords,
	double 		sign,
	COMPONENT 	*compold,
	double 		*rho0,
	double 		*rhold,
	double 		*corold,
	double 		*riv)
{
	INTERFACE 	*intfc = fr->interf;
	SESAME_EOS	*seos = Ses_front_seos(fr);
	double 		rmin = fr->rect_grid->L[0];
	double 		rmax = fr->rect_grid->U[0];
	int 		nx = fr->rect_grid->gmax[0];
	double 		drho = (rmax-rmin)/((nx-1)*NINT);
	CURVE 		**cur, *phsbdry;
	COMPONENT 	comp;
	double 		S, rph, var[NUM_SES_VAR], c;
	double 		lrho0 = *rho0;
	double 		lrhold = *rhold;
	double 		lcoords[MAXD];
	double 		dr;

	set_ses_intrp_flag_all(SESAME_RHO_ENTROPY);

	S = coords[1];

	comp = nearest_interior_comp(YES,COMP_PURE_PHASE,coords,intfc);
	if (is_exterior_comp(*compold,intfc))
	{
	    lcoords[0] = lrhold;
	    lcoords[1] = S;
	    *compold = nearest_interior_comp(YES,COMP_PURE_PHASE,lcoords,intfc);
	}
	if(comp != *compold)
	{
	    /* Crossing phase boundary */
	    for(cur = intfc->curves; cur && *cur;cur++)
	    {
	    	if(wave_type(*cur) == PHASE_BOUNDARY)
	    		phsbdry = *cur;
	    }
	
	    get_rho_on_phase_bdry_at_entropy(S,lrho0,lrhold,&rph,phsbdry);
	    lcoords[0] = rph;
	    lcoords[1] = S;
	    if ((comp == COMP_PURE_PHASE) || (*compold == COMP_MIXED_PHASE))
	    	ses_solution(lcoords,*compold,Hyper_surf(phsbdry),POSITIVE_SIDE,
			     fr,wv,var);	
	    else
	    	ses_solution(lcoords,*compold,Hyper_surf(phsbdry),NEGATIVE_SIDE,
			     fr,wv,var);	
	    c = var[RS_RP]*ses_rs_rho_from_grid(rph,seos)*
	    	ses_rs_temp_from_var(var[RS_T],seos) + var[RS_CP];
	    c = var[RS_AG]*c/ses_rs_rho_from_grid(rph,seos);
	    c = pow(c,0.5);
	    c = c/ses_rs_rho_from_grid(rph,seos);
	    dr = ses_rs_rho_from_grid(rph,seos) -
		 ses_rs_rho_from_grid(*rhold,seos);
	    /* TODO: Would Simpson's Rule be better?*/
	    *riv = *riv + 0.5*(c + *corold)*dr;
	    if(comp == COMP_PURE_PHASE)
	    	ses_solution(lcoords,comp,Hyper_surf(phsbdry),NEGATIVE_SIDE,
			     fr,wv,var);	
	    else
	    	ses_solution(lcoords,comp,Hyper_surf(phsbdry),POSITIVE_SIDE,
			     fr,wv,var);	
	    c = var[RS_RP]*ses_rs_rho_from_grid(rph,seos)*
	    	ses_rs_temp_from_var(var[RS_T],seos) + var[RS_CP];
	    c = var[RS_AG]*c/ses_rs_rho_from_grid(rph,seos);
	    c = pow(c,0.5);
	    c = c/ses_rs_rho_from_grid(rph,seos);
	    *corold = c;
	    *rhold = rph;
	    *compold = comp;
	    *rho0 = *rho0 - sign*drho;
	    return;
	}
	else
	{
	    ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);	
	    c = var[RS_RP]*ses_rs_rho_from_grid(*rho0,seos)*
	    	ses_rs_temp_from_var(var[RS_T],seos) + var[RS_CP];
	    c = var[RS_AG]*c/ses_rs_rho_from_grid(*rho0,seos);
	    c = pow(c,0.5);
	    c = c/ses_rs_rho_from_grid(*rho0,seos);
	    dr = ses_rs_rho_from_grid(*rho0,seos) -
		ses_rs_rho_from_grid(*rhold,seos);
	    /* TODO: Would Simpson's Rule be better?*/
	    *riv = *riv + 0.5*(c + *corold)*dr;
	    *corold = c;
	    *rhold = *rho0;
	    *compold = comp;
	    return;
	}
}		/*end get_riv_int*/

/*
*		get_rho_on_phase_bdry_at_entropy():
*
*	Finds the density on the phase boundary at the point
*	with entropy S and density between rmin and rmax.
*/

LOCAL	void get_rho_on_phase_bdry_at_entropy(
	double 		S,
	double 		rmin,
	double 		rmax,
	double 		*rho,
	CURVE 		*phsbdry)
{
	BOND 		*b;
	double 		slope;

	for (b = phsbdry ->first; b != NULL; b = b->next) 
	{
		if(!Between(S,Coords(b->start)[1],Coords(b->end)[1]))
			continue;
		slope = (Coords(b->end)[0] - Coords(b->start)[0]) / 
			(Coords(b->end)[1] - Coords(b->start)[1]);
		*rho = slope*(S - Coords(b->start)[1]) + Coords(b->start)[0];
		if (Between(*rho,rmin,rmax))
			return;
	}
}		/*end get_rho_on_phase_bdry_at_entropy*/

/*
*
*			get_state_on_ref_curve()
*
*	Determine the density on the reference curve for computing the
*	Riemann invariant. The curve is given by x = y, where
*	x = (rho - rmin)/(rmax - rmin) and y = (Smax - S)/(Smax - Smin).
*
*/

LOCAL	void get_state_on_ref_curve(
	Front 		*fr,
	Wave 		*wv,
	double 		S,
	double 		*rho0,
	double 		*cor0,
	COMPONENT 	*comp)
{
	INTERFACE 	*intfc = fr->interf;
	SESAME_EOS	*seos = Ses_front_seos(fr);
	double 		Smax = fr->rect_grid->U[1];
	double 		Smin = fr->rect_grid->L[1];
	double 		rmin = fr->rect_grid->L[0];
	double 		rmax = fr->rect_grid->U[0];
	double 		y, var[NUM_SES_VAR];
	double 		coords[MAXD];


	y = (Smax - S)/(Smax - Smin);
	*rho0 = y*(rmax - rmin) + rmin;

	coords[0] = *rho0;	coords[1] = S;
	*comp = nearest_interior_comp(YES,COMP_PURE_PHASE,coords,intfc);
	set_ses_intrp_flag_all(SESAME_RHO_ENTROPY);
	ses_solution(coords,*comp,NULL,POSITIVE_SIDE,fr,wv,var);
	*cor0 = var[RS_RP]*ses_rs_rho_from_grid(*rho0,seos)*
		ses_rs_temp_from_var(var[RS_T],seos) + var[RS_CP];
	*cor0 = var[RS_AG]*(*cor0)/ses_rs_rho_from_grid(*rho0,seos);
	*cor0 = pow(*cor0,0.5);
	*cor0 = *cor0/ses_rs_rho_from_grid(*rho0,seos);
}		/*end get_state_on_ref_curve*/


LOCAL	void set_smax(
	SESAME_EOS 	*seos)
{
	BOND 		*b;
	CURVE 		**cur;
	INTERFACE 	*intfc = seos->fr[SESAME_RHO_ENTROPY]->interf;
	double 		s;

	DEBUG_ENTER(set_smax)

	for(cur = intfc->curves; cur && *cur; cur++)
	{
		if(wave_type(*cur) == PHASE_BOUNDARY)
		{
			s = Coords((*cur)->first->start)[1];
			for (b = (*cur)->first; b != NULL; b = b->next) 
			{
				if(Coords(b->end)[1] > s)
					s = Coords(b->end)[1];
			}
		}
	}
	seos->Smax = s;
	if (DEBUG) (void) printf("%g\n",seos->Smax);
	DEBUG_LEAVE(set_smax)
}		/*end set_smax*/
#endif /* defined(PHASE_CODE) */
#endif /* defined(SESAME_CODE) && defined(TWOD) */
