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
*				gprcur.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains routines for printing states on curves.
*
*/

#if defined(TWOD)

#include <gdecs/gdecs.h>


typedef struct {
	double rho, mom[MAXD], vel[MAXD];
	double vnorm, vtan, p, ent, ieng, temp, c, ener;
	double n_fast, n_slow;
	double t_fast, t_slow;
#if defined(COMBUSTION_CODE)
	double prod;
#endif /* defined(COMBUSTION_CODE) */
} PRINT_STATE;

	/* LOCAL Function Declarations */
LOCAL	void	load_state_info_for_printing(double*,double*,int,
	        			     PRINT_STATE*,Locstate);
LOCAL	void	print_wave_type_and_start_end_status(CURVE*);
LOCAL	void	sprint_position(char*,POINT*,int);

EXPORT void g_fprint_header_for_graph_curve_states(
	FILE		*file,
	Front		*front,
	const char	*message)
{
	int	   i, dim = front->rect_grid->dim;
	const char **s;
	static const char *Cname[3] = { "X", "Y", "Z" };
	static const char *LMname[3] = { "L_X_MOM", "L_Y_MOM", "L_Z_MOM" };
	static const char *RMname[3] = { "R_X_MOM", "R_Y_MOM", "R_Z_MOM" };
	static const char *LVname[3] = { "L_X_VEL", "L_Y_VEL", "L_Z_VEL" };
	static const char *RVname[3] = { "R_X_VEL", "R_Y_VEL", "R_Z_VEL" };
	static const char *othervars[] = {	"L_PRESSURE",	"R_PRESSURE",
	        			"L_ENTROPY",	"R_ENTROPY",
	        			"L_INT_ENG",	"R_INT_ENG",
	        			"L_TEMP",	"R_TEMP",
	        			"L_SOUND_SP",	"R_SOUND_SP",
	        			"L_NORM_VEL",	"R_NORM_VEL",
	        			"L_TANG_VEL",	"R_TANG_VEL",
	        			"VEL_SHEAR",
	        			"L_N_FASTCH",	"R_N_FASTCH",
	        			"L_T_FASTCH",	"R_T_FASTCH",
	        			"L_N_SLOWCH",	"R_N_SLOWCH",
	        			"L_T_SLOWCH",	"R_T_SLOWCH",
	        			"L_MACH_NUM",	"R_MACH_NUM",
	        			NULL};


	(void) foutput(file);
	(void) fprintf(file,"%-10s %-10s ","TIME","TIME_STEP");
	for (i = 0; i < dim; i++)
	    (void) fprintf(file,"%-10s ",Cname[i]);
	(void) fprintf(file,"%-10s ","ARCLENGTH");
	(void) fprintf(file,"%-10s %-10s %-10s %-10s ",
	               "L_DENS","R_DENS","L_ENERGY","R_ENERGY");
	for (i = 0; i < dim; i++)
	    (void) fprintf(file,"%-10s %-10s ",LMname[i],RMname[i]);
	for (i = 0; i < dim; i++)
	    (void) fprintf(file,"%-10s %-10s ",LVname[i],RVname[i]);
	for (s = othervars; *s != NULL; s++)
	    (void) fprintf(file,"%-10s ",*s);
#if defined(COMBUSTION_CODE)
	if (g_composition_type() == ZND)
	    (void) fprintf(file,"%-10s %-10s ","L_PROD_DENS","R_PROD_DENS");
#endif /* defined(COMBUSTION_CODE) */
	(void) fprintf(file,"%s\n",message);
}		/*end g_fprint_header_for_graph_curve_states*/

/*
*			g_fgraph_curve_states():
*
*	Prints coordinates, arclength, left and right density, left and right
*	energy density, left and right velocity and momentum, left and right 
*	normal and tangential velocity in a columnar format suitable
*	for graphs.
*
*	The starting arclength is *length; the arclength at the end
*	of the curve is returned through ft_assignment to *length.
*	Also note that a terminating newline should be printed by the caller.
*	(These two quirks are kludges so that successive curves may
*	effectively be joined onto one graph.)
*/

EXPORT void g_fgraph_curve_states(
	FILE		*file,
	CURVE		*c,
	Front		*front,
	double		*length)
{
	double		vars[50]; /*Enough to hold all data*/
	double		arclength,nor[MAXD],tngt[MAXD];
	PRINT_STATE	Psl, Psr;
	int             wtype = wave_type(c);
	int		i, dim = front->interf->dim;
	int		first_call_to_loop;
	int		n_vars;
	Locstate	sl,sr;
	BOND		*b;

	b = c->first;
	first_call_to_loop = YES;
	while (b != NULL) 
	{
	    BOND *next_b;
	    POINT *p;

	    if (first_call_to_loop)
	    {
	    	next_b = b;
	    	first_call_to_loop = NO;
	    	p = b->start;
	    	arclength = *length;
	    }
	    else
	    {
	    	next_b = b->next;
	    	p = b->end;
	    	arclength += bond_length(b);
	    }
	    normal(p,Hyper_surf_element(b),Hyper_surf(c),nor,front);
	    tngt[0] = -nor[1];	tngt[1] = nor[0];
	    slsr(p,Hyper_surf_element(b),Hyper_surf(c),&sl,&sr);
	    load_state_info_for_printing(nor,tngt,dim,&Psl,sl);
	    load_state_info_for_printing(nor,tngt,dim,&Psr,sr);
	    n_vars = 0;
	    vars[n_vars++] = front->time;
	    vars[n_vars++] = front->step;
	    for (i = 0; i < dim; i++)
	    	vars[n_vars++] = Coords(p)[i];
	    vars[n_vars++] = arclength;
	    vars[n_vars++] = Psl.rho;
	    vars[n_vars++] = Psr.rho;
	    vars[n_vars++] = Psl.ener;
	    vars[n_vars++] = Psr.ener;
	    for (i = 0; i < dim; i++)
	    {
	    	vars[n_vars++] = Psl.mom[i];
	    	vars[n_vars++] = Psr.mom[i];
	    }
	    for (i = 0; i < dim; i++)
	    {
	    	vars[n_vars++] = Psl.vel[i];
	    	vars[n_vars++] = Psr.vel[i];
	    }
	    vars[n_vars++] = Psl.p;
	    vars[n_vars++] = Psr.p;
	    vars[n_vars++] = Psl.ent;
	    vars[n_vars++] = Psr.ent;
	    vars[n_vars++] = Psl.ieng;
	    vars[n_vars++] = Psr.ieng;
	    vars[n_vars++] = Psl.temp;
	    vars[n_vars++] = Psr.temp;
	    vars[n_vars++] = Psl.c;
	    vars[n_vars++] = Psr.c;
	    vars[n_vars++] = Psl.vnorm;
	    vars[n_vars++] = Psr.vnorm;
	    vars[n_vars++] = Psl.vtan;
	    vars[n_vars++] = Psr.vtan;
	    vars[n_vars++] = Psl.vtan-Psr.vtan;
	    vars[n_vars++] = Psl.n_fast;
	    vars[n_vars++] = Psr.n_fast;
	    vars[n_vars++] = Psl.t_fast;
	    vars[n_vars++] = Psr.t_fast;
	    vars[n_vars++] = Psl.n_slow;
	    vars[n_vars++] = Psr.n_slow;
	    vars[n_vars++] = Psl.t_slow;
	    vars[n_vars++] = Psr.t_slow;
	    if (is_shock_wave(wtype))
	    {
	    	double dvn = fabs(Psl.vnorm-Psr.vnorm);
	    	double drho = fabs(Psl.rho - Psr.rho);
	    	vars[n_vars++] = (Psr.rho/drho)*(dvn/Psl.c);
	    	vars[n_vars++] = (Psl.rho/drho)*(dvn/Psr.c);
	    }
	    else if (is_rarefaction_wave(wtype))
	    {
	    	vars[n_vars++] = 1.0;
	    	vars[n_vars++] = 1.0;
	    }
	    else if (is_scalar_wave(wtype))
	    {
	    	vars[n_vars++] = 0.0;
	    	vars[n_vars++] = 0.0;
	    }
	    else
	    {
	    	vars[n_vars++] = -HUGE_VAL;
	    	vars[n_vars++] = -HUGE_VAL;
	    }
#if defined(COMBUSTION_CODE)
	    if (g_composition_type() == ZND)
	    {
	    	vars[n_vars++] = Psl.prod;
	    	vars[n_vars++] = Psr.prod;
	    }
#endif /* defined(COMBUSTION_CODE) */
	    if (is_binary_output() == YES)
	    {
	    	(void) fprintf(file,"\f%c",n_vars);
	    	(void) fwrite((const void *)vars,sizeof(double),n_vars,file);
	    }
	    else
	    {
	    	(void) fprintf(file,"%"FFMT,vars[0]);
	    	for (i = 1; i < n_vars; i++)
	    	    (void) fprintf(file," %"FFMT,vars[i]);
	        (void) fprintf(file,"\n");
	    }
	    b = next_b;
	}
	*length = arclength;
}		/*end g_fgraph_curve_states*/

LOCAL void load_state_info_for_printing(
	double		*nor,
	double		*tngt,
	int		dim,
	PRINT_STATE	*pst,
	Locstate	st)
{
	static Locstate state = NULL;
	int		i;


	if (is_obstacle_state(st)) 
	{
	    zero_scalar(pst,sizeof(PRINT_STATE));
	    return;
	}

	if (state == NULL) 
	    (*Params(st)->_alloc_state)(&state,sizeof(VGas));

	set_state(state,VGAS_STATE,st);
	pst->vnorm = pst->vtan = 0.0;
	for (i = 0; i < dim; i++)
	{
	    pst->vel[i] = vel(i,state);
	    pst->mom[i] = Mom(st)[i];
	    pst->vnorm += nor[i] * pst->vel[i];
	    pst->vtan += tngt[i] * pst->vel[i];
	}
	pst->rho = Dens(state);
	pst->p = pressure(state);
	pst->ieng = specific_internal_energy(state);
	pst->ener = Energy(st);
	pst->ent = entropy(state);
	pst->temp = temperature(state);
	pst->c = sound_speed(state);
	pst->n_fast = pst->vnorm + pst->c;
	pst->n_slow = pst->vnorm - pst->c;
	pst->t_fast = pst->vtan + pst->c;
	pst->t_slow = pst->vtan - pst->c;
#if defined(COMBUSTION_CODE)
	pst->prod = (Composition_type(st) == ZND) ? Prod(st) : 0.0;
#endif /* defined(COMBUSTION_CODE) */
}		/*end load_state_info_for_printing*/

EXPORT void print_header_for_self_similar_states_along_curve(
	FILE		*file,
	char		*title,
	int		dim)
{
	static char	Cname[3][3] = { "X", "Y", "Z" };
	static char	LMname[3][8] = { "L_X_MOM", "L_Y_MOM", "L_Z_MOM" };
	static char	RMname[3][8] = { "R_X_MOM", "R_Y_MOM", "R_Z_MOM" };
	int		i;

	(void) foutput(file);
	(void) fprintf(file,"%-10s %-10s ","TIME","TIME_STEP");
	for (i = 0; i < dim; i++)
	    (void) fprintf(file,"%-10s ",Cname[i]);
	(void) fprintf(file,"%-10s ","ARCLENGTH");
	(void) fprintf(file,"%-10s %-10s %-10s %-10s ",
	               "L_DENS","R_DENS","L_ENERGY","R_ENERGY");
	for (i = 0; i < dim; i++)
	    (void) fprintf(file,"%-10s %-10s ",LMname[i],RMname[i]);
	(void) fprintf(file,"%-10s %-10s %-10s %-10s ",
	               "L_PRESSURE","R_PRESSURE","L_ENTROPY","R_ENTROPY");
	(void) fprintf(file,"%-10s %-10s %-10s %-10s ",
	               "L_TEMP","R_TEMP","L_SOUND_SP","R_SOUND_SP");
	(void) fprintf(file,"%-10s %-10s %-10s %-10s ",
	               "L_N_VEL","R_N_VEL","L_T_VEL","R_T_VEL");
	(void) fprintf(file,"%-10s %-10s %-10s %-10s ",
	               "L_N_SS_VEL","R_N_SSVEL","L_T_SS_VEL","R_T_SS_VEL");
	(void) fprintf(file,"%-10s %-10s %-15s %-15s ",
	               "L_MACH_NUM","R_MACH_NUM",
	               "L_SS_MACH_NUM","R_SS_MACH_NUM");
	(void) fprintf(file,"%s\n",title);

	/* TODO: code needed for combustion */
}		/*end print_header_for_self_similar_states_along_curve*/

/*
*		print_self_similar_front_states_along_curve():
*
*	Prints coordinates, arclength, left and right density, left and right
*	energy density, left and right momentuma,
*	left and right normal and tangential regular and self-similar velocity 
*	along a curve in a columnar format suitable for "graphs".
*
*	The starting arclength is *length; the arclength at the end
*	of the curve is returned through ft_assignment to *length.
*	Also note that a terminating newline should be printed by the caller.
*	(These two quirks are kludges so that successive curves may
*	effectively be joined onto one graph.)
*
*	Currently, this function is only used for RAMP_REFLECTION problems.
*	For future applications, beware of code specific to above.
*/

EXPORT void print_self_similar_front_states_along_curve(
	FILE		*file,
	CURVE		*c,
	Front		*front,
	double		*ss_origin,
	double		time_elapsed,
	double		*length)
{
	double		arclength;
	double		nor[MAXD];
	double		tngt[MAXD];
	double		tmp[MAXD];
	PRINT_STATE	Psl, Psr;
	double		ssvnorm, ssvtan;
	double		Mnumr, Mnuml;
	double		ssMnumr, ssMnuml;
	double		ssvnorml, ssvnormr;
	double		ssvtanl, ssvtanr;
	int		i, dim = front->interf->dim;
	int		first_call_to_loop = YES;
	Locstate	sl,sr;
	BOND		*b;

	    /* TODO: code needed for combustion */

	b = c->first;
	while (b != NULL)
	{
	    BOND *next_b;
	    POINT *p;
	    if (first_call_to_loop)
	    {
	    	next_b = b;
	    	first_call_to_loop = NO;
	    	p = b->start;
	    	arclength = *length;
	    }
	    else
	    {
	    	next_b = b->next;
	    	p = b->end;
	    	arclength += bond_length(b);
	    }
	    normal(p,Hyper_surf_element(b),Hyper_surf(c),nor,front);
	    tngt[0] = -nor[1];	tngt[1] = nor[0];
	    slsr(p,Hyper_surf_element(b),Hyper_surf(c),&sl,&sr);
	    ssvnorm = ssvtan = 0.0;
	    for (i = 0; i < dim; i++)
	    {
	    	ssvnorm += nor[i] * (Coords(b->start)[i]-ss_origin[i]);
	    	ssvtan  += tngt[i] * (Coords(b->start)[i]-ss_origin[i]);
	    }
	    ssvnorm /= time_elapsed;
	    ssvtan  /= time_elapsed;

	    load_state_info_for_printing(nor,tngt,dim,&Psl,sl);
	    if (is_obstacle_state(sl))
	    {
	    	ssvnorml = 0.;	ssvtanl = 0.;
	    	ssMnuml = Mnuml = 0.;
	    }
	    else
	    {
	    	ssvnorml = Psl.vnorm - ssvnorm;
	    	ssvtanl = Psl.vtan - ssvtan;

	    	tmp[0] = Psl.vnorm;	tmp[1] = Psl.vtan;
	    	Mnuml = mag_vector(tmp,dim) / Psl.c;

	    	tmp[0] = ssvnorml; tmp[1] = ssvtanl;
	    	ssMnuml = mag_vector(tmp,dim) / Psl.c;
	    }

	    load_state_info_for_printing(nor,tngt,dim,&Psr,sr);
	    if (is_obstacle_state(sr))
	    {
	    	ssvnormr = 0.;	ssvtanr = 0.;
	    	ssMnumr = Mnumr = 0.;
	    }
	    else
	    {
	    	ssvnormr = Psr.vnorm - ssvnorm;
	    	ssvtanr = Psr.vtan - ssvtan;

	    	tmp[0] = Psr.vnorm;	tmp[1] = Psr.vtan;
	    	Mnumr = mag_vector(tmp,dim) / Psr.c;

	    	tmp[0] = ssvnormr; tmp[1] = ssvtanr;
	    	ssMnumr = mag_vector(tmp,dim) / Psr.c;
	    }

	    (void) fprintf(file,"%g %d ",front->time,front->step);
	    for(i = 0; i < dim; i++)
	    	(void) fprintf(file,"%g ",Coords(p)[i]);
	    (void) fprintf(file,"%g ",arclength);
	    (void) fprintf(file,"%g %g %g %g ",
	    	           Psl.rho,Psr.rho,Psl.ener,Psr.ener);
	    for(i = 0; i < dim; i++)
	    	(void) fprintf(file,"%g %g ",Psl.mom[i],Psr.mom[i]);
	    fprint_line_of_floats(file,20,
	        	          Psl.p,Psr.p,
	        		  Psl.ent,Psr.ent,
	        		  Psl.temp,Psr.temp,
	        	          Psl.c,Psr.c,
	        		  Psl.vnorm,Psr.vnorm,
	        		  Psl.vtan,Psr.vtan,
	        	          ssvnorml,ssvnormr,
	        		  ssvtanl,ssvtanr,
	        	          Mnuml,Mnumr,
	        		  ssMnuml,ssMnumr);
	    b = next_b;
	}
	*length = arclength;

	    /* mark endpoint (for corner point in reflection problems) */

	b = c->last;
	for(i = 0; i < dim; i++)
	    (void) fprintf(file,"%g ",Coords(b->end)[i]);
	(void) fprintf(file,"%g ",arclength);
	(void) fprintf(file,"%g %g %g %g ",
	               1.01*Psl.rho,1.01*Psr.rho,
	               1.01*Psl.ener,1.01*Psr.ener);
	for(i = 0; i < dim; i++)
	    (void) fprintf(file,"%g %g ",1.01*Psl.mom[i],1.01*Psr.mom[i]);
	fprint_line_of_floats(file,20,
	                      1.01*Psl.p,1.01*Psr.p,
	        	      1.01*Psl.ent,1.01*Psr.ent,
	                      1.01*Psl.temp,1.01*Psr.temp,
	        	      1.01*Psl.c,1.01*Psr.c,
	                      1.01*Psl.vnorm,1.01*Psr.vnorm,
	        	      1.01*Psl.vtan,1.01*Psr.vtan,
	                      1.01*ssvnorml,1.01*ssvnormr,
	        	      1.01*ssvtanl,1.01*ssvtanr,
	                      1.01*Mnuml,1.01*Mnumr,
	        	      1.01*ssMnuml,1.01*ssMnumr);
	for(i = 0; i < dim; i++)
	    (void) fprintf(file,"%g ",Coords(b->end)[i]);
	(void) fprintf(file,"%g ",arclength);
	(void) fprintf(file,"%g %g %g %g ",
	               0.99*Psl.rho,0.99*Psr.rho,
	               0.99*Psl.ener,0.99*Psr.ener);
	for(i = 0; i < dim; i++)
	    (void) fprintf(file,"%g %g ",0.99*Psl.mom[i],0.99*Psr.mom[i]);
	fprint_line_of_floats(file,20,
	                      0.99*Psl.p,0.99*Psr.p,
	        	      0.99*Psl.ent,0.99*Psr.ent,
	                      0.99*Psl.temp,0.99*Psr.temp,
	        	      0.99*Psl.c,0.99*Psr.c,
	                      0.99*Psl.vnorm,0.99*Psr.vnorm,
	        	      0.99*Psl.vtan,0.99*Psr.vtan,
	                      0.99*ssvnorml,0.99*ssvnormr,
	        	      0.99*ssvtanl,0.99*ssvtanr,
	                      0.99*Mnuml,0.99*Mnumr,
	        	      0.99*ssMnuml,0.99*ssMnumr);

}		/*end print_self_similar_front_states_along_curve*/

EXPORT void verbose_print_curve_states(
	CURVE		*curve)
{
	char		mesg[120];
	char		posn[120];
	int		dim = curve->interface->dim;
	BOND		*bb;

	if (curve == NULL) return;
	(void) printf("Points and states on curve %llu\n",curve_number(curve));
	print_wave_type_and_start_end_status(curve);
	sprint_position(posn,curve->start->posn,dim);
	(void) sprintf(mesg,"Left state at %s\n",posn);
	verbose_print_state(mesg,left_start_state(curve));
	(void) sprintf(mesg,"Right state at %s\n",posn);
	verbose_print_state(mesg,right_start_state(curve));
	for (bb = curve->first; bb != NULL; bb = bb->next) 
	{
	    sprint_position(posn,bb->end,dim);
	    (void) sprintf(mesg,"Left state at %s\n",posn);
	    verbose_print_state(mesg,
	        left_state_at_point_on_curve(bb->end,bb,curve));
	    (void) sprintf(mesg,"Right state at %s\n",posn);
	    verbose_print_state(mesg,
	        right_state_at_point_on_curve(bb->end,bb,curve));
	}
	(void) printf("End of states on curve %llu\n",curve_number(curve));
}		/*end verbose_print_curve_states*/

LOCAL void sprint_position(
	char		*posn,
	POINT		*p,
	int		dim)
{
	char		coords[80];
	int		i;

	(void) strcpy(posn,"(");
	for (i = 0; i < dim; i++)
	{
	    (void) sprintf(coords,"%g",Coords(p)[i]);
	    (void) strcat(posn,coords);
	    if (i < (dim - 1)) (void) strcat(posn,", ");
	}
	(void) strcat(posn,")");
}		/*end sprint_position*/

#if defined(DEBUG_NODE_PROPAGATE)
EXPORT void verbose_print_bond_states(
	const char	*message,
	BOND		*b,
	CURVE		*c)
{
	int		i, dim = c->interface->dim;

	if (b == NULL) 
	{
	    (void) printf("bond = NULL\n");
	    return;
	}
	(void) printf("%s bond (%llu): ",message,bond_number(b,c->interface));
	for (i = 0; i < dim; i++)
	    (void) printf("%g ",Coords(b->start)[i]);
	(void) printf("-> ");
	for (i = 0; i < dim; i++)
	    (void) printf("%g ",Coords(b->end)[i]);
	(void) printf("prev = %llu next = %llu\n",
	              bond_number(b->prev,c->interface),
	              bond_number(b->next,c->interface));
	verbose_print_state("left state b->start",
	    left_state_at_point_on_curve(b->start,b,c));
	verbose_print_state("right state b->start",
	    right_state_at_point_on_curve(b->start,b,c));
	verbose_print_state("left state b->end",
	    left_state_at_point_on_curve(b->end,b,c));
	verbose_print_state("right state b->end",
	    right_state_at_point_on_curve(b->end,b,c));
}		/*end verbose_print_bond_states*/
#endif /* defined(DEBUG_NODE_PROPAGATE) */

LOCAL void print_wave_type_and_start_end_status(
	CURVE		*curve)
{
	print_wave_type("\n\tcurve->wave_type = ",wave_type(curve),
	                "\n",curve->interface);
	print_curve_status("\tcurve->start_status = ",start_status(curve));
	print_curve_status("\tcurve->end_status = ",end_status(curve));
	(void) printf("\n");
}		/*end print_wave_type_and_start_end_status*/
#endif /* defined(TWOD) */
