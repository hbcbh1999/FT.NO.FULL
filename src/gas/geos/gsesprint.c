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
*				gsesprint.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains the interpolators for Sesame hyp solutions.
*
*
*	Contains functions for the printout of the eos tri solution
*	states.
*
*	See gsesintrp.c for a description of the sesame hyp solution
*	states.
*/

#if defined(SESAME_CODE) && defined(TWOD)
#include <geos/sesame.h>
#include <sys/types.h>
#include <time.h>

	/* LOCAL Function Prototypes */
LOCAL	void	print_ses_tri_solution_data(FILE*,Front*,Wave*,const char*,
					    const char*,size_t,
					    double (**)(double*,Front*,POINTER,
						       COMPONENT,Locstate),
					    SESAME_EOS*);

LOCAL	double	re_entropy(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	re_log_pressure(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	re_riv(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	re_log_temperature(double*,Front*,POINTER,COMPONENT,Locstate);

LOCAL	double	rt_entropy(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	rt_log_energy(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	rt_log_pressure(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	rt_riv(double*,Front*,POINTER,COMPONENT,Locstate);

LOCAL	double	rs_log_pressure(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	rs_log_energy(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	rs_log_temperature(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	rs_riv(double*,Front*,POINTER,COMPONENT,Locstate);

LOCAL	double	ps_log_density(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	ps_log_energy(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	ps_log_temperature(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	ps_riv(double*,Front*,POINTER,COMPONENT,Locstate);

LOCAL	void	verbose_ses_print_curve_states(CURVE*,SESAME_TABLE_TYPE,
					       SESAME_EOS*);
LOCAL	void	verbose_ses_print_state(const char*,double*,Locstate,
					SESAME_TABLE_TYPE,SESAME_EOS*);

LOCAL	double	vp_log_energy(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	vp_entropy(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	vp_log_temperature(double*,Front*,POINTER,COMPONENT,Locstate);
LOCAL	double	vp_riv(double*,Front*,POINTER,COMPONENT,Locstate);

EXPORT	void	print_title_for_sesame(
	FILE		*file,
	SESAME_EOS	*seos)
{
	char   title[1024];

	(void) sprintf(title,"\n\t\tWELCOME TO SESAME EOS DATA FILE\n"
	                    "\n\t\tMATERIAL ID %d\n",seos->ids2);
	print_title(file,title);

}		/*end print_title_for_sesame*/

EXPORT	void	fprint_SESAME_params(
	FILE		*file,
	SESAME_EOS	*seos)
{
	boolean	bio = is_binary_output();
	(void) fprintf(file,"\tEquation of state = %d SESAME\n",SESAME);
	(void) fprintf(file,"\tsesame library = %s\n",seos->seslib);
	(void) fprintf(file,"\tSesame material index = %d\n",seos->ids2);

	fwrite_float(file,"\tDiscreteness = ",seos->eps,bio,"%"FFMT,"\n");

	(void) fprintf(file,"\tDensity-Temperature window parameters\n");

	fwrite_float(file,"\t    minimum density = ",Rho_min(seos),bio,
		     "%"FFMT" g/cc",", ");
	fwrite_float(file,"maximum density = ",Rho_max(seos),bio,
		     "%"FFMT" g/cc",", ");
	fwrite_float(file,"reference density = ",Rho_ref(seos),bio,
		     "%"FFMT" g/cc","\n");

	fwrite_float(file,"\t    minimum temperature = ",
		     Temp_min(seos),bio,"%"FFMT" K"," ");
	fwrite_float(file,"maximum temperature = ",Temp_max(seos),
		     bio,"%"FFMT" K",", ");
	fwrite_float(file,"\t    reference temperature = ",
		     Temp_ref(seos),bio,"%"FFMT" K","\n");

	fwrite_float(file,"\t    reference sound speed = ",
		     Reference_sound_speed(seos),bio,"%"FFMT" km/sec","\n");

	(void) fprintf(file,"    \tdensity mesh = %d temperature mesh = %d\n",
		Nrho_hyp(seos),Ntemp_hyp(seos));

	fwrite_float(file,"\tAbsolute error = ",seos->abser0,bio,"%"FFMT," ");
	fwrite_float(file,"relative error = ",seos->reler0,bio,"%"FFMT,"\n");

	(void) fprintf(file,"\tPressure Limits\n");
	fwrite_float(file,"\t    minimum pressure = ",
		     Pressure_min(seos),bio,"%"FFMT" GPa",", ");
	fwrite_float(file,"maximum pressure = ",
		     Pressure_max(seos),bio,"%"FFMT" GPa",", ");
	fwrite_float(file,"reference pressure = ",
		     Pressure_ref(seos),bio,"%"FFMT" GPa","\n");

	(void) fprintf(file,"\tEnergy Limits\n");
	fwrite_float(file,"\t    minimum energy = ",
		     Energy_min(seos),bio,"%"FFMT" kJ/g",", ");
	fwrite_float(file,"maximum energy = ",
		     Energy_max(seos),bio,"%"FFMT" kJ/g",", ");
	fwrite_float(file,"reference energy = ",
		     Energy_ref(seos),bio,"%"FFMT" kJ/g","\n");

	(void) fprintf(file,"\tEntropy Limits\n");
	fwrite_float(file,"\t    minimum entropy = ",
		     Entropy_min(seos),bio,"%"FFMT," ");
	fwrite_float(file,"maximum entropy = ",Entropy_max(seos),bio,"%"FFMT,
		     "\n");
	(void) fprintf(file,"\tEntropy grid mapping parameters\n");
	fwrite_float(file,"\t    density-entropy entropy scale = ",
		     RS_entropy_scale(seos),bio,"%"FFMT,"\n");
	fwrite_float(file,"\t    density-entropy entropy shift = ",
		     RS_entropy_shift(seos),bio,"%"FFMT,"\n");
	fwrite_float(file,"\t    pressure-entropy entropy scale = ",
		     PS_entropy_scale(seos),bio,"%"FFMT,"\n");
	fwrite_float(file,"\t    pressure-entropy entropy shift = ",
		     PS_entropy_shift(seos),bio,"%"FFMT,"\n");

	(void) fprintf(file,"\tMultiphase eos = %s\n",
		       (multiphase_eos(seos) == YES) ? "YES" : "NO");
}		/*end fprint_SESAME_params*/

EXPORT	void	read_print_SESAME_params(
	SESAME_EOS    *seos,
	const IO_TYPE *io_type)
{
	char	s[10];
	FILE	*file = io_type->file;
	(void) fgetstring(file,"sesame library = ");
	(void) fscanf(file,"%s",seos->seslib);
	(void) fgetstring(file,"Sesame material index = ");
	(void) fscanf(file,"%d",&seos->ids2);

	seos->eps = fread_float("Discreteness = ",io_type);

	Rho_min(seos) = fread_float("minimum density = ",io_type);
	Rho_max(seos) = fread_float("maximum density = ",io_type);
	Rho_ref(seos) = fread_float("reference density = ",io_type);
	Temp_min(seos) = fread_float("minimum temperature = ",io_type);
	Temp_max(seos) = fread_float("maximum temperature = ",io_type);
	Temp_ref(seos) = fread_float("reference temperature = ",io_type);
	Reference_sound_speed(seos) =
		fread_float("reference sound speed = ",io_type);
	(void) fgetstring(file,"density mesh = ");
	(void) fscanf(file,"%d",&Nrho_hyp(seos));
	(void) fgetstring(file,"temperature mesh = ");
	(void) fscanf(file,"%d",&Ntemp_hyp(seos));
	seos->abser0 = fread_float("Absolute error = ",io_type);
	seos->reler0 = fread_float("relative error = ",io_type);

	Pressure_min(seos) = fread_float("minimum pressure = ",io_type);
	Pressure_max(seos) = fread_float("maximum pressure = ",io_type);
	Pressure_ref(seos) = fread_float("reference pressure = ",io_type);

	Energy_min(seos) = fread_float("minimum energy = ",io_type);
	Energy_max(seos) = fread_float("maximum energy = ",io_type);
	Energy_ref(seos) = fread_float("reference energy = ",io_type);

	Entropy_min(seos) = fread_float("minimum entropy = ",io_type);
	Entropy_max(seos) = fread_float("maximum entropy = ",io_type);
	RS_entropy_scale(seos) =
	    fread_float("density-entropy entropy scale = ",io_type);
	RS_entropy_shift(seos) =
	    fread_float("density-entropy entropy shift = ",io_type);
	PS_entropy_scale(seos) =
	    fread_float("pressure-entropy entropy scale = ",io_type);
	PS_entropy_shift(seos) =
	    fread_float("pressure-entropy entropy shift = ",io_type);

	(void) fgetstring(file,"Multiphase eos = ");
	(void) fscanf(file,"%s",s);
	multiphase_eos(seos) = (s[0] == 'Y') ? YES : NO;
}		/*end read_print_SESAME_params*/

EXPORT	void verbose_ses_show_intfc_states(
	INTERFACE		*intfc,
	SESAME_TABLE_TYPE	eos_type,
	SESAME_EOS		*seos)
{
	const char	*s;
	CURVE		**c;

	s = sesame_table_names()[eos_type];
	(void) printf("\t\tEOS STATES ON THE %s INTERFACE %llu\n\n",s,
		interface_number(intfc));
	for( c = intfc->curves;	 *c != NULL;  c++ )
		verbose_ses_print_curve_states(*c,eos_type,seos);
	(void) printf("\n\n");
	(void) printf("\t\tEND OF EOS STATES ON THE %s INTERFACE %llu\n\n",
	       s,interface_number(intfc));
}		/*end verbose_ses_show_intfc_states*/

LOCAL	void verbose_ses_print_state(
	const char 	  *message,
	double 		  *coords,
	Locstate 	  state,
	SESAME_TABLE_TYPE eos_type,
	SESAME_EOS	  *seos)
{
	double 		log_rho, log_T, log_p, log_e;
	double		red_p, red_e, cold_p, cold_e;
	double 		rho, T, p, e, adb_gam, gru_gam;
	double		S, c, stf_gam, ply_gam, R, p0, S0;

	switch (eos_type)
	{
	case SESAME_RHO_TEMP:
		rho = ses_rt_rho_from_grid(coords[0],seos);
		T = ses_rt_temp_from_grid(coords[1],seos);
		log_rho = log(rho);
		log_T = log(T);
		red_p = ses_rt_redp(state);
		red_e = ses_rt_rede(state);
		cold_p = ses_rt_coldp(state);
		cold_e = ses_rt_colde(state);
		p = rho*T*red_p + cold_p;
		e = T*red_e + cold_e;
		log_p = log(p);
		log_e = log(e);
		S = ses_rt_S(state);
		adb_gam = ses_rt_adb_gam(state);
		gru_gam = ses_rt_gru_gam(state);
		break;

	case SESAME_RHO_ENERGY:
		rho = ses_re_rho_from_grid(coords[0],seos);
		e = ses_re_engy_from_grid(coords[1],seos);
		log_rho = log(rho);
		log_e = log(e);
		red_p = ses_re_redp(state);
		cold_p = ses_re_coldp(state);
		T = ses_re_temp_from_var(ses_re_Tvar(state),seos);
		p = rho*T*red_p + cold_p;
		log_T = log(T);
		log_p = log(p);
		S = ses_re_S(state);
		adb_gam = ses_re_adb_gam(state);
		gru_gam = ses_re_gru_gam(state);
		break;

	case SESAME_PRESS_ENTROPY:
		p = ses_ps_press_from_grid(coords[0],seos);
		log_p = log(p);
		S = ses_ps_entpy_from_grid(coords[1],seos);
		red_e = ses_ps_rede(state);
		cold_e = ses_ps_colde(state);
		rho = ses_ps_rho_from_var(ses_ps_rho_var(state),seos);
		log_rho = log(rho);
		T = ses_ps_temp_from_var(ses_ps_Tvar(state),seos);
		e = T*red_e + cold_e;
		log_T = log(T);
		log_e = log(e);
		adb_gam = ses_ps_adb_gam(state);
		gru_gam = ses_ps_gru_gam(state);
		break;

	case SESAME_RHO_ENTROPY:
		rho = ses_rs_rho_from_grid(coords[0],seos);
		S = ses_rs_rho_from_grid(coords[1],seos);
		log_rho = log(rho);
		red_p = ses_rs_redp(state);
		red_e = ses_rs_rede(state);
		cold_p = ses_rs_coldp(state);
		cold_e = ses_rs_colde(state);
		T = ses_rs_temp_from_var(ses_rs_Tvar(state),seos);
		log_T = log(T);
		p = rho*T*red_p + cold_p;
		e = T*red_e + cold_e;
		log_p = log(p);
		log_e = log(e);
		adb_gam = ses_rs_adb_gam(state);
		gru_gam = ses_rs_gru_gam(state);
		break;

	default:
		screen("Unknown equation of state type ");
		screen("in ses_print_state()\n");
		clean_up(ERROR);
	}
	c = (adb_gam*p/rho);
	(void) printf("%s\n",message);
	(void) printf("\tdensity          = %-11g  log(density)       = %-11g\n",
		rho,log_rho);
	(void) printf("\tTemperature      = %-11g  log(Temperature)   = %-11g\n",
		T,log_T);
	(void) printf("\tpressure         = %-11g  log(pressure)      = %-11g\n",
		p,log_p);
	(void) printf("\treduced pressure = %-11g  cold pressure      = %-11g\n",
		red_p,cold_p);
	(void) printf("\tenergy           = %-11g  log(energy)        = %-11g\n",
		e,log_e);
	(void) printf("\treduced energy   = %-11g  cold energy        = %-11g\n",
		red_e,cold_e);
	(void) printf("\texp(entropy)     = %-11g  entropy            = %-11g\n",
		exp(S),S);
	(void) printf("\tadiabatic gamma  = %-11g  Gruneisen gamma = %-11g\n",
			adb_gam,gru_gam);
	(void) printf("\t\tPolytropic Parameters\n");
	ply_gam = 1.0 + p/(e*rho);
	R = p/(rho*T);
	S0 = S - R*(log_p - ply_gam*log_rho)/(ply_gam - 1.0);
	(void) printf("\tgamma            = %-11g  R                  = %-11g\n",
		ply_gam,R);
	(void) printf("\tS0               = %-11g\n",S0);
	(void) printf("\t\tStiffened Polytropic Parameters\n");
	p0 = p*((adb_gam/ply_gam) - 1.0)/(1.0 + c*c/(e*ply_gam));
	stf_gam = (p + e*rho)/(e*rho - p0);
	R = (p + p0)/(rho*T);
	S0 = S - R*(log(p + p0) - stf_gam*log_rho)/(stf_gam - 1.0);
	(void) printf("\tgamma            = %-11g  R                  = %-11g\n",
		stf_gam,R);
	(void) printf("\tp0               = %-11g  S0                 = %-11g\n",
		p0,S0);
	(void) printf("\n");
}		/*end verbose_ses_print_state*/


/*ARGSUSED*/
EXPORT	void	ses_rt_fprint_state_data(
	FILE		*file,
	Locstate	state,
	INTERFACE	*intfc)
{
	(void) fprintf(file,"cold_p    = %-14g, cold_e    = %-14g\n",
				ses_rt_coldp(state),ses_rt_colde(state));
	(void) fprintf(file,"reduced_p = %-14g, reduced_p = %-14g\n",
				ses_rt_redp(state),ses_rt_rede(state));
	(void) fprintf(file,"adb_gam   = %-14g, gru_gam   = %-14g\n",
				ses_rt_adb_gam(state),ses_rt_gru_gam(state));
	(void) fprintf(file,"S         = %-14g, riv       = %-14g\n",
				ses_rt_S(state),ses_rt_riv(state));
	(void) fprintf(file,"\n");
}		/*end ses_rt_fprint_state_data*/

/*ARGSUSED*/
EXPORT	void	ses_re_fprint_state_data(
	FILE		*file,
	Locstate	state,
	INTERFACE	*intfc)
{
	(void) fprintf(file,"cold_p    = %-14g, reduced_p = %-14g\n",
				ses_re_coldp(state),ses_re_redp(state));
	(void) fprintf(file,"T_var     = %-14g\n",ses_re_Tvar(state));
	(void) fprintf(file,"adb_gam   = %-14g, gru_gam   = %-14g\n",
				ses_re_adb_gam(state),ses_re_gru_gam(state));
	(void) fprintf(file,"S         = %-14g, riv       = %-14g\n",
				ses_re_S(state),ses_re_riv(state));
	(void) fprintf(file,"\n");
}		/*end ses_re_fprint_state_data*/

/*ARGSUSED*/
EXPORT	void	ses_rs_fprint_state_data(
	FILE		*file,
	Locstate	state,
	INTERFACE	*intfc)
{
	(void) fprintf(file,"cold_p    = %-14g, cold_e    = %-14g\n",
				ses_rs_coldp(state),ses_rs_colde(state));
	(void) fprintf(file,"reduced_p = %-14g, reduced_p = %-14g\n",
				ses_rs_redp(state),ses_rs_rede(state));
	(void) fprintf(file,"adb_gam   = %-14g, gru_gam   = %-14g\n",
				ses_rs_adb_gam(state),ses_rs_gru_gam(state));
	(void) fprintf(file,"T_var     = %-14g, riv       = %-14g\n",
				ses_rs_Tvar(state),ses_rs_riv(state));
	(void) fprintf(file,"\n");
}		/*end ses_rs_fprint_state_data*/

/*ARGSUSED*/
EXPORT	void	ses_ps_fprint_state_data(
	FILE		*file,
	Locstate	state,
	INTERFACE	*intfc)
{
	(void) fprintf(file,"cold_e    = %-14g, reduced_e = %-14g\n",
				ses_ps_colde(state),ses_ps_rede(state));
	(void) fprintf(file,"T_var     = %-14g, rho_var   = %-14g\n",
				ses_ps_Tvar(state),ses_ps_rho_var(state));
	(void) fprintf(file,"adb_gam   = %-14g, gru_gam   = %-14g\n",
				ses_rs_adb_gam(state),ses_rs_gru_gam(state));
	(void) fprintf(file,"riv       = %-14g\n",ses_rs_riv(state));
	(void) fprintf(file,"\n");
}		/*end ses_ps_fprint_state_data*/

/*ARGSUSED*/
EXPORT	void	ses_vp_fprint_state_data(
	FILE		*file,
	Locstate	state,
	INTERFACE	*intfc)
{
	(void) fprintf(file,"cold_e    = %-14g, reduced_e = %-14g\n",
				ses_vp_colde(state),ses_vp_rede(state));
	(void) fprintf(file,"T_var     = %-14g, S         = %-14g\n",
				ses_vp_Tvar(state),ses_vp_S(state));
	(void) fprintf(file,"adb_gam   = %-14g, gru_gam   = %-14g\n",
				ses_vp_adb_gam(state),ses_vp_gru_gam(state));
	(void) fprintf(file,"riv       = %-14g\n",ses_vp_riv(state));
	(void) fprintf(file,"\n");
}		/*end ses_ps_fprint_state_data*/

LOCAL void verbose_ses_print_curve_states(
	CURVE 		  *curve,
	SESAME_TABLE_TYPE eos_type,
	SESAME_EOS	  *seos)
{
	BOND 	*bb;

	(void) printf("EOS states on curve %llu\n",curve_number(curve));
	verbose_ses_print_state("Left state",Coords(curve->start->posn),
	      left_start_state(curve),eos_type,seos);
	verbose_ses_print_state("Right state",Coords(curve->start->posn),
	      right_start_state(curve),eos_type,seos);
	for (bb = curve->first; bb != NULL; bb = bb->next)
	{
	    verbose_ses_print_state("Left state",Coords(bb->end),
			            left_state_at_point_on_curve(bb->end,bb,
								 curve),
			            eos_type,seos);
	    verbose_ses_print_state("Right state",Coords(bb->end),
			            right_state_at_point_on_curve(bb->end,bb,
								  curve),
			            eos_type,seos);
	}
	(void) printf("End of EOS states on curve %llu\n",curve_number(curve));
}		/*end verbose_ses_print_curve_states*/

LOCAL	void	print_ses_tri_solution_data(
	FILE		*file,
	Front		*fr,
	Wave		*wv,
	const char	*tri_header,
	const char	*table_name,
	size_t		num_plots,
	double		(**plot_fns)(double*,Front*,POINTER,COMPONENT,Locstate),
	SESAME_EOS	*seos)
{
	RECT_GRID	*rgr;
	double		area;
	static const char	*FORMAT =
		"\n      stop_time = %-10g               stop_step = %-10d\n";

	if (file != stdout)
	{
	    record_print_version(file);
	    print_title_for_sesame(file,seos);
	}
	rgr = fr->rect_grid;
	(void) foutput(file);
	(void) fprintf(file,"\t\t\tINITIAL DATA:\n\n\n");
	fprint_rectangular_grid(file,rgr);
	(void) fprintf(file,FORMAT,0.0,0);
	area = (rgr->U[0] - rgr->L[0])*(rgr->U[1] - rgr->L[1]);
	(void) fprintf(file,"\n\t\tComputational Area = %g\n",area);
	(void) fprintf(file,"\n\t\tRemap Geometry:  IDENTITY_REMAP\n");
	(void) fprintf(file,"\n\t\tPrinting Interval:  1 mesh units\n\n\n\n\n");
	(void) fprintf(file,"\n\n\n\n");
	(void) foutput(file);
	(void) fprintf(file,"%s INTERFACE\n",table_name);
	(void) fprintf(file,"\n\n\n\n");
	(void) foutput(file);
	(void) fprintf(file,"\t\t\tFRONT DATA:\n");
	(void) fprintf(file," \n\t\t\tFront Rectangular Grid:\n\n");
	fprint_rectangular_grid(file,rgr);
	(void) fprintf(file,"\n\t\t\tInterface Topological Grid:\n\n");
	fprint_rectangular_grid(file,&topological_grid(fr->interf));
	fprint_interface(file,fr->interf);
	(void) fprintf(file,"\n\n\n\n");
	(void) foutput(file);
	(void) fprintf(file,"\t\t\tEND OF FRONT DATA:\n");

	(void) foutput(file);
	(void) fprintf(file,"\t\t\tSTATE DATA:\n");
	(void) foutput(file);
	(void) fprintf(file,"%s",tri_header);
	(void) fprintf(file,"\n");
	(void) fprintf(file,"#Point Source Data\n#0\n");
	print_tri_soln(file,fr,wv,wave_tri_soln(wv),num_plots,plot_fns);
	(void) fprintf(file,"END OF TRI_SOLN\n");
	(void) fprintf(file,"\n\n\n\n");
	(void) foutput(file);
	(void) fprintf(file,"\t\t\tEND OF STATE DATA\n");
}		/*end print_ses_tri_solution_data*/

EXPORT	void	print_rt_tri_soln(
	FILE		*file,
	SESAME_EOS	*seos)
{
	Front	   *fr = seos->fr[SESAME_RHO_TEMP];
	Wave	   *wv = seos->wave[SESAME_RHO_TEMP];
	double	   (*plot_fns[4])(double*,Front*,POINTER,COMPONENT,Locstate);
	static const char *tri_header =
		"TRI_SOLN: RT_LOG_PRESSURE RT_LOG_ENERGY RT_ENTROPY RT_RIV";

	plot_fns[0] = rt_log_pressure;
	plot_fns[1] = rt_log_energy;
	plot_fns[2] = rt_entropy;
	plot_fns[3] = rt_riv;

	print_ses_tri_solution_data(file,fr,wv,tri_header,
				    ses_table_name(seos,SESAME_RHO_TEMP),
				    4,plot_fns,seos);
}		/*end print_rt_tri_soln*/

EXPORT	void	print_re_tri_soln(
	FILE		*file,
	SESAME_EOS	*seos)
{
	Front	   *fr = seos->fr[SESAME_RHO_ENERGY];
	Wave	   *wv = seos->wave[SESAME_RHO_ENERGY];
	double	   (*plot_fns[4])(double*,Front*,POINTER,COMPONENT,Locstate);
	static const char *tri_header =
	    "TRI_SOLN: RE_LOG_PRESSURE RE_LOG_TEMPERATURE RE_ENTROPY RE_RIV";

	plot_fns[0] = re_log_pressure;
	plot_fns[1] = re_log_temperature;
	plot_fns[2] = re_entropy;
	plot_fns[3] = re_riv;

	print_ses_tri_solution_data(file,fr,wv,tri_header,
				    ses_table_name(seos,SESAME_RHO_ENERGY),
				    4,plot_fns,seos);
}		/*end print_re_tri_soln*/


EXPORT	void	print_rs_tri_soln(
	FILE		*file,
	SESAME_EOS	*seos)
{
	Front	   *fr = seos->fr[SESAME_RHO_ENTROPY];
	Wave	   *wv = seos->wave[SESAME_RHO_ENTROPY];
	double	   (*plot_fns[4])(double*,Front*,POINTER,COMPONENT,Locstate);
	static const char *tri_header =
	    "TRI_SOLN: RS_LOG_PRESSURE RS_LOG_TEMPERATURE RS_LOG_ENERGY RS_RIV";

	plot_fns[0] = rs_log_pressure;
	plot_fns[1] = rs_log_temperature;
	plot_fns[2] = rs_log_energy;
	plot_fns[3] = rs_riv;

	print_ses_tri_solution_data(file,fr,wv,tri_header,
				    ses_table_name(seos,SESAME_RHO_ENTROPY),
				    4,plot_fns,seos);
}		/*end print_rs_tri_soln*/


EXPORT	void	print_ps_tri_soln(
	FILE		*file,
	SESAME_EOS	*seos)
{
	Front	   *fr = seos->fr[SESAME_PRESS_ENTROPY];
	Wave	   *wv = seos->wave[SESAME_PRESS_ENTROPY];
	double	   (*plot_fns[4])(double*,Front*,POINTER,COMPONENT,Locstate);
	static const char *tri_header =
	    "TRI_SOLN: PS_LOG_ENERGY PS_LOG_TEMPERATURE PS_LOG_DENSITY PS_RIV";

	plot_fns[0] = ps_log_energy;
	plot_fns[1] = ps_log_temperature;
	plot_fns[2] = ps_log_density;
	plot_fns[3] = ps_riv;

	print_ses_tri_solution_data(file,fr,wv,tri_header,
				    ses_table_name(seos,SESAME_PRESS_ENTROPY),
				    4,plot_fns,seos);
}		/*end print_ps_tri_soln*/


EXPORT	void	print_vp_tri_soln(
	FILE		*file,
	SESAME_EOS	*seos)
{
	Front	   *fr = seos->fr[SESAME_VOLUME_PRESSURE];
	Wave	   *wv = seos->wave[SESAME_VOLUME_PRESSURE];
	double	   (*plot_fns[4])(double*,Front*,POINTER,COMPONENT,Locstate);
	static const char *tri_header =
	    "TRI_SOLN: VP_ENTROPY VP_LOG_TEMPERATURE VP_LOG_ENERGY VP_RIV";

	plot_fns[0] = vp_entropy;
	plot_fns[1] = vp_log_temperature;
	plot_fns[2] = vp_log_energy;
	plot_fns[3] = vp_riv;

	print_ses_tri_solution_data(file,fr,wv,tri_header,
				    ses_table_name(seos,SESAME_VOLUME_PRESSURE),
				    4,plot_fns,seos);
}		/*end print_vp_tri_soln*/



/*ARGSUSED*/
LOCAL double rt_log_pressure(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	return plog(ses_rt_redp(state)*
		ses_rt_temp_from_grid(coords[1],Ses_front_seos(front))*
		ses_rt_rho_from_grid(coords[0],Ses_front_seos(front)) +
		ses_rt_coldp(state));
}		/*end rt_log_pressure*/

/*ARGSUSED*/
LOCAL double rt_log_energy(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	double	T = ses_rt_temp_from_grid(coords[1],Ses_front_seos(front));
	return plog(ses_rt_rede(state)*T+ses_rt_colde(state));
}		/*end rt_log_energy*/

/*ARGSUSED*/
LOCAL double rt_entropy(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	return ses_rt_S(state);
}		/*end rt_entropy*/

/*ARGSUSED*/
LOCAL double rt_riv(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	SESAME_EOS	*seos = Ses_front_seos(front);
	double		cref = Reference_sound_speed(seos);
	double	rho = ses_rt_rho_from_grid(coords[0],seos);
	double	T = ses_rt_temp_from_grid(coords[1],seos);
	double	p = ses_rt_coldp(state) + ses_rt_redp(state)*rho*T;
	double	g = ses_rt_adb_gam(state);
	return 0.5*(sqrt(g*p/rho)+cref)*ses_rt_riv(state);
}		/*end rt_riv*/

/*ARGSUSED*/
LOCAL double re_log_pressure(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	SESAME_EOS	*seos = Ses_front_seos(front);
	double	T = ses_re_temp_from_var(ses_re_Tvar(state),seos);
	double	rho = ses_re_rho_from_grid(coords[0],seos);
	return plog(ses_re_redp(state)*T*rho + ses_re_coldp(state));
}		/*end re_log_pressure*/

/*ARGSUSED*/
LOCAL double re_log_temperature(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	SESAME_EOS	*seos = Ses_front_seos(front);
	return plog(ses_re_temp_from_var(ses_re_Tvar(state),seos));
}		/*end re_log_temperature*/

/*ARGSUSED*/
LOCAL double re_entropy(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	return ses_re_S(state);
}		/*end re_entropy*/

/*ARGSUSED*/
LOCAL double re_riv(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	SESAME_EOS	*seos = Ses_front_seos(front);
	double		cref = Reference_sound_speed(seos);
	double	rho = ses_re_rho_from_grid(coords[0],seos);
	double	T = ses_re_temp_from_var(ses_re_Tvar(state),seos);
	double	p = ses_re_coldp(state) + ses_re_redp(state)*rho*T;
	double	g = ses_re_adb_gam(state);
	return 0.5*(sqrt(g*p/rho)+cref)*ses_re_riv(state);
}		/*end re_riv*/

/*ARGSUSED*/
LOCAL double rs_log_pressure(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	SESAME_EOS	*seos = Ses_front_seos(front);
	double	T = ses_rs_temp_from_var(ses_rs_Tvar(state),seos);
	double	rho = ses_rs_rho_from_grid(coords[0],seos);
	return plog(ses_rs_redp(state)*T*rho + ses_rs_coldp(state));
}		/*end rs_log_pressure*/

/*ARGSUSED*/
LOCAL double rs_log_temperature(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	SESAME_EOS	*seos = Ses_front_seos(front);
	return plog(ses_rs_temp_from_var(ses_rs_Tvar(state),seos));
}		/*end rs_log_temperature*/

/*ARGSUSED*/
LOCAL double rs_log_energy(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	SESAME_EOS	*seos = Ses_front_seos(front);
	double	T = ses_rs_temp_from_var(ses_rs_Tvar(state),seos);

	return plog(ses_rs_rede(state)*T+ses_rs_colde(state));
}		/*end rs_log_energy*/

/*ARGSUSED*/
LOCAL double rs_riv(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	SESAME_EOS	*seos = Ses_front_seos(front);
	double		cref = Reference_sound_speed(seos);
	double	rho = ses_rs_rho_from_grid(coords[0],seos);
	double	T = ses_rs_temp_from_var(ses_rs_Tvar(state),seos);
	double	p = ses_rs_coldp(state) + ses_rs_redp(state)*rho*T;
	double	g = ses_rs_adb_gam(state);
	return 0.5*(sqrt(g*p/rho)+cref)*ses_rs_riv(state);
}		/*end rs_riv*/

/*ARGSUSED*/
LOCAL double ps_log_energy(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	SESAME_EOS	*seos = Ses_front_seos(front);
	double	T = ses_ps_temp_from_var(ses_ps_Tvar(state),seos);

	return plog(ses_ps_rede(state)*T+ses_ps_colde(state));
}		/*end ps_log_energy*/

/*ARGSUSED*/
LOCAL double ps_log_temperature(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	SESAME_EOS	*seos = Ses_front_seos(front);
	return plog(ses_ps_temp_from_var(ses_ps_Tvar(state),seos));
}		/*end ps_log_temperature*/

/*ARGSUSED*/
LOCAL double ps_log_density(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	SESAME_EOS	*seos = Ses_front_seos(front);
	return plog(ses_ps_rho_from_var(ses_ps_rho_var(state),seos));
}		/*end ps_log_density*/

/*ARGSUSED*/
LOCAL double ps_riv(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	SESAME_EOS	*seos = Ses_front_seos(front);
	double		cref = Reference_sound_speed(seos);
	double	rho = ses_ps_rho_from_var(ses_ps_rho_var(state),seos);
	double	p = ses_ps_press_from_grid(coords[0],seos);
	double	g = ses_ps_adb_gam(state);
	return 0.5*(sqrt(g*p/rho)+cref)*ses_ps_riv(state);
}		/*end ps_riv*/

/*ARGSUSED*/
LOCAL double vp_entropy(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	return ses_vp_S(state);
}		/*end vp_entropy*/

/*ARGSUSED*/
LOCAL double vp_log_temperature(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	SESAME_EOS	*seos = Ses_front_seos(front);
	return plog(ses_vp_temp_from_var(ses_vp_Tvar(state),seos));
}		/*end vp_log_temperature*/

/*ARGSUSED*/
LOCAL double vp_log_energy(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	SESAME_EOS	*seos = Ses_front_seos(front);
	double	T = ses_vp_temp_from_var(ses_vp_Tvar(state),seos);

	return plog(ses_vp_rede(state)*T+ses_vp_colde(state));
}		/*end vp_log_energy*/

/*ARGSUSED*/
LOCAL double vp_riv(
	double		*coords,
	Front		*front,
	POINTER		wave,
	COMPONENT	comp,
	Locstate	state)
{
	SESAME_EOS	*seos = Ses_front_seos(front);
	double		cref = Reference_sound_speed(seos);
	double	rho = 1.0/ses_vp_vol_from_grid(coords[0],seos);
	double	p = ses_vp_press_from_grid(coords[1],seos);
	double	g = ses_vp_adb_gam(state);
	return 0.5*(sqrt(g*p/rho)+cref)*ses_vp_riv(state);
}		/*end ps_riv*/


#endif /* defined(SESAME_CODE) && defined(TWOD) */
