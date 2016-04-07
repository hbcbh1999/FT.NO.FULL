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
*				gprstate.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains printing routines for gas dynamics.
*
*	g_print_state() is set to a function pointer in the Front structure.
*
*	The following routines are (also) available throughout the gas code:
*
*	g_print_state(), verbose_print_state().
*/

#include <gdecs/gdecs.h>

	/* LOCAL Function Declarations */
LOCAL	void	g_fprint_Estate(FILE*,Locstate);
LOCAL	void	g_fprint_Fstate(FILE*,Locstate);
LOCAL	void	g_fprint_Tstate(FILE*,Locstate);
LOCAL	void	g_fprint_raw_Estate(FILE*,Locstate,int);
LOCAL	void	g_fprint_raw_Fstate(FILE*,Locstate,int);
LOCAL	void	g_fprint_raw_Tstate(FILE*,Locstate,int);
LOCAL	void	g_fprint_raw_state(FILE*,Locstate,int);
LOCAL	void	g_verbose_fprint_state_wrapper(FILE*,Locstate);
LOCAL	void	set_prt_node_sts_params(RP_DATA*,int,int*,int*,
					const char**,int*,int*);
LOCAL	void 	fft_output2d(const char*,const char*,int,RECT_GRID*,COMPLEX**);

EXPORT	void g_fprint_hsbdry_type(
	FILE		*file,
	const char	*mesg1,
	int		hsbdry_type,
	const char	*mesg2,
	INTERFACE	*intfc)
{
	const char *suffix;

	switch (intfc->dim)
	{
	case 3:
	    suffix = "HSBDRY";
	    break;
	case 2:
	    suffix = "NODE";
	    break;
	case 1:
	    screen("ERROR in g_fprint_hsbdry_type(), "
		   "invalid dimension %d\n",intfc->dim);
	    clean_up(ERROR);
	}
	if (hsbdry_type < FIRST_PHYSICS_HSBDRY_TYPE)
	{
	    f_fprint_hsbdry_type(file,mesg1,hsbdry_type,mesg2,intfc);
	    return;
	}
	if (mesg1 != NULL)
	    (void) fprintf(file,"%s",mesg1);
	switch(hsbdry_type) 
	{
	case B_REFLECT_HSBDRY:
	    (void) fprintf(file,"B_REFLECT_%s",suffix);
	    break;
	case MACH_HSBDRY:
	    (void) fprintf(file,"MACH_%s",suffix);
	    break;
	case ATTACHED_B_HSBDRY:
	    (void) fprintf(file,"ATTACHED_B_%s",suffix);
	    break;
	case CROSS_HSBDRY:
	    (void) fprintf(file,"CROSS_%s",suffix);
	    break;
	case OVERTAKE_HSBDRY:
	    (void) fprintf(file,"OVERTAKE_%s",suffix);
	    break;
	case DIFFRACTION_HSBDRY:
	    (void) fprintf(file,"DIFFRACTION_%s",suffix);
	    break;
	case TRANSMISSION_HSBDRY:
	    (void) fprintf(file,"TRANSMISSION_%s",suffix);
	    break;
	case CC_HSBDRY:
	    (void) fprintf(file,"CC_%s",suffix);
	    break;
	case WAVE_END_HSBDRY:
	    (void) fprintf(file,"WAVE_END_%s",suffix);
	    break;
	case TOT_INT_REFL_HSBDRY:
	    (void) fprintf(file,"TOT_INT_REFL_%s",suffix);
	    break;
	/*#bjet2 */
	case NEUMANN_CURVE:
	    (void) fprintf(file,"NEUMANN_CURVE");
	    break;
	default:
	    (void) fprintf(file,"%d -- ** UNKNOWN %s TYPE **",
			   hsbdry_type,suffix);
	    break;
	}
	if (mesg2 != NULL)
	    (void) fprintf(file,"%s",mesg2);
}		/*end g_fprint_hsbdry_type*/

EXPORT	int g_read_hsbdry_type_from_string(
	const char	*type,
	INTERFACE	*intfc)
{
	int		hsbdry_type = UNKNOWN_HSBDRY_TYPE;

	switch(type[0]) 
	{
	case 'B':
	    hsbdry_type = B_REFLECT_HSBDRY;
	    break;
	case 'M':
	    hsbdry_type = MACH_HSBDRY;
	    break;
	case 'A':
	    hsbdry_type = ATTACHED_B_HSBDRY;
	    break;
#if defined(FULL_PHYSICS)
	case 'C':
	    if (type[1] == 'C')
	    	hsbdry_type = CC_HSBDRY;
	    else if (type[1] == 'R')
	    	hsbdry_type = CROSS_HSBDRY;
	    break;
	case 'D':
	    if (type[2] == 'F')
	    	hsbdry_type = DIFFRACTION_HSBDRY;
	    break;
	case 'O':
	    hsbdry_type = OVERTAKE_HSBDRY;
	    break;
	case 'T':
	    if (type[1] == 'R')
	    	hsbdry_type = TRANSMISSION_HSBDRY;
	    else if (type[1] == 'O')
	    	hsbdry_type = TOT_INT_REFL_HSBDRY;
	    break;
	case 'W':
	    hsbdry_type = WAVE_END_HSBDRY;
	    break;
#endif /* defined(FULL_PHYSICS) */
	default:
	    break;
	}
	if (hsbdry_type != UNKNOWN_HSBDRY_TYPE)
	    return hsbdry_type;
	return f_read_hsbdry_type_from_string(type,intfc);
}		/*end g_read_hsbdry_type_from_string*/


#if defined(TWOD)
#if defined(FULL_PHYSICS) && defined(DEBUG_NODE_PROPAGATE)
EXPORT	void print_diffraction_status(
	const char	*message,
	int		status)
{
	switch (status)
	{
	case ERROR_DIFFRACTION:
	    (void) printf("%sERROR_DIFFRACTION\n",message);
	    break;
	case REGULAR_DIFFRACTION:
	    (void) printf("%sREGULAR_DIFFRACTION\n",message);
	    break;
	case ANOMALOUS_REFLECTION:
	    (void) printf("%sANOMALOUS_REFLECTION\n",message);
	    break;
	case REGULAR_TO_MACH_DIFFRACTION:
	    (void) printf("%sREGULAR_TO_MACH_DIFFRACTION\n",message);
	    break;
	case PRECURSOR_WITH_REFLECTED_RAREFACTION:
	    (void) printf("%sPRECURSOR_WITH_REFLECTED_RAREFACTION\n",message);
	    break;
	case PRECURSOR_WITH_REFLECTED_SHOCK:
	    (void) printf("%sPRECURSOR_WITH_REFLECTED_SHOCK\n",message);
	    break;
	default:
	    (void) printf("%s%d -- ** UNKNOWN DIFFRACTION STATUS **\n",
			  message,status);
	    break;
	}
}		/*end print_diffraction_status*/
#endif /* defined(FULL_PHYSICS) && defined(DEBUG_NODE_PROPAGATE) */


EXPORT	void fprint_curve_status(
	FILE		*file,
	const char	*message,
	int		status)
{
	switch(status) 
	{
	case PASSIVE:
	    (void) fprintf(file,"%sPASSIVE\n",message);
	    break;
	case FIXED:
	    (void) fprintf(file,"%sFIXED\n",message);
	    break;
	case INCIDENT:
	    (void) fprintf(file,"%sINCIDENT\n",message);
	    break;
	case TRANSMITTED:
	    (void) fprintf(file,"%sTRANSMITTED\n",message);
	    break;
	case REFLECTED:
	    (void) fprintf(file,"%sREFLECTED\n",message);
	    break;
	case MACH_STEM:
	    (void) fprintf(file,"%sMACH_STEM\n",message);
	    break;
	case SLIP:
	    (void) fprintf(file,"%sSLIP\n",message);
	    break;
	case CONTACT_TARGET:
	    (void) fprintf(file,"%sCONTACT_TARGET\n",message);
	    break;
	case OVERTOOK:
	    (void) fprintf(file,"%sOVERTOOK\n",message);
	    break;
	case VIRTUAL:
	    (void) fprintf(file,"%sVIRTUAL\n",message);
	    break;
	case UNKNOWN_CURVE_STATUS:
	    (void) fprintf(file,"%sUNKNOWN_CURVE_STATUS\n",message);
	    break;
	default:
	    (void) fprintf(file,"%s%d -- ** UNKNOWN CURVE STATUS **\n",
			   message,status);
	    break;
	}
}		/*end fprint_curve_status*/

EXPORT	int read_curve_status_from_string(
	const char *type)
{
	int status = UNKNOWN_CURVE_STATUS;

	switch(type[0]) 
	{
	case 'P':
	case 'p':
	    status = PASSIVE;
	    break;
	case 'F':
	case 'f':
	    status = FIXED;
	    break;
	case 'I':
	case 'i':
	    status = INCIDENT;
	    break;
	case 'T':
	case 't':
	    status = TRANSMITTED;
	    break;
	case 'R':
	case 'r':
	    status = REFLECTED;
	    break;
	case 'M':
	case 'm':
	    status = MACH_STEM;
	    break;
	case 'S':
	case 's':
	    status = SLIP;
	    break;
	case 'O':
	case 'o':
	    status = OVERTOOK;
	    break;
	case 'C':
	case 'c':
	    status = CONTACT_TARGET;
	    break;
	case 'V':
	case 'v':
	    status = VIRTUAL;
	    break;
	default:
	    screen("ERROR in read_curve_status_from_string(), "
		   "unknown status\n");
	    clean_up(ERROR);
	    break;
	}
	return status;
}		/*end read_curve_status_from_string*/
#endif /* defined(TWOD) */

EXPORT	const char *g_wave_type_as_string(
	int w_type)
{
	switch(w_type) 
	{
	case BACKWARD_SHOCK_WAVE:
	    return "BACKWARD_SHOCK_WAVE";
	case BACKWARD_SOUND_WAVE_LE:
	    return "BACKWARD_SOUND_WAVE_LE";
	case BACKWARD_SOUND_WAVE_TE:
	    return "BACKWARD_SOUND_WAVE_TE";
	case CONTACT:
	    return "CONTACT";
	case THIN_FLAME:
	    return "THIN_FLAME";
	case FORWARD_SHOCK_WAVE:
	    return "FORWARD_SHOCK_WAVE";
	case FORWARD_SOUND_WAVE_LE:
	    return "FORWARD_SOUND_WAVE_LE";
	case FORWARD_SOUND_WAVE_TE:
	    return "FORWARD_SOUND_WAVE_TE";
	case TIME_DIRICHLET_BOUNDARY:
	    return "TIME_DIRICHLET_BOUNDARY";
	case MOVABLE_BODY_BOUNDARY:
	    return "MOVABLE_BODY_BOUNDARY";
	case VELOCITY_SPECIFIED:
	    return "VELOCITY_SPECIFIED";
	case RIEMANN_PROBLEM_WAVE:
	    return "RIEMANN_PROBLEM_WAVE";
	case UNKNOWN_WAVE_TYPE:
	    return "UNKNOWN_WAVE_TYPE";
	default:
	    return f_wave_type_as_string(w_type);
	}
}		/*end g_wave_type_as_string*/

EXPORT	int g_read_wave_type_from_string(
	const char	*type)
{
	int i;
	int w_type = UNKNOWN_WAVE_TYPE;
	static struct { const char *name; int type; } wave_type_map[] = {
	    {"BACKWARD_SOUND_WAVE_LE",  BACKWARD_SOUND_WAVE_LE},
	    {"BL",                      BACKWARD_SOUND_WAVE_LE},
	    {"BACKWARD_SOUND_WAVE_TE",  BACKWARD_SOUND_WAVE_TE},
	    {"BT",                      BACKWARD_SOUND_WAVE_TE},
	    {"BACKWARD_SHOCK_WAVE",     BACKWARD_SHOCK_WAVE},
	    {"B",                       BACKWARD_SHOCK_WAVE},
	    {"CONTACT",                 CONTACT},
	    {"C",                       CONTACT},
	    {"THIN_FLAME",              THIN_FLAME},
	    {"TF",              	THIN_FLAME},
	    {"FORWARD_SOUND_WAVE_LE",   FORWARD_SOUND_WAVE_LE},
	    {"FL",                      FORWARD_SOUND_WAVE_LE},
	    {"FORWARD_SOUND_WAVE_TE",   FORWARD_SOUND_WAVE_TE},
	    {"FT",                      FORWARD_SOUND_WAVE_TE},
	    {"FORWARD_SHOCK_WAVE",      FORWARD_SHOCK_WAVE},
	    {"F",                       FORWARD_SHOCK_WAVE},
	    {"TIME_DIRICHLET_BOUNDARY", TIME_DIRICHLET_BOUNDARY},
	    {"MOVABLE_BODY_BOUNDARY",MOVABLE_BODY_BOUNDARY},
	    {"MN",			MOVABLE_BODY_BOUNDARY},
	    {"VELOCITY_SPECIFIED",      VELOCITY_SPECIFIED},
	    {"RIEMANN_PROBLEM_WAVE",    RIEMANN_PROBLEM_WAVE},
	    {NULL,UNKNOWN_WAVE_TYPE}
	};
	for (i = 0; wave_type_map[i].name != NULL; ++i)
	{
	    if (strcasecmp(type,wave_type_map[i].name) == 0)
	    {
	        w_type = wave_type_map[i].type;
		break;
	    }
	}

	if (w_type != UNKNOWN_WAVE_TYPE)
	    return w_type;
	return f_read_wave_type_from_string(type);
}		/*end g_read_wave_type_from_string*/


EXPORT	const char *rsoln_wave_name(
	RIEMANN_SOLVER_WAVE_TYPE wave)
{
	static char s[128];
	switch(wave)
	{
	case SHOCK:
	    return "SHOCK";
	case RAREFACTION:
	    return "RAREFACTION";
#if defined(COMBUSTION_CODE)
	case STRONG_DET:
	    return "STRONG_DET";
	case CJ_DET:
	    return "CJ_DET";
#endif /* defined(COMBUSTION_CODE) */
	case UNKNOWN_WAVE_TYPE:
	    return "UNKNOWN_WAVE_TYPE";
	default:
	    (void) sprintf(s,"%d UNDEFINED WAVE TYPE",wave);
	    return s;
	}
}		/*end rsoln_wave_name*/

EXPORT	void print_rsoln_wave(
	const char	         *message,
	RIEMANN_SOLVER_WAVE_TYPE wave,
	const char	         *end)
{
	(void) printf("%s%s%s",message,rsoln_wave_name(wave),end);
}		/*end print_rsoln_wave*/

EXPORT	void print_state_type(
	const char	*message,
	int		state_type)
{
	fprint_state_type(stdout,message,state_type);
}		/*end print_state_type*/

EXPORT	const char *state_type_name(
	int state_type)
{
        static char s[180];
	switch (state_type)
	{
	case GAS_STATE:
	    return "GAS_STATE";
	case EGAS_STATE:
	    return "EGAS_STATE";
	case TGAS_STATE:
	    return "TGAS_STATE";
	case FGAS_STATE:
	    return "FGAS_STATE";
	case VGAS_STATE:
	    return "VGAS_STATE";
#if defined(COMBUSTION_CODE)
	case ZGAS_STATE:
	    return "ZGAS_STATE";
	case CGAS_STATE:
	    return "CGAS_STATE";
#endif /* defined(COMBUSTION_CODE) */
	default:
	    (void) sprintf(s,"%d***UNKNOWN_STATE***",state_type);
	    return s;
	}
}		/*end state_type_name*/

EXPORT	void fprint_state_type(
	FILE		*file,
	const char	*message,
	int		state_type)
{
        (void) fprintf(file,"%s%s\n",message,state_type_name(state_type));
}		/*end fprint_state_type*/

EXPORT	int g_read_state_type_from_string(
	const char *stype)
{
	if (strcmp(stype,"GAS_STATE") == 0)
	    return GAS_STATE;
	if (strcmp(stype,"EGAS_STATE") == 0)
	    return EGAS_STATE;
	if (strcmp(stype,"TGAS_STATE") == 0)
	    return TGAS_STATE;
	if (strcmp(stype,"FGAS_STATE") == 0)
	    return FGAS_STATE;
	if (strcmp(stype,"VGAS_STATE") == 0)
	    return VGAS_STATE;
#if defined(COMBUSTION_CODE)
	if (strcmp(stype,"ZGAS_STATE") == 0)
	    return ZGAS_STATE;
	if (strcmp(stype,"CGAS_STATE") == 0)
	    return CGAS_STATE;
#endif /* defined(COMBUSTION_CODE) */
	screen("ERROR in read_state_type_from_string(), "
	       "unknown state type = %s\n",stype);
	clean_up(ERROR);
	return UNKNOWN_STATE;
}		/*end read_state_type_from_string*/

/*ARGSUSED*/
EXPORT	void g_fprint_intfc_state(
	FILE		*file,
	Locstate	state,
	INTERFACE	*intfc)
{
	g_fprint_state(file,state);
}		/*end g_fprint_intfc_state*/

/*ARGSUSED*/
EXPORT	void g_verbose_fprint_intfc_state(
	FILE		*file,
	Locstate	state,
	INTERFACE	*intfc)
{
	verbose_fprint_state(file,"",state);
}		/*end g_verbose_fprint_intfc_state*/

/*ARGSUSED*/
EXPORT	void g_fprint_state_data(
	FILE		*file,
	Locstate	state,
	INTERFACE	*intfc)
{
	fprint_gas_data(file,state);
}		/*end g_fprint_state_data*/

EXPORT	void fprint_gas_data(
	FILE		*file,
	Locstate	state)
{
	(void) fprintf(file,"State information for the ");
	if (state == NULL)
	{
	    (void) fprintf(file,"NULL state 0x%p\n\n",(POINTER)state);
	    return;
	}
	if (is_obstacle_state(state)) 
	{
	    (void) fprintf(file,"OBSTACLE state 0x%p\n\n",(POINTER)state);
	    return;
	}
	(void) fprintf(file,"state 0x%p\n",(POINTER)state);
	(void) fprintf(file,"\tState Data ");
	if (is_binary_output() == YES)
	{
	    uint64_t prms;
	    double    *x;
	    int      stype = state_type(state);
	    int      failed = material_failure(state);

	    (void) fprintf(file,"\f%c",(char)Params(state)->sizest);
	    x = &Dens(state);
	    (void) fwrite((const void *)x,FLOAT,2+SMAXD,file);
	    prms = gas_param_number(Params(state));
	    (void) fwrite((const void *)&prms,sizeof(uint64_t),1,file);
	    (void) fwrite((const void *)&stype,sizeof(int),1,file);
	    (void) fwrite((const void *)&failed,sizeof(int),1,file);
#if defined(COMBUSTION_CODE)
	    switch (Composition_type(state))
	    {
	    case ZND:
	    case PTFLAME:
	        (void) fwrite((const void *)pdens(state),FLOAT,1,file);
		break;
	    case TWO_CONSTITUENT_REACTIVE:
	        (void) fwrite((const void *)pdens(state),FLOAT,2,file);
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
                        (void) fwrite((const void *)&(pdens(state)[i]),FLOAT,1,file);
                }
            }
	    (void) fprintf(file,"\n");
	    return;
	}
	(void) fprintf(file,"\n");
	switch (state_type(state))
	{
	case GAS_STATE:
	    g_fprint_state(file,state);
	    break;

	case EGAS_STATE:
	    g_fprint_Estate(file,state);
	    break;

	case TGAS_STATE:
	    g_fprint_Tstate(file,state);
	    break;

	case FGAS_STATE:
	    g_fprint_Fstate(file,state);
	    break;

	case VGAS_STATE:
	    verbose_fprint_state(file,"",state);
	    break;

	default:
	    screen("ERROR in fprint_gas_data(), "
	           "unknown state type %d\n",state_type(state));
	    clean_up(ERROR);
	}
}		/*end fprint_gas_data*/

EXPORT	void	fprint_raw_gas_data(
	FILE		*file,
	Locstate	state,
	int		dim)
{
	(void) fprintf(file,"state 0x%p\n",(POINTER)state);
	if (state == NULL)
	{
	    (void) printf("NULL state\n");
	    return;
	}
	(void) fprintf(file,"\tState Data ");
	(void) fprintf(file,"\n");
	switch (state_type(state))
	{
	case GAS_STATE:
	    g_fprint_raw_state(file,state,dim);
	    break;

	case EGAS_STATE:
	    g_fprint_raw_Estate(file,state,dim);
	    break;

	case TGAS_STATE:
	    g_fprint_raw_Tstate(file,state,dim);
	    break;

	case FGAS_STATE:
	    g_fprint_raw_Fstate(file,state,dim);
	    break;

	case OBSTACLE_STATE:
	    (void) fprintf(file,"Obstacle state type,  "
				"printing as GAS_STATE\n");
	    g_fprint_raw_state(file,state,dim);
	    break;

	case VGAS_STATE:
	    g_fprint_raw_Tstate(file,state,dim);
	    (void) printf("Specific internal energy = %"FFMT"\n",Int_en(state));
	    (void) printf("Entropy = %"FFMT"\n",Entropy(state));
	    (void) printf("Sound_speed = %"FFMT"\n",sound_speed(state));
#if defined(VERBOSE_GAS_PLUS)
	    (void) printf("Enthalpy = %"FFMT"\n",Enthalpy(state));
	    (void) printf("Temperature = %"FFMT"\n",Temp(state));
#endif /* defined(VERBOSE_GAS_PLUS) */
	    break;

	case UNKNOWN_STATE:
	    (void) fprintf(file,"Unknown state type,  printing as GAS_STATE\n");
	    g_fprint_raw_state(file,state,dim);
	    break;

	default:
	    screen("ERROR in fprint_raw_gas_data(), "
	           "unknown state type %d\n",state_type(state));
	    clean_up(ERROR);
	}
}		/*end fprint_raw_gas_data*/


/*
*			g_print_state():
*
*       Print the density, energy density, and x,y,z-momentum density
*	of a given state.
*/

EXPORT void g_print_state(
	Locstate	state)
{
	g_fprint_state(stdout,state);
}		/*end g_print_state*/

EXPORT void g_fprint_state(
	FILE		*file,
	Locstate	state)
{
	int dim;
	if (is_obstacle_state(state))
	{
	    (void) fprintf(file,"state %p (OBSTACLE STATE)\n\n",state);
	    return;
	}
	if (current_interface())
	    dim = current_interface()->dim;
	else
	    dim = Params(state)->dim;
	g_fprint_raw_state(file,state,dim);

#if !defined(COMBUSTION_CODE)
	(void) fprintf(file,"\n");
#else /* !defined(COMBUSTION_CODE) */
	if (Composition_type(state) == PURE_NON_REACTIVE)
	{
	    (void) fprintf(file,"\n");
	    return;
	}
	(void) fprintf(file,"burned = %s   q = %"FFMT"   t_crit = %"FFMT"\n",
			    Burned(state)? "BURNED" : "UNBURNED",
			    Params(state)->q,
			    Params(state)->critical_temperature);

	if ((Composition_type(state) == PTFLAME) ||
	    (Composition_type(state) == THINFLAME))
	{
	    (void) fprintf(file,"\n");
	    return;
	}

	(void) fprintf(file,"product density = %"FFMT"\n",Prod(state));

	if (Composition_type(state) == ZND)
	{
	    (void) fprintf(file,"\n");
	    return;
	}

	(void) fprintf(file," rho1 = %"FFMT"\n\n",Dens1(state));
#endif /* !defined(COMBUSTION_CODE) */
	(void) fprintf(file,"%-24s = %-24s %-24s = %-"FFMT"\n\n","gamma_set",
                      Local_gamma_set(state)?"YES" : "NO","local_gamma",
                      Local_gamma(state));
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(state) != NULL &&
               Params(state)->n_comps != 1)
            {
                int i;
                for(i = 0; i < Params(state)->n_comps; i++)
                    (void) fprintf(file,"partial density[%2d] = %"FFMT"\n",
                       i,pdens(state)[i]);
                (void) fprintf(file,"\n");
            }
        }
}		/*end g_fprint_state*/

LOCAL void g_fprint_raw_state(
	FILE		*file,
	Locstate	state,
	int		dim)
{
	int		i;
	static	char	mname[3][3] = { "mx", "my", "mz"};

	(void) fprintf(file,"state %p\n",state);
	/*(void) fprintf(file,"\t%-7s = %"FFMT" %-7s = %"FFMT"\n\t","density", */
	/*	            Dens(state),"energy",Energy(state)); */
	(void) fprintf(file,"\t%-7s = %22.16g %-7s = %22.16g\n\t","density",
		            Dens(state),"energy",Energy(state));
	for (i = 0; i < dim; i++)
	    (void) fprintf(file,"%-7s = %"FFMT" ",mname[i],Mom(state)[i]);
	(void) fprintf(file,"\n");
	(void) fprintf(file,"\ttype = %u, failed = %u\n",
	                    state_type(state),material_failure(state));
	if (debugging("prt_params"))
	    fprint_Gas_param(file,Params(state));
	else
	    (void) fprintf(file,"Gas_param = %llu\n",
		                gas_param_number(Params(state)));
	if (debugging("local_gamma"))
	{
	    (void) fprintf(file,"Local gamma set: %s\n",
			   Local_gamma_set(state) ? "YES" : "NO");
	    if (Local_gamma_set(state))
		(void) fprintf(file,"Local gamma = %"FFMT"\n",
			       Local_gamma(state));
	}
}		/*end g_fprint_raw_state*/

/*
*			g_fprint_Estate():
*
*       Prints the density, specific internal energy, x,y,z-velocity, and 
*	of a given state, as specified by an EGas variable.
*/

LOCAL void g_fprint_Estate(
	FILE		*file,
	Locstate	state)
{
	if (is_obstacle_state(state))
	{
	    (void) fprintf(file,"(OBSTACLE STATE)\n\n");
	    return;
	}
	g_fprint_raw_Estate(file,state,Params(state)->dim);

#if !defined(COMBUSTION_CODE)
	(void) fprintf(file,"\n");
#else /* !defined(COMBUSTION_CODE) */
	if (Composition_type(state) == PURE_NON_REACTIVE)
	{
	    (void) fprintf(file,"\n");
	    return;
	}

	(void) fprintf(file,"burned = %s   q = %"FFMT"   t_crit = %"FFMT"\n",
	                    Burned(state) ? "BURNED" : "UNBURNED",
	                    Params(state)->q,
			    Params(state)->critical_temperature);

	if (Composition_type(state) == PTFLAME)
	{
	    (void) fprintf(file,"\n");
	    return;
	}

	(void) fprintf(file,"reaction progress = %"FFMT"\n",React(state));

	if (Composition_type(state) == ZND)
	{
	    (void) fprintf(file,"\n");
	    return;
	}

	(void) fprintf(file," rho1 = %"FFMT"\n\n",Dens1(state));
#endif /* !defined(COMBUSTION_CODE) */
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(state) != NULL &&
               Params(state)->n_comps != 1)
            {
                int i;
                for(i = 0; i < Params(state)->n_comps; i++)
                    (void) fprintf(file,"partial density[%2d] = %"FFMT"\n",
                       i,pdens(state)[i]);
                (void) fprintf(file,"\n");
            }
        }
}		/*end g_fprint_Estate*/

LOCAL void g_fprint_raw_Estate(
	FILE		*file,
	Locstate	state,
	int		dim)
{
	int		i;
	static	char	vname[3][3] = { "vx", "vy", "vz"};

	(void) fprintf(file,"\tdensity = %"FFMT" energy = %"FFMT" ",
		            Dens(state),Energy(state));
	for (i = 0; i < dim; i++)
	    (void) fprintf(file,"%-8s = %"FFMT" ",vname[i],Vel(state)[i]);
	if (debugging("prt_params"))
	    fprint_Gas_param(file,Params(state));
	else
	    (void) fprintf(file,"Gas_param = %llu\n",
		                gas_param_number(Params(state)));

}		/*end g_fprint_raw_Estate*/

/*
*			g_fprint_Tstate():
*
*       Prints the density, pressure, x,y,z-velocity of a
*	given state, as specified by a thermodynamic (TGas) variable.
*
*			uses TGas
*/

LOCAL void g_fprint_Tstate(
	FILE		*file,
	Locstate	state)
{
	if (is_obstacle_state(state))
	{
		(void) fprintf(file,"(OBSTACLE STATE)\n\n");
		return;
	}
	g_fprint_raw_Tstate(file,state,Params(state)->dim);

#if !defined(COMBUSTION_CODE)
	(void) fprintf(file,"\n");
#else /* !defined(COMBUSTION_CODE) */
	if (Composition_type(state) == PURE_NON_REACTIVE)
	{
		(void) fprintf(file,"\n");
		return;
	}

	(void) fprintf(file,"burned = %s   q = %"FFMT"   t_crit = %"FFMT"\n",
		Burned(state) ? "BURNED" : "UNBURNED",
		Params(state)->q,Params(state)->critical_temperature);

	if ((Composition_type(state) == PTFLAME) ||
	    (Composition_type(state) == THINFLAME))
	{
		(void) fprintf(file,"\n");
		return;
	}

	(void) fprintf(file,"reaction progress = %"FFMT"\n",React(state));

	if (Composition_type(state) == ZND)
	{
		(void) fprintf(file,"\n");
		return;
	}

	(void) fprintf(file," rho1 = %"FFMT"\n\n",Dens1(state));
#endif /* !defined(COMBUSTION_CODE) */
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(state) != NULL &&
               Params(state)->n_comps != 1)
            {
                int i;
                for(i = 0; i < Params(state)->n_comps; i++)
                    (void) fprintf(file,"partial density[%2d] = %"FFMT"\n",
                       i,pdens(state)[i]);
                (void) fprintf(file,"\n");
            }
        }
}		/*end g_fprint_Tstate*/

LOCAL void g_fprint_raw_Tstate(
	FILE		*file,
	Locstate	state,
	int		dim)
{
	int		i;
	static	char	vname[3][3] = { "vx", "vy", "vz"};

	(void) fprintf(file,"\tdensity = %"FFMT" pressure = %"FFMT" ",
		            Dens(state),Press(state));
	for (i = 0; i < dim; i++)
	    (void) fprintf(file,"%-8s = %"FFMT" ",vname[i],Vel(state)[i]);
	if (debugging("prt_params"))
	    fprint_Gas_param(file,Params(state));
	else
	    (void) fprintf(file,"Gas_param = %llu\n",
		                gas_param_number(Params(state)));


}		/*end g_fprint_raw_Tstate*/

/*
*			g_fprint_Fstate():
*
*       Prints the density, temperature, x,y,z-velocity of a
*	given state, as specified by a thermodynamic (FGas) variable.
*
*			uses FGas
*/

LOCAL void g_fprint_Fstate(
	FILE		*file,
	Locstate	state)
{
	if (is_obstacle_state(state))
	{
	    (void) fprintf(file,"(OBSTACLE STATE)\n\n");
	    return;
	}
	g_fprint_raw_Fstate(file,state,Params(state)->dim);

#if !defined(COMBUSTION_CODE)
	(void) fprintf(file,"\n");
#else /* !defined(COMBUSTION_CODE) */
	if (Composition_type(state) == PURE_NON_REACTIVE)
	{
	    (void) fprintf(file,"\n");
	    return;
	}

	(void) fprintf(file,"burned = %s   q = %"FFMT"   t_crit = %"FFMT"\n",
		Burned(state) ? "BURNED" : "UNBURNED",
		Params(state)->q,Params(state)->critical_temperature);

	if (Composition_type(state) == PTFLAME)
	{
	    (void) fprintf(file,"\n");
	    return;
	}

	(void) fprintf(file,"reaction progress = %"FFMT"\n",React(state));

	if (Composition_type(state) == ZND)
	{
	    (void) fprintf(file,"\n");
	    return;
	}

	(void) fprintf(file," rho1 = %"FFMT"\n\n",Dens1(state));
#endif /* !defined(COMBUSTION_CODE) */
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(state) != NULL &&
               Params(state)->n_comps != 1)
            {
                int i;
                for(i = 0; i < Params(state)->n_comps; i++)
                    (void) fprintf(file,"partial density[%2d] = %"FFMT"\n",
                       i,pdens(state)[i]);
                (void) fprintf(file,"\n");
            }
        }
}		/*end g_fprint_Fstate*/

LOCAL void g_fprint_raw_Fstate(
	FILE		*file,
	Locstate	state,
	int		dim)
{
	int		i;
	static	char	vname[3][3] = { "vx", "vy", "vz"};

	(void) fprintf(file,"\tdensity = %"FFMT" temperature = %"FFMT" ",
		            Dens(state),Temperature(state));
	for (i = 0; i < dim; i++)
	    (void) fprintf(file,"%-8s = %"FFMT" ",vname[i],Vel(state)[i]);
	if (debugging("prt_params"))
	    fprint_Gas_param(file,Params(state));
	else
	    (void) fprintf(file,"Gas_param = %llu\n",
		                gas_param_number(Params(state)));


}		/*end g_fprint_raw_Fstate*/


/*
*			g_verbose_fprint_state():
*
*       Prints the density, energy density, x y z momentum density, reaction 
*	progress, product density, pressure, entropy, temperature, sound speed,
*	velocity angle, x,y,z-velocity of a given gas state.
*	If flag equals GAS_STATE the state is assumed to be in Gas form, while
*	if flag equals TGAS_STATE the state is assumed to be in TGas form.
*/

EXPORT void g_verbose_fprint_state(
	FILE	   *file,
	const char *name,
	Locstate   state)
{
	boolean		bin_out;
	int		i, dim;
	double		p, c, S, v[SMAXD], speed;
	static	char	vname[3][3] = { "vx", "vy", "vz"};
	static	char	mname[3][3] = { "mx", "my", "mz"};

	if (name == NULL)
	    name = "";
	(void) fprintf(file,"\n%s:\n",name);
	(void) fprintf(file,"address %p\n",state);
	if (state == NULL || is_obstacle_state(state))
	{
	    (void) fprintf(file,"(OBSTACLE STATE)\n\n");
	    return;
	}
	dim = Params(state)->dim;

	p = pressure(state);
	c = sound_speed(state);
	S = entropy(state);

	(void) fprintf(file,"%-24s = %-"FFMT" %-24s = %-"FFMT"\n",
		       "density",Dens(state),
		       "specific internal energy",
		       specific_internal_energy(state));
	(void) fprintf(file,"%-24s = %-"FFMT" %-24s = %-"FFMT"\n","pressure",p,
		       "sound speed",c);
	(void) fprintf(file,"%-24s = %-"FFMT" %-24s = %-"FFMT"\n","temperature",
		       temperature(state),"specific entropy",S);

	speed = 0.0;
	for (i = 0; i < dim; i++)
	{
	    v[i] = vel(i,state);	speed += sqr(v[i]);
	    (void) fprintf(file,"%-24s = %-"FFMT" %-24s = %-"FFMT"\n",
			   mname[i],mom(i,state),vname[i],v[i]);
	}
	speed = sqrt(speed);

	(void) fprintf(file,"%-24s = %-"FFMT"","total energy",energy(state));
	if (c > 0. && Dens(state) > 0.)
	   (void) fprintf(file," %-24s = %-"FFMT"\n","Mach number",speed / c);
	else
	   (void) fprintf(file,"\n");

#if defined(TWOD)
	if (dim == 2)
	    (void) fprintf(file,"%-24s = %-"FFMT"\n","velocity angle",
			   degrees(angle(v[0],v[1])));
#endif /* defined(TWOD) */


	fprint_state_type(file,"State type = ",state_type(state));
	(void) fprintf(file,"Params state = %llu\n",
		       gas_param_number(Params(state)));

	bin_out = is_binary_output();
	set_binary_output(NO);
	if (debugging("prt_params"))
	    fprint_Gas_param(file,Params(state));
	else
	    (void) fprintf(file,"Gas_param = %llu\n",
		                gas_param_number(Params(state)));
	set_binary_output(bin_out);


#if !defined(COMBUSTION_CODE)
	(void) fprintf(file,"\n");
#else /* !defined(COMBUSTION_CODE) */
	if (Composition_type(state) == PURE_NON_REACTIVE)
	{
	    (void) fprintf(file,"\n");
	    return;
	}

	(void) fprintf(file,"%-24s = %-12s   %-24s = %-"FFMT"\n","burned",
		       Burned(state) ? "BURNED" : "UNBURNED",
		       "q",Params(state)->q);
	(void) fprintf(file,"%-24s = %-"FFMT"\n","t_crit",
		       Params(state)->critical_temperature);

	if (Composition_type(state) == PTFLAME)
	{
	    (void) fprintf(file,"\n");
	    return;
	}

	(void) fprintf(file,"%-24s = %-"FFMT"\n","product density",Prod(state));
	(void) fprintf(file,"%-24s = %-"FFMT"\n",
		       "reaction progress",React(state));

	if (Composition_type(state) == ZND)
	{
	    (void) fprintf(file,"\n");
	    return;
	}

	(void) fprintf(file,"%-24s = %-"FFMT"\n","rho1",Dens1(state));
#endif /* !defined(COMBUSTION_CODE) */
	(void) fprintf(file,"%-24s = %-24s %-24s = %-"FFMT"\n\n","gamma_set",
		      Local_gamma_set(state)?"YES" : "NO","local_gamma",
		      Local_gamma(state));
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(state) != NULL &&
               Params(state)->n_comps != 1)
            {
                for(i = 0; i < Params(state)->n_comps; i++)
                    (void) fprintf(file,"partial density[%2d] = %"FFMT"\n",
                       i,pdens(state)[i]);
                (void) fprintf(file,"\n");
                 /* TMP the following print is for the debuging display purpose */
                for(i = 0; i < Params(state)->n_comps; i++)
                    (void) fprintf(file,"[%d]Mass fraction = %-"FFMT"\n",
                             i, pdens(state)[i]/Dens(state));
                (void) fprintf(file,"\n");
            }
        }
}		/*end g_verbose_fprint_state*/

LOCAL	void g_verbose_fprint_state_wrapper(
	FILE		*file,
	Locstate	state)
{
	verbose_fprint_state(file,"",state);
}		/*end g_verbose_fprint_state_wrapper*/

EXPORT	void g_verbose_print_state(
	Locstate	state)
{
	g_verbose_fprint_state_wrapper(stdout,state);
}		/*end g_verbose_print_state*/

#if defined(FULL_PHYSICS)
EXPORT	void print_RP_node_states(
	const char	*message,
	double		*nod_v,
	RP_DATA		*RP,
	int		n_type)
{
	double	   M;
	double	   v[MAX_N_CURVES][MAXD];
	double	   q[MAX_N_CURVES], theta[MAX_N_CURVES];
	double	   dtheta;
	int	   i, j;
	int	   dim = Params(RP->state[0])->dim;
	int	   nangs;
	int	   iang[MAX_N_CURVES];
	int	   ist[MAX_N_CURVES];
	int	   num_sts;
	char	   mesg[80];
	const char *ang_names[7];

	set_prt_node_sts_params(RP,n_type,&num_sts,&nangs,ang_names,ist,iang);

	if (message != NULL)
	    (void) printf("%s\n",message);

	print_angle_direction("RP->ang_dir =",RP->ang_dir,"\n");
	(void) printf("Node velocity = <%"FFMT", %"FFMT">\n",nod_v[0],nod_v[1]);
	for (i = 0; i < num_sts; i++)
	{
	    M = mach_number(RP->state[ist[i]],nod_v);
	    (void) sprintf(mesg,"state%d is %s, M%d = %"FFMT,ist[i],
	    	(M >= 1.0) ? "supersonic" : "subsonic",ist[i],M);
	    verbose_print_state(mesg,RP->state[ist[i]]);
	}


	(void) printf("\nWave Angles:\n");
	for (i = 0; i < nangs; i++)
	{
	    (void) printf("RP->ang[%d],\n\t%s = %"FFMT" (%g degrees)\n",iang[i],
	    	          ang_names[iang[i]],RP->ang[iang[i]],
	    	          degrees(RP->ang[iang[i]]));
	    dtheta = RP->ang[iang[(i+1)%nangs]] - RP->ang[iang[i]];
	    (void) printf("\tRP->ang[%d] - RP->ang[%d] = ",
	    	          iang[(i+1)%nangs],iang[i]);
	    (void) printf("%"FFMT" (%g degrees)\n",dtheta,degrees(dtheta));
	}

	(void) printf("\nSteady state velocities\n");
	for (i = 0; i < num_sts; i++)
	{
	    q[i] = 0;
	    for (j = 0; j < dim; j++)
	    {
	    	v[i][j] = vel(j,RP->state[ist[i]]) - nod_v[j];
	    	q[i] += sqr(v[i][j]);
	    }
	    theta[i] = angle(v[i][0],v[i][1]);
	    q[i] = sqrt(q[i]);
	}
	for (i = 0; i < num_sts; i++)
	{
	    dtheta = theta[(i+1)%num_sts] - theta[i];
	    (void) printf("vel_ang(RP->state[%d]) - vel_ang(RP->state[%d]) = ",
			(ist[i]+1)%num_sts,ist[i]);
	    (void) printf("%"FFMT" (%g degrees)\n",dtheta,degrees(dtheta));
	    (void) printf("RP->state[%d],\n\t",ist[i]);
	    (void) printf("v = <");
	    for (j = 0; j < dim; j++)
	    	(void) printf("%"FFMT"%s",v[i][j],(j==(dim-1))?">, ":",");
	    (void) printf("q = %"FFMT", theta = %"FFMT" (%g degrees)\n",
			  q[i],theta[i],degrees(theta[i]));
	}
	(void) printf("\n");
}		/*end print_RP_node_states*/

LOCAL	void set_prt_node_sts_params(
	RP_DATA		*RP,
	int		n_type,
	int		*num_sts,
	int		*nangs,
	const char	**ang_names,
	int		*ist,
	int		*iang)
{
	double		p[4];
	int		i;

	switch (n_type)
	{
	case DIFFRACTION_NODE:
	case TOT_INT_REFL_NODE:
	    *num_sts = 7;
	    for (i = 0; i < *num_sts; i++)
	        ist[i] = i;
	    ang_names[0] = "Incident shock angle";
	    ang_names[1] = "Reflected rarefaction leading edge angle";
	    ang_names[2] = "Reflected shock angle";
	    ang_names[3] = "Reflected rarefaction trailing edge angle";
	    ang_names[4] = "Downstream contact angle";
	    ang_names[5] = "Transmitted shock angle";
	    ang_names[6] = "Upstream contact angle";
	    iang[0] = 0;
	    if (pressure(RP->state[4]) < pressure(RP->state[1]))
	    {
	    	*nangs = 6;
	    	iang[1] = 1; iang[2] = 3; iang[3] = 4;
	    	iang[4] = 5; iang[5] = 6;
	    }
	    else
	    {
	    	*nangs = 5;
	    	iang[1] = 2; iang[2] = 4; iang[3] = 5; iang[4] = 6;
	    }
	    break;
	case TRANSMISSION_NODE:
	    *num_sts = 4;
	    ist[0] = 0; ist[1] = 1; ist[2] = 3; ist[3] = 4;
	    *nangs = 4;
	    iang[0] = 0; iang[1] = 1; iang[2] = 3; iang[3] = 4;
	    ang_names[0] = "Upstream contact angle";
	    ang_names[1] = "Incident shock angle";
	    ang_names[2] = "Downstream contact angle";
	    ang_names[3] = "Downstream contact angle";
	    ang_names[4] = "Transmitted shock angle";
	    break;
	case CROSS_NODE:
	    *num_sts = 5;
	    for (i = 0; i < 5; i++) ist[i] = i;
	    for (i = 0; i < 4; i++)
	    	p[i] = pressure(RP->state[i+1]);
	    ang_names[0] = "Incident shock angle";
	    iang[0] = 0;
	    if (p[1] < p[0])
	    {
	    	*nangs = 6;
	    	iang[1] = 5;
	    	ang_names[1] =
	    		"Reflected rarefaction leading edge angle";
	    	iang[2] = 6;
	    	ang_names[2] =
	    		"Reflected rarefaction trailing edge angle";
	    	iang[3] = 2;
	    	ang_names[3] = "Contact angle";
	    	iang[4] = 3;
	    	ang_names[4] = "Reflected shock angle";
	    	iang[5] = 4;
	    	ang_names[5] = "Incident shock angle";
	    }
	    else if (p[2] < p[3])
	    {
	    	*nangs = 6;
	    	iang[1] = 1;
	    	ang_names[1] = "Reflected shock angle";
	    	iang[2] = 2;
	    	ang_names[2] = "Contact angle";
	    	iang[3] = 5;
	    	ang_names[3] =
	    		"Reflected rarefaction leading edge angle";
	    	iang[4] = 6;
	    	ang_names[4] =
	    		"Reflected rarefaction trailing edge angle";
	    	iang[5] = 4;
	    	ang_names[5] = "Incident shock angle";
	    }
	    else
	    {
	    	*nangs = 5;
	    	for (i = 1; i < *nangs; i++) iang[i] = i;
	    	ang_names[1] = "Reflected shock angle";
	    	ang_names[2] = "Contact angle";
	    	ang_names[3] = "Reflected shock angle";
	    	ang_names[4] = "Incident shock angle";
	    }
	    break;
	case OVERTAKE_NODE:
	    *num_sts = 7;
	    for (i = 0; i < *num_sts; i++) ist[i] = i;
	    iang[0] = 0;
	    ang_names[0] = "Ahead incident shock angle";
	    iang[1] = 1;
	    ang_names[1] = "Overtaking incident shock angle";
	    iang[0] = 0;
	    if (pressure(RP->state[5]) < pressure(RP->state[2]))
	    {
	    	*nangs = 6;
	    	iang[2] = 2;
	    	ang_names[2] =
	    		"Reflected rarefaction leading edge angle";
	    	iang[3] = 4;
	    	ang_names[3] =
	    		"Reflected rarefaction trailing edge angle";
	    	iang[4] = 5;
	    	ang_names[4] = "Contact angle";
	    	iang[5] = 6;
	    	ang_names[5] = "Transmitted shock angle";
	    }
	    else
	    {
	    	*nangs = 5;
	    	iang[2] = 2;
	    	ang_names[2] = "Reflected shock angle";
	    	iang[3] = 5;
	    	ang_names[4] = "Contact angle";
	    	iang[4] = 6;
	    	ang_names[5] = "Transmitted shock angle";
	    }
	    break;
	case MACH_NODE:
	    *num_sts = 4;
	    for (i = 0; i < *num_sts; i++) ist[i] = i;
	    *nangs = 4;
	    for (i = 0; i < *nangs; i++) iang[i] = i;
	    ang_names[0] = "Incident shock angle";
	    ang_names[1] = "Reflected shock angle";
	    ang_names[2] = "Contact angle";
	    ang_names[3] = "Mach stem angle";
	    break;
	case B_REFLECT_NODE:
	    *num_sts = 3;
	    for (i = 0; i < *num_sts; i++) ist[i] = i;
	    *nangs = 4;
	    for (i = 0; i < *nangs; i++) iang[i] = i;
	    ang_names[0] = "Upstream wall angle";
	    ang_names[1] = "Incident shock angle";
	    ang_names[2] = "Reflected shock angle";
	    ang_names[3] = "Downstream wall angle";
	    break;
	default:
	    return;
	}
}		/*end print_node_sts_params*/

/*
 *	Print time boundary state data with time dependent pressure
*/

EXPORT	void	g_fprint_tdp_boundary_state_data(
        FILE            *file,
        INTERFACE       *intfc,
        BOUNDARY_STATE  *bstate)
{
	FD_DATA *fd_data = (FD_DATA*)bstate->_boundary_state_data;
	(void) f_fprint_boundary_state_data(file,intfc,bstate);
	(void) fprintf(file,"Rise time = %g\n",fd_data->tr);
	(void) fprintf(file,"Peak time = %g\n",fd_data->tp);
	(void) fprintf(file,"Shut-off time = %g\n",fd_data->ts);
	(void) fprintf(file,"Ambient pressure = %g\n",fd_data->pr_a);
	(void) fprintf(file,"Peak pressure = %g\n",fd_data->pr_p);
        
	/* Feb 19 2004: print reference state */
        if (fd_data->state != NULL)
        {
            (void) fprintf(file,"Reference state:\n");
            Energy(fd_data->state)=energy(fd_data->state);
            set_type_of_state(fd_data->state,GAS_STATE);
            fprint_state_data(file,fd_data->state,intfc);
        }

}	/* end g_fprint_tdp_boundary_state_data */

EXPORT	boolean output_spectral_analysis(
	char		*basename,
	Wave		*wave,
	Front		*front)
{
	COMPONENT	comp;
	Locstate	state;
	double		*coords;
	double		*L = wave->rect_grid->L;
	double		*U = wave->rect_grid->U;
	double		*h = wave->rect_grid->h;
	int		icoords[MAXD];
	int		dim = wave->rect_grid->dim;
	int		status;
	int		step = front->step;
	char		energy_name[100],vorticity_name[100],
			enstrophy_name[100],dens_name[100],pres_name[100];
	FILE		*energy_file,*vorticity_file,
			*enstrophy_file,*dens_file,*pres_file;

	debug_print("fft","Entered fft_energy_spectral()\n");

	(void) sprintf(energy_name,"%s.energy%s.dat",basename,
                         right_flush(step,TSTEP_FIELD_WIDTH));
	(void) sprintf(vorticity_name,"%s.vorticity%s.dat",basename,
                         right_flush(step,TSTEP_FIELD_WIDTH));
	(void) sprintf(enstrophy_name,"%s.enstrophy%s.dat",basename,
                         right_flush(step,TSTEP_FIELD_WIDTH));
	(void) sprintf(dens_name,"%s.density%s.dat",basename,
                         right_flush(step,TSTEP_FIELD_WIDTH));
	(void) sprintf(pres_name,"%s.pressure%s.dat",basename,
                         right_flush(step,TSTEP_FIELD_WIDTH));
	energy_file    = fopen(energy_name,"w");
	vorticity_file = fopen(vorticity_name,"w");
	enstrophy_file = fopen(enstrophy_name,"w");
	dens_file = fopen(dens_name,"w");
	pres_file = fopen(pres_name,"w");

	if (wave->sizest == 0)
	{
	    debug_print("fft","Left fft_energy_spectral()\n");
	    return FUNCTION_FAILED;
	}

	switch (dim)
	{
#if defined(ONED)
	case 1:
	{
	    int		ix;
	    int		xmax;
	    COMPLEX	*mesh_energy;

	    xmax = wave->rect_grid->gmax[0];
	    uni_array(&mesh_energy,xmax,sizeof(COMPLEX));
	    for (ix = 0; ix < xmax; ++ix)
	    {
	    	icoords[0] = ix;
	    	coords = Rect_coords(icoords,wave);
	    	comp = Rect_comp(icoords,wave);
	    	state = Rect_state(icoords,wave);
		mesh_energy[ix].real = Energy(state);
		mesh_energy[ix].imag = 0.0;
	    }
	    break;
	}
#endif /* defined(ONED) */
#if defined(TWOD)
	case 2:
	{
	    int		ix, iy;
	    int		xmax, ymax, mx, my, dummy;
	    COMPLEX	**mesh_energy,**mesh_vorticity,**mesh_enstrophy;
	    Locstate	lstate,rstate,bstate,tstate;
	    double	kk,kx,ky,dk;

	    xmax = wave->rect_grid->gmax[0];
	    ymax = wave->rect_grid->gmax[1];
	    if (!Powerof2(xmax,&mx,&dummy) || !Powerof2(ymax,&my,&dummy))
	    {
		screen("fft_energy_spectral() cannot analyze "
				"mesh not power of 2\n");
		screen("xmax = %d  ymax = %d\n",xmax,ymax);
		return FUNCTION_FAILED;
	    }
	    bi_array(&mesh_energy,xmax,ymax,sizeof(COMPLEX));
	    bi_array(&mesh_vorticity,xmax,ymax,sizeof(COMPLEX));
	    fprintf(energy_file,"zone  i=%d, j=%d\n",xmax,ymax);
	    fprintf(vorticity_file,"zone  i=%d, j=%d\n",xmax,ymax);
	    fprintf(dens_file,"zone  i=%d, j=%d\n",xmax,ymax);
	    fprintf(pres_file,"zone  i=%d, j=%d\n",xmax,ymax);
	    fprintf(enstrophy_file,"zone  i=%d, j=%d\n",xmax,ymax);
	    for (iy = 0; iy < ymax; ++iy)
	    {
	    	for (ix = 0; ix < xmax; ++ix)
	    	{
	    	    icoords[0] = ix; icoords[1] = iy;
	    	    coords = Rect_coords(icoords,wave);
	    	    comp = Rect_comp(icoords,wave);
	    	    state = Rect_state(icoords,wave);
		    mesh_energy[ix][iy].real = kinetic_energy(state);
		    mesh_energy[ix][iy].imag = 0.0;

		    if (ix != 0) icoords[0] = ix - 1;
		    else icoords[0] = ix;
	    	    lstate = Rect_state(icoords,wave);
		    if (ix != xmax-1) icoords[0] = ix + 1;
		    else icoords[0] = ix;
	    	    rstate = Rect_state(icoords,wave);

		    icoords[0] = ix;
		    if (iy != 0) icoords[1] = iy - 1;
		    else icoords[1] = iy;
	    	    bstate = Rect_state(icoords,wave);
		    if (iy != ymax-1) icoords[1] = iy + 1;
		    else icoords[1] = iy;
	    	    tstate = Rect_state(icoords,wave);
		    mesh_vorticity[ix][iy].real = (Mom(rstate)[1]/Dens(rstate) 
			    		- Mom(lstate)[1]/Dens(lstate))/h[0]
					- (Mom(tstate)[0]/Dens(tstate)
					- Mom(bstate)[0]/Dens(bstate))/h[1];
		    mesh_vorticity[ix][iy].imag = 0.0;
		    fprintf(energy_file,"%lf\n",kinetic_energy(state));
		    fprintf(vorticity_file,"%lf\n",mesh_vorticity[ix][iy].real);
		    fprintf(enstrophy_file,"%lf\n",
		    		sqr(mesh_vorticity[ix][iy].real));
		    fprintf(dens_file,"%lf\n",Dens(state));
		    fprintf(pres_file,"%lf\n",pressure(state));
	    	}
	    }
	    fft_output2d(basename,"energy",step,wave->rect_grid,
	    				mesh_energy);
	    fft_output2d(basename,"vorticity",step,wave->rect_grid,
	    				mesh_vorticity);
	    free_these(2,mesh_energy,mesh_vorticity);
	    break;
	}
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:
	{
	    int		ix, iy, iz;
	    int		xmax, ymax, zmax;
	    COMPLEX	***mesh_energy;

	    xmax = wave->rect_grid->gmax[0];
	    ymax = wave->rect_grid->gmax[1];
	    zmax = wave->rect_grid->gmax[2];
	    tri_array(&mesh_energy,xmax,ymax,zmax,sizeof(COMPLEX));
	    for (iz = 0; iz < zmax; ++iz)
	    {
	    	icoords[2] = iz;
	    	for (iy = 0; iy < ymax; ++iy)
	    	{
	    	    icoords[1] = iy;
	    	    for (ix = 0; ix < xmax; ++ix)
	    	    {
	    	    	icoords[0] = ix;
	    	    	coords = Rect_coords(icoords,wave);
	    	    	comp = Rect_comp(icoords,wave);
	    	    	state = Rect_state(icoords,wave);
			mesh_energy[ix][iy][iz].real = Energy(state);
		    	mesh_energy[ix][iy][iz].imag = 0.0;
	    	    }
	    	}
	    }
	    break;
	}
#endif /* defined(THREED) */
	}
	fclose(energy_file);
	fclose(vorticity_file);
	fclose(enstrophy_file);
	fclose(dens_file);
	fclose(pres_file);

	debug_print("fft","Left fft_energy_spectral()\n");
	return FUNCTION_SUCCEEDED;
}		/*end fft_energy_spectral*/

LOCAL	void fft_output2d(
	const char *basename,
	const char *var_name,
	int step,
	RECT_GRID *grid,
	COMPLEX **cc)
{
	int i, ix, iy;
	int kmax, xmax, ymax;
	double kx, ky, kk, dk, var_k;
	double *h, *L, *U;
	char fft_name[120], mfft_name[120];
	FILE *fft_file, *mfft_file;

	xmax = grid->gmax[0];
	ymax = grid->gmax[1];
	h = grid->h;
	L = grid->L;
        U = grid->U;

	(void) sprintf(fft_name,"%s.%s.fft%s.dat",basename,var_name,
                         right_flush(step,TSTEP_FIELD_WIDTH));
	(void) sprintf(mfft_name,"%s.%s.mfft%s.dat",basename,var_name,
                         right_flush(step,TSTEP_FIELD_WIDTH));
	fft_file  = fopen(fft_name,"w");
	mfft_file = fopen(mfft_name,"w");

	fft2d(cc,xmax,ymax,1);

	for (iy = 0; iy < ymax/2; ++iy)
	{
	    for (ix = 0; ix < xmax/2; ++ix)
	    {
	    	var_k = 2.0*sqrt(sqr(cc[ix][iy].real) + sqr(cc[ix][iy].imag));
		fprintf(fft_file,"%f\n",var_k);
	    }
	}
	kmax = xmax/2;
        dk = 2.0*PI/(U[0] - L[0]);
        fprintf(mfft_file,"VARIABLES=k,'%s(k)'\n",var_name);
	for (i = 0; i < kmax; ++i)
	{
	    var_k = 0.0;
	    kk = 2.0*i*PI/(U[0] - L[0]);
	    for (iy = 0; iy < ymax/2; ++iy)
	    {
	    	for (ix = 0; ix < xmax/2; ++ix)
		{
		    kx = 2.0*ix*PI/(U[0] - L[0]);
		    ky = 2.0*iy*PI/(U[1] - L[1]);
		    if (sqrt(sqr(kx)+sqr(ky)) >= kk-0.5*dk &&
		        sqrt(sqr(kx)+sqr(ky)) < kk+0.5*dk)
		    {
		    	var_k += 2.0*sqrt(sqr(cc[ix][iy].real) +
				          sqr(cc[ix][iy].imag));
		    }
		}
	    }
	    fprintf(mfft_file,"%lf\t%lf\n",kk,var_k);
	}

	fclose(fft_file);
	fclose(mfft_file);
}	/* end fft_output */


EXPORT	boolean output_spectral_in_time(
	char		*basename,
	Wave		*wave,
	Front		*front)
{
	COMPONENT	comp;
	Locstate	state;
	double		*coords;
	double		*L = wave->rect_grid->L;
	double		*U = wave->rect_grid->U;
	double		*h = wave->rect_grid->h;
	double		kk,kx,ky,dk,Ek,Vk;
	int		icoords[MAXD];
	int		dim = wave->rect_grid->dim;
	int		status;
	char		eng_name[100],vor_name[100];
	static	FILE	*eng_file,*vor_file;
	static  int	first = YES;
	int		i, ix, iy;
	int		xmax, ymax, kmax, mx, my, dummy;
	COMPLEX		**mesh_eng,**mesh_vor;
	Locstate	lstate,rstate,bstate,tstate;

	debug_print("fft","Entered output_spectral_in_time()\n");

	if (first)
	{
	    first = NO;
	    (void) sprintf(eng_name,"%s.energy-time.dat",basename);
	    (void) sprintf(vor_name,"%s.vorticity-time.dat",basename);
	    eng_file    = fopen(eng_name,"w");
	    vor_file = fopen(vor_name,"w");
	    fprintf(eng_file,"VARIABLES=k,E(k)\n",xmax,ymax);
	    fprintf(vor_file,"VARIABLES=k,V(k)\n",xmax,ymax);
	}


 	xmax = wave->rect_grid->gmax[0];
 	ymax = wave->rect_grid->gmax[1];
	if (!Powerof2(xmax,&mx,&dummy) || !Powerof2(ymax,&my,&dummy))
	{
	    screen("output_spectral_in_time() cannot analyze "
				"mesh not power of 2\n");
	    screen("xmax = %d  ymax = %d\n",xmax,ymax);
	    return FUNCTION_FAILED;
	}
	bi_array(&mesh_eng,xmax,ymax,sizeof(COMPLEX));
	bi_array(&mesh_vor,xmax,ymax,sizeof(COMPLEX));
	for (iy = 0; iy < ymax; ++iy)
	{
	    for (ix = 0; ix < xmax; ++ix)
	    {
	    	icoords[0] = ix; icoords[1] = iy;
	    	coords = Rect_coords(icoords,wave);
	    	comp = Rect_comp(icoords,wave);
	    	state = Rect_state(icoords,wave);

		if (ix != 0) icoords[0] = ix - 1;
		else icoords[0] = ix;
	    	lstate = Rect_state(icoords,wave);
		if (ix != xmax-1) icoords[0] = ix + 1;
		else icoords[0] = ix;
	    	rstate = Rect_state(icoords,wave);

		icoords[0] = ix;
		if (iy != 0) icoords[1] = iy - 1;
		else icoords[1] = iy;
	    	bstate = Rect_state(icoords,wave);
		if (iy != ymax-1) icoords[1] = iy + 1;
		else icoords[1] = iy;
	    	tstate = Rect_state(icoords,wave);

		mesh_eng[ix][iy].real = kinetic_energy(state);
		mesh_eng[ix][iy].imag = 0.0;
		mesh_vor[ix][iy].real = (Mom(rstate)[1]/Dens(rstate) 
			    	- Mom(lstate)[1]/Dens(lstate))/h[0]
				- (Mom(tstate)[0]/Dens(tstate)
				- Mom(bstate)[0]/Dens(bstate))/h[1];
		mesh_vor[ix][iy].imag = 0.0;
	    }
	}
	fft2d(mesh_eng,xmax,ymax,1);
	fft2d(mesh_vor,xmax,ymax,1);

	kmax = xmax/2;
        dk = 2.0*PI/(U[0] - L[0]);
	fprintf(eng_file,"ZONE\n",kk,Ek);
	fprintf(vor_file,"ZONE\n",kk,Vk);
	for (i = 0; i < kmax; ++i)
	{
	    Ek = 0.0;
	    Vk = 0.0;
	    kk = 2.0*i*PI/(U[0] - L[0]);
	    for (iy = 0; iy < ymax/2; ++iy)
	    {
	    	for (ix = 0; ix < xmax/2; ++ix)
		{
		    kx = 2.0*ix*PI/(U[0] - L[0]);
		    ky = 2.0*iy*PI/(U[1] - L[1]);
		    if (sqrt(sqr(kx)+sqr(ky)) >= kk-0.5*dk &&
		        sqrt(sqr(kx)+sqr(ky)) < kk+0.5*dk)
		    {
		    	Ek += 2.0*sqrt(sqr(mesh_eng[ix][iy].real) +
				       sqr(mesh_eng[ix][iy].imag));
		    	Vk += 2.0*sqrt(sqr(mesh_vor[ix][iy].real) +
				       sqr(mesh_vor[ix][iy].imag));
		    }
		}
	    }
	    fprintf(eng_file,"%lf\t%lf\n",kk,Ek);
	    fprintf(vor_file,"%lf\t%lf\n",kk,Vk);
	}
	fflush(eng_file);
	fflush(vor_file);

	free_these(2,mesh_eng,mesh_vor);

	debug_print("fft","Left output_spectral_in_time()\n");
	return FUNCTION_SUCCEEDED;
}		/*end output_spectral_in_time*/

#endif /* defined(FULL_PHYSICS) */
