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
*
*				giniteos.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains equation of state dependent functions needed for 
*	initialization and restart.
*
*	Functions for initialization
*
*		prompt_for_equation_of_state
*
*/

#define	DEBUG_STRING	"init_eos"
#include <geos/geosdecs.h>

	/* LOCAL Function Declarations */
LOCAL	EOS	*set_equation_of_state(INIT_DATA*,int);

EXPORT	Gas_param *read_print_EOS_data(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	Gas_param     *params)
{
	FILE    *file = io_type->file;
	int	equation_of_state;

	DEBUG_ENTER(read_print_EOS_data)
	if (!fgetstring(file,"Equation of state = "))
	{
	    screen("ERROR in read_print_EOS_data(), "
	           "can't read equation of state type\n");
	    clean_up(ERROR);
	}
	(void) fscanf(file,"%d",&equation_of_state);
	switch (equation_of_state)
	{
	case	OBSTACLE_EOS:
	    DEBUG_LEAVE(read_print_EOS_data)
	    return Params(return_obst_state());
	default:
	    if (params == NULL)
	    	params = alloc_Gas_param(sizeof(Gas_param));
	    params->eos = set_equation_of_state(init,equation_of_state);
	    read_print_EOS_params(init,io_type,params);
	    break;
	}

	DEBUG_LEAVE(read_print_EOS_data)
	return params;
}		/*end read_print_EOS_data*/

EXPORT void g_prompt_for_equation_of_state(
	INIT_DATA  *init,
	Gas_param  **params,
	const char *msg1,
	const char *msg2,
	INIT_PHYSICS    *ip)
{
	F_USER_INTERFACE *fuh;
	char	*message1 = NULL, *message2 = NULL;
	int	equation_of_state;

        static          Prompt_type Promtype[] =
	{
            {"Obstacle (behind reflecting wall)", "O", 1,OBSTACLE_EOS},
	    {"Polytropic (gamma law) gas", "P", 1,POLYTROPIC},
	    {"Stiffened polytropic gas", "SP", 2,STIFFENED_POLYTROPIC},
	    {"Stellar (helmholtz eos) gas", "ST", 2,STELLAR},
#if defined(MULTI_COMPONENT)
	    {"Multiple component polytropic gas","MP",2,MULTI_COMP_POLYTROPIC},
#endif /* defined(MULTI_COMPONENT) */
#if defined(SESAME_CODE) && defined(TWOD)
            {"Sesame table lookup", "SE", 2, SESAME},
#endif /* defined(SESAME_CODE) && defined(TWOD) */
	    {"JWL Equation of state", "J", 1,JWL},
	    {"Mie Gruneisen", "M", 1,MIE_GRUNEISEN},
	    {"Isentropic two phase eos", "S2PH", 4,ISENTROPIC_TWO_PHASE},
	    {"Generic Test (minimal SPOLY for GENERIC testing)","GT",2,GENTEST},
            {NULL, NULL, 0, UNKNOWN_EOS}
        };
	
        DEBUG_ENTER(g_prompt_for_equation_of_state)
	equation_of_state = prompt_for_equation_of_state_type(msg1,msg2,
			                                      &message1,
							      &message2,
							      &Promtype[0],
							      POLYTROPIC);

	switch (equation_of_state)
	{
	case	OBSTACLE_EOS:
	    *params = Params(return_obst_state());
	    break;
	default:
	    *params = alloc_Gas_param(sizeof(Gas_param));
	    (*params)->eos = set_equation_of_state(init,equation_of_state);
	    (*params)->dim = ip->root->front->rect_grid->dim;
	    fuh = f_user_hook(ip->root->front->rect_grid->dim);
	    (*params)->_alloc_state = fuh->_alloc_state;
	    (*params)->_alloc_intfc_state = fuh->_alloc_intfc_state;
	    (*params)->_clear_state = fuh->_clear_state;
	    (*params)->sizest = ip->root->front->sizest;
	    prompt_for_EOS_params(init,*params,msg1,msg2);
	    break;
	}
	if (message1 != NULL)
	    free(message1);
	if (message2 != NULL)
	    free(message2);
	DEBUG_LEAVE(g_prompt_for_equation_of_state)
}		/*end g_prompt_for_equation_of_state*/

EXPORT  int    		prompt_for_equation_of_state_type(
        const char  *msg1,
        const char  *msg2,
        char        **message1,
        char        **message2,
	Prompt_type *Promtype,
	int	    default_equation_of_state)
{
        char            s[Gets_BUF_SIZE];
        int             equation_of_state;
        int             i;
        size_t		len1, len2;
 
        DEBUG_ENTER(prompt_for_equation_of_state_type)
        if (msg1 == NULL)
	    msg1 = "";
        len1 = strlen(msg1);
	if (*message1 == NULL)
            scalar(message1,len1+2);
        if (len1 != 0)
        {
            if (msg1[0] == ' ')
		(void) strcpy(*message1,msg1);
            else
		(void) sprintf(*message1," %s",msg1);
        }
        if (msg2 == NULL)
	    msg2 = "";
        len2 = strlen(msg2);
	if (*message2 == NULL)
            scalar(message2,len2+2);
        if (len2 != 0)
        {
            if (msg2[0] == ' ')
		(void) strcpy(*message2,msg2);
            else
		(void) sprintf(*message2," %s",msg2);
        }
        
        screen("Enter the equation of state type for the%s material%s.\n",
                *message1,*message2);
        screen("Current choices are\n");
        for (i = 0; Promtype[i].prompt != NULL; i++)
        {
            screen("\t\t");
            if (Promtype[i+1].prompt == NULL)
		screen("or ");
            screen("%s (%s)",Promtype[i].prompt,Promtype[i].select);
            screen("%s\n",(Promtype[i+1].prompt == NULL) ? "." : ",");
        }
        screen("\tEnter choice here (dflt = %s): ",Promtype[1].select);
        (void) Gets(s);

        equation_of_state = default_equation_of_state;
        for (i = 0; Promtype[i].prompt != NULL; i++)
        {
            if (strncasecmp(s,Promtype[i].select,Promtype[i].ncmp)==0)
            {
                equation_of_state = Promtype[i].type.itype;
                break;
            }
        }
        if (equation_of_state == UNKNOWN_EOS)
        {
            screen("ERROR in prompt_for_equation_of_state_type(), "
		   "unrecognized equation_of_state type\n");
            clean_up(ERROR);
        }
	DEBUG_LEAVE(prompt_for_equation_of_state_type)
        return	equation_of_state;
}               /*end prompt_for_equation_of_state_type*/

LOCAL	EOS	*set_equation_of_state(
	INIT_DATA *init,
	int	  equation_of_state)
{
	EOS	*eos;
	char	pathstore[256],*path;

        DEBUG_ENTER(set_equation_of_state)
	strcpy(pathstore,output_filename(init));
	path = get_dirname(dirname(pathstore));
	switch (equation_of_state)
	{
	case	POLYTROPIC:
	    eos = set_POLY_eos(NULL);
	    break;

	case	STIFFENED_POLYTROPIC:
	    eos = set_SPOLY_eos(NULL);
	    break;

#if defined(COMBUSTION_CODE)
	case	STELLAR:
	    eos = set_STELLAR_eos(NULL,path);
	    break;
#endif /* defined(COMBUSTION_CODE) */

#if defined(MULTI_COMPONENT)
	case MULTI_COMP_POLYTROPIC:
	    eos = set_MPOLY_eos(NULL);
	    break;
#endif /* defined(MULTI_COMPONENT) */
	
#if defined(SESAME_CODE) && defined(TWOD)
	case	SESAME:
	    eos = set_SESAME_eos(init,NULL);
	    break;
#endif /* defined(SESAME_CODE) && defined(TWOD) */

	case	JWL:
	    eos = set_JWL_eos(NULL);
	    break;

	case	MIE_GRUNEISEN:
	    eos = set_MG_eos(NULL);
	    break;

	case	GENTEST:
	    eos = set_GENTEST_eos(NULL);
	    break;

	case	ISENTROPIC_TWO_PHASE:
	    eos = set_S2PHASE_eos(NULL);
	    break;

	default:
	    eos = NULL;
	    screen("ERROR in set_equation_of_state() "
	           "Unknown equation of state\n");
	    clean_up(ERROR);
	    break;
	}
	DEBUG_LEAVE(set_equation_of_state)
	return eos;
}	/*end set_equation_of_state*/
