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
*				gibifur.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Prompts for dynamic tracking options.
*
*	Exported function:
*		g_prompt_for_tracked_bifurcations
*/

#if defined(FULL_PHYSICS)
#include <ginit/ginit.h>

	/* LOCAL Function Declarations */
LOCAL	void	prompt_for_tracking_at_node(int,int,int,SCAT_WV_TOL*);
LOCAL	void	set_tracking_defaults(SCAT_WV_TOL*,int,UNTRACK_NODE_FLAG*,
				      int,int);
LOCAL	void	print_tracking_choices(SCAT_WV_TOL*,int,UNTRACK_NODE_FLAG*,int);
#if defined(TWOD)
LOCAL	void	prompt_for_node_untrack(int,UNTRACK_NODE_FLAG*);
#endif /* defined(TWOD) */

EXPORT	void	g_prompt_for_tracked_bifurcations(
	INIT_DATA	*init)
{
	char		s[Gets_BUF_SIZE];
	SCAT_WV_TOL		*swt;
	UNTRACK_NODE_FLAG	*unf;
	int                     dim = i_intfc(init)->dim;
	int			n_swtl, n_unfl;

	scattered_wave_tolerance_list(init) = NULL;
	untrack_node_options_list(init) = NULL;
	if (strcmp(ex_name(init),"spolars") == 0)
	    return;
	n_swtl = NUM_PHYS_NODE_TYPES*NUM_WAVE_CLASSES*NUM_NODE_STATUS;
	n_unfl = NUM_PHYS_NODE_TYPES;
	uni_array(&swt,n_swtl,sizeof(SCAT_WV_TOL));
	scattered_wave_tolerance_list(init) = swt;
	uni_array(&unf,n_unfl,sizeof(UNTRACK_NODE_FLAG));
	untrack_node_options_list(init) = unf;
	set_tracking_defaults(swt,n_swtl,unf,n_unfl,dim);
	print_tracking_choices(swt,n_swtl,unf,n_unfl);

	screen("Use default settings for dynamic tracking (default = y): ");
	(void) Gets(s);
	if (s[0] != 'n' && s[0] != 'N')
	    return;

	switch (dim)
	{
#if defined(ONED)
	case 1:
	    prompt_for_tracking_at_node(ONED_INTERACTION,SHOCK_WAVE,
					ONED_WAVE,swt);
	    prompt_for_tracking_at_node(ONED_INTERACTION,RAREF_LEADING_EDGE,
					ONED_WAVE,swt);
	    prompt_for_tracking_at_node(ONED_INTERACTION,RAREF_TRAILING_EDGE,
					ONED_WAVE,swt);
	    prompt_for_tracking_at_node(ONED_INTERACTION,CONTACT_WAVE,
					ONED_WAVE,swt);
	    break;
#endif /*defined(ONED)*/
#if defined(TWOD)
	case 2:
	    prompt_for_tracking_at_node(B_REFLECT_NODE,SHOCK_WAVE,
					REFLECTED,swt);
	    prompt_for_tracking_at_node(ATTACHED_B_NODE,SHOCK_WAVE,
					REFLECTED,swt);
	    prompt_for_tracking_at_node(MACH_NODE,SHOCK_WAVE,MACH_STEM,swt);
	    prompt_for_tracking_at_node(MACH_NODE,SHOCK_WAVE,REFLECTED,swt);
	    prompt_for_tracking_at_node(MACH_NODE,CONTACT_WAVE,SLIP,swt);
	    prompt_for_tracking_at_node(CROSS_NODE,SHOCK_WAVE,REFLECTED,swt);
	    prompt_for_tracking_at_node(CROSS_NODE,CONTACT_WAVE,SLIP,swt);
	    prompt_for_tracking_at_node(OVERTAKE_NODE,SHOCK_WAVE,REFLECTED,swt);
	    prompt_for_tracking_at_node(OVERTAKE_NODE,
					RAREF_LEADING_EDGE,REFLECTED,swt);
	    prompt_for_tracking_at_node(OVERTAKE_NODE,
					RAREF_TRAILING_EDGE,REFLECTED,swt);
	    prompt_for_tracking_at_node(OVERTAKE_NODE,SHOCK_WAVE,
					TRANSMITTED,swt);
	    prompt_for_tracking_at_node(OVERTAKE_NODE,CONTACT_WAVE,SLIP,swt);
	    prompt_for_tracking_at_node(DIFFRACTION_NODE,SHOCK_WAVE,
					REFLECTED,swt);
	    prompt_for_tracking_at_node(DIFFRACTION_NODE,
				        RAREF_LEADING_EDGE,REFLECTED,swt);
	    prompt_for_tracking_at_node(DIFFRACTION_NODE,
				        RAREF_TRAILING_EDGE,REFLECTED,swt);
	    prompt_for_tracking_at_node(DIFFRACTION_NODE,SHOCK_WAVE,TRANSMITTED,
				        swt);
	    prompt_for_tracking_at_node(DIFFRACTION_NODE,CONTACT_WAVE,SLIP,swt);
	    prompt_for_tracking_at_node(TRANSMISSION_NODE,
				        SHOCK_WAVE,TRANSMITTED,swt);
	    prompt_for_tracking_at_node(TRANSMISSION_NODE,CONTACT_WAVE,
					SLIP,swt);
	    screen("\n");

	    /*
	     * Currently, these are the only nodes with untracking algorithms.
	     */

	    prompt_for_node_untrack(B_REFLECT_NODE,unf);
	    prompt_for_node_untrack(MACH_NODE,unf);
	    prompt_for_node_untrack(OVERTAKE_NODE,unf);
	    prompt_for_node_untrack(PRECURSOR_RR_DIFFRACTION,unf);
	    break;
#endif /* defined(TWOD) */

#if defined(THREED)
	case 3: /*Not implemented*/
	    break;
#endif /* defined(THREED) */
	}

	print_tracking_choices(swt,n_swtl,unf,n_unfl);
	screen("This completes the prompting for dynamical "
	       "tracking decisions.\n");
}		/*end g_prompt_for_tracked_bifurcations*/

LOCAL	void	print_tracking_choices(
	SCAT_WV_TOL	  *swt,
	int		  n_swtl,
	UNTRACK_NODE_FLAG *unf,
	int		  n_unfl)
{
	char	   s[80], line[1024];
	const char *mesg;
	int	   i;

	screen("\nDynamic tracking decision variables\n");
	screen_print_long_string(
	    "Tracking decisions on dynamically produced waves are based on a "
	    "floating point cutoff on the wave strength. "
	    "Strengths are normalized to zero for weak waves, "
	    "so a tolerance of 0.0 will always signal tracking, and a very "
	    "large  tolerance will always signal not to track. "
	    "For each interaction type, you will be asked to enter the "
	    "cutoff tolerance and a tolerance type for determining whether "
	    "scattered waves of the indicated type should be tracked when "
	    "produced by a specific bifurcation type.\n"
	);
	screen("The currently supported tolerance types are\n\t"
	       "Never track (Never)\n\t"
	       "Always track (Always)\n\t"
	       "Pressure ratio minus one across the wave (Pressure)\n\t"
	       "Absolute value of the Atwood number across the "
	           "wave (Atwood)\n\t"
	       "Mach number minus one for the state "
	           "ahead of the wave (Mach)\n");

	for (i = 0; i < n_swtl; i++)
	{
	    if (swt[i].wave_name == NULL)
		continue;
	    mesg = swt[i].wave_name;
	    if ((swt[i].wave_tol == 0.0) || (swt[i].tol_type == ALWAYS_TRACK))
		(void) sprintf(s,"always track");
	    else if ((swt[i].wave_tol == HUGE_VAL) ||
		     (swt[i].tol_type == NEVER_TRACK))
		(void) sprintf(s,"never track");
	    else
		(void) sprintf(s,"%g %s",swt[i].wave_tol,swt[i].tol_type_name);
	    (void) sprintf(line,"Wave strength tolerance for "
				"tracking %s = %s\n",mesg,s);
	    screen_print_long_string(line);
	}

	for (i = 0; i < n_unfl; i++)
	{
	    const char *mesg = unf[i]._node_name;
	    if (mesg == NULL)
	        continue;
	    if (unf[i]._untrack == YES)
	    {
	        (void) sprintf(line,"Turn off tracking of %s "
				    "if node propagation fails",mesg);
	    }
	    else if (unf[i]._untrack == NO)
	    {
	        (void) sprintf(line,"Don't Turn off tracking of %s "
				    "if node propagation fails\n",mesg);
	    }
	    screen_print_long_string(line);
	}
	screen("End Dynamic tracking decision variables\n\n");
}		/*end print_tracking_choices*/

LOCAL	void	set_tracking_defaults(
	SCAT_WV_TOL	  *swt,
	int		  n_swtl,
	UNTRACK_NODE_FLAG *unf,
	int		  n_unfl,
	int               dim)
{
	int	i;

	/* Set defaults */

	for (i = 0; i < n_swtl; i++)
	{
	    swt[i].wave_name = NULL;
	    swt[i].wave_tol = HUGE_VAL;
	    swt[i].tol_type = NEVER_TRACK;
	}
	for (i = 0; i < n_unfl; i++)
	{
	    unf[i]._node_name = NULL;
	    unf[i]._untrack = NO;
	}

#if defined(TWOD)
	if (dim == 2)
	{
	    unf[untrack_node_index(B_REFLECT_NODE)]._node_name =
	        "regular reflection node";
	    unf[untrack_node_index(MACH_NODE)]._node_name = "Mach node";
	    unf[untrack_node_index(OVERTAKE_NODE)]._node_name = "overtake node";
	    unf[untrack_node_index(PRECURSOR_RR_DIFFRACTION)]._node_name =
	        "precursor rr diffraction (cluster)";
	}
#endif /* defined(TWOD) */


	switch (dim)
	{
#if defined(ONED)
	case 1:
	    /* one dimensional Riemann problems */
	    i = scat_wv_index(ONED_INTERACTION,SHOCK_WAVE,ONED_WAVE);
	    swt[i].wave_tol = 0.1;
	    swt[i].tol_type = MACH_NUMBER_AHEAD;
	    swt[i].tol_type_name = "shock ahead Mach number minus one";
	    swt[i].wave_name = "shocks at one dimensional wave interactions";

	    i = scat_wv_index(ONED_INTERACTION,RAREF_LEADING_EDGE,ONED_WAVE);
	    swt[i].wave_tol = HUGE_VAL;
	    swt[i].tol_type = NEVER_TRACK;
	    swt[i].tol_type_name = "never track";
	    swt[i].wave_name = "rarefaction leading edges at one dimensional "
			       "wave interactions";

	    i = scat_wv_index(ONED_INTERACTION,RAREF_TRAILING_EDGE,ONED_WAVE);
	    swt[i].wave_tol = HUGE_VAL;
	    swt[i].tol_type = NEVER_TRACK;
	    swt[i].tol_type_name = "never track";
	    swt[i].wave_name = "rarefaction trailing edges at one dimensional "
			       "wave interactions";

	    i = scat_wv_index(ONED_INTERACTION,CONTACT_WAVE,ONED_WAVE);
	    swt[i].wave_tol = 0.1;
	    swt[i].tol_type = ATWOOD_NUMBER;
	    swt[i].tol_type_name = "Atwood number across wave";
	    swt[i].wave_name = "contacts at one dimensional wave interactions";
	    break;
#endif /* defined(ONED) */
#if defined(TWOD)
	case 2:
	    /* Boundary reflect nodes */
	    i = scat_wv_index(B_REFLECT_NODE,SHOCK_WAVE,REFLECTED);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = PRESSURE_RATIO;
	    swt[i].tol_type_name = "pressure ratio minus one across wave";
	    swt[i].wave_name = "reflected shocks at regular reflections";

	    /* Attached boundary nodes */
	    i = scat_wv_index(ATTACHED_B_NODE,SHOCK_WAVE,REFLECTED);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = PRESSURE_RATIO;
	    swt[i].tol_type_name = "pressure ratio minus one across wave";
	    swt[i].wave_name =
	        "reflected shocks at attached boundary reflection nodes";

	    /* Mach nodes */
	    i = scat_wv_index(MACH_NODE,SHOCK_WAVE,REFLECTED);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = PRESSURE_RATIO;
	    swt[i].tol_type_name = "pressure ratio minus one across wave";
	    swt[i].wave_name = "reflected shocks at Mach reflections";

	    i = scat_wv_index(MACH_NODE,SHOCK_WAVE,MACH_STEM);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = PRESSURE_RATIO;
	    swt[i].tol_type_name = "pressure ratio minus one across wave";
	    swt[i].wave_name = "the Mach stem at Mach reflections";

	    i = scat_wv_index(MACH_NODE,CONTACT_WAVE,SLIP);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = ATWOOD_NUMBER;
	    swt[i].tol_type_name = "Atwood number across wave";
	    swt[i].wave_name = "the slip line at Mach reflections";

	    /* Cross nodes */
	    i = scat_wv_index(CROSS_NODE,CONTACT_WAVE,SLIP);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = ATWOOD_NUMBER;
	    swt[i].tol_type_name = "Atwood number wave";
	    swt[i].wave_name = "slip lines produced by shock crossings";

	    i = scat_wv_index(CROSS_NODE,SHOCK_WAVE,REFLECTED);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = PRESSURE_RATIO;
	    swt[i].tol_type_name = "pressure ratio minus one across wave";
	    swt[i].wave_name = "reflected shocks at shock crossings";

	    /* Overtake nodes */
	    i = scat_wv_index(OVERTAKE_NODE,CONTACT_WAVE,SLIP);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = ATWOOD_NUMBER;
	    swt[i].tol_type_name = "Atwood number across wave";
	    swt[i].wave_name = "slip lines at shock overtakes";

	    i = scat_wv_index(OVERTAKE_NODE,SHOCK_WAVE,REFLECTED);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = PRESSURE_RATIO;
	    swt[i].tol_type_name = "pressure ratio minus one across wave";
	    swt[i].wave_name = "reflected shocks at shock overtakes";

	    i = scat_wv_index(OVERTAKE_NODE,SHOCK_WAVE,TRANSMITTED);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = PRESSURE_RATIO;
	    swt[i].tol_type_name = "pressure ratio minus one across wave";
	    swt[i].wave_name = "transmitted shocks at shock overtakes";

	    i = scat_wv_index(OVERTAKE_NODE,RAREF_LEADING_EDGE,REFLECTED);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = PRESSURE_RATIO;
	    swt[i].tol_type_name = "pressure ratio minus one across wave";
	    swt[i].wave_name =
	        "reflected rarefaction leading edges at shock overtakes";

	    i = scat_wv_index(OVERTAKE_NODE,RAREF_TRAILING_EDGE,REFLECTED);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = PRESSURE_RATIO;
	    swt[i].tol_type_name = "pressure ratio minus one across wave";
	    swt[i].wave_name = "reflected rarefaction trailing edges "
	    		   "at shock overtakes";

	    /* diffraction nodes */
	    i = scat_wv_index(DIFFRACTION_NODE,CONTACT_WAVE,SLIP);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = ATWOOD_NUMBER;
	    swt[i].tol_type_name = "Atwood number across wave";
	    swt[i].wave_name = "material interfaces at shock-contact "
			       "diffractions";

	    i = scat_wv_index(DIFFRACTION_NODE,SHOCK_WAVE,REFLECTED);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = PRESSURE_RATIO;
	    swt[i].tol_type_name = "pressure ratio minus one across wave";
	    swt[i].wave_name = "reflected shocks at shock-contact diffractions";

	    i = scat_wv_index(DIFFRACTION_NODE,SHOCK_WAVE,TRANSMITTED);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = PRESSURE_RATIO;
	    swt[i].tol_type_name = "pressure ratio minus one across wave";
	    swt[i].wave_name = "transmitted shocks at shock-contact "
			       "diffractions";

	    i = scat_wv_index(DIFFRACTION_NODE,RAREF_LEADING_EDGE,REFLECTED);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = PRESSURE_RATIO;
	    swt[i].tol_type_name = "pressure ratio minus one across wave";
	    swt[i].wave_name = "reflected rarefaction leading edges "
	    		   "at shock-contact diffractions";

	    i = scat_wv_index(DIFFRACTION_NODE,RAREF_TRAILING_EDGE,REFLECTED);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = PRESSURE_RATIO;
	    swt[i].tol_type_name = "pressure ratio minus one across wave";
	    swt[i].wave_name = "reflected rarefaction trailing edges at "
	                       "shock-contact diffractions";

	    /* transmission nodes */
	    i = scat_wv_index(TRANSMISSION_NODE,CONTACT_WAVE,SLIP);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = ATWOOD_NUMBER;
	    swt[i].tol_type_name = "Atwood number across wave";
	    swt[i].wave_name = "material interfaces at shock-contact "
	    		   "transmission nodes";

	    i = scat_wv_index(TRANSMISSION_NODE,SHOCK_WAVE,TRANSMITTED);
	    swt[i].wave_tol = 0.0;
	    swt[i].tol_type = PRESSURE_RATIO;
	    swt[i].tol_type_name = "pressure ratio minus one across wave";
	    swt[i].wave_name = "transmitted shocks at shock-contact "
	    		   "transmission nodes";
	    break;
#endif /* defined(TWOD) */
	}
}		/*end set_tracking_defaults*/


LOCAL	void prompt_for_tracking_at_node(
	int		n_type,
	int		w_class,
	int		status,
	SCAT_WV_TOL	*swtl)
{
	int	    i = scat_wv_index(n_type,w_class,status);
	SCAT_WV_TOL *swt = swtl+i;
	const char  *mesg = swt->wave_name;
	char	    s[Gets_BUF_SIZE];
	char	    line[1024];

	if ((swt->wave_tol == 0.0) || (swt->tol_type == ALWAYS_TRACK))
	    (void) sprintf(s,"yes");
	else if ((swt->wave_tol == HUGE_VAL) || (swt->tol_type == NEVER_TRACK))
	    (void) sprintf(s,"no");
	else
	{
	    (void) sprintf(s,"%g %s",swt->wave_tol,swt->tol_type_name);
	}

	(void) sprintf(line,
	     "Wave strength tolerance for tracking %s, "
	     "yes = always track, no = never track, "
	     "otherwise enter the tolerance and tolerance type "
	     "(default = %s)",mesg,s);
	screen_print_long_string(line);
	screen(": ");
	(void) Gets(s);
	if (s[0] == 'n' || s[0] == 'N')
	{
	    swt->wave_tol = HUGE_VAL;
	    swt->tol_type = NEVER_TRACK;
	    swt->tol_type_name = "never track";
	}
	else if (s[0] == 'y' || s[0] == 'Y' || s[0] == 'a' || s[0] == 'A')
	{
	    swt->wave_tol = 0.0;
	    swt->tol_type = ALWAYS_TRACK;
	    swt->tol_type_name = "always track";
	}
	else if (s[0] != '\0')
	{
	    char *c;

	    for (c = s; *c != '\0'; ++c)
	        *c = tolower(*c);

	    if (strstr(s,"never track") != NULL)
	    {
	        swt->wave_tol = HUGE_VAL;
	        swt->tol_type = NEVER_TRACK;
	        swt->tol_type_name = "never track";
	    }
	    else if (strstr(s,"always track") != NULL)
	    {
	        swt->wave_tol = 0.0;
	        swt->tol_type = ALWAYS_TRACK;
	        swt->tol_type_name = "always track";
	    }
	    else
	    {
	        if (sscan_float(s,&swt->wave_tol) == 1)
		{
	            if (strstr(s,"pressure ratio") != NULL)
		    {
	                swt->tol_type = PRESSURE_RATIO;
	                swt->tol_type_name =
			    "pressure ratio minus one across wave";
		    }
	            else if ((strstr(s,"mach number") != NULL) &&
		             (w_class == SHOCK_WAVE))
		    {
	                swt->tol_type = MACH_NUMBER_AHEAD;
	                swt->tol_type_name =
			    "shock ahead Mach number minus one";
		    }
	            else if ((strstr(s,"atwood number") != NULL) ||
		             (strstr(s,"density ratio") != NULL))
		    {
	                swt->tol_type = ATWOOD_NUMBER;
	                swt->tol_type_name = "Atwood number across wave";
		    }
	        }
		else
	        {
	            swt->wave_tol = HUGE_VAL;
	            swt->tol_type = NEVER_TRACK;
	            swt->tol_type_name = "never track";
	        }
	    }
	}
}		/*end prompt_for_tracking_at_node*/


#if defined(TWOD)
LOCAL void prompt_for_node_untrack(
	int		  n_type,
	UNTRACK_NODE_FLAG *unfl)
{
	int               i = untrack_node_index(n_type);
	UNTRACK_NODE_FLAG *unf = unfl + i;
	const char	  *mesg = unf->_node_name;
	char	          s[Gets_BUF_SIZE];

	if (mesg == NULL)
	    return;

	screen("Type 'y' to use untracking if propagation of\n");
	screen("\t%s fails (default = %s): ",mesg,
	       (unf->_untrack == NO) ? "'n'" : "'y'");
	(void) Gets(s);

	if (s[0] == 'n' || s[0] == 'N')
	    unf->_untrack = NO;
	else if (s[0] == 'y' || s[0] == 'Y')
	    unf->_untrack = YES;
}		/*end prompt_for_node_untrack*/
#endif /* defined(TWOD) */

#endif /* defined(FULL_PHYSICS) */
