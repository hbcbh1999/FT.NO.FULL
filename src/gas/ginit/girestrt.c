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
*				girestrt.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains restart initialization routines for gas dynamics.
*
*	g_restart_initializer(), and g_restart_intfc_initializer(), called
*	from init_states() in dinit.c.
*/

#include <ginit/ginit.h>

	/* LOCAL Function Declarations */
LOCAL	void	set_restart_params_for_comp(Locstate,COMPONENT);

/*ARGSUSED*/
EXPORT void g_restart_initializer(
	int		*icoords,
	int		n_vars,
	COMPONENT	comp,
	INPUT_SOLN	**restart_soln,
	Locstate	state,
	INIT_DATA	*init)
{
	int		dim = restart_soln[0]->grid.dim;
	int		i, var;

	debug_print("init_states","Entering g_restart_initializer()\n");
	if (debugging("init_states"))
	{	
	    double coords[MAXD];
	    RECT_GRID *gr = &restart_soln[0]->grid;

	    print_int_vector("icoords = ",icoords,dim," ");
	    for (i = 0; i < dim; i++)
	    	coords[i] = cell_edge(icoords[i],i,gr);
	    print_general_vector(", coords = ",coords,dim,"\n");
	}

	set_restart_params_for_comp(state,comp);

	var = 0;
	set_type_of_state(state,GAS_STATE);
	Dens(state) = is_state(restart_soln[var++],icoords);
	Energy(state) = is_state(restart_soln[var++],icoords);
	for (i = 0; i < dim; i++)
	    Mom(state)[i] = is_state(restart_soln[var++],icoords);
	reset_gamma(state);
#if defined(COMBUSTION_CODE)
	if (Params(state) && Composition_type(state) == ZND)
	    Prod(state) = is_state(restart_soln[var++],icoords);
#endif /* defined(COMBUSTION_CODE) */
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            int num_comps;
            num_comps = n_vars - var;

            for(i = 0; i < num_comps; i++)
            {
                pdens(state)[i] = is_state(restart_soln[var++],icoords);
                /* New 051005 */
                if(fabs(pdens(state)[i]) < 10.0*MACH_EPS && pdens(state)[i] < 0.0)
                    pdens(state)[i] = 0.0;
                /*
                else if(fabs(pdens(state)[i]) > 10.0*MACH_EPS && pdens(state)[i] < 0.0)
                {
                    printf("ERROR in g_restart_initializer()\n");
                    printf("partial density %20.18g < 0.0\n",pdens(state)[i]);
                    clean_up(ERROR);
                }
                */
                /* End of New 051005 */
            }
        }

	if (debugging("init_states"))
	{
	    g_print_state(state);
	}
	if (debugging("check_restart"))
	{
	    if (is_bad_state(state,YES,"g_restart_initializer"))
	    {
	        screen("ERROR in g_restart_initializer(), bad state "
		       "initialized\n");
	        g_print_state(state);
		clean_up(ERROR);
	    }
	}
	debug_print("init_states","Leaving g_restart_initializer()\n");
}		/*end g_restart_initializer*/




/*ARGSUSED*/
EXPORT void g_restart_intfc_initializer(
	POINT		   *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF	   *hs,
	Locstate	   lstate,
	Locstate	   rstate,
	INIT_DATA	   *init)
{
	INTERFACE	*intfc = hs->interface;
	debug_print("init_states","Entering g_restart_intfc_initializer()\n");
	if (debugging("init_states"))
	    print_general_vector("Coords(p) = ",Coords(p),intfc->dim,"\n");

	switch (wave_type(hs))
	{
	case PASSIVE_BOUNDARY:
	case SUBDOMAIN_BOUNDARY:
	    obstacle_state(intfc,lstate,size_of_state(intfc));
	    obstacle_state(intfc,rstate,size_of_state(intfc));
	    return;

	default:
	    break;
	}

	set_restart_params_for_comp(lstate,negative_component(hs));
	set_restart_params_for_comp(rstate,positive_component(hs));

	if ((comp_type(negative_component(hs))->type == EXTERIOR) ||
	    (comp_type(negative_component(hs))->type == OBSTACLE))
	{
	    switch (wave_type(hs)) 
	    {
	    case DIRICHLET_BOUNDARY:
	    case NEUMANN_BOUNDARY:
	    case PASSIVE_BOUNDARY:
	    case SUBDOMAIN_BOUNDARY:
	    case MOVABLE_BODY_BOUNDARY:
	    	obstacle_state(intfc,lstate,size_of_state(intfc));
	    	break;

	    default:
	    	screen("ERROR: unknown bdry type "
	    	       "in g_restart_intfc_initializer\n");
	    	clean_up(ERROR);
	    }
	}

	if ((comp_type(positive_component(hs))->type == EXTERIOR) ||
	    (comp_type(positive_component(hs))->type == OBSTACLE))
	{
	    switch (wave_type(hs)) 
	    {
	    case DIRICHLET_BOUNDARY:
	    case NEUMANN_BOUNDARY:
	    case PASSIVE_BOUNDARY:
	    case SUBDOMAIN_BOUNDARY:
	    case MOVABLE_BODY_BOUNDARY:
	    	obstacle_state(intfc,rstate,size_of_state(intfc));
	    	break;

	    default:
	    	screen("ERROR: unknown bdry type "
	    	       "in g_restart_intfc_initializer\n");
	    	clean_up(ERROR);
	    }
	}

	if (debugging("init_states")) 
	{
	    (void) printf("left state\n");
	    g_print_state(lstate);
	    (void) printf("right state\n");
	    g_print_state(rstate);
	}
	if (debugging("check_restart"))
	{
	    boolean lbad = is_bad_state(lstate,NO,"g_restart_intfc_initializar");
	    boolean rbad = is_bad_state(rstate,NO,"g_restart_intfc_initializar");
	    if (lbad || rbad)
	    {
	        screen("ERROR in g_restart_intfc_initializer(), bad %s state "
		       "initialized\n",(lbad && !rbad) ? "left" :
		                       (rbad && !lbad) ? "right" :
				       "left && right");
		if (lbad)
	            g_print_state(lstate);
		if (rbad)
	            g_print_state(rstate);
		clean_up(ERROR);
	    }
	}
	debug_print("init_states","Leaving g_restart_intfc_initializer()\n");
}	    /*end g_restart_intfc_initializar*/


/*ARGSUSED*/
EXPORT void g_restart_set_intfc_states(
	double		   *sl,
	double		   *sr,
	int		   var,
	POINT		   *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF	   *hs)
{
	int		dim = hs->interface->dim;
	Locstate	lstate, rstate;

	slsr(p,hse,hs,&lstate,&rstate);
	set_type_of_state(lstate,GAS_STATE);
	set_type_of_state(rstate,GAS_STATE);
	if (var == 0)
	{
	    Dens(lstate) = *sl;
	    Dens(rstate) = *sr;
	}
	else if (var == 1)
	{
	    Energy(lstate) = *sl;
	    Energy(rstate) = *sr;
	}
	else if ((2 <= var) && (var < (2+dim)))
	{
	    Mom(lstate)[var-2] = *sl;
	    Mom(rstate)[var-2] = *sr;
	}
#if defined(COMBUSTION_CODE)
	else if (var == (2+dim))
	{
	   if (comp_type(positive_component(hs))->params->composition_type==ZND)
		Prod(rstate) = *sr;
	   if (comp_type(negative_component(hs))->params->composition_type==ZND)
		Prod(lstate) = *sl;
	}
#endif /* defined(COMBUSTION_CODE) */
	reset_gamma(lstate);
	reset_gamma(rstate);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(comp_type(positive_component(hs))->params != NULL &&
               comp_type(positive_component(hs))->params->n_comps != 1)
            {
                pdens(rstate)[var-(2+dim)] = *sr;
                /* New 051005 */
                if(fabs(pdens(rstate)[var-(2+dim)]) < 10.0*MACH_EPS && 
			pdens(rstate)[var-(2+dim)] < 0.0)
                    pdens(rstate)[var-(2+dim)] = 0.0;
                /* End of New 051005 */
            }
            if(comp_type(negative_component(hs))->params != NULL &&
               comp_type(negative_component(hs))->params->n_comps != 1)
            {
                pdens(lstate)[var-(2+dim)] = *sl;
                if(fabs(pdens(lstate)[var-(2+dim)]) < 10.0*MACH_EPS && 
			pdens(lstate)[var-(2+dim)] < 0.0)
                    pdens(lstate)[var-(2+dim)] = 0.0;
            }
        }
}		/*end g_restart_set_intfc_states*/


LOCAL void set_restart_params_for_comp(
	Locstate	state,
	COMPONENT	comp)
{
	Init_params(state,comp_type(comp)->params);
}		/*end set_restart_params_for_comp*/
