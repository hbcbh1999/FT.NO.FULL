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
*				gisc.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains initialization routines for shock-contact
*	interactions.
*
*/


#if defined(FULL_PHYSICS)
#include <ginit/ginit.h>

enum _INIT_TAYLOR_WAVE_TYPE {
	DIRECT	 = 0,
	INDIRECT = 1
};
typedef enum _INIT_TAYLOR_WAVE_TYPE INIT_TAYLOR_WAVE_TYPE;

	/* LOCAL Function Declarations */
LOCAL	void	get_untracked_shock_state(double*,Locstate,
	                                  COMP_TYPE*,HYPER_SURF*,INTERFACE*,
					  INIT_DATA*,int);
LOCAL	void	free_untracked_shock_wave_comp_type(COMP_TYPE*);

#if defined(TWOD)
LOCAL	_TAYLOR_WAVE	*set_taylor_wave_comp_type(COMP_TYPE*,_TAYLOR_WAVE*,
	                                           Front*);
LOCAL	FOURIER_POLY	*prompt_for_fourier_poly(RECT_GRID*,double);
LOCAL	boolean	idmscs_func(double,double*,POINTER);
LOCAL	double	init_meshkov_front_geometry(FOURIER_POLY**,_TAYLOR_WAVE**,
	                                    int*,char*,int*,double*,COMPONENT*,
	                                    int*,int*,int*,
	                                    double*,double*,Front*);
LOCAL	int	init_trans_node_param_choice(void);
LOCAL	void	adjust_to_off_grid_line(double*,int,double,Front*);
LOCAL	void	free_Prandtl_Meyer_wave_comp_type(COMP_TYPE*);
LOCAL	void	free_nwave_comp_type(COMP_TYPE*);
LOCAL	void	free_taylor_wave_comp_type(COMP_TYPE*);
LOCAL	void	get_state_in_Prandtl_Meyer_wave(double*,Locstate,COMP_TYPE*,
	                                        HYPER_SURF*,INTERFACE*,
						INIT_DATA*,int);
LOCAL	void	get_state_nwave(double*,Locstate,COMP_TYPE*,
	                        HYPER_SURF*,INTERFACE*,INIT_DATA*,int);
LOCAL	void	get_state_taylor_wave(double*,Locstate,COMP_TYPE*,
	                              HYPER_SURF*,INTERFACE*,INIT_DATA*,int);
LOCAL	void	init_centered_Prandtl_Meyer_wave(Locstate,Locstate,double,
	                                         double,double*,int,
	                                         double*,COMP_TYPE*);
LOCAL	void	init_diffraction_node(Gas_param*,Gas_param*,Front*,INIT_DATA*);
LOCAL	INIT_TAYLOR_WAVE_TYPE init_meshkov_second_contact_states(Locstate*,
								 Gas_param**,
	                                                         double,
								 _TAYLOR_WAVE*,
								 Front*);
LOCAL	void	init_meshkov_single_contact_states(Locstate*,Gas_param**,
	                                           double,_TAYLOR_WAVE*,Front*);
LOCAL	void	init_pre_diffraction_node(Gas_param*,Gas_param*,
	                                  Front*,boolean,INIT_DATA*);
LOCAL	void	init_taylor_wave_component(COMP_TYPE*,_TAYLOR_WAVE*,int,double,
	                                   Locstate,Locstate,Locstate,Front*);
LOCAL	void	init_transmission_node(Gas_param*,Gas_param*,RP_DATA**,Front*);
LOCAL	double	meshkov_amb_vert_vel(double,double,double,double,Gas_param*,
	                             Gas_param*);
LOCAL	void	print_TAYLOR_WAVE_structure(_TAYLOR_WAVE*);
LOCAL	void	set_Prandtl_Meyer_wave_comp_type(COMP_TYPE*,Front*);
LOCAL	void	set_nwave_comp_type(COMP_TYPE*,Front*);
#endif /* defined(TWOD) */


EXPORT	void	set_untracked_shock_wave_comp_type(
	COMP_TYPE	*comp_type,
	UT_SHOCK	*utsw,
	Front		*front)
{
	if (comp_type->type == UNTRACKED_SHOCK) /*ALREADY SET*/
	    return;

	comp_type->type = UNTRACKED_SHOCK;
	if (utsw == NULL)
	{
	    scalar(&utsw,sizeof(UT_SHOCK));
	    alloc_state(front->interf,&utsw->state0,front->sizest);
	    alloc_state(front->interf,&utsw->state1,front->sizest);
	    utsw->free_with_comp_type = YES;
	}
	else
	    utsw->free_with_comp_type = NO;
	comp_type->extra = (POINTER)utsw;
	comp_type->_get_state = get_untracked_shock_state;
	comp_type->free_comp_type_extra = free_untracked_shock_wave_comp_type;
}		/*end untracked_shock_wave_comp_type*/

LOCAL	void	free_untracked_shock_wave_comp_type(
	COMP_TYPE	*comp_type)
{
	UT_SHOCK *utsw;

	if (comp_type->type != UNTRACKED_SHOCK)
	    return;

	utsw = (UT_SHOCK*)comp_type->extra;
	if ((utsw == NULL) || (utsw->free_with_comp_type == NO))
	    return;
	free(utsw->state0);
	free(utsw->state1);
	free(utsw);
	comp_type->extra = NULL;
}		/*end free_untracked_shock_wave_comp_type*/


/*
*		get_untracked_shock_state():
*
*	Initializes the states in a region containing an untracked shock.
*/

LOCAL	void 	get_untracked_shock_state(
	double		*coords,
	Locstate	state,
	COMP_TYPE	*comp_type,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		stype)
{
	UT_SHOCK 	*utsw = Untracked_shock(comp_type);
	Locstate	state0 = utsw->state0;
	Locstate	state1 = utsw->state1;
	double		coords0[MAXD], coords1[MAXD];
	double		dposn[MAXD];
	double		*posn = utsw->posn;
	double		*nor = utsw->nor;
	double		width = utsw->width;
	double		sp, alpha;
	int		i, dim = Params(state0)->dim;

	if ( comp_type->type != UNTRACKED_SHOCK )
	{
	    screen("ERROR: in get_untracked_shock_state(), "
	           "inconsistent comp_type->type\n");
	    clean_up(ERROR);
	}
	for (i = 0; i < dim; i++)
	{
	    coords0[i] = coords[i] + 0.5*width*nor[i];
	    coords1[i] = coords[i] - 0.5*width*nor[i];
	}

	for (i = 0; i < dim; i++)
	    dposn[i] = coords[i] - posn[i];
	sp = scalar_product(dposn,nor,dim);
	alpha = 0.5 - sp/width;
	if (alpha < 0.0)
	{
	    if (utsw->ctype0 != NULL)
	    	Get_state(coords0,state0,utsw->ctype0,hs,intfc,init,stype);
	    copy_state(state,state0);
	}
	else if (alpha > 1.0)
	{
	    if (utsw->ctype1 != NULL)
	    	Get_state(coords1,state1,utsw->ctype1,hs,intfc,init,stype);
	    copy_state(state,state1);
	}
	else
	{
	    if (utsw->ctype0 != NULL)
	    	Get_state(coords0,state0,utsw->ctype0,hs,intfc,init,stype);
	    if (utsw->ctype1 != NULL)
	    	Get_state(coords1,state1,utsw->ctype1,hs,intfc,init,stype);
	    bi_interpolate_intfc_states(current_interface(),(1.0-alpha),
	                                alpha,coords,state0,coords,
	                                state1,state);
	}
}		/*end get_untracked_shock_state*/


#if defined(TWOD)

/*
*			init_meshkov():
*
*	Initializes several versions or the Richtmyer-Meshkov
*	instability.  The most general case consists of seven
*	states separated by various waves.
*
*	st[0] - State on downstream side of the target contact
*	st[1] - State on upstream side of the target contact
*	st[2] - State immediately behind the incident shock
*	st[3] - State on downstream side of the second contact
*	st[4] - State on upstream side of the second contact
*	st[5] - State at the leading edge of the Taylor wave
*	st[6] - State at the trailing edge of the Taylor wave
*/

EXPORT	void init_meshkov(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front		*front = ip->root->front;
	INTERFACE	*intfc = front->interf;
	FOURIER_POLY	*fpoly;
	NODE		*ns, *ne;
	CURVE		*cur;
	Locstate	st[7];
	Gas_param	*params[4];
	COMPONENT	comp[4];
	COMPONENT	lcomp, rcomp;
	_TAYLOR_WAVE	*twave = NULL, *twave1 = NULL;
	double		XL = front->rect_grid->L[0];
	double		XU = front->rect_grid->U[0];
	double		y_shock;
	double		coords[MAXD];
	double		sign;
	double		st1, st2;
	int		num_points;
	int		read_contact_from_file;
	int		track[7];
	int		second_contact = NO;
	int		tbar_given = NO;
	int		i;
	INIT_TAYLOR_WAVE_TYPE	twavests;
	size_t		sizest = front->sizest;
	char		s[Gets_BUF_SIZE], contact_file[Gets_BUF_SIZE];

	debug_print("init_meshkov","Entered init_meshkov()\n");
	comp[0] = FIRST_DYNAMIC_COMPONENT;
	comp[1] = comp[2] = comp[3] = comp[0] + 1;
	params[0] = params[1] = params[2] = params[3] = NULL;
	track[0] = track[6] = YES;
	for (i = 0; i < 7; i++)
	    alloc_state(front->interf,&st[i],sizest);

	screen("\nInitialization of Richtmyer-Meshkov Interaction\n");

	st1 = st2 = 0.0;
	sign = init_meshkov_front_geometry(&fpoly,&twave,
	                                   &read_contact_from_file,
	                                   contact_file,&num_points,&y_shock,
	                                   comp,&track[0],&tbar_given,
	                                   &second_contact,
	                                   &st1,&st2,front);

	/* Initialize Gas Params */

	screen("\nInitialize Equation of State Information\n");
	params[1] = init_eos_params(init,ip,
	                            (sign==1.0) ? " above the target contact" :
	                                          " below the target contact",
	                            YES);
	params[2] = params[3] = params[1];
	screen("Do you wish to enter another equation of state\n\t");
	screen("for the material on the other ");
	screen("side of the target contact [y]: ");
	(void) Gets(s);
	if (s[0] == 'n' || s[0] == 'N')
	    params[0] = params[1];
	else
	{
	    params[0] = init_eos_params(init,ip,
	                                (sign==1.0) ?
	                                        " below the target contact" :
	                                        " above the target contact",
	                                YES);
	}
	if (second_contact)
	{
	    screen("Do you wish to enter a different ");
	    screen("equation of state for ");
	    screen("the material in the expanding region ");
	    screen("inside the rarefaction (yes(dflt) or no): ");
	    (void) Gets(s);
	    if (s[0] != 'n' && s[0] != 'N')
	    {
	    	params[3] = init_eos_params(init,ip," in the Taylor wave",YES);
	    }
	}

	/* Initialize State Information */

	if (second_contact)
	{
	    twavests = init_meshkov_second_contact_states(st,params,
	                                                  sign,twave,front);
	}
	else
	{
	    twavests = DIRECT;
	    init_meshkov_single_contact_states(st,params,sign,twave,front);
	}
	if (debugging("init_meshkov"))
	{
	    (void) printf("Initialized states for Meshkov Instability\n");
	    for (i = 0; i < 7; i++)
	    {
	    	(void) printf("\n");
	    	(void) sprintf(s,"st[%d]",i);
	    	verbose_print_state(s,st[i]);
	    }
	    (void) printf("\n");
	}

	/* Set comp_type information */

	/* comp[0] */
	set_ambient_comp_type(comp_type(comp[0]),front);
	set_state(Ambient(comp_type(comp[0])),GAS_STATE,st[0]);

	screen("Type 'y' if you wish to track the reflected wave");
	screen(" (shock or rarefaction): ");
	(void) Gets(s);
	track[1] = track[2] = track[3] =
	    ((s[0] == 'y') || (s[0] == 'Y')) ? YES : NO;

	screen("Type 'y' if you wish to track the transmitted shock: ");
	(void) Gets(s);
	track[5] = ((s[0] == 'y') || (s[0] == 'Y')) ? YES : NO;

	if (track[0])
	{
	    /* comp[1] */
	    if (track[5])
	    {
	    	screen("Is the flow on the opposite side "
	    	       "of the target contact\n\tfrom the "
	    	       "incident shock constant? (dflt = no): ");
	    	(void) Gets(s);
	    	if (s[0] == 'y' || s[0] == 'Y')
	    	    (void)SetConstantFlowRegion(comp[0],
	                                        Ambient(comp_type(comp[0])),
	                                        intfc);
	    }
	    set_ambient_comp_type(comp_type(comp[1]),front);
	    set_state(Ambient(comp_type(comp[1])),GAS_STATE,st[1]);
	    screen("Is the flow on the side of the target contact\n\t"
	           "adjacent to the incident shock constant? (dflt = no): ");
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y')
	    	(void)SetConstantFlowRegion(comp[1],Ambient(comp_type(comp[1])),
	                                    intfc);

	    /* comp[2] */
	    if (twave == NULL)
	    {
	        set_ambient_comp_type(comp_type(comp[2]),front);
	    	set_state(Ambient(comp_type(comp[2])),GAS_STATE,st[2]);
	    	if (track[1] && track[2] && track[3])
	    	{
	    	    screen("Is the flow behind the incident "
	    	           "shock constant? (dflt = no): ");
	    	    (void) Gets(s);
	    	    if (s[0] == 'y' || s[0] == 'Y')
	    	    	(void)SetConstantFlowRegion(comp[2],
	    	    				    Ambient(comp_type(comp[2])),
	                                            intfc);
	    	}
	    }
	    else if (second_contact == NO)
	    {
	    	init_taylor_wave_component(comp_type(comp[2]),
	                                   twave,tbar_given,sign,st[1],
	                                   st[4],st[6],front);
	    }
	    else
	    {
		switch (twavests)
	    	{
		case DIRECT:
	    	    scalar(&twave1,sizeof(_TAYLOR_WAVE));
	    	    twave1->z0 = twave1->zs = y_shock;
	    	    twave1->z1 = twave->z0;
	    	    init_taylor_wave_component(comp_type(comp[2]),twave1,NO,
	                                       sign,st[1],st[2],st[3],front);
		    break;
		case INDIRECT:
	            set_ambient_comp_type(comp_type(comp[2]),front);
	    	    set_state(Ambient(comp_type(comp[2])),GAS_STATE,st[2]);
		    break;
	    	}

	    	twave->zs = twave->z0;
	    	init_taylor_wave_component(comp_type(comp[3]),twave,tbar_given,
	                                   sign,st[4],st[5],st[6],front);

	    }
	}
	else if (twave != NULL)
	{
	    if (second_contact == NO)
	    {
	    	init_taylor_wave_component(comp_type(comp[1]),twave,tbar_given,
	                                   sign,st[1],st[4],st[6],front);
	    }
	    else
	    {
	    	switch (twavests)
	    	{
		case DIRECT:
	    	    scalar(&twave1,sizeof(_TAYLOR_WAVE));
	    	    twave1->z0 = twave1->zs = y_shock;
	    	    twave1->z1 = twave->z0;
	    	    init_taylor_wave_component(comp_type(comp[1]),twave1,NO,
	                                       sign,st[1],st[2],st[3],front);
		    break;

		case INDIRECT:
	            set_ambient_comp_type(comp_type(comp[1]),front);
	    	    set_state(Ambient(comp_type(comp[1])),GAS_STATE,st[2]);
		    break;
	    	}

	    	twave->zs = twave->z0;
	    	init_taylor_wave_component(comp_type(comp[3]),twave,tbar_given,
	                                   sign,st[4],st[5],st[6],front);

	    }
	}
	else
	{
	    UT_SHOCK *utsw;

	    set_untracked_shock_wave_comp_type(comp_type(comp[1]),NULL,front);
	    utsw = (UT_SHOCK*)comp_type(comp[1])->extra;
	    set_state(utsw->state0,GAS_STATE,st[1]);
	    set_state(utsw->state1,GAS_STATE,st[2]);
	    utsw->nor[0] = 0.0;
	    utsw->nor[1] = (sign > 0.0) ? -1.0 : 1.0;
	    utsw->posn[1] = y_shock;
	    utsw->posn[0] = 0.5*(XU + XL);
	    utsw->width = 10.0*front->rect_grid->h[1];
	    utsw->ctype0 = utsw->ctype1 = NULL;
	    utsw->_wave_type = (sign > 0.0) ?
	    			BACKWARD_SHOCK_WAVE : FORWARD_SHOCK_WAVE;
	}


	/* Make target contact */

	if (sign > 0.0)
	{
	    lcomp = comp[0];
	    rcomp = comp[1];
	}
	else
	{
	    lcomp = comp[1];
	    rcomp = comp[0];
	}
	if (read_contact_from_file)
	{
	    coords[0] = 0.0;	coords[1] = 0.0;
	    ns = make_node(Point(coords));
	    ne = make_node(Point(coords));
	    cur = read_curve_from_file(lcomp,rcomp,ns,ne,contact_file);
	    wave_type(cur) = CONTACT;
	    surface_tension(cur) = st1;
	    start_status(cur) = end_status(cur) = INCIDENT;
	}
	else
	{
	    g_make_fourier_curve(CONTACT,num_points,XL,XU,
	    		         fpoly,lcomp,rcomp,st1);
	    free(fpoly);
	    fpoly = NULL;
	}
	
	/* Make incident shock */

	if (track[0])
	{
	    if (sign > 0.0)
	    {
	    	lcomp = comp[1];
	    	rcomp = comp[2];
	    }
	    else
	    {
	    	lcomp = comp[2];
	    	rcomp = comp[1];
	    }
	    coords[0] = XU;	coords[1] = y_shock;
	    ns = make_node(Point(coords));
	    coords[0] = XL;	coords[1] = y_shock;
	    ne = make_node(Point(coords));
	    cur = make_curve(lcomp,rcomp,ns,ne);
	    wave_type(cur) = (sign > 0.0) ? BACKWARD_SHOCK_WAVE :
	                                    FORWARD_SHOCK_WAVE;
	    start_status(cur) = end_status(cur) = INCIDENT;
	}

	/* Make second contact */

	if (second_contact)
	{
	    if (sign > 0.0)
	    {
	    	lcomp = comp[2];
	    	rcomp = comp[3];
	    }
	    else
	    {
	    	lcomp = comp[3];
	    	rcomp = comp[2];
	    }
	    coords[0] = XU;	coords[1] = twave->z0;
	    ns = make_node(Point(coords));
	    coords[0] = XL;	coords[1] = twave->z0;
	    ne = make_node(Point(coords));
	    cur = make_curve(lcomp,rcomp,ns,ne);
	    wave_type(cur) = CONTACT;
	    surface_tension(cur) = st2;
	    start_status(cur) = end_status(cur) = INCIDENT;
	}
	for (i = 0; i < 7; i++) free(st[i]);

	if (debugging("init_meshkov"))
	{
	    (void) printf("Initialized states by component\n");
	    for (i = 0; i < 4; i++)
	    {
	    	(void) printf("\n");
	    	(void) sprintf(s,"comp[%d] = %d",i,comp[i]);
	    	verbose_print_state(s,Ambient(comp_type(comp[i])));
	    }
	    (void) printf("\n");
	}
	debug_print("init_meshkov","Left init_meshkov()\n");
}		/*end init_meshkov*/


EXPORT	void init_shock_diffraction(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front		*front = ip->root->front;
	INTERFACE	*intfc = front->interf;
	char		s[Gets_BUF_SIZE];
	Gas_param	*params, *params2;
	boolean		constant_ahead_of_incident;
	boolean		constant_ahead_of_transmitted;

	/* Input gas parameters */

	params = init_eos_params(init,ip," with the incident shock ",YES);
	params2 = init_eos_params(init,ip," with the transmitted shock ",YES);

	screen("Is the flow ahead of the incident shock "
	       "constant? (dflt = no): ");
	(void) Gets(s);
	constant_ahead_of_incident = (s[0] == 'y' || s[0] == 'Y') ? YES : NO;
	screen("Is the flow ahead of the transmitted shock "
	       "constant? (dflt = no): ");
	(void) Gets(s);
	constant_ahead_of_transmitted = (s[0] == 'y' || s[0] == 'Y') ? YES : NO;

	screen("\nThere are two versions of the initialization "
	       "for the shock contact interaction.\n\t"
	       "Version 1 assumes that the interaction has already "
	       "occurred and\n\t"
	       "initializes a diffraction node "
	       "configuration.  Version 2 assumes that\n\t"
	       "no collision "
	       "between the shock and contact has occurred yet and\n\t"
	       "initializes the two curves prior to "
	       "interaction.\n"
	       "\tVersion 2 allows initialization of a tail shock.\n"
	       "Type 2 for version 2, otherwise version 1 will be assumed: "
	       "will be assumed: ");
	(void) Gets(s);
	if (s[0] != '2') 
	    init_diffraction_node(params,params2,front,init);
	else 
	{
	    screen("Enter 'y' to initialize a trailing shock: ");
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y')
	        init_pre_diffraction_node(params,params2,front,YES,init);
	    else
	        init_pre_diffraction_node(params,params2,front,NO,init);
	}
	if (constant_ahead_of_incident == YES)
	   (void)SetConstantFlowRegion(COMPIF,Ambient(comp_type(COMPIF)),intfc);
	if (constant_ahead_of_transmitted == YES)
	   (void)SetConstantFlowRegion(COMPTF,Ambient(comp_type(COMPTF)),intfc);

}		/*end init_shock_diffraction*/

EXPORT	void init_shock_transmission(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front		*front = ip->root->front;
	double		*L = front->rect_grid->L, *U = front->rect_grid->U;
	double		*h = front->rect_grid->h;
	double		coords[MAXD];
	double		surf_ten;
	Gas_param	*params, *params2;
	NODE		*tran_node;
	NODE		*ns, *ne;
	CURVE		*cur;
	RP_DATA		*RP = NULL;

	/* Input gas parameters */
	params = init_eos_params(init,ip," for the gas above",YES);
	params2 = init_eos_params(init,ip," for the gas below",YES);

	set_ambient_comp_type(comp_type(COMPAR),front);
	set_ambient_comp_type(comp_type(COMPBR),front);
	set_ambient_comp_type(comp_type(COMPAL),front);
	set_ambient_comp_type(comp_type(COMPBL),front);

	init_transmission_node(params,params2,&RP,front);
	
	surf_ten = prompt_for_surface_tension(CONTACT,"for the contact ");

	copy_state(Ambient(comp_type(COMPAR)),RP->state[0]);
	copy_state(Ambient(comp_type(COMPAL)),RP->state[1]);
	copy_state(Ambient(comp_type(COMPBL)),RP->state[2]);
	copy_state(Ambient(comp_type(COMPBR)),RP->state[3]);


	coords[0] = L[0]+.25*(U[0]-L[0]);	coords[1] = .5*(L[1]+U[1]);
	ns = tran_node = make_node(Point(coords));
	copy_RP_DATA_structure(Rp_data(tran_node),RP);
	coords[0] = U[0];	coords[1] = .5*(L[1]+U[1]);
	ne = make_node(Point(coords));
	node_type(ns) = TRANSMISSION_NODE;
	cur = make_curve(COMPAR,COMPBR,ns,ne);
	wave_type(cur) = CONTACT;
	start_status(cur) = CONTACT_TARGET;
	end_status(cur) = INCIDENT;

	make_rotated_curve(h,RP->ang[4],L,U,NULL,NULL,NULL,tran_node,
	                   COMPBL,COMPBR,&cur,NEGATIVE_ORIENTATION,
	                   FORWARD_SHOCK_WAVE,TRANSMITTED,INCIDENT,0.0,front);

	make_rotated_curve(h,RP->ang[1],L,U,NULL,NULL,NULL,tran_node,
	                   COMPAL,COMPAR,&cur,POSITIVE_ORIENTATION,
	                   FORWARD_SHOCK_WAVE,INCIDENT,INCIDENT,0.0,front);

	make_rotated_curve(h,RP->ang[3],L,U,NULL,NULL,NULL,tran_node,
	                   COMPAL,COMPBL,&cur,NEGATIVE_ORIENTATION,
	                   CONTACT,SLIP,INCIDENT,surf_ten,front);
}		/*end init_shock_trasmission*/

LOCAL	double init_meshkov_front_geometry(
	FOURIER_POLY	**fpoly,
	_TAYLOR_WAVE	**twave,
	int		*read_contact_from_file,
	char		*contact_file,
	int		*num_points,
	double		*y_shock,
	COMPONENT	*comp,
	int		*track_incident_shock,
	int		*tbar_given,
	int		*second_contact,
	double		*st1,
	double		*st2,
	Front		*front)
{
	char		s[Gets_BUF_SIZE];
	int		dim = front->rect_grid->dim;
	double		*L = front->rect_grid->L;
	double		*U = front->rect_grid->U;
	double		*h = front->rect_grid->h;
	double		sign;

	*tbar_given = NO;
	*second_contact = NO;
	screen("Shock is incident from above (a, default) or below (b): ");
	(void) Gets(s);
	sign = (s[0] == 'b' || s[0] == 'B') ? -1.0 : 1.0;
	screen("\nThe initial shock is incident from %s.  ",
	        (sign == 1.0) ? "above" : "below");
	screen("The shock wave and\n\tcontact discontinuity curves "
	       "can be initialized in one\n\tof two different ways.  "
	       "A flat shock colliding with a sinusoidal\n\t"
	       "contact (choice a or default), or a flat shock "
	       "colliding with an \n\tarbitrarily shaped contact, "
	       "whose coordinates are read from some\n\tfile ");
	screen("(choice b).\n");
	screen("Enter choice (a or b): ");
	(void) Gets(s);
	if (s[0] == 'b' || s[0] == 'B')
	{
	    *read_contact_from_file = YES;
	    screen("Enter the file name containing the\n\t");
	    screen("coordinates of the contact (stdin default): ");
	    (void) Gets(contact_file);
	    if (contact_file[0] == '\0')
	    {
	    	(void) strcpy(contact_file,"stdin");
	    	screen("You will be prompted later for "
	    	       "the contact discontinuity coordinates\n");
	    }
	}
	else
	{
	    *read_contact_from_file = NO;
	    *fpoly = prompt_for_fourier_poly(front->rect_grid,sign);
	    *num_points = 2*front->rect_grid->gmax[0];
	}
	*st1 = prompt_for_surface_tension(CONTACT,"for the target contact ");


	screen("Enter the y coordinate of the shock: ");
	(void) Scanf("%f\n",y_shock);
	adjust_to_off_grid_line(y_shock,1,-sign,front);
	screen("\nDo you wish to track the incident shock?[y(dflt), n]: ");
	(void) Gets(s);
	*track_incident_shock = YES;
	if (s[0] == 'n' || s[0] == 'N') *track_incident_shock = NO;
	else				comp[2] = comp[3] = comp[1] + 1;

	*twave = NULL;
	screen("Do you wish to initialize a Taylor wave behind ");
	screen("the incident shock?(dflt=n): ");
	(void) Gets(s);
	if (s[0] == 'y' || s[0] == 'Y')
	{
	    scalar(twave,sizeof(_TAYLOR_WAVE));
	    (*twave)->zs = *y_shock;
	    screen("The Taylor wave can be initialized in one of ");
	    screen("two ways,  an untracked\n\t");
	    screen("rarefaction wave overtaking the incident shock ");
	    screen("(choice a, default),  or\n\tan expanding ");
	    screen("material behind the incident shock (choice b).\n");
	    screen("Enter choice (a or b): ");
	    (void) Gets(s);
	    if (s[0] == 'b' || s[0] == 'B')
	    {
	        *second_contact = YES;
	        *st2 = prompt_for_surface_tension(CONTACT,
	                                         "for the additional contact ");
	        comp[3] = comp[2] + 1;
	        screen("Enter the distance from the leading edge of ");
	        screen("the rarefaction wave\n\t");
	        screen("to the shock (default = 0.5*h[%d]): ",dim-1);
	        (void) Gets(s);
	        if (s[0] == '\0')
	            (*twave)->z0 = *y_shock + sign*0.5*h[dim-1];
	        else
	        {
	            (void) sscan_float(s,&(*twave)->z0);
	            (*twave)->z0 += sign*(*y_shock);
	        }
	        adjust_to_off_grid_line(&(*twave)->z0,1,-sign,front);
	        screen("The trailing edge position can be entered "
	               "in one of two ways,\n\t");
	        screen("direct entry of the distance of the trailing "
	               "edge from\n\t");
	        screen("the %s boundary, (choice a and default), ",
	               (sign == 1.0) ? "upper" : "lower");
	        screen("or by entering\n\t");
	        screen("the time elapsed since the formation of ");
	        screen("the rarefaction\n\twave (choice b).");
	        screen("Enter choice (a or b): ");
	        (void) Gets(s);
	        if (s[0] == 'b' || s[0] == 'B')
	        {
	            *tbar_given = YES;
	            screen("Enter the time elapsed since the formation "
	                   "of the rarefaction wave: ");
	            (void) Scanf("%f\n",&(*twave)->tbar);
	        }
	        else
	        {
	            screen("Enter the distance from the trailing edge of "
	                   "the rarefaction wave\n\t");
	            screen("to the %s boundary (default = 0): ",
	                   (sign == 1.0) ? "upper" : "lower");
	            (void) Gets(s);
	            if (s[0] == '\0')
	            {
	                (*twave)->z1 = (sign == 1.0) ? U[dim-1] : L[dim-1];
	            }
	            else
	            {
	                (void) sscan_float(s,&(*twave)->z1);
	                (*twave)->z1 = (sign == 1.0) ? U[dim-1] - (*twave)->z1:
	                        	               L[dim-1] + (*twave)->z1;
	            }
	            if (sign*(*twave)->z0 > sign*(*twave)->z1)
	            {
	                screen("ERROR in init_meshkov_front_geometry(), "
	                       "Invalid rarefaction edges\n");
	                screen("Low pressure edge adjacent to shock\n");
	                clean_up(ERROR);
	            }
	        }
	    }
	    else
	    {
	        screen("Enter the distance from the leading edge of ");
	        screen("the rarefaction wave\n\t");
	        screen("to the shock (default = 0): ");
	        (void) Gets(s);
	        if (s[0] == '\0')
	            (*twave)->z0 = *y_shock;
	        else
	        {
	            (void) sscan_float(s,&(*twave)->z0);
	            (*twave)->z0 += sign*(*y_shock);
	        }
	        screen("Enter the distance from the trailing edge of "
	               "the rarefaction wave\n\t");
	        screen("to the %s boundary (default = 0): ",
	               (sign == 1.0) ? "upper" : "lower");
	        (void) Gets(s);
	        if (s[0] == '\0')
	            (*twave)->z1 = (sign == 1.0) ? U[dim-1] : L[dim-1];
	        else
	        {
	            (void) sscan_float(s,&(*twave)->z1);
	            (*twave)->z1 = (sign == 1.0) ? U[dim-1] - (*twave)->z1:
	                                           L[dim-1] + (*twave)->z1;
	        }
	        if (sign*(*twave)->z0 > sign*(*twave)->z1)
	        {
	            screen("ERROR in init_meshkov_front_geometry(), "
	                   "Invalid rarefaction edges\n");
	            screen("Low pressure edge adjacent to shock\n");
	            clean_up(ERROR);
	        }
	    }
	}
	return sign;
}		/*end init_meshkov_front_geometry*/

LOCAL	FOURIER_POLY *prompt_for_fourier_poly(
	RECT_GRID	*rgr,
	double		sign)
{
	FOURIER_POLY	*fpoly;
	char		s[Gets_BUF_SIZE];
	int		i, j;
	int		dim = rgr->dim;
	int		nmodes;
	double		*L = rgr->L;
	double		*U = rgr->U;
	double		nu[MAXD-1];
	double		z0;

	screen("Enter the mean position of front above L[%d]: ",dim-1);
	(void) Scanf("%f\n",&z0);
	z0 += L[dim-1];
	screen("Enter the number of Fourier modes ");
	screen("in the contact (default = 1): ");
	(void) Gets(s);
	if (s[0] != '\0')	(void) sscanf(s,"%d",&nmodes);
	else			nmodes = 1;
	fpoly = allocate_fourier_poly(nmodes,dim,NULL);
	fpoly->z0 = z0;
	fpoly->num_modes = nmodes;
	fpoly->dim = rgr->dim;

	for (i = 0; i < nmodes; i++)
	{
	    screen("\tEnter amplitude %d: ",i);
	    (void) Scanf("%f\n",&fpoly->A[i]);
	    for (j = 0; j < dim-1; j++)
	    {
	    	screen("\tEnter frequency %d ",i);
	    	screen("in coordinate direction %d\n\t\t", j);
	    	screen("(integer for periodic): ");
	    	(void) Scanf("%f\n",&nu[j]);
	    }
	    fpoly->phase[i] = (sign == 1.0) ? 1.5*PI : 0.5*PI;
	    screen("\tEnter phase %d in degrees ",i);
	    screen("(default = %s gives cosine curve): ",
	                           (sign == 1.0) ? "270" : "90");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	    	(void) sscan_float(s,&fpoly->phase[i]);
	    	fpoly->phase[i] = radians(fpoly->phase[i]);
	    }
	    for (j = 0; j < dim-1; j++)
	    {
	    	fpoly->nu[i][j] = 2.0*PI*nu[j]/(U[j] - L[j]);
	    	fpoly->phase[i] -= 2.0*PI*nu[j]*L[j]/(U[j] - L[j]);
	    }
	}
	return fpoly;
}		/*end prompt_for_fourier_poly*/

typedef struct {
	Locstate *st;
	double pstar, ustar, p6, v6;
	double sign;
} IDMSCS;

LOCAL	INIT_TAYLOR_WAVE_TYPE init_meshkov_second_contact_states(
	Locstate	*st,
	Gas_param	**params,
	double		sign,
	_TAYLOR_WAVE	*twave,
	Front		*front)
{
	char                     s[Gets_BUF_SIZE];
	double                    v[2];
	double                    p1, r1, r2;
	INIT_TAYLOR_WAVE_TYPE    status;
	RIEMANN_SOLVER_WAVE_TYPE l_wave, r_wave;
	int                      i;

	screen("Current reference frame assumes target contact is at rest\n"
	       "X component of velocity will be set to zero\n");

	/* Thermodynamics of State on upstream side of target contact */

	(void) sprintf(s," %s the target contact",
	               (sign == 1.0) ? "below" : "above");
	prompt_for_thermodynamics(st[0],params[0],s);

	/* Thermodynamics of State on downstream side of target contact */

	set_type_of_state(st[1],TGAS_STATE);
	Init_params(st[1],params[1]);
	screen("Enter the density %s the target contact: ",
	        (sign == 1.0) ? "above" : "below");
	(void) Scanf("%f\n",&Dens(st[1]));
	Press(st[1]) = pressure(st[0]);
        reset_gamma(st[1]);
	zero_state_velocity(st[1],params[1]->dim);

	screen("The states behind the incident shock ");
	screen("and in the Taylor wave can be entered\n");
	screen("\tin two different ways.  The direct ");
	screen("input [d, or default] gives the\n");
	screen("\tstate of the flow behind the ");
	screen("incident shock in the material\n");
	screen("\t%s the second contact just prior ",
	    (sign == 1.0) ? "above" : "below");
	screen("to the collision of that wave\n");
	screen("\twith the second contact.  A ");
	screen("Riemann problem is solved for the\n");
	screen("\tinteraction between this state ");
	screen("and the state on the %s side\n",
	    (sign == 1.0) ? "upper" : "lower");
	screen("\tof the target contact.  The ");
	screen("solution to this Riemann problem gives\n");
	screen("\tthe scattered wave information ");
	screen("and hence the strength of the\n");
	screen("\tincident shock after it is ");
	screen("diffracted by the second contact.\n");
	screen("\tSubsequently you will enter the ");
	screen("pressure behind the Taylor wave\n");
	screen("\tto complete the initialization.\n");
	screen("\tThe second method of initialization ");
	screen("[i] is an indirect input of\n");
	screen("\tdata.  In this case you enter ");
	screen("the pressure behind the incident\n");
	screen("\tshock after it has been scattered ");
	screen("by the second contact.\n");
	screen("\tYou then enter the full state ");
	screen("of the flow behind the Taylor\n");
	screen("\twave (ie the state at the bottom ");
	screen("of the computational rectangle).\n");
	screen("\tAn inverse Riemann problem is ");
	screen("solved with this data to give\n");
	screen("\tthe configuration about the second contact.\n");
	screen("Enter choice [d or i]: ");
	(void) Gets(s);
	if (s[0] == 'i' || s[0] == 'I')
	{
	    IDMSCS    Idmscs;
	    double    v6, p5, p6, pu;
	    double    vdiff, vell, velu;
	    double    pstar, ustar;
	    int i;

	    status = INDIRECT;
	    screen("Enter the pressure behind the incident shock after\n");
	    screen("\tit has diffracted through the second contact: ");
	    (void) Scanf("%f\n",&pstar);
	    state_w_pr_on_Hugoniot(st[1],pstar,st[2],TGAS_STATE);
	    p1 = pressure(st[1]);
	    r1 = Dens(st[1]);
	    r2 = Dens(st[2]);
	    Vel(st[2])[0] = 0.0;
	    Vel(st[2])[1] = ustar = vel(1,st[1]) -
	            sign*sqrt((pstar - p1)*(1.0/r1 - 1.0/r2));
	    ft_assign(st[3],st[2],front->sizest);
	    set_type_of_state(st[4],TGAS_STATE);
	    set_type_of_state(st[5],TGAS_STATE);
	    Vel(st[4])[0] = Vel(st[5])[0] = 0.0;

	    screen("There are two ways to enter the state at the ");
	    screen("trailing edge of the rarefaction.\n");
	    screen("\tYou way either enter this state directly [direct], ");
	    screen("or indirectly\n");
	    screen("\tby entering the CJ state at the leading edge of ");
	    screen("the rarefaction\n");
	    screen("\tand the pressure at the trailing edge.\n");
	    screen("Enter choice [d or i(default)]: ");
	    (void) Gets(s);
	    if (s[0] == 'd' || s[0] == 'D')
	    {
	        prompt_for_state(st[6],TGAS_STATE,params[3],
	            "\n\tat the trailing edge of the rarefaction");
	    }
	    else
	    {
	        double p_trail;
	        Locstate CJ;

	        alloc_state(front->interf,&CJ,front->sizest);
	        (void) sprintf(s,"\n\t%s%s","of the CJ state ",
	                       "at the leading edge of the rarefaction");
	        prompt_for_state(CJ,TGAS_STATE,params[3],s);
	        screen("Enter the pressure at the trailing edge ");
	        screen("of the rarefaction: ");
	        (void) Scanf("%f\n",&p_trail);
	        state_on_adiabat_with_pr(CJ,p_trail,st[6],TGAS_STATE);
	        Vel(st[6])[0] = 0.0;
	        Vel(st[6])[1] = vel(1,CJ) - sign*riemann_wave_curve(CJ,p_trail);
	        free(CJ);
	    }

	    Idmscs.pstar = pstar;
	    Idmscs.ustar = ustar;
	    Idmscs.st = st;
	    Idmscs.sign = sign;
	    Idmscs.p6 = p6 = pressure(st[6]);
	    Idmscs.v6 = v6 = vel(1,st[6]);

	    vdiff = sign*(ustar - v6);
	    (void) idmscs_func(p6,&vell,(POINTER)&Idmscs);
	    pu = max(p6,pstar);
	    (void) idmscs_func(pu,&velu,(POINTER)&Idmscs);
	    for (i = 0; i < 20; i++)
	    {
	        if ((vdiff-vell)*(vdiff-velu) <= 0.0)
	            break;
	        pu *= 2.0;
	        (void) idmscs_func(pu,&velu,(POINTER)&Idmscs);
	    }

	    if (find_root(idmscs_func,(POINTER)&Idmscs,
	              vdiff,&p5,p6,pu,EPS,EPS) == FUNCTION_FAILED)
	    {
	        screen("ERROR in ");
	        screen("init_meshkov_second_contact_states(), ");
	        screen("indirect initialization failed\n");
	        clean_up(ERROR);
	    }
	}
	else
	{
	    double    pstarl, pstarr, ustarl, ustarr, ml, mr;

	    status = DIRECT;

	    prompt_for_state(st[5],TGAS_STATE,params[3],
	            "\n\tat the leading edge of the rarefaction");
	    Vel(st[5])[0] = 0.0;
	    if (sign == 1.0)
	    {
	        (void) find_mid_state(st[1],st[5],0.0,&pstarl,&pstarr,
	                              &ustarl,&ustarr,&ml,&mr,&l_wave,&r_wave);
	        if (l_wave == RAREFACTION)
	        {
	            set_state(st[2],state_type(st[1]),st[1]);
	            state_on_adiabat_with_pr(st[1],pstarl,st[3],TGAS_STATE);
	            Vel(st[3])[1] = ustarl;
	            Vel(st[3])[0] = 0.0;
	        }
	        else
	        {
	            state_w_pr_on_Hugoniot(st[1],pstarl,st[2],TGAS_STATE);
	            Vel(st[2])[1] = ustarl;
	            Vel(st[2])[0] = 0.0;
	            set_state(st[3],state_type(st[2]),st[2]);
	        }
	        if (r_wave == RAREFACTION)
	        {
	            state_on_adiabat_with_pr(st[5],pstarr,st[4],TGAS_STATE);
	            Vel(st[4])[1] = ustarr;
	            Vel(st[4])[0] = 0.0;
	        }
	        else
	        {
	            state_w_pr_on_Hugoniot(st[5],pstarr,st[4],TGAS_STATE);
	            Vel(st[4])[1] = ustarr;
	            Vel(st[4])[0] = 0.0;
	        }
	    }
	    else
	    {
	        (void) find_mid_state(st[5],st[1],0.0,&pstarl,&pstarr,
	                              &ustarl,&ustarr,&ml,&mr,&l_wave,&r_wave);
	        if (r_wave == RAREFACTION)
	        {
	            set_state(st[2],state_type(st[1]),st[1]);
	            state_on_adiabat_with_pr(st[1],pstarr,st[3],TGAS_STATE);
	            Vel(st[3])[1] = ustarr;
	            Vel(st[3])[0] = 0.0;
	        }
	        else
	        {
	            state_w_pr_on_Hugoniot(st[1],pstarr,st[2],TGAS_STATE);
	            Vel(st[2])[1] = ustarr;
	            Vel(st[2])[0] = 0.0;
	            set_state(st[3],state_type(st[2]),st[2]);
	        }
	        if (l_wave == RAREFACTION)
	        {
	            state_on_adiabat_with_pr(st[5],pstarl,st[4],TGAS_STATE);
	            Vel(st[4])[1] = ustarl;
	            Vel(st[4])[0] = 0.0;
	        }
	        else
	        {
	            state_w_pr_on_Hugoniot(st[5],pstarl,st[4],TGAS_STATE);
	            Vel(st[4])[1] = ustarl;
	            Vel(st[4])[0] = 0.0;
	        }
	    }

	    if (twave != NULL)
	    {
	        screen("Enter the pressure at the trailing edge\n\t");
	        screen("of the rarefaction in the Taylor wave: ");
	        (void) Scanf("%f\n",&p1);
	        state_on_adiabat_with_pr(st[5],p1,st[6],TGAS_STATE);
	        Vel(st[6])[0] = 0.0;
	        Vel(st[6])[1] = vel(1,st[5])
	                          - sign*riemann_wave_curve(st[5],p1);
	    }
	    else
	    {
	        set_state(st[6],state_type(st[5]),st[5]);
	    }
	}

	screen("A particular frame of reference may be ");
	screen("specified by entering the vertical\n");
	screen("\tcomponent of the gas between the incident ");
	screen("shock and\n\tcontact discontinuity.  ");
	screen("Otherwise a reference frame will be taken\n");
	screen("\tso that the zero amplitude approximation ");
	screen("gives the contact\n\tdiscontinuity at rest ");
	screen("after the interaction.\n");
	screen("Enter velocity here: ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    (void) sscan_float(s,&v[1]);
	}
	else
	{
	    v[1] = sign*meshkov_amb_vert_vel(Dens(st[1]),Dens(st[0]),
	                                     pressure(st[0]),pressure(st[2]),
	                                     params[1],params[0]);
	}
	v[0] = 0.0;
	for (i = 0; i < 7; i++)
	    add_velocity_to_state(st[i],v);

	return status;

}	    /*end init_meshkov_second_contact_states*/

LOCAL	boolean idmscs_func(
	double		p,
	double		*vel,
	POINTER		parameters)
{
	double		p6 = ((IDMSCS *)parameters)->p6;
	double		v6 = ((IDMSCS *)parameters)->v6;
	double		pstar = ((IDMSCS *)parameters)->pstar;
	double		sign = ((IDMSCS *)parameters)->sign;
	double		v1, v2;
	Locstate	*st = ((IDMSCS *)parameters)->st;

	state_on_adiabat_with_pr(st[6],p,st[5],state_type(st[5]));
	v1 = riemann_wave_curve(st[5],p6);
	Vel(st[5])[1] = v6 + sign*v1;
	v2 = riemann_wave_curve(st[5],pstar);
	if (pstar > p)
	    state_w_pr_on_Hugoniot(st[5],pstar,st[4],state_type(st[4]));
	else
	    state_on_adiabat_with_pr(st[5],pstar,st[4],state_type(st[4]));
	Vel(st[4])[1] = Vel(st[5])[1] + sign*v2;

	*vel = v1 + v2;
	return FUNCTION_SUCCEEDED;
}		/*end idmscs_func*/


LOCAL void init_meshkov_single_contact_states(
	Locstate	*st,
	Gas_param	**params,
	double		sign,
	_TAYLOR_WAVE	*twave,
	Front		*front)
{
	char  s[Gets_BUF_SIZE];
	double p2, v[MAXD], s_n[MAXD];
	double shock_speed;
	int   i;

	/* Thermodynamics of State on upstream side of target contact */

	(void) sprintf(s," %s the target contact",
	               (sign == 1.0) ? "below" : "above");
	prompt_for_thermodynamics(st[0],params[0],s);

	/* Thermodynamics of State on downstream side of target contact */

	set_type_of_state(st[1],TGAS_STATE);
	Init_params(st[1],params[1]);
	screen("Enter the density %s the target contact: ",
	        (sign == 1.0) ? "above" : "below");
	(void) Scanf("%f\n",&Dens(st[1]));
	Press(st[1]) = pressure(st[0]);
        reset_gamma(st[1]);
	zero_state_velocity(st[1],params[1]->dim);

	screen("Enter the pressure behind the incident shock: ");
	(void) Scanf("%f\n",&p2);

	screen("A particular frame of reference may be ");
	screen("specified by entering the vertical\n");
	screen("\tcomponent of the gas between the incident ");
	screen("shock and\n\tcontact discontinuity.  ");
	screen("Otherwise a reference frame will be taken\n");
	screen("\tso that the zero amplitude approximation ");
	screen("gives the contact\n\tdiscontinuity at rest ");
	screen("after the interaction.\n");
	screen("Enter velocity here: ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    (void) sscan_float(s,&v[1]);
	}
	else
	{
	    v[1] =  sign*meshkov_amb_vert_vel(Dens(st[1]),Dens(st[0]),
	                                      pressure(st[0]),p2,params[1],
	                                      params[0]);
	}
	v[0] = 0.0;

	add_velocity_to_state(st[0],v);
	add_velocity_to_state(st[1],v);

	set_state(st[1],GAS_STATE,st[1]);
	s_n[0] = 0.0;	s_n[1] = -sign;
	(void) s_polar_4(BEHIND_PRESSURE,p2,&shock_speed,s_n,st[1],
	                 st[2],GAS_STATE);

	for (i = 3; i < 6; i++)
	{
	    ft_assign(st[i],st[2],front->sizest);
	}

	if (twave != NULL)
	{
	    double p1;

	    screen("Enter the pressure at the trailing edge\n\t");
	    screen("of the rarefaction in the Taylor wave: ");
	    (void) Scanf("%f\n",&p1);
	    state_on_adiabat_with_pr(st[5],p1,st[6],TGAS_STATE);
	    Vel(st[6])[0] = 0.0;
	    Vel(st[6])[1] = vel(1,st[5]) - sign*riemann_wave_curve(st[5],p1);
	}
	else
	{
	    ft_assign(st[6],st[5],front->sizest);
	}
}		/*end init_meshkov_single_contact_states*/



/* 
*			init_diffraction_node():
*
*	This routine initializes the states about the diffraction node.
*/


LOCAL void init_diffraction_node(
	Gas_param	*params,
	Gas_param	*params2,
	Front		*front,
	INIT_DATA	*init)
{
	size_t		sizest = front->sizest;
	RECT_GRID	*rect_grid = front->rect_grid;
	Locstate	st_0, st_1, st_4, st_5, st_6;
	NODE_FLAG	flag;
	double		*L = rect_grid->L, *U = rect_grid->U;
	double		*h = rect_grid->h;
	double		frame_ang, incident_ang;
	double		mf, M0;
	double		surf_ten;
	double		x;
	double		delta_ang = 0.0;
	int		i, dim = params->dim;
	int		w_type;
	boolean		sonic_incident_shock = NO;
	int		is_plus_orientation, track_reflected_wave;
	boolean		is_reflected_shock;
	int		change_node_vel;
	char		s[Gets_BUF_SIZE];
	double		rho0,rho6,pr0,pr1;
	double		v0[MAXD], q0,q6;
	double		c0sq, c1sq, V0, V1;
	double 		slip, abs_v[SMAXD];
	double		**t;
	double		v_fac = 1.0;
	COMPONENT	l_comp, r_comp;
	RP_DATA		*RP;
	CURVE		*cur;
	NODE		*diff_node;
	POINT		*p;
	SIN_SQR_PERT_PARAMS	pert_params;
	static	double	nor[] = {1.0, 0.0, 0.0};
	static const char *fmt = "%lf %lf";

	debug_print("init","Entering init_diffraction_node()\n");

	set_to_next_node_only(flag);

	        /* Allocate storage for states around node */

	set_ambient_comp_type(comp_type(COMPIF),front);
	st_0 = Ambient(comp_type(COMPIF));
	set_ambient_comp_type(comp_type(COMPTB),front);
	st_5 = Ambient(comp_type(COMPTB));
	set_ambient_comp_type(comp_type(COMPTF),front);
	st_6 = Ambient(comp_type(COMPTF));

	RP = allocate_RP_DATA_structure(sizest,YES,GAS_STATE);

	/* Input frame reference angle */

	frame_ang = 0.;
	screen("\nEnter the angle (in degrees) which the initial velocity "
	       "(in the steady frame) of\n\t"
	       "the diffraction node makes with the positive x axis ");
	print_angle("(default equals ",frame_ang,").\n");
	screen("Enter choice here: ");
	(void) Gets(s);
	if (s[0] != '\0') 
	{
	    (void) sscan_float(s,&frame_ang);
	    frame_ang = normalized_angle(radians(frame_ang));
	}
	bi_array(&t,2,SMAXD,FLOAT);
	t[1][0] = cos(frame_ang);
	t[1][1] = sin(frame_ang);

	/* Is a reflected shock to be tracked? */

	screen("\nType n if the reflected wave is not to be tracked: ");
	(void) Gets(s);
	track_reflected_wave = (s[0] == 'n' || s[0] == 'N') ? NO : YES;

	/* Input orientation type */

	screen("\nThere are two possible orientations for the diffraction ");
	screen("node.\n\t");
	screen("Positive orientation, in which the incident shock to ");
	screen("front contact\n\t");
	screen("direction is clockwise, and negative ");
	screen("orientation in which this\n\tdirection is counter clockwise. ");
	screen("Type n for negative orientation\n\t");
	screen("(shocks will be ");
	screen("backward sound waves) otherwise positive orientation\n\t");
	screen("(shocks will be forward sound waves) will be used.\n");
	screen("Enter choice here: ");
	(void) Gets(s);
	is_plus_orientation = (s[0] == 'n' || s[0] == 'N') ? NO : YES;
	RP->ang_dir = (is_plus_orientation) ? COUNTER_CLOCK : CLOCKWISE;

	change_node_vel = NO;
	screen("\nThe possible choices for the reference frame are:\n");
	screen("\t(1, default) Ahead state at rest\n");
	screen("\t(2) Node at rest\n");
	screen("\t(3) Reduce node velocity by fraction\n");
	screen("\t(4) User defined ahead velocity\n");
	screen("Enter choice here: ");
	(void) Gets(s);
	switch (s[0])
	{
	case '4':
	    screen("Enter the x component of velocity (dflt = 0): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscan_float(s,v0);
	    screen("Enter the y component of velocity (dflt = 0): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscan_float(s,v0+1);
	    break;
	case '3':
	    change_node_vel = YES;
	    screen("Enter the velocity change factor (dflt = 1): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscan_float(s,&v_fac);
	    v0[0] = 0.0;
	    v0[1] = 0.0;
	    break;
	case '2':
	    change_node_vel = YES;
	    v0[0] = 0.0;
	    v0[1] = 0.0;
	    break;
	case '1':
	default:
	    v0[0] = 0.0;
	    v0[1] = 0.0;
	    break;
	}


	/* Input densities ahead */

	screen("\nEnter the density rho0 in front of the incident shock "
	       "and the density rho6\n\t"
	       "in front of the transmitted shock.\n");
	screen("Enter densities: ");
	(void) Scanf("%f %f\n",&rho0,&rho6);

	/* Input pressures ahead */

	screen("\nEnter the pressure pr0 in front of the incident shock: ");
	(void) Scanf("%f\n",&pr0);
	Dens(RP->state[0]) = rho0;
	Press(RP->state[0]) = pr0;
	Init_params(RP->state[0],params);
	set_type_of_state(RP->state[0],TGAS_STATE);
	reset_gamma(RP->state[0]);
	zero_state_velocity(RP->state[0],dim);

	prompt_for_behind_shock_state(RP->state[0],RP->state[1],YES,nor,
	                              TGAS_STATE,YES,init);
	pr1 = pressure(RP->state[1]);

	screen("\nGeometric configuration is specified by either\n"
	       "\tthe incident angle (I, default),\n"
	       "\tthe incident Mach number (M),\n"
	       "\tor sonic incident shock (S)\n");
	screen("Enter choice: ");
	(void) Gets(s);
	switch (s[0])
	{
	case 'S':
	case 's':
	    sonic_incident_shock = YES;
	    Dens(RP->state[0]) = rho0;
	    Press(RP->state[0]) = pr0;
	    Init_params(RP->state[0],params);
	    set_type_of_state(RP->state[0],TGAS_STATE);
	    reset_gamma(RP->state[0]);
	    zero_state_velocity(RP->state[0],dim);
	    state_w_pr_on_Hugoniot(RP->state[0],pr1,RP->state[1],TGAS_STATE);
	    c0sq = sound_speed_squared(RP->state[0]);
	    c1sq = sound_speed_squared(RP->state[1]);
	    V0 = 1.0/rho0;
	    V1 = 1.0/Dens(RP->state[1]);
	    M0 = sqrt((c1sq + (pr1 - pr0)*(V0 + V1))/c0sq);
	    x = mass_flux(pr1,RP->state[0]) / acoustic_impedance(RP->state[0]);
	    incident_ang = asin(x/M0);
	    debug_print("init","Incident angle = %g\n",incident_ang);
	    break;

	case 'M':
	case 'm':

	    /* Input incident Mach number */

	    screen("Enter the incident Mach number: ");
	    (void) Scanf("%f\n",&M0);
	    Dens(RP->state[0]) = rho0;
	    Press(RP->state[0]) = pr0;
	    Init_params(RP->state[0],params);
	    set_type_of_state(RP->state[0],TGAS_STATE);
	    reset_gamma(RP->state[0]);
	    zero_state_velocity(RP->state[0],dim);
	    x = mass_flux(pr1,RP->state[0])/acoustic_impedance(RP->state[0]);
	    incident_ang = asin(x/M0);
	    debug_print("init","Incident angle = %g\n",incident_ang);
	    break;

	case 'I':
	case 'i':
	default:
	    /* Input incident angle */

	    screen("\nEnter the incident angle: ");
	    (void) Scanf("%f\n",&incident_ang);
	    incident_ang = radians(incident_ang);
	}

	t[0][0] = cos(frame_ang + incident_ang);
	t[0][1] = sin(frame_ang + incident_ang);

	/* Input the slip */

	screen("\nEnter the slip (dflt = 0): ");
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&slip);

	p = Point(NULL);
	Coords(p)[0] = L[0] + 0.25*(U[0] - L[0]);
	Coords(p)[1] = L[1] + 0.5*(U[1] - L[1]);
	screen("\nTo specify an initial position for the diffraction node,\n"
	       "\tenter the x and y coordinates (dflt = %g %g): ",
	        Coords(p)[0],Coords(p)[1]);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if (sscanf(s,fmt,Coords(p),Coords(p)+1) != 2)
	    {
	    	screen("ERROR in init_diffraction_node(), "
	    	       "Improper input of node position\n");
	    	clean_up(ERROR);
	    }
	}

	pert_params.nu = pert_params.epsilon = 0.;
	screen("\nIn order to have a sinusoidally perturbed"
	       " front contact, enter the frequency\n\t"
	       "nu and the amplitude epsilon of the perturbation: "); 
	(void) Gets(s);
	if (s[0] != '\0') 
	{
	    if (sscanf(s,fmt,&pert_params.nu,&pert_params.epsilon) != 2)
	    {
	    	screen("ERROR in init_diffraction_node(), "
	    	       "Improper input of perturbation data\n");
	    	clean_up(ERROR);
	    }
	}

	surf_ten = prompt_for_surface_tension(CONTACT,"for the contact ");


	/* Set state in front of incident shock */

	Init_params(RP->state[0],params);
	Dens(RP->state[0]) = rho0;
	Vel(RP->state[0])[0] = v0[0];
	Vel(RP->state[0])[1] = v0[1];
	Press(RP->state[0]) = pr0;
	set_type_of_state(RP->state[0],TGAS_STATE);
	set_state(RP->state[0],GAS_STATE,RP->state[0]);


	mf = mass_flux(pr1,RP->state[0]);
	q0 = mf / fabs(rho0*sin(incident_ang));
	debug_print("init","q0 = %g\nIncident_angle = %g (recalculated = %g)\n",
	      q0,degrees(incident_ang),
	      degrees(asin(mf/(rho0*q0))));

	abs_v[0] = v0[0] + q0*t[1][0];
	abs_v[1] = v0[1] + q0*t[1][1];

	q6 = q0 - slip;


	    /* Set state in front of transmitted shock */

	Init_params(RP->state[6],params2);
	Dens(RP->state[6]) = rho6;
	Vel(RP->state[6])[0] = (abs_v[0] - q6*t[1][0]);
	Vel(RP->state[6])[1] = (abs_v[1] - q6*t[1][1]);
	Press(RP->state[6]) = pr0;
	set_type_of_state(RP->state[6],TGAS_STATE);
	set_state(RP->state[6],GAS_STATE,RP->state[6]);

	    /* Set state behind incident shock */

	Init_params(RP->state[1],params);
	Dens(RP->state[1]) = 1.;	/* Artifical values, used to pass */
	Vel(RP->state[1])[0] = 0.;	/* pressure via Gas structure     */
	Vel(RP->state[1])[1] = 0.;
	Press(RP->state[1]) = pr1;
	set_type_of_state(RP->state[1],TGAS_STATE);
	set_state(RP->state[1],GAS_STATE,RP->state[1]);



	switch(is_regular_diffraction_node(Coords(p),abs_v,NULL,t,RP,
	                                   NULL,&is_reflected_shock,
	                                   front,DIFFRACTION_NODE,flag)) 
	{
	case REGULAR_DIFFRACTION:
	    break;


	case ANOMALOUS_REFLECTION:
	case REGULAR_TO_MACH_DIFFRACTION:
	case ERROR_DIFFRACTION:
	default:
	    screen("ERROR in init_diffraction_node(), "
	           "is_regular_diffraction_node() failed\n"
	           "Possible bifurcation,  CODE NEEDED\n");
	    clean_up(ERROR);
	}
	free(t);

	delta_ang = 0.0;
	if (sonic_incident_shock == YES)
	{
	    screen("Enter the perturbation angle of the ahead contact ");
	    screen("(in degrees) (dflt = 0): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscan_float(s,&delta_ang);
	    delta_ang = radians(delta_ang);
	}

	if (change_node_vel == YES)
	{
	    double rv[MAXD];

	    for (i = 0; i < dim; i++)
	        rv[i] = -v_fac*abs_v[i];
	    for (i = 0; i < 7; i++)
	        add_velocity_to_state(RP->state[i],rv);
	    for (i = 0; i < dim; i++)
	        abs_v[i] += rv[i];
	}

	ft_assign(st_0,RP->state[0],sizest);
	ft_assign(st_5,RP->state[5],sizest);
	ft_assign(st_6,RP->state[6],sizest);

	if (track_reflected_wave) 
	{
	    set_ambient_comp_type(comp_type(COMPIB),front);
	    st_1 = Ambient(comp_type(COMPIB));
	    ft_assign(st_1,RP->state[1],sizest);
	    set_ambient_comp_type(comp_type(COMPRB),front);
	    st_4 = Ambient(comp_type(COMPRB));
	    ft_assign(st_4,RP->state[4],sizest);
	    if (!is_reflected_shock)
	    {
	        set_Prandtl_Meyer_wave_comp_type(comp_type(COMPRM),front);
	        init_centered_Prandtl_Meyer_wave(st_1,st_4,RP->ang[1],
	                                         RP->ang[3],Coords(p),
	                                         is_plus_orientation,abs_v,
	                                         comp_type(COMPRM));
	    }
	}
	else
	{
	    if (is_reflected_shock)
	    {
	        UT_SHOCK *utsw;

	       	set_untracked_shock_wave_comp_type(comp_type(COMPIB),
	                                           NULL,front);
	        utsw = (UT_SHOCK*)comp_type(COMPIB)->extra;
	        set_state(utsw->state0,GAS_STATE,RP->state[1]);
	        set_state(utsw->state1,GAS_STATE,RP->state[4]);
	        if (is_plus_orientation)
	        {
	            utsw->nor[0] =  sin(RP->ang[2]);
	            utsw->nor[1] = -cos(RP->ang[2]);
	        }
	        else
	        {
	            utsw->nor[0] = -sin(RP->ang[2]);
	            utsw->nor[1] =  cos(RP->ang[2]);
	        }
	        for (i = 0; i < dim; i++)
	            utsw->posn[i] = Coords(p)[i];
	        utsw->width = 10.0*front->rect_grid->h[1];
	        utsw->_wave_type = (is_plus_orientation) ? FORWARD_SHOCK_WAVE :
	            			                   BACKWARD_SHOCK_WAVE;
	        utsw->ctype0 = utsw->ctype1 = NULL;
	    }
	    else
	    {
	        set_Prandtl_Meyer_wave_comp_type(comp_type(COMPIB),front);
	        init_centered_Prandtl_Meyer_wave(RP->state[1],RP->state[4],
	    		                         RP->ang[1],RP->ang[3],
	                                         Coords(p),is_plus_orientation,
	                                         abs_v,comp_type(COMPIB));
	    }
	}
	RP->ang[6] += delta_ang;

	diff_node = make_node(p);
	node_type(diff_node) = DIFFRACTION_NODE;
	copy_RP_DATA_structure(Rp_data(diff_node),RP);
	Node_vel(diff_node)[0] = abs_v[0];
	Node_vel(diff_node)[1] = abs_v[1];

	        /* Make ahead contact */

	if (is_plus_orientation)
	{
	    l_comp = COMPIF;
	    r_comp = COMPTF;
	}
	else
	{
	    l_comp = COMPTF;
	    r_comp = COMPIF;
	}
	make_rotated_curve(h,RP->ang[6],L,U,sin_sqr_pert,sin_sqr_pert_prime,
	                   (POINTER) &pert_params,diff_node,l_comp,r_comp,
	                   &cur,POSITIVE_ORIENTATION,CONTACT,
	                   CONTACT_TARGET,INCIDENT,
	                   surf_ten,front);



	    /* Make incident shock wave */

	if (is_plus_orientation) 
	{
	    l_comp = COMPIB; r_comp = COMPIF;
	    w_type = FORWARD_SHOCK_WAVE;
	}
	else 
	{
	    l_comp = COMPIF; r_comp = COMPIB;
	    w_type = BACKWARD_SHOCK_WAVE;
	}
	make_rotated_curve(h,RP->ang[0],L,U,NULL,NULL,NULL,diff_node,
	                   l_comp,r_comp,&cur,POSITIVE_ORIENTATION,
	                   w_type,INCIDENT,INCIDENT,0.0,front);

	/* Make transmitted shock wave */

	if (is_plus_orientation) 
	{
	    l_comp = COMPTB; r_comp = COMPTF;
	    w_type = FORWARD_SHOCK_WAVE;
	}
	else 
	{
	    l_comp = COMPTF; r_comp = COMPTB;
	    w_type = BACKWARD_SHOCK_WAVE;
	}
	make_rotated_curve(h,RP->ang[5],L,U,NULL,NULL,NULL,diff_node,
	                   l_comp,r_comp,&cur,NEGATIVE_ORIENTATION,w_type,
	                   TRANSMITTED,INCIDENT,0.0,front);

	if (track_reflected_wave) 
	{
	    /* Make deflected contact */

	    if (is_plus_orientation) 
	    {
	    	l_comp = COMPRB;
	    	r_comp = COMPTB;
	    }
	    else
	    {
	    	l_comp = COMPTB;
	    	r_comp = COMPRB;
	    }
	    make_rotated_curve(h,RP->ang[4],L,U,NULL,NULL,NULL,diff_node,
	                       l_comp,r_comp,&cur,NEGATIVE_ORIENTATION,
	                       CONTACT,SLIP,INCIDENT,surf_ten,front);

	    if (is_reflected_shock)
	    {
	    	/* Make reflected shock */

	    	if (is_plus_orientation) 
	    	{
	    	    l_comp = COMPRB; r_comp = COMPIB;
	            w_type = FORWARD_SHOCK_WAVE;
	        }
	        else 
	        {
	            l_comp = COMPIB; r_comp = COMPRB;
	            w_type = BACKWARD_SHOCK_WAVE;
	        }
	        make_rotated_curve(h,RP->ang[2],L,U,NULL,NULL,NULL,diff_node,
	                           l_comp,r_comp,&cur,POSITIVE_ORIENTATION,
	                           w_type,REFLECTED,INCIDENT,0.0,front);
	    }
	    else
	    {
	    	if (is_plus_orientation) 
	    	{
	    	    l_comp = COMPRM; r_comp = COMPIB;
	            w_type = FORWARD_SOUND_WAVE_LE;
	        }
	        else 
	        {
	            l_comp = COMPIB; r_comp = COMPRM;
	            w_type = BACKWARD_SOUND_WAVE_LE;
	        }
	        make_rotated_curve(h,RP->ang[1],L,U,NULL,NULL,NULL,diff_node,
	                           l_comp,r_comp,&cur,POSITIVE_ORIENTATION,
	                           w_type,REFLECTED,INCIDENT,0.0,front);

	        /* Make reflected rarefaction trailing edge */

	        if (is_plus_orientation) 
	        {
	            l_comp = COMPRB; r_comp = COMPRM;
	            w_type = FORWARD_SOUND_WAVE_TE;
	        }
	        else 
	        {
	            l_comp = COMPRM; r_comp = COMPRB;
	            w_type = BACKWARD_SOUND_WAVE_TE;
	        }
	        make_rotated_curve(h,RP->ang[3],L,U,NULL,NULL,NULL,diff_node,
	                           l_comp,r_comp,&cur,POSITIVE_ORIENTATION,
	                           w_type,REFLECTED,INCIDENT,0.0,front);
	    }
	}
	else 
	{
	    /* Make deflected contact */

	    if (is_plus_orientation) 
	    {
	    	l_comp = COMPIB;
	    	r_comp = COMPTB;
	    }
	    else
	    {
	    	l_comp = COMPTB;
	    	r_comp = COMPIB;
	    }
	    make_rotated_curve(h,RP->ang[4],L,U,NULL,NULL,NULL,diff_node,
	                       l_comp,r_comp,&cur,NEGATIVE_ORIENTATION,
	                       CONTACT,SLIP,INCIDENT,surf_ten,front);
	}
	debug_print("init","Leaving init_diffraction_node()\n");
}		/*end init_diffraction_nodes*/

/*
*		centered_Prandtl_Meyer_wave_state():
*	
*	Initializes the states in a region which contains a
*	centered rarefaction wave.
*/

	/* Initialization info for centered rarefaction waves */

typedef struct {
	Locstate	state0, state1;
	double		ang0, ang1;
	double		theta0, theta1;
	double		q0[MAXD], q1[MAXD];
	double		A0, A1;
	double		abs_v[MAXD];
	double		qhat2;
	double		center[MAXD];
	int		is_C_minus_wave;
} CENTRD_RAREF;

/*ARGSUSED*/
LOCAL	void	get_state_in_Prandtl_Meyer_wave(
	double		*coords,
	Locstate	state,
	COMP_TYPE	*ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		stype)
{
	CENTRD_RAREF	*cntrd_raref = (CENTRD_RAREF *) ct->extra;
	Locstate	state0 = cntrd_raref->state0;
	Locstate	state1 = cntrd_raref->state1;
	double		*center = cntrd_raref->center;
	double		A0 = cntrd_raref->A0;
	double		d[MAXD], x, y;
	double		q, H, A, turn_ang;
	double		*q0 = cntrd_raref->q0;
	double		w, v[MAXD];
	double		ang;
	double		Psi;
	int		is_C_minus_wave = cntrd_raref->is_C_minus_wave;

	debug_print("init_states","Entered get_state_in_Prandtl_Meyer_wave()\n");
	d[0] = coords[0] - center[0];	d[1] = coords[1] - center[1];
	ang = angle(d[0],d[1]);
	if ((ang <= cntrd_raref->ang0 && is_C_minus_wave) || 
	    (ang >= cntrd_raref->ang0 && !is_C_minus_wave)) 
	{
	    set_state(state,stype,state0);
	    return;
	}
	else if ((ang >= cntrd_raref->ang1 && is_C_minus_wave) ||
	         (ang <= cntrd_raref->ang1 && !is_C_minus_wave)) 
	{
	    set_state(state,stype,state1);
	    return;
	}
	x = q0[0]*d[0] + q0[1]*d[1];
	y = q0[0]*d[1] - q0[0]*d[0];
	Psi = angle(x,y);
	w = (is_C_minus_wave) ? -A0 - Psi : Psi - A0;
	A = state_in_prandtl_meyer_wave(w,A0,state0,cntrd_raref->A1,state1,
	                                state,TGAS_STATE);
	turn_ang = (is_C_minus_wave) ? (A - A0) - w : w - (A - A0);
	H = specific_enthalpy(state);
	q = sqrt(cntrd_raref->qhat2 - 2.0*H);
	v[0] = q * cos(cntrd_raref->theta0 + turn_ang);
	v[1] = q * sin(cntrd_raref->theta0 + turn_ang);
	Vel(state)[0] = v[0] + cntrd_raref->abs_v[0];
	Vel(state)[1] = v[1] + cntrd_raref->abs_v[1];
	set_state(state,stype,state);
}		/*end get_state_in_Prandtl_Meyer_wave*/

LOCAL	void init_centered_Prandtl_Meyer_wave(
	Locstate	state0,
	Locstate	state1,
	double		angle0,
	double		angle1,
	double		*cen,
	int		is_C_minus_wave,
	double		*abs_v,
	COMP_TYPE	*comp_type)
{
	CENTRD_RAREF	*cntrd_raref;
	double		c0sq, c1sq, q0sq, q1sq, theta0, theta1;
	double		M0sq, M1sq;
	double		*v0, *v1;
	double		H0, H1;

	if (Different_params(state0,state1))
	{
	    screen("ERROR in init_centered_Prandtl_Meyer_wave(), "
	           "different params\n");
	    clean_up(ERROR);
	}
	cntrd_raref = (CENTRD_RAREF *) comp_type->extra;
	copy_state(cntrd_raref->state0,state0);
	cntrd_raref->ang0 = angle0;
	copy_state(cntrd_raref->state1,state1);
	cntrd_raref->ang1 = angle1;
	cntrd_raref->abs_v[0] = abs_v[0];
	cntrd_raref->abs_v[1] = abs_v[1];
	cntrd_raref->center[0] = cen[0];
	cntrd_raref->center[1] = cen[1];
	cntrd_raref->is_C_minus_wave = is_C_minus_wave;

	v0 = cntrd_raref->q0;	v1 = cntrd_raref->q1;
	c0sq = sound_speed_squared(state0);
	c1sq = sound_speed_squared(state1);
	H0 = specific_enthalpy(state0);
	H1 = specific_enthalpy(state1);

	v0[0] = vel(0,state0) - abs_v[0];
	v0[1] = vel(1,state0) - abs_v[1];
	q0sq = sqr(v0[0]) + sqr(v0[1]);
	theta0 = angle(v0[0],v0[1]);

	v1[0] = vel(0,state1) - abs_v[0];
	v1[1] = vel(1,state1) - abs_v[1];
	q1sq = sqr(v1[0]) + sqr(v1[1]);
	theta1 = angle(v1[0],v1[1]);

	cntrd_raref->qhat2 = 0.5*(q0sq + 2.0*H0 + q1sq + 2.0*H1);
	M0sq = q0sq / c0sq;
	cntrd_raref->A0 = asin(1.0/M0sq);
	M1sq = q1sq / c1sq;
	cntrd_raref->A1 = asin(1.0/M1sq);
	if (M0sq < SONIC_MINUS_SQR || M1sq < SONIC_MINUS_SQR)
	{
	    screen("ERROR in init_centered_Prandtl_Meyer_wave(), "
	           "Rarefaction wave in subsonic state\n");
	    (void) printf("M0sq = %g, M1sq = %g\n",M0sq,M1sq);
	    clean_up(ERROR);
	}
	cntrd_raref->theta0 = theta0;
	cntrd_raref->theta1 = theta1;

}		/*end init_centered_Prandtl_Meyer_wave*/

        /* Initialization info for two pre-diffraction shocks. */

typedef struct {
    Locstate st_at_front;	/* the states just behind the leading shock   */
    Locstate st_at_tail;	/*     and just ahead of the trailing shock   */
    double x_edge;               /* the x coord of the left edge of the contact*/
    double y_front;              /* the height of the leading shock above the  */
                                /*     left edge of the contact               */
    double y_tail;               /* the height of the trailing shock above the */
                                /*     left edge of the contact               */
    double theta;                /* the angle of the shocks wrt the horizontal */
} NWAVE_PARAMS;

LOCAL void init_pre_diffraction_node(
	Gas_param	*params,
	Gas_param	*params2,
	Front		*front,
	boolean		init_trailing_shock,
	INIT_DATA	*init)
{
	RECT_GRID	*rg = front->rect_grid;
	const int	dim = rg->dim;
	double		*h = rg->h;
	double		*L = rg->L;
	double		*U = rg->U;

	FOURIER_POLY	*fpoly = NULL;
	NWAVE_PARAMS	*nwp = NULL;
	NODE		*ns, *ne;
	CURVE		*cur;
	SIN_SQR_PERT_PARAMS pert_params;

	Locstate	st_0, st_1, st_2, st_3, st_4;

	double		surf_ten;
	double		p0, p2, rho0, rho4, v0, v_shear = 0.0;
	double		contact_ang, contact_height, sh_to_ct_d;
	double		incident_angle;
	double		coords[MAXD], newcoords[MAXD];
	char		s[Gets_BUF_SIZE];
	double		c_t[MAXD];
	double		s_n[MAXD], s_t[MAXD];
	static const char	*fmt = "%lf %lf";

	set_ambient_comp_type(comp_type(COMPIF),front);
	st_0 = Ambient(comp_type(COMPIF));
	set_ambient_comp_type(comp_type(COMPTF),front);
	st_4 = Ambient(comp_type(COMPTF));

	if (init_trailing_shock == YES)
	{
	    set_nwave_comp_type(comp_type(COMPIB),front);
	    nwp = (NWAVE_PARAMS*)comp_type(COMPIB)->extra;
	    st_1 = nwp->st_at_front;
	    st_2 = nwp->st_at_tail;
	    set_ambient_comp_type(comp_type(COMPBB),front);
	    st_3 = Ambient(comp_type(COMPBB));
	}
	else
	{
	    set_ambient_comp_type(comp_type(COMPIB),front);
	    st_1 = Ambient(comp_type(COMPIB));
	}

	    /* Initialize front contact */

	screen("Do you wish a perturbed material interface (dflt = no): ");
	(void) Gets(s);
	if (s[0] == 'y' || s[0] == 'Y')
	{
	    fpoly = prompt_for_fourier_poly(rg,1.0);
	    contact_height = fpoly->z0;
	    contact_ang = 0.0;
	}
	else
	{
	    screen("Enter the angle (in degrees, between 0 and 90) ");
	    screen("which the\n");
	    screen("\tcontact makes with the x axis: ");
	    (void) Scanf("%f\n",&contact_ang);
	    contact_ang = radians(contact_ang);

	    screen("Enter the height of the material interface above L[1]\n");
	    screen("\tat the left boundary: ");
	    (void) Scanf("%f\n",&contact_height);

	    screen("\nIn order to have a sinusoidally perturbed");
	    screen(" material interface,\nenter the frequency nu and");
	    screen(" the amplitude epsilon of the perturbation.\n"); 
	    screen("\tEnter choice here: ");
	    pert_params.nu = pert_params.epsilon = 0.0;
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscanf(s,fmt,&pert_params.nu,&pert_params.epsilon);
	}
	c_t[0] = cos(contact_ang);	c_t[1] = sin(contact_ang);

	screen("The%s shock is assumed to be incident on the material ",
	        (nwp != NULL) ? " leading" : "");
	screen("interface\n\tfrom the left. ");
	screen("Enter the height of the shock above the\n");
	screen("\tleft hand edge of the material interface: ");
	(void) Scanf("%f\n",&sh_to_ct_d);

	screen("Enter the angle (in degrees, between 0 and 90) ");
	screen("of this shock\n");
	screen("\twith respect to the x axis: ");
	(void) Scanf("%f\n",&incident_angle);
	incident_angle = radians(incident_angle);
	s_t[0] = cos(incident_angle);	s_t[1] = sin(incident_angle);
	s_n[0] = s_t[1];		s_n[1] = -s_t[0];

	/* Input the state information ahead */

	screen("\nEnter the density and pressure in front of the %s shock\n",
	        (nwp != NULL) ? "leading" : "incident");
	screen("\tAND the density on the other side of the\n");
	screen("\tmaterial interface: ");
	(void) Scanf("%f %f %f\n",&rho0,&p0,&rho4);

	v0 = 0.0;			/* default */
	screen("Enter the velocity in front of the incident shock, ");
	screen("(dflt = %g): ",v0);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&v0);
	Vel(st_0)[0] = v0*c_t[0];	Vel(st_0)[1] = v0*c_t[1];

	Init_params(st_0,params);
	Dens(st_0) = rho0;
	Press(st_0) = p0;
	set_type_of_state(st_0,TGAS_STATE);
	set_state(st_0,GAS_STATE,st_0);

	prompt_for_behind_shock_state(st_0,st_1,YES,s_n,TGAS_STATE,YES,init);

	if (nwp != NULL)
	{
	    nwp->x_edge = L[0];
	    nwp->theta = incident_angle;
	    nwp->y_front = sh_to_ct_d;

	    screen("Enter the height of the trailing shock above the\n");
	    screen("\tleft hand edge of the material interface: ");
	    (void) Scanf("%f\n",&nwp->y_tail);

	    screen("The velocity just ahead of this shock is assumed to be ");
	    screen("that which exists\n\tjust behind the leading shock, ");
	    screen("while the region between\n\tthese shocks is isentropic.\n");
	    screen("Enter the pressure in front of the trailing shock: ");
	    (void) Scanf("%f\n",&p2);

	    Vel(st_2)[0] = vel(0,st_0); Vel(st_2)[1] = vel(1,st_0);

	    Init_params(st_2,params);
	    Init_params(st_3,params);
	    Dens(st_2) = dens_rarefaction(p2,st_0);
	    Press(st_2) = p2;
	    set_type_of_state(st_2,TGAS_STATE);
	    set_state(st_2,GAS_STATE,st_2);

	    prompt_for_behind_shock_state(st_2,st_3,YES,s_n,
	                                  TGAS_STATE,YES,init);
	}

	surf_ten = prompt_for_surface_tension(CONTACT,"for the contact ");

	/* Store data for front states */

	v_shear = 0.0;			/* default */
	screen("Enter a shear velocity across the material interface, ");
	screen("(dflt = %g): ",v_shear);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&v_shear);

	Init_params(st_4,params2);
	Dens(st_4) = rho4;
	Vel(st_4)[0] = Vel(st_0)[0] + v_shear*c_t[0];
	Vel(st_4)[1] = Vel(st_0)[1] + v_shear*c_t[1];
	Press(st_4) = p0;
	set_type_of_state(st_4,TGAS_STATE);
	set_state(st_4,GAS_STATE,st_4);

	    /* Make leading incident shock */

	coords[0] = L[0];
	coords[1] = L[1] + contact_height + sh_to_ct_d;
	ns = make_node(Point(coords));
	if (!intersect_ray_with_boundary(coords,s_t,L,U,newcoords,dim))
	{
	    screen("ERROR in init_pre_diffraction_node(), "
	           "can't find opposite point of incident shock\n");
	    clean_up(ERROR);
	}
	ne = make_node(Point(newcoords));
	cur = make_curve(COMPIB,COMPIF,ns,ne);
	wave_type(cur) = FORWARD_SHOCK_WAVE;
	start_status(cur) = end_status(cur) = INCIDENT;
	screen("Type 'y' to turn off tracking of the %s shock: ",
	        (nwp != NULL) ? "leading" : "incident");
	(void) Gets(s);
	if ((s[0] == 'y') || (s[0] == 'Y'))
	    untracked_hyper_surf(cur) = YES;

	    /* Make trailing incident shock (if any) */
	
	if (nwp != NULL)
	{
	    coords[0] = L[0];
	    coords[1] = L[1] + contact_height + nwp->y_tail;
	    ns = make_node(Point(coords));
	    if (!intersect_ray_with_boundary(coords,s_t,L,U,newcoords,dim))
	    {
	        screen("ERROR in init_pre_diffraction_node(), ");
	        screen("can't find opposite point of incident shock\n");
	        clean_up(ERROR);
	    }
	    ne = make_node(Point(newcoords));
	    cur = make_curve(COMPBB,COMPIB,ns,ne);
	    wave_type(cur) = FORWARD_SHOCK_WAVE;
	    start_status(cur) = end_status(cur) = INCIDENT;
	    screen("Type 'y' to turn off tracking of the trailing shock: ");
	    (void) Gets(s);
	    if ((s[0] == 'y') || (s[0] == 'Y'))
	        untracked_hyper_surf(cur) = YES;
	}
	    /* Make front contact */

	if (fpoly != NULL)
	{
	    int num_points = rg->gmax[0];
	    g_make_fourier_curve(CONTACT,num_points,L[0],U[0],
	    		         fpoly,COMPIF,COMPTF,surf_ten);
	    free(fpoly);
	    fpoly = NULL;
	}
	else
	{
	    coords[0] = L[0];
	    coords[1] = L[1] + contact_height;
	    ns = make_node(Point(coords));
	    make_rotated_curve(h,contact_ang,L,U,sin_sqr_pert,
	                       sin_sqr_pert_prime,(POINTER) &pert_params,ns,
	                       COMPIF,COMPTF,&cur,POSITIVE_ORIENTATION,CONTACT,
	                       INCIDENT,INCIDENT,surf_ten,front);
	}
}		/*end init_pre_diffraction_node*/

/*
*			init_transmission_node(): 
*
*	This routine initializes the states about the transmission node.
*	There are two possible choices for one of the fundamental parameters
*	used in find_transmission_node_states to calculate the states about
*	the transmission node, the slip or the incident angle.  If the 
*	incident angle is used, the velocity below right is first set equal to 
*	the velocity above right.  The correct value is then ft_assigned
*	in find_transmission_node_states.
*/

LOCAL void init_transmission_node(
	Gas_param	*params,
	Gas_param	*params2,
	RP_DATA		**pRP,
	Front		*front)
{
	RP_DATA	  *RP;
	NODE_FLAG flag;
	size_t	  sizest = front->sizest;
	double	  rho_a,rho_b,pr_r,pr_l;
	double	  **t;
	double	  q[5];
	double	  slip,vx_av,vy_av;
	double	  inc_ang;

	debug_print("init","Entering init_transmission_node()\n");

	clear_node_flag(flag);
	bi_array(&t,2,MAXD,FLOAT);


	/* Input densities ahead */

	screen("Enter the densities above and below: ");
	(void) Scanf("%f %f\n",&rho_a,&rho_b);

	/* Input pressures ahead */

	screen("Enter the pressures left (behind) and right (ahead): ");
	(void) Scanf("%f %f\n",&pr_l,&pr_r);

	switch (init_trans_node_param_choice())
	{
	case USE_SLIP:
	    t[0][0] = t[0][1] = ERROR_FLOAT;
	    screen("\nEnter the slip: ");
	    (void) Scanf("%f\n",&slip);

	/* Input the average of the velocities above right and below right */

	    screen("Enter x component of mean velocity ahead: ");
	    (void) Scanf("%f\n",&vx_av);
	    screen("Enter y component of mean velocity ahead: ");
	    (void) Scanf("%f\n",&vy_av);
	    break;

	case USE_INCIDENT_ANGLE:

	    screen("\nEnter the incident angle in degrees: ");
	    (void) Scanf("%f\n",&inc_ang);
	    inc_ang = radians(inc_ang);
	    t[0][0] = cos(inc_ang);	t[0][1] = sin(inc_ang);
	    slip = 0.0;	/* For convenience */

	/* Input the velocity above right */

	    screen("Enter x component of velocity above right: ");
	    (void) Scanf("%f\n",&vx_av);
	    screen("Enter y component of velocity above right: ");
	    (void) Scanf("%f\n",&vy_av);
	    break;

	default:
	    screen("ERROR: in init_transmission_node(), "
	           "Unknown parameter choice\n");
	    clean_up(ERROR);
	}

	*pRP = RP = allocate_RP_DATA_structure(sizest,YES,GAS_STATE);
	RP->ang_dir = COUNTER_CLOCK;

	    /* Set state above to the right */

	Init_params(RP->state[0],params);
	Dens(RP->state[0]) = rho_a;
	Vel(RP->state[0])[0] = (vx_av - 0.5*slip);
	Vel(RP->state[0])[1] = vy_av;
	Press(RP->state[0]) = pr_r;
	set_type_of_state(RP->state[0],TGAS_STATE);
	set_state(RP->state[0],GAS_STATE,RP->state[0]);

	    /* Set state below to the right */

	Init_params(RP->state[4],params2);
	Dens(RP->state[4]) = rho_b;
	Vel(RP->state[4])[0] = (vx_av + .5*slip);
	Vel(RP->state[4])[1] = vy_av;
	Press(RP->state[4]) = pr_r;
	set_type_of_state(RP->state[4],TGAS_STATE);
	set_state(RP->state[4],GAS_STATE,RP->state[4]);


	    /* Set state behind, above contact */

	Init_params(RP->state[1],params);
	Dens(RP->state[1]) = 1.0;	/* Artifical values, used to pass */
	Vel(RP->state[1])[0] = 0.0;	/* pressure via Gas structure     */
	Vel(RP->state[1])[1] = 0.0;
	Press(RP->state[1]) = pr_l;
	set_type_of_state(RP->state[1],TGAS_STATE);
	set_state(RP->state[1],GAS_STATE,RP->state[1]);

	    /* Set state behind, below contact */

	Init_params(RP->state[3],params2);


	t[1][0] = 1.0;		t[1][1] = 0.0;
	(void) find_transmission_node_states(q,t,RP,
	                                     trans_node_parameter(),WEAK,flag);
	free(t);
	t = NULL;


	if (debugging("init")) 
	{
	    (void) printf("States after initialization.\n");
	    verbose_print_state("above right",RP->state[0]);
	    verbose_print_state("below right",RP->state[4]);
	    verbose_print_state("above left",RP->state[1]);
	    verbose_print_state("below left",RP->state[3]);
	    print_angle("Incident angle =",RP->ang[1],"\n");
	    print_angle("Transmitted angle =",RP->ang[4],"\n");
	    print_angle("Contact angle =",RP->ang[3],"\n");
	}
	debug_print("init","Leaving init_transmission_node()\n");
}		/*end init_transmission_node*/


LOCAL double meshkov_amb_vert_vel(
	double		rho_a,
	double		rho_b,
	double		p0,
	double		p1,
	Gas_param	*params,
	Gas_param	*params2)
{
	double		         ustarl, ustarr;
	double		         pstarl, pstarr;
	double		         shock_speed, ml, mr;
	RIEMANN_SOLVER_WAVE_TYPE l_wave, r_wave;
	size_t		         sizest = params->sizest;
	int		         dim = params->dim;
	Locstate	         st_l, st_r;
	static double	         s_n[MAXD] = {-1.0, 0.0};
	
	(*params->_alloc_state)(&st_l,sizest);
	(*params->_alloc_state)(&st_r,sizest);
	Init_params(st_l,params);
	zero_state_velocity(st_l,dim);
	Dens(st_l) = rho_a;
	Press(st_l) = p0;
	set_type_of_state(st_l,TGAS_STATE);
	set_state(st_l,GAS_STATE,st_l);
	(void) s_polar_4(BEHIND_PRESSURE,p1,&shock_speed,s_n,st_l,st_r,
	                 GAS_STATE);
	set_state(st_r,TGAS_STATE,st_r);
	Init_params(st_l,params2);
	zero_state_velocity(st_l,dim);
	Press(st_l) = p0;
	Dens(st_l) = rho_b;
	reset_gamma(st_l);
	(void) find_mid_state(st_l,st_r,0.0,&pstarl,&pstarr,&ustarl,&ustarr,
	                      &ml,&mr,&l_wave,&r_wave);

	free_these(2,st_l,st_r);
	return -0.5*(ustarl + ustarr);
}		/*end meshkov_amb_vert_vel*/

LOCAL	void init_taylor_wave_component(
	COMP_TYPE	*comp_type,
	_TAYLOR_WAVE	*twave,
	int		tbar_given,
	double		sign,
	Locstate	stas,
	Locstate	st0,
	Locstate	st1,
	Front		*front)
{
	double		c0, c1;
	double		l0, l1;
	int		i, dim = front->rect_grid->dim;

	twave = set_taylor_wave_comp_type(comp_type,twave,front);

	set_state(twave->st0,GAS_STATE,st0);
	set_state(twave->st1,GAS_STATE,st1);
	set_state(twave->stas,GAS_STATE,stas);

	twave->sizest = front->sizest;
	twave->sgn = sign;

	for (i = 0; i < dim; i++)
	    twave->va[i] = vel(i,st0);
	twave->vz0 = vel(1,st0);

	twave->l_or_r = (sign == 1.0) ? LEFT_FAMILY : RIGHT_FAMILY;

	c0 = sound_speed(twave->st0);
	l0 = vel(1,twave->st0) - sign*c0;
	c1 = sound_speed(twave->st1);
	l1 = vel(1,twave->st1) - sign*c1;

	if (tbar_given)
	    twave->z1 = twave->z0 + twave->tbar * (l1 - l0);
	else
	    twave->tbar = (twave->z1 - twave->z0)/(l1 - l0);

	twave->zbar = 0.5*(twave->z0 + twave->z1 - (l0 + l1)*twave->tbar);

	if (debugging("taylor_wave")) print_TAYLOR_WAVE_structure(twave);
}		/*end init_taylor_wave_component*/


LOCAL	void print_TAYLOR_WAVE_structure(
	_TAYLOR_WAVE	*twave)
{
	int		dim = Params(twave->st0)->dim;

	(void) printf("Printout of _TAYLOR_WAVE structure %p\n",(POINTER)twave);
	(void) printf("z0 = %g, z1 = %g, zs = %g\n",twave->z0,twave->z1,
	              twave->zs);
	(void) printf("sign = %g\n",twave->sgn);
	(void) printf("vz0 = %g, vz1 = %g, ",twave->vz0,twave->vz1);
	print_general_vector("va = ",twave->va,dim,"\n");
	(void) printf("zbar = %g, tbar = %g\n",twave->zbar,twave->tbar);
	(void) printf("sizest = %d\n",(int)twave->sizest);
	switch(twave->l_or_r)
	{
	case LEFT_FAMILY:
	    (void) printf("l_or_r = LEFT_FAMILY (%d)\n",twave->l_or_r);
	    break;
	case RIGHT_FAMILY:
	    (void) printf("l_or_r = RIGHT_FAMILY (%d)\n",twave->l_or_r);
	    break;
	default:
	    (void) printf("l_or_r = UNKNOWN WAVE FAMILY (%d)\n",twave->l_or_r);
	}
	verbose_print_state("stas",twave->stas);
	verbose_print_state("st0",twave->st0);
	verbose_print_state("st1",twave->st1);
	(void) printf("End of printout of _TAYLOR_WAVE structure %p\n",
	       (POINTER)twave);
}		/*end print_TAYLOR_WAVE_structure*/

LOCAL void adjust_to_off_grid_line(
	double		*posn,
	int		idir,
	double		sign,
	Front		*front)
{
	double	    L = front->rect_grid->L[idir];
	double	    h = front->rect_grid->h[idir];
	double	    np, p, newposn;
	static const double MIN_GRID_SEP = 0.01; /*TOLERANCE*/

	np = rint((*posn - L)/h - 0.5);
	p = L + (np + 0.5)*h;
	if (fabs(*posn - p) >= MIN_GRID_SEP*h) return;

	newposn = p + sign*MIN_GRID_SEP*h;
	if (debugging("adjust"))
	    (void) printf("Adjusting position %g to %g\n",*posn,newposn);
	*posn = newposn;

}		/*end adjust_to_off_grid_line*/

LOCAL int init_trans_node_param_choice(void)
{
	char		s[Gets_BUF_SIZE];

	screen("Type (S) or (I) to use the slip or incident ");
	screen("angle to compute a\n");
	screen("\ttransmission node configuration.\n");
	screen("Enter choice here (default = S): ");

	(void) Gets(s);
	if (s[0] == 'i' || s[0] == 'I')
	    return set_trans_node_parameter(USE_INCIDENT_ANGLE);
	return set_trans_node_parameter(USE_SLIP);
}		/*end init_trans_node_param_choice*/


LOCAL	void	set_Prandtl_Meyer_wave_comp_type(
	COMP_TYPE	*comp_type,
	Front		*front)
{
	CENTRD_RAREF	*crf;

	if (comp_type->type == PRANDTL_MEYER_WAVE) /*ALREADY SET*/
	    return;

	comp_type->type = PRANDTL_MEYER_WAVE;
	scalar(&crf,sizeof(CENTRD_RAREF));
	alloc_state(front->interf,&crf->state0,front->sizest);
	alloc_state(front->interf,&crf->state1,front->sizest);
	comp_type->extra = (POINTER)crf;
	comp_type->_get_state = get_state_in_Prandtl_Meyer_wave;
	comp_type->free_comp_type_extra = free_Prandtl_Meyer_wave_comp_type;
}		/*end set_Prandtl_Meyer_wave_comp_type*/

LOCAL	void	free_Prandtl_Meyer_wave_comp_type(
	COMP_TYPE	*comp_type)
{
	CENTRD_RAREF *crf;

	if (comp_type->type != PRANDTL_MEYER_WAVE)
	    return;

	crf = (CENTRD_RAREF*)comp_type->extra;
	if (crf == NULL)
	    return;

	free(crf->state0);
	free(crf->state1);
	free(crf);
	comp_type->extra = NULL;
}		/*end free_Prandtl_Meyer_wave_comp_type*/
	
/*
*			get_state_taylor_wave():
*
*	Returns the state inside a one dimensional rarefaction wave
*	at the location coords.
*/

/*ARGSUSED*/
LOCAL	void 	get_state_taylor_wave(
	double		*coords,
	Locstate	state,
	COMP_TYPE	*comp_type,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		stype)
{
	_TAYLOR_WAVE	*twave = Taylor_wave(comp_type);
	double		speed;
	double		vy;
	double		spdnew;
	int		dim = Params(twave->st0)->dim;
	int		i, k = dim-1;

	debug_print("init_states","Entered get_state_taylor_wave()\n");
	if ( comp_type->type != TAYLOR_WAVE )
	{
	    screen("ERROR: in get_state_taylor_wave()\n"
	           "inconsistent comp_type->type\n");
	    clean_up(ERROR);
	}

	debug_print("taylor_wave","Entered get_state_taylor_wave()\n");
	if (debugging("taylor_wave"))
	{
	    print_general_vector("coords = ",coords,dim,"\n");
	}
	if (twave->stas != NULL)
	{ 
	    if ((twave->l_or_r == LEFT_FAMILY && coords[k] < twave->zs) 
	     || (twave->l_or_r == RIGHT_FAMILY && coords[k] > twave->zs) )
	    {
	        set_state(state,stype,twave->stas); 
	        if (debugging("taylor_wave"))
	            (void) printf("State set to ahead shock state\n");
	        debug_print("taylor_wave","Leaving get_state_taylor_wave()\n");
	        return; 
	    }
	}
	if ((twave->l_or_r == LEFT_FAMILY && coords[k] < twave->z0) 
	 || (twave->l_or_r == RIGHT_FAMILY && coords[k] > twave->z0) )
	{
	    set_state(state,stype,twave->st0); 
	    if (debugging("taylor_wave"))
	    	(void) printf("State set to state 0\n");
	    debug_print("taylor_wave","Leaving get_state_taylor_wave()\n");
	    return; 
	}
	if ((twave->l_or_r == LEFT_FAMILY && coords[k] > twave->z1) 
	 || (twave->l_or_r == RIGHT_FAMILY && coords[k] < twave->z1) )
	{
	    set_state(state,stype,twave->st1); 
	    if (debugging("taylor_wave"))
	    	(void) printf("State set to state 1\n");
	    debug_print("taylor_wave","Leaving get_state_taylor_wave()\n");
	    return; 
	}
	speed = (coords[k] - twave->zbar)/twave->tbar;
	(void) oned_state_in_rarefaction_fan(speed,twave->vz0,twave->st0,
	                                     twave->st1,state,TGAS_STATE,
	                                     &spdnew,twave->l_or_r);
	vy = Vel(state)[0];
	for (i = 0; i < k; i++)
	    Vel(state)[i] = twave->va[i];
	Vel(state)[k] = vy;
	set_state(state,stype,state);
	if (debugging("taylor_wave"))
	{
	    (void) printf("State set inside rarefaction fan\n");
	}
	debug_print("taylor_wave","Leaving taylor_wave_state()\n");
}		/*end get_state_taylor_wave*/

LOCAL	_TAYLOR_WAVE	*set_taylor_wave_comp_type(
	COMP_TYPE	*comp_type,
	_TAYLOR_WAVE    *twave,
	Front		*front)
{
	if (comp_type->type == TAYLOR_WAVE) /*ALREADY SET*/
	    return Taylor_wave(comp_type);

	comp_type->type = TAYLOR_WAVE;
	if (twave == NULL)
	    scalar(&twave,sizeof(_TAYLOR_WAVE));
	comp_type->extra = (POINTER) twave;
	if (twave->st0 == NULL)
	    alloc_state(front->interf,&twave->st0,front->sizest);
	if (twave->st1 == NULL)
	    alloc_state(front->interf,&twave->st1,front->sizest);
	if (twave->stas == NULL)
	    alloc_state(front->interf,&twave->stas,front->sizest);
	comp_type->_get_state = get_state_taylor_wave;
	comp_type->free_comp_type_extra = free_taylor_wave_comp_type;
	return twave;
}		/*end set_taylor_wave_comp_type*/

LOCAL	void	free_taylor_wave_comp_type(
	COMP_TYPE	*comp_type)
{
	_TAYLOR_WAVE    *twave;
	if (comp_type->type != TAYLOR_WAVE)
	    return;
	
	twave = (_TAYLOR_WAVE*)comp_type->extra;
	if (twave == NULL)
	    return;
	if (twave->st0 != NULL)
	    free(twave->st0);
	if (twave->st1 != NULL)
	    free(twave->st1);
	if (twave->stas != NULL)
	    free(twave->stas);
	comp_type->extra = NULL;
}		/*end free_taylor_wave_comp_type*/

LOCAL	void	set_nwave_comp_type(
	COMP_TYPE	*comp_type,
	Front		*front)
{
	NWAVE_PARAMS	*nwp;
	if (comp_type->type == NWAVE) /*ALREADY SET*/
	    return;

	comp_type->type = NWAVE;
	scalar(&nwp,sizeof(NWAVE));
	alloc_state(front->interf,&nwp->st_at_front,front->sizest);
	alloc_state(front->interf,&nwp->st_at_tail,front->sizest);
	comp_type->extra = (POINTER) nwp;
	comp_type->_get_state = get_state_nwave;
	comp_type->free_comp_type_extra = free_nwave_comp_type;
}		/*end set_nwave_comp_type*/

LOCAL	void	free_nwave_comp_type(
	COMP_TYPE	*comp_type)
{
	NWAVE_PARAMS	*nwp;

	if (comp_type->type != NWAVE)
	    return;
	
        nwp = (NWAVE_PARAMS*)comp_type->extra;
	if (nwp == NULL)
	    return;
	if (nwp->st_at_front != NULL)
	    free(nwp->st_at_front);
	if (nwp->st_at_tail != NULL)
	    free(nwp->st_at_tail);
	free(nwp);
	comp_type->extra = NULL;
}		/*end free_nwave_comp_type*/

/* 			get_state_nwave():
*
*	Compute the state at a position s between two parallel shocks incident
*	on a contact surface, assuming that the pressure varies linearly along
*	the normal distance between the shocks, while the velocity is unchanged.
*/

/*ARGSUSED*/
LOCAL	void	get_state_nwave(
	double		*coords,
	Locstate	s,
	COMP_TYPE	*comp_type,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		stype)
{
        const NWAVE_PARAMS *nwp = (NWAVE_PARAMS*)comp_type->extra;

	/* the y-coordinate where the segment connecting coords to the
	   horizontal axis intersects the leading shock: */

	double y_int = (coords[0] - nwp->x_edge)*tan(nwp->theta) + nwp->y_front;

	/* the ratio of the normal distance from coords to the leading shock
	   to the distance between the shocks: */

	double r = (coords[1] - y_int)/(nwp->y_tail - nwp->y_front);

	double pf = pressure(nwp->st_at_front);
	double pt = pressure(nwp->st_at_tail);

	set_state(s,TGAS_STATE,nwp->st_at_front);
	Press(s) = r*pt + (1.0-r)*pf;
	Dens(s) = dens_rarefaction(Press(s),nwp->st_at_front);
	set_state(s,stype,s);

	if (debugging("nwave"))
	{
	    (void) printf("In set_state_nwave(), "
			  "coords[0] = %g; coords[1] = %g\n",
	      coords[0],coords[1]);
	    (void) printf("\tpre at front = %g; pre at tail = %g; r = %g\n",
	                  pf,pt,r);
	    (void) printf("\tden = %g; pre = %g; v_x = %g; v_y = %g\n",
	                  Dens(s),pressure(s),vel(0,s),vel(1,s));
	}

	return;
}		/*end set_state_nwave*/

#endif /* defined(TWOD) */
#endif /* defined(FULL_PHYSICS) */
