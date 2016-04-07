gbifur/gbifur.c:LOCAL	void	normal_to_regular_reflection(float,float*,BOND*,O_CURVE*,
gbifur/gbifur.c:LOCAL	void	normal_to_mach_reflection(float,BOND*,O_CURVE*,O_CURVE*,
gbifur/gbifur.c:LOCAL	int	untrack_B_reflect_node(Front*,Wave*,O_CURVE*,O_CURVE*,
gbifur/gbifur.c:*	A B_node can bifurcate into a regular reflection, Mach stem or
gbifur/gbifur.c:*       then incident angle < PI/2 means a transition to a regular reflection
gbifur/gbifur.c:	debug_print("gBnode","Calling is_regular_reflection\n");
gbifur/gbifur.c:	if (is_regular_reflection(node_v,fr,RP))
gbifur/gbifur.c:	    if (!track_scattered_wave(B_REFLECT_NODE,SHOCK_WAVE,
gbifur/gbifur.c:					 REFLECTED,RP->state[1],
gbifur/gbifur.c:	    normal_to_regular_reflection(initial_time_elapsed,node_v,
gbifur/gbifur.c:	    normal_to_mach_reflection(initial_time_elapsed,bcorner,oldcphys,
gbifur/gbifur.c:*			normal_to_regular_reflection():
gbifur/gbifur.c:*	of topology from a normal shock (B_NODE) to a regular reflection
gbifur/gbifur.c:*	configuration (B_REFLECT_NODE).
gbifur/gbifur.c:LOCAL void normal_to_regular_reflection(
gbifur/gbifur.c:	debug_print("gBnode","Creating normal_to_regular_reflection bifurcation\n");
gbifur/gbifur.c:	    screen("ERROR in normal_to_regular_reflection(), "
gbifur/gbifur.c:	    (void) printf("real reflected wall length = %g\n",
gbifur/gbifur.c:	node_type(rnode) = B_REFLECT_NODE;
gbifur/gbifur.c:		 track_scattered_wave(B_REFLECT_NODE,SHOCK_WAVE,REFLECTED,
gbifur/gbifur.c:	    screen("ERROR in normal_to_regular_reflection(), "
gbifur/gbifur.c:}		/*end normal_to_regular_reflection*/
gbifur/gbifur.c:*			normal_to_mach_reflection():
gbifur/gbifur.c:*	of topology from a normal shock (B_NODE) to a Mach reflection
gbifur/gbifur.c:LOCAL void normal_to_mach_reflection(
gbifur/gbifur.c:	debug_print("gBnode","Creating Mach reflection bifurcation\n");
gbifur/gbifur.c:				 bubble,NORMAL_TO_MACH_REFLECTION,fr)
gbifur/gbifur.c:	    screen("ERROR in normal_to_mach_reflection(), "
gbifur/gbifur.c:	if (track_scattered_wave(MACH_NODE,SHOCK_WAVE,REFLECTED,
gbifur/gbifur.c:}		/*end normal_to_mach_reflection*/
gbifur/gbifur.c:*			g_reflect_node_bifurcation():
gbifur/gbifur.c:*	This routine gives the bifurcation of a regular reflection node to
gbifur/gbifur.c:EXPORT void g_reflect_node_bifurcation(
gbifur/gbifur.c:	debug_print("grBnode","Entered g_reflect_node_bifurcation\n");
gbifur/gbifur.c:	    screen("ERROR in g_reflect_node_bifurcation(), "
gbifur/gbifur.c:	    	(void) printf("reflected wave untracked\n");
gbifur/gbifur.c:				 bubble,REGULAR_TO_MACH_REFLECTION,fr))
gbifur/gbifur.c:	    if (untrack_node(fr,B_REFLECT_NODE))
gbifur/gbifur.c:	    	if (!untrack_B_reflect_node(fr,wave,oldcinc,cinc,oldcrefl,
gbifur/gbifur.c:		    screen("ERROR in g_reflect_node_bifurcation(), "
gbifur/gbifur.c:		debug_print("grBnode","Left g_reflect_node_bifurcation\n");
gbifur/gbifur.c:	    	screen("ERROR in g_reflect_node_bifurcation(), "
gbifur/gbifur.c:	    screen("ERROR in g_reflect_node_bifurcation(), "
gbifur/gbifur.c:	    /* Transfer reflected curve to Mach triple point */
gbifur/gbifur.c:	debug_print("grBnode","Left g_reflect_node_bifurcation\n");
gbifur/gbifur.c:}		/*end g_reflect_node_bifurcation*/
gbifur/gbifur.c:LOCAL int untrack_B_reflect_node(
gbifur/gbifur.c:	    screen("ERROR in untrack_B_reflect_node(), "
gbifur/gbifur.c:}		/*end untrack_B_reflect_node*/
gbifur/gbifur.c:*	This function is used to compute the node velocity as a reflected
gbifur/gbifur.c:	    else /* REGULAR_REFLECTION */
gbifur/gbifur.c:	        debug_print("gBnode","Regular reflection node case\n");
gbifur/gbifur.c:	    case B_REFLECT_NODE:
gbifur/gbifur.c:	            if(is_regular_reflection(node_v,fr,RP))
gbifur/gbifur.c:		case REFLECTION_BOUNDARY:
gbifur/gbifur.c:		case REFLECTION_BOUNDARY:
gbifur/gbifur.c:*	This function creates the reflected or bow shock when initializing
gbifur/gbifur.c:*	a reflection (reg or mach).
gbifur/gbifur.c:	    set_status_at_node(bow_shock->curve,bow_shock->orient,REFLECTED);
gbifur/gbifurprotos.h:IMPORT	void	g_reflect_node_bifurcation(Front*,Wave*,O_CURVE*,O_CURVE*,
gbifur/grefl.c:*	reflection or refraction problems.
gbifur/grefl.c:*	This routine initializes the states for a regular reflection,
gbifur/grefl.c:*	to fix the approximate shape of the reflection bubble, for 
gbifur/grefl.c:*		reflection state
gbifur/grefl.c:*		reflected shock angle
gbifur/grefl.c:*			reflection wall length
gbifur/grefl.c:*		reflection wall length
gbifur/grefl.c:	    print_angle("reflected shock angle relative to ahead wall =",
gbifur/grefl.c:	    verbose_print_state("reflection",refl);
gbifur/grefl.c:*	This routine initializes the states for a Mach reflection
gbifur/grefl.c:*	reflection where the incident angle has shrunk/grown resp.
gbifur/grefl.c:*	parameters to fix approximately the shape of the reflection
gbifur/grefl.c:*	reflection.
gbifur/grefl.c:*			reflection wall length
gbifur/grefl.c:*		reflected state - behind refl shock at triple point
gbifur/grefl.c:*		reflection wall length
gbifur/grefl.c:*		refl angle - pos x axis to reflected shock
gbifur/grefl.c:	if (flag == NORMAL_TO_MACH_REFLECTION)
gbifur/grefl.c:	else if (flag == REGULAR_TO_MACH_REFLECTION)
gbifur/grefl.c:			screen("for regular to mach reflection\n");
gbifur/grefl.c:	 * the reflected shock and the contact.
gbifur/grefl.c:	if ((flag == NORMAL_TO_MACH_REFLECTION) && !bubble->is_attached)
gbifur/grefl.c:	else if (flag == REGULAR_TO_MACH_REFLECTION)
gbifur/grefl.c:		print_angle("reflected shock angle wrt to ahead wall =",
gbifur/grefl.c:		if (flag == NORMAL_TO_MACH_REFLECTION)
gbifur/grefl.c:			(void) printf("In NORMAL_TO_MACH_REFLECTION case.\n");
gbifur/grefl.c:			(void) printf("In REGULAR_TO_MACH_REFLECTION case.\n");
gbifur/grefl.c:*	are bifurcating from regular to mach reflection and there is
gbifur/grefl.c:*	of the reflected shock bubble for regular and Mach reflection.
gbifur/grefl.c:*	For regular reflection, wall_to_node_ang = 0.
gbifur/grefl.c:*		refl_wall_length - corner to reflection point (reg)
gbifur/grefl.c:*	mach reflection), and the time elapsed since bifurcation.
gbifur/grefl.c:*	the reflection (RP->state[1] and RP->state[2]).
gbifur/grefl.c:		print_general_vector("Reflected shock tangent = ",
gbifur/grefl.c:		    (void) printf("Reflected wall length = %g\n",*refl_length);
gbifur/grefl.c:	if (start_status(c) == REFLECTED)
gbifur/grefl.c:	else if (end_status(c) == REFLECTED)
gbifur/grefl.c:	    {		/* mach reflection */
gbifur/grefl.c:	    {		/* regular reflection */
gbifur/grefl.c:	    {		/* mach reflection */
gbifur/grefl.c:	    {		/* regular reflection */
gbifur/grefl.c:	else if (node_type(c->start) == B_REFLECT_NODE)
gbifur/grefl.c:	{		/* refl wall of reg reflection */
gbifur/grefl.c:	else if (node_type(c->end) == B_REFLECT_NODE)
gbifur/grefl.c:	{		/* refl wall of reg reflection */
gbifur/grefl.c:	    else if (status_at_node(cps,cps_or) == REFLECTED)
gbifur/grefl.c:	    else if (status_at_node(cpe,cpe_or) == REFLECTED)
gbifur/grefl.c:	    case REFLECTED:		/* wall -- bow shock to corner */
gbifur/grefl.c:	    case REFLECTED:		/* wall -- corner to bow shock */
gbifur/grp.c:	    free_o_curve_family(rpn_reflected1(rp_node)); 
gbifur/grp.c:	    rpn_reflected1(rp_node) = NULL;
gbifur/grp.c:	    free_o_curve_family(rpn_reflected2(rp_node)); 
gbifur/grp.c:	    rpn_reflected2(rp_node) = NULL;
gbifur/grp.c:	    			    &c2,&orient2,REFLECTED);
gbifur/grp.c:		init_cfamily(&rpn_reflected1(rp_node),c1,orient1);
gbifur/grp.c:		init_cfamily(&rpn_reflected2(rp_node),c2,orient2);
gbifur/grp.c:	    replace_null_curves_in_family(&rpn_reflected1(rp_node),rp);
gbifur/grp.c:	    replace_null_curves_in_family(&rpn_reflected2(rp_node),rp);
gbifur/grp.c:	    relocate_null_pointer((POINTER *)&rpn_reflected1(rp_node),
gbifur/grp.c:			          (POINTER *)&rpn_reflected2(rp_node));
gbifur/grp.c:*	The continuation of a reflected wave through a node is not defined,
gbifur/grp.c:	case B_REFLECT_NODE:
gbifur/grp.c:	case B_REFLECT_NODE:
gbifur/grp.c:	case B_REFLECT_NODE:
gbifur/grp.c:*	possible to have reflected waves at A and B (on a Neumann boundary),
gbifur/grp.c:*	are no reflected waves (case 0 below).  In both cases, the null
gbifur/grp.c:	    if (find_curves_at_rp_with_status(oldoc,newoc,2,rp,REFLECTED) ==
gbifur/grp.c:	                     "can't find reflected curves\n");
gbifur/grp.c:	                     "code needed for other than 2 reflected curves\n");
gbifur/grp.c:	find_curve_with_status(rpn->node,&cp->curve,&cp->orient,REFLECTED);
gbifur/grp.c:	if (c != NULL) set_status_at_node(c,orient,REFLECTED);
gbifur/grp.c:	identify_curves_with_status(cnode,oc[1],oc[3],REFLECTED);
gbifur/grp.c:*	only the reflected, contact and transmitted waves (if tracked).  We
gbifur/grp.c:	    case REFLECTED:
gbifur/grp.c:	    if (rpn_reflected1(rp_node))
gbifur/grp.c:	    if (rpn_reflected2(rp_node))
gbifur/grp.c:	(void) printf("number of reflected curves = %d\n",rp_num_refl(rp));
gbifur/grp.c:	    rpn_reflected1(rpn) = rpn_reflected2(rpn) = NULL;
gbifur/grp.c:	    delete_curve_from_o_curve_family(curve,&rpn_reflected1(rpn));
gbifur/grp.c:	    delete_curve_from_o_curve_family(curve,&rpn_reflected2(rpn));
gbifur/grp.c:	free_o_curve_family(rpn_reflected1(rpn));
gbifur/grp.c:	free_o_curve_family(rpn_reflected2(rpn));
gbifur/grp.c:	(void) printf("Reflected curve families:\n");
gbifur/grp.c:	print_o_curve_family(rpn_reflected1(rpn));
gbifur/grp.c:	print_o_curve_family(rpn_reflected2(rpn));
gbifur/grp.c:	(void) printf("Number of reflected curves = %d\n",rp_num_refl(rp));
gbifur/grp.c:	rpn_reflected1(rpn) = NULL;
gbifur/grp.c:	rpn_reflected2(rpn) = NULL;
gbifur/gsc1.c:*	1.  Leading edge of reflected rarefaction wave (if it exists)
gbifur/gsc1.c:*	2.  Reflected shock (if it exists)
gbifur/gsc1.c:*	3.  Trailing edge of reflected rarefaction wave (if it exists)
gbifur/gsc1.c:*	4.  Contact behind the reflected wave
gbifur/gsc1.c:*	2.  Reflected shock (if it exists)
gbifur/gsc1.c:*	3.  Contact behind the reflected wave
gbifur/gsc1.c:	bool		is_reflected_shock;
gbifur/gsc1.c:				      REFLECTED,
gbifur/gsc1.c:				      REFLECTED,
gbifur/gsc1.c:				      REFLECTED,
gbifur/gsc1.c:					    &is_reflected_shock,front,
gbifur/gsc1.c:	case ANOMALOUS_REFLECTION:
gbifur/gsc1.c:	if (!angles_consistent_at_diff_node(RP,is_reflected_shock))
gbifur/gsc1.c:							 RP,is_reflected_shock,
gbifur/gsc1.c:	*  Make reflected and transmitted curves, set status at node
gbifur/gsc1.c:	int		is_reflected_shock)
gbifur/gsc1.c:	    (void) printf("is_reflected_shock = %s, ",
gbifur/gsc1.c:			 (is_reflected_shock) ? "YES" : "NO");
gbifur/gsc1.c:	    if (is_reflected_shock)
gbifur/gsc1.c:	if (is_reflected_shock)
gbifur/gsc1.c:	bool		is_reflected_shock,
gbifur/gsc1.c:	    (void) printf("is_reflected_shock = %s\n",
gbifur/gsc1.c:	                  y_or_n(is_reflected_shock));
gbifur/gsc1.c:	if (is_reflected_shock)
gbifur/gsc1.c:	        track_scattered_wave(DIFFRACTION_NODE,SHOCK_WAVE,REFLECTED,
gbifur/gsc1.c:				     REFLECTED,RP->state[1],RP->state[4],fr))
gbifur/gsc1.c:	    (void) printf("Reflected side boundary, ");
gbifur/gsc1.c:	/* Realign reflected boundary */
gbifur/gsc2.c:	    	    (void) printf("Reflected leading edge rarefaction\n");
gbifur/gsc2.c:	    	    (void) printf("Reflected shock\n");
gbifur/gsc2.c:	    	    (void) printf("Reflected trailing edge rarefaction\n");
gbifur/gsc3.c:	    case REFLECTED:
gbifur/gsc3.c:	    	    (void) printf("Reflected curve found\n");
gbifur/gsc3.c:	                  "\nReflected Side Boundary\n");
gbifur/gsc3.c:	*	the shocks or the bifurcation to mach reflection.
gbifur/gsc3.c:*	B_REFLECT_NODE:
gbifur/gsc3.c:	case B_REFLECT_NODE:
gbifur/guntan1d.c:LOCAL	void	reflect_point_at_neumann_boundary(POINT*,POINT*,SIDE,Front*);
gbifur/guntan1d.c:		rect_boundary_type(intfc,0,0) == REFLECTION_BOUNDARY)
gbifur/guntan1d.c:	    reflect_state(sl,intfc,Coords(ptr),gr->L,nor);
gbifur/guntan1d.c:		rect_boundary_type(intfc,0,1) == REFLECTION_BOUNDARY)
gbifur/guntan1d.c:	    reflect_state(sr,intfc,Coords(ptl),gr->U,nor);
gbifur/guntan1d.c:	case REFLECTION_BOUNDARY:
gbifur/guntan1d.c:	    reflect_point_at_neumann_boundary(bdry_pt,interior_pt,
gbifur/guntan1d.c:*		reflect_point_at_neumann_boundary():
gbifur/guntan1d.c:LOCAL	void reflect_point_at_neumann_boundary(
gbifur/guntan1d.c:	debug_print("guntan","Entered reflect_point_at_neumann_boundary()\n");
gbifur/guntan1d.c:	    reflect_state(sl,intfc,coords,Coords(bdry_pt),nor);
gbifur/guntan1d.c:	    reflect_state(sr,intfc,coords,Coords(bdry_pt),nor);
gbifur/guntan1d.c:	    screen("ERROR in reflect_point_at_neumann_boundary(), "
gbifur/guntan1d.c:	debug_print("guntan","Left reflect_point_at_neumann_boundary()\n");
gbifur/guntan1d.c:}		/*end reflect_point_at_neumann_boundary*/
gbifur/guntan2d.c:	bool		is_reflected_shock[2];
gbifur/guntan2d.c:	static	int	status[7] = { INCIDENT, REFLECTED, REFLECTED,
gbifur/guntan2d.c:				      REFLECTED, SLIP, TRANSMITTED,
gbifur/guntan2d.c:						   &is_reflected_shock[i],
gbifur/guntan2d.c:	        case ANOMALOUS_REFLECTION:
gbifur/guntan2d.c:	                                               &is_reflected_shock[i],
gbifur/guntan2d.c:							is_reflected_shock[i],
gbifur/guntan2d.c:	if ((cross1!=NULL) && (is_reflected_shock[0]!=is_reflected_shock[1]))
gbifur/guntan2d.c:	*	Assign components and start-end states of reflected, transmitted
gbifur/gvecuntan.c:	static int	status[7] = { OVERTOOK, INCIDENT, REFLECTED, REFLECTED,
gbifur/gvecuntan.c:					REFLECTED, SLIP, TRANSMITTED };
gbifur/gvecuntan.c:	    /* Set tracking of reflected waves */
gbifur/gvecuntan.c:					  REFLECTED,RP->state[2],RP->state[5],
gbifur/gvecuntan.c:	    				         REFLECTED,RP->state[2],
gbifur/gvecuntan.c:				      REFLECTED,
gbifur/gvecuntan.c:				      REFLECTED,
gbifur/gvecuntan.c:	    /* Find the incident and reflected curves 
gbifur/gvecuntan.c:	    		 	  1 reflected (4)	
gbifur/gvecuntan.c:	    			  3 reflected (0) */
gbifur/gvecuntan.c:	/* Join components together if there are no tracked reflected or
gbifur/gvecuntan.c:		trk[i] =  (node_type(n[i]) == B_REFLECT_NODE) ?
gbifur/gvecuntan.c:				track_scattered_wave(B_REFLECT_NODE,SHOCK_WAVE,
gbifur/gvecuntan.c:					REFLECTED,Rp_data(n[i])->state[1],
gbifur/gvecuntan.c:			(void) printf("deleting reflected curve, trk = (%s, %s), ",
gbifur/gvecuntan.c:		i_to_f_dir,B_REFLECT_NODE,Node_vel(n),front))
gbifur/gvecuntan.c:	if (!is_regular_reflection(Node_vel(n),front,RP))
gbifur/gvecuntan.c:		(void) printf("irregular boundary reflection\n");
gbifur/gvecuntan.c:	node_type(n) = B_REFLECT_NODE;
gbifur/gvecuntan.c:	set_status_at_node(cext,Opposite_orient(cp_orient),REFLECTED);
gbifur/gvecuntan.c:*	Join is set true if there are no tracked reflected or 
gbifur/gvecuntan.c:	track[1] = track_scattered_wave(CROSS_NODE,SHOCK_WAVE,REFLECTED,
gbifur/gvecuntan.c:	track[3] = track_scattered_wave(CROSS_NODE,SHOCK_WAVE,REFLECTED,
gdecs/gdecs.h:	RAMP_REFLECTION,
gdecs/gdecs.h:		/* Values for flag in bifurcation of reflections */
gdecs/gdecs.h:	NORMAL_TO_MACH_REFLECTION    = 1,
gdecs/gdecs.h:	NORMAL_TO_REGULAR_REFLECTION = 2,
gdecs/gdecs.h:	REGULAR_TO_MACH_REFLECTION   = 3
gdecs/gdecs.h:	/* Initialization info for the bubble in mach or regular reflections.
gdecs/gdecs.h:	   The states and angles around the reflection point are assumed to be
gdecs/gdecs.h:						/* reflection point */
gdecs/grdecs.h:	ANOMALOUS_REFLECTION,
gdecs/grdecs.h:	PRECURSOR_WITH_REFLECTED_RAREFACTION,
gdecs/grdecs.h:	PRECURSOR_WITH_REFLECTED_SHOCK
gdecs/guserint.h:	B_REFLECT_HSBDRY  = FIRST_PHYSICS_HSBDRY_TYPE,
gdecs/guserint.h:	B_REFLECT_NODE	  = B_REFLECT_HSBDRY,
gdecs/guserint.h:				/* 1 incident, 1 reflected, 2 boundaries */
gdecs/guserint.h:				/* 1 incident, 1 reflected, */
gdecs/guserint.h:				/* total internal reflection node */
gdecs/guserint.h:				   * total internal reflection, and possibly
gdecs/guserint.h:	NUM_PHYS_NODE_TYPES	= ONED_INTERACTION - B_REFLECT_NODE + 1,
gdecs/guserrp.h:	O_CURVE_FAMILY *_reflected1;
gdecs/guserrp.h:	O_CURVE_FAMILY *_reflected2;
geos/giniteos.c:            {"Obstacle (behind reflecting wall)", "O", 1,OBSTACLE_EOS},
ghyp/gvisc.c:	/* coords_ref is the location of the reflection of the point coords
ghyp/gvisc.c:	    || !reflect_pt_about_Nbdry(coords,coords_ref,nor,
ginit/gcauchy.c: * return NO for default reflection boundary state. */
ginit/gibifur.c:	    prompt_for_tracking_at_node(B_REFLECT_NODE,SHOCK_WAVE,
ginit/gibifur.c:					REFLECTED,swt);
ginit/gibifur.c:					REFLECTED,swt);
ginit/gibifur.c:	    prompt_for_tracking_at_node(MACH_NODE,SHOCK_WAVE,REFLECTED,swt);
ginit/gibifur.c:	    prompt_for_tracking_at_node(CROSS_NODE,SHOCK_WAVE,REFLECTED,swt);
ginit/gibifur.c:	    prompt_for_tracking_at_node(OVERTAKE_NODE,SHOCK_WAVE,REFLECTED,swt);
ginit/gibifur.c:					RAREF_LEADING_EDGE,REFLECTED,swt);
ginit/gibifur.c:					RAREF_TRAILING_EDGE,REFLECTED,swt);
ginit/gibifur.c:					REFLECTED,swt);
ginit/gibifur.c:				        RAREF_LEADING_EDGE,REFLECTED,swt);
ginit/gibifur.c:				        RAREF_TRAILING_EDGE,REFLECTED,swt);
ginit/gibifur.c:	    prompt_for_node_untrack(B_REFLECT_NODE,unf);
ginit/gibifur.c:	    unf[untrack_node_index(B_REFLECT_NODE)]._node_name =
ginit/gibifur.c:	        "regular reflection node";
ginit/gibifur.c:	    /* Boundary reflect nodes */
ginit/gibifur.c:	    i = scat_wv_index(B_REFLECT_NODE,SHOCK_WAVE,REFLECTED);
ginit/gibifur.c:	    swt[i].wave_name = "reflected shocks at regular reflections";
ginit/gibifur.c:	    i = scat_wv_index(ATTACHED_B_NODE,SHOCK_WAVE,REFLECTED);
ginit/gibifur.c:	        "reflected shocks at attached boundary reflection nodes";
ginit/gibifur.c:	    i = scat_wv_index(MACH_NODE,SHOCK_WAVE,REFLECTED);
ginit/gibifur.c:	    swt[i].wave_name = "reflected shocks at Mach reflections";
ginit/gibifur.c:	    swt[i].wave_name = "the Mach stem at Mach reflections";
ginit/gibifur.c:	    swt[i].wave_name = "the slip line at Mach reflections";
ginit/gibifur.c:	    i = scat_wv_index(CROSS_NODE,SHOCK_WAVE,REFLECTED);
ginit/gibifur.c:	    swt[i].wave_name = "reflected shocks at shock crossings";
ginit/gibifur.c:	    i = scat_wv_index(OVERTAKE_NODE,SHOCK_WAVE,REFLECTED);
ginit/gibifur.c:	    swt[i].wave_name = "reflected shocks at shock overtakes";
ginit/gibifur.c:	    i = scat_wv_index(OVERTAKE_NODE,RAREF_LEADING_EDGE,REFLECTED);
ginit/gibifur.c:	        "reflected rarefaction leading edges at shock overtakes";
ginit/gibifur.c:	    i = scat_wv_index(OVERTAKE_NODE,RAREF_TRAILING_EDGE,REFLECTED);
ginit/gibifur.c:	    swt[i].wave_name = "reflected rarefaction trailing edges "
ginit/gibifur.c:	    i = scat_wv_index(DIFFRACTION_NODE,SHOCK_WAVE,REFLECTED);
ginit/gibifur.c:	    swt[i].wave_name = "reflected shocks at shock-contact diffractions";
ginit/gibifur.c:	    i = scat_wv_index(DIFFRACTION_NODE,RAREF_LEADING_EDGE,REFLECTED);
ginit/gibifur.c:	    swt[i].wave_name = "reflected rarefaction leading edges "
ginit/gibifur.c:	    i = scat_wv_index(DIFFRACTION_NODE,RAREF_TRAILING_EDGE,REFLECTED);
ginit/gibifur.c:	    swt[i].wave_name = "reflected rarefaction trailing edges at "
ginit/gimkbub.c:            if(rect_boundary_type(intfc,0,0) == REFLECTION_BOUNDARY)
ginit/gimkbub.c:            if(rect_boundary_type(intfc,1,0) == REFLECTION_BOUNDARY)
ginit/gimkbub.c:            if(rect_boundary_type(intfc,0,0) == REFLECTION_BOUNDARY)
ginit/gimkbub.c:            if(rect_boundary_type(intfc,1,0) == REFLECTION_BOUNDARY)
ginit/gimkbub.c:                /** For conservative insertion, do not allow across reflection line **/
ginit/gimkbub.c:                if(rect_boundary_type(intfc,i,0) == REFLECTION_BOUNDARY)
ginit/gimkbub.c:                if(rect_boundary_type(intfc,i,1) == REFLECTION_BOUNDARY)
ginit/gimkbub.c:                if(rect_boundary_type(intfc,i,0) == REFLECTION_BOUNDARY)
ginit/gimkbub.c:                if(rect_boundary_type(intfc,i,1) == REFLECTION_BOUNDARY)
ginit/gimkbub.c:                /** For conservative insertion, do not allow across reflection line **/
ginit/gimkbub.c:                if(rect_boundary_type(intfc,i,0) == REFLECTION_BOUNDARY)
ginit/gimkbub.c:                if(rect_boundary_type(intfc,i,1) == REFLECTION_BOUNDARY)
ginit/ginit.h:	bool	_reflect_small_loop_shocks;
ginit/ginit.h:#define	reflect_small_loop_shocks(init)					\
ginit/ginit.h:	g_init_data(init)->_reflect_small_loop_shocks
ginit/ginitintfc.c:	    {{"a ramp reflection problem", "RR", 2,
ginit/ginitintfc.c:		{RAMP_REFLECTION}}, init_ramp_reflection},
ginit/ginitintfc.c:		case REFLECTION_BOUNDARY:
ginit/ginitintfc.c:		case REFLECTION_BOUNDARY:
ginit/ginitintfc.c:	        wave_type(*s) = (w_type == REFLECTION_BOUNDARY) ?
ginit/ginitintfc.c:		case REFLECTION_BOUNDARY:
ginit/ginitintfc.c:	case REFLECTION_BOUNDARY:
ginit/ginitintfc.c:	case REFLECTION_BOUNDARY:
ginit/ginitintfc.c:	    case REFLECTION_BOUNDARY:
ginit/ginitintfc.c:	    case REFLECTION_BOUNDARY:
ginit/ginitintfc.c:	rect_boundary_type(intfc,0,rside) = REFLECTION_BOUNDARY;
ginit/ginitintfc.c:	screen("Use Neumman (N) or reflecting (R) "
ginit/ginitintfc.c:	       (rect_boundary_type(intfc,0,1) == REFLECTION_BOUNDARY) ?
ginit/ginitintfc.c:	    rect_boundary_type(intfc,0,rside) = REFLECTION_BOUNDARY;
ginit/ginitintfc.c:	    /* Reflected shock from wall collides with contact */
ginit/ginitintfc.c:	    /* Reflected shock from wall collides with contact */
ginit/ginitintfc.c:	(void) printf("Reflected wave = %s, velocity = %g\n",
ginit/ginitintfc.c:	verbose_print_state("reflected mid state",state[4]);
ginit/ginitintfc.c:	(void) printf("Reflected wave at wall = %s, velocity = %g\n",
ginit/ginitintfc.c:	(void) printf("Reflected wave = %s, velocity = %g\n",
ginit/ginitintfc.c:	verbose_print_state("reflected mid state",state[5]);
ginit/ginitphys.c:	    {"Reflecting", "RE", 1, {REFLECTION_BOUNDARY}   },
ginit/ginitphys.c:	fuh->_reflect_state = g_reflect_state;
ginit/ginitphys.c:	if (reflect_small_loop_shocks(init) == YES)
ginit/ginitprotos.h:IMPORT	void	init_ramp_reflection(INIT_DATA*,INIT_PHYSICS*);
ginit/gipert.c:		rect_boundary_type(intfc,j,0) = REFLECTION_BOUNDARY;
ginit/gipert.c:		rect_boundary_type(intfc,j,1) = REFLECTION_BOUNDARY;
ginit/gipert.c:	       "\treflecting boundary state (F, default).\n"
ginit/giphysprompt.c:	    screen("Reflect small loop shocks (dflt = %s): ",
ginit/giphysprompt.c:		   y_or_n(reflect_small_loop_shocks(init)));
ginit/giphysprompt.c:		reflect_small_loop_shocks(init) = YES;        
ginit/giphysprompt.c:	reflect_small_loop_shocks(init) = NO;        
ginit/giphysprompt.c:	    reflect_small_loop_shocks(init) = YES;        
ginit/giphysprompt.c:	    reflect_small_loop_shocks(init) = NO;        
ginit/giprt.c:	case RAMP_REFLECTION:
ginit/giprt.c:	    stat_prompt = "the ramp reflection statistics";
ginit/giprt.c:	    fr_func = show_front_states_for_ramp_reflections;
ginit/giprt.c:*	to the corner position in a ramp reflection problem type.
ginit/giprt.c:*	WARNING: this diagnostic is specific to RAMP_REFLECTION problems.
ginit/giprt.c:	ramp_reflection_corner_posn(ss_origin,NO,dim);
ginit/girefl.c:*	reflection or refraction problems.
ginit/girefl.c:LOCAL	void	init_mach_reflection(int,CURVE*,Bubble*,Front*);
ginit/girefl.c:LOCAL	void	init_regular_reflection(float*,int,CURVE*,Bubble*,Front*);
ginit/girefl.c:*			init_ramp_reflection():
ginit/girefl.c:*	This function drives the initialization of reflections created by
ginit/girefl.c:*	result can be either a regular reflection, or a mach reflection,
ginit/girefl.c:EXPORT void init_ramp_reflection(
ginit/girefl.c:	    screen("ERROR in init_ramp_reflection(), "
ginit/girefl.c:	    screen("ERROR in init_ramp_reflection(), Invalid choice %s\n",s);
ginit/girefl.c:	ramp_reflection_corner_posn(corner_posn,YES,dim);
ginit/girefl.c:	    if (is_regular_reflection(node_v,front,bubble->RP))
ginit/girefl.c:	    	init_regular_reflection(node_v,prob_type,ramp,bubble,front);
ginit/girefl.c:	    	 * reg vs mach nodes.  Compare B_reflect_node_propagate()
ginit/girefl.c:	    	init_mach_reflection(prob_type,ramp,bubble,front);
ginit/girefl.c:	        (void) printf("not past corner in init_ramp_reflection()\n");
ginit/girefl.c:	*  a corner at the end opposite from the reflection. For example
ginit/girefl.c:}		/*end init_ramp_reflection*/
ginit/girefl.c:LOCAL	void init_regular_reflection(
ginit/girefl.c:	node_type(ns) = B_REFLECT_NODE;
ginit/girefl.c:		/* Make reflected shock */
ginit/girefl.c:	start_status(bow) = REFLECTED;
ginit/girefl.c:}		/*end init_regular_reflection*/
ginit/girefl.c:LOCAL void init_mach_reflection(
ginit/girefl.c:	screen("In Mach reflection case,\n\t");
ginit/girefl.c:				 bubble,NORMAL_TO_MACH_REFLECTION,front))
ginit/girefl.c:		screen("ERROR in init_mach_reflection(), ");
ginit/girefl.c:		/* Make reflected shock */
ginit/girefl.c:	start_status(bow) = REFLECTED;
ginit/girefl.c:}		/*end init_mach_reflection*/
ginit/girefl.c:*	anom_ang -- the transition point to anomolous reflection.  p1 has
ginit/girefl.c:*	detach_ang -- the theoretical limit of regular reflection.  At this
ginit/girefl.c:*		point, the reflected polar is tangent to the vertical axis
ginit/girefl.c:*		reflected polar intersects the vertical axis at the max
ginit/girefl.c:*	sonic_ang -- this where the reflected polar intersects the vertical
ginit/girefl.c:**	th01, th12, th13 -- the turning angles across the incident, reflecte
ginit/girefl.c:*	beta1 -- the angle between the reflected and its behind streamline
ginit/girefl.c:*	omega1 -- the angle between the reflected shock and its incoming
ginit/girefl.c:*	Omega1 -- the angle between the incident and reflected shocks
ginit/girefl.c:*	omega1p -- the angle between the reflected shock and the incoming
ginit/girefl.c:*	point to anomolous reflection, and corresponds to the point where
ginit/girefl.c:*	reflection is possible.  The point where it begins to fail is 
ginit/girefl.c:*	the point where the reflected polar is tangent to the verticle
ginit/girefl.c:*	at which the reflected polar crosses the verticle axis at the
ginit/girefl.c:*	and the pressure at the sonic point on the reflected polar.  The
ginit/girefl.c:*	regular reflection solution (cannot compute p2), the pressure at
ginit/girefl.c:*	regular reflection solution (cannot compute p2), the pressure at
ginit/girefl.c:*	of incidence corresponding to the onset of von Neumann reflection.
ginit/girefl.c:*	at the interesection of the incident and reflected polars) is 
ginit/girefl.c:*	greater (less) than the max pressure on the reflected polar.  The
ginit/girm_linear.c:	(void) printf("\nThe reflected wave is a %s wave.\n",
ginit/girm_linear.c:	           "Reflected wave is not a shock or a rarefaction!\n");
ginit/girm_linear.c:	    (void) printf("Reflected rarefaction strength = %g\n",strength);
ginit/girm_linear.c:	    (void) printf("Reflected shock strength = %g\n",strength);
ginit/girm_linear.c:	       "interval between the\n\ttransmitted and reflected waves "
ginit/girm_linear.c:	/* between cont. surf. and reflected shock or TE of rarefaction wave */
ginit/girm_linear.c:	/* on the reflected shock or trailing edge of RR */
ginit/girm_linear.c:*	on the transmitted and reflected waves don't intersect with the
ginit/girm_linear.c:	    /* Print pressure between contact and reflected shock 
ginit/girm_linear.c:	    /* Finally, print pressure ahead of reflected wave */
ginit/girm_linear.c:*	For the region between the transmitted shock and either a reflected
ginit/girm_linear.c:	/* Create reflected wave(s) */
ginit/gisc.c:	screen("Type 'y' if you wish to track the reflected wave");
ginit/gisc.c:	int		is_plus_orientation, track_reflected_wave;
ginit/gisc.c:	bool		is_reflected_shock;
ginit/gisc.c:	/* Is a reflected shock to be tracked? */
ginit/gisc.c:	screen("\nType n if the reflected wave is not to be tracked: ");
ginit/gisc.c:	track_reflected_wave = (s[0] == 'n' || s[0] == 'N') ? NO : YES;
ginit/gisc.c:	                                   NULL,&is_reflected_shock,
ginit/gisc.c:	case ANOMALOUS_REFLECTION:
ginit/gisc.c:	if (track_reflected_wave) 
ginit/gisc.c:	    if (!is_reflected_shock)
ginit/gisc.c:	    if (is_reflected_shock)
ginit/gisc.c:	if (track_reflected_wave) 
ginit/gisc.c:	    if (is_reflected_shock)
ginit/gisc.c:	    	/* Make reflected shock */
ginit/gisc.c:	                           w_type,REFLECTED,INCIDENT,0.0,front);
ginit/gisc.c:	                           w_type,REFLECTED,INCIDENT,0.0,front);
ginit/gisc.c:	        /* Make reflected rarefaction trailing edge */
ginit/gisc.c:	                           w_type,REFLECTED,INCIDENT,0.0,front);
ginit/spolars.c:LOCAL	int	init_plot_total_reflection(INIT_DATA*,INIT_PHYSICS*,float*,
ginit/spolars.c:LOCAL	int	init_plot_wall_reflection(INIT_DATA*,INIT_PHYSICS*,float*,
ginit/spolars.c:#pragma	noinline	init_plot_total_reflection
ginit/spolars.c:#pragma	noinline	init_plot_wall_reflection
ginit/spolars.c:	screen("\tRegular reflection (R, or %d)\n",i);
ginit/spolars.c:	screen("\tMach reflection (M, or %d)\n",i);
ginit/spolars.c:	screen("\tTotal internal reflection (TIR, or %d)\n",i);
ginit/spolars.c:	    num_polars = init_plot_wall_reflection(init,ip,theta0,state,
ginit/spolars.c:	    num_polars = init_plot_wall_reflection(init,ip,theta0,state,
ginit/spolars.c:	    num_polars = init_plot_total_reflection(init,ip,theta0,state,npts,
ginit/spolars.c:LOCAL int init_plot_wall_reflection(
ginit/spolars.c:	    init_plot_wall_reflection) 
ginit/spolars.c:	set_plot_control("reflected ",&plot_control[1],&npts[1],state[1]);
ginit/spolars.c:}		/*end init_plot_wall_reflection*/
ginit/spolars.c:	bool		         is_reflected_shock;
ginit/spolars.c:	    REFLECTED_ANGLE_1,
ginit/spolars.c:	    REFLECTED_ANGLE_2,
ginit/spolars.c:	    REFLECTED_ANGLE_3,
ginit/spolars.c:	    REFLECTED_TURN_ANGLE,
ginit/spolars.c:	    BEHIND_REFLECTED_MACH_NUMBER,
ginit/spolars.c:	header[REFLECTED_ANGLE_1] = "REFLECTED_ANGLE_1";
ginit/spolars.c:	header[REFLECTED_ANGLE_2] = "REFLECTED_ANGLE_2";
ginit/spolars.c:	header[REFLECTED_ANGLE_3] = "REFLECTED_ANGLE_3";
ginit/spolars.c:	header[REFLECTED_TURN_ANGLE] = "REFLECTED_TURN_ANGLE";
ginit/spolars.c:	header[BEHIND_REFLECTED_MACH_NUMBER] = "BEHIND_REFLECTED_MACH_NUMBER";
ginit/spolars.c:	    vars[0][REFLECTED_ANGLE_1]		= 180.0;
ginit/spolars.c:	    vars[0][REFLECTED_ANGLE_2]		= 180.0;
ginit/spolars.c:	    vars[0][REFLECTED_ANGLE_3]		= 180.0;
ginit/spolars.c:	    vars[0][REFLECTED_TURN_ANGLE]	= 0.0;
ginit/spolars.c:				                &is_reflected_shock,NULL,
ginit/spolars.c:	    vars[1][REFLECTED_ANGLE_1] = degrees(RP->ang[1]);
ginit/spolars.c:	    vars[1][REFLECTED_ANGLE_2] = degrees(RP->ang[2]);
ginit/spolars.c:	    vars[1][REFLECTED_ANGLE_3] = degrees(RP->ang[3]);
ginit/spolars.c:	    vars[1][REFLECTED_TURN_ANGLE] = degrees(RP->theta[2]);
ginit/spolars.c:	    			 vars[1][REFLECTED_TURN_ANGLE] -
ginit/spolars.c:	    vars[0][BEHIND_REFLECTED_MACH_NUMBER]	=
ginit/spolars.c:	    vars[1][BEHIND_REFLECTED_MACH_NUMBER]	= RP->M[4];
ginit/spolars.c:				            &is_reflected_shock,NULL,
ginit/spolars.c:	    vars[nlines][REFLECTED_ANGLE_1] = degrees(RP->ang[1]);
ginit/spolars.c:	    vars[nlines][REFLECTED_ANGLE_2] = degrees(RP->ang[2]);
ginit/spolars.c:	    vars[nlines][REFLECTED_ANGLE_3] = degrees(RP->ang[3]);
ginit/spolars.c:	    vars[nlines][REFLECTED_TURN_ANGLE] = degrees(RP->theta[2]);
ginit/spolars.c:	    			      vars[nlines][REFLECTED_TURN_ANGLE] -
ginit/spolars.c:	    vars[nlines][BEHIND_REFLECTED_MACH_NUMBER]	 = RP->M[4];
ginit/spolars.c:	set_plot_control("reflected ",&plot_control[2],&npts[2],state[2]);
ginit/spolars.c:	    bool is_reflected_shock;
ginit/spolars.c:	                                    &is_reflected_shock,NULL,
ginit/spolars.c:LOCAL	int init_plot_total_reflection(
ginit/spolars.c:	    init_plot_total_reflection) 
ginit/spolars.c:	set_plot_control("reflected ",&plot_control[1],&npts[1],state[1]);
ginit/spolars.c:}		/*end init_plot_total_reflection*/
ginit/spolars.c:	set_plot_control("reflected ",&plot_control[2],&npts[2],state[2]);
ginit/spolars.c:	bool           is_reflected_shock;
ginit/spolars.c:	                                &is_reflected_shock,NULL,
ginit/spolars.c:	                                    &is_reflected_shock,NULL,
ginit/spolars.c:	bool  is_reflected_shock;
ginit/spolars.c:	                                 &is_reflected_shock,NULL,
ginit/spolars.c:	                                    &is_reflected_shock,NULL,
ginit/spolars.c:	set_plot_control("first reflected ",&plot_control[1],&npts[1],state[1]);
ginit/spolars.c:	set_plot_control("second reflected ",&plot_control[2],&npts[2],state[2]);
ginit/spolars.c:	verbose_print_state("\nState behind reflected wave",sml);
ginit/spolars.c:	(void) printf("Reflected wave mass flux = %g\n",ml);
ginit/spolars.c:	(void) printf("Reflected wave pressure ratio = %g\n",pml/pressure(sl));
ginit/spolars.c:	(void) printf("Reflected wave density ratio = %g\n",Dens(sml)/Dens(sl));
ginit/spolars.c:	    (void) printf("Reflected shock mass flux = %g\n",ml);
ginit/spolars.c:	    (void) printf("Reflected shock velocity = %g\n",s);
ginit/spolars.c:	    (void) printf("Reflected shock ahead Mach number = %g\n",
ginit/spolars.c:	    (void) printf("Reflected shock behind Mach number = %g\n",
ginit/spolars.c:	    (void) printf("Reflected rarefaction ");
ginit/spolars.c:	    (void) printf("Reflected rarefaction ");
ginit/spolars.c:	    (void) printf("Reflected rarefaction ");
ginit/spolars.c:	    (void) printf("Reflected rarefaction ");
gintfc/gintfcprotos.h:IMPORT	void	g_reflect_point(POINT*,float*,float*,INTERFACE*);
gintfc/gintfcprotos.h:IMPORT	void	g_reflect_node2d(NODE*,float*,float*);
gintfc/gtop.c:EXPORT	void	g_reflect_point(
gintfc/gtop.c:	POINT		*point,/* point being reflected */
gintfc/gtop.c:	float		*p,	/* point on reflection plane */
gintfc/gtop.c:	INTERFACE	*intfc)	/* interface being reflected */
gintfc/gtop.c:	f_reflect_point(point,p,n,intfc);
gintfc/gtop.c:}		/*end f_reflect_point*/
gintfc/gtop.c:EXPORT  void    g_reflect_node2d(
gintfc/gtop.c:	NODE		*node,/* node being reflected */
gintfc/gtop.c:	float		*p,     /* point on reflection plane */
gintfc/gtop.c:	    	reflect_state(RP->state[i],intfc,pt,p,n);
gintfc/gtop.c:	f_reflect_node2d(node,p,n);
gintfc/gtop.c:}		/*end g_reflect_node2d*/
gintfc/guserintfc.c:	    iuh->_reflect_point = g_reflect_point;
gintfc/guserintfc.c:	    iuh->_reflect_node = g_reflect_node2d;
gintfc/guserintfc.c:	* The following code is designed to enforce reflection symmetry
gintfc/guserintfc.c:	* along a curve crossing a reflection boundary.  First points
gintfc/guserintfc.c:	        if ((rect_boundary_type(intfc,i,d) == REFLECTION_BOUNDARY))
gintfc/guserintfc.c:        * The following code is designed to enforce reflection symmetry
gintfc/guserintfc.c:        * along a curve crossing a reflection boundary.  First points
gintfc/guserintfc.c:                if ((rect_boundary_type(intfc,i,d) == REFLECTION_BOUNDARY))
gnode/ganom.c:*	Transforms a regular shock diffraction into anomalous reflection
gnode/ganom.c:EXPORT	int anomalous_reflection_propagate(
gnode/ganom.c:	debug_print("diffraction","Entered anomalous_reflection_propagate()\n");
gnode/ganom.c:	    (void) printf("WARNING in anomalous_reflection_propagate(), "
gnode/ganom.c:		/*     Propagate leading edge of reflected	*/
gnode/ganom.c:	    	    (void) printf("Unable to propagate anomalous reflection\n");
gnode/ganom.c:	    	    (void) printf("Unable to propagate anomalous reflection\n");
gnode/ganom.c:	    (void) printf("WARNING in anomalous_reflection_propagate(), "
gnode/ganom.c:	    (void) printf("WARNING in anomalous_reflection_propagate(), "
gnode/ganom.c:	debug_print("diffraction","Left anomalous_reflection_propagate(), ");
gnode/ganom.c:}		/*end anomalous_reflection_propagate*/
gnode/ganom.c:			(void) printf("anomalous reflection\n");
gnode/gbnode.c:LOCAL	void	modify_regular_reflection_node(POINT*,BOND*,BOND*,O_CURVE*,
gnode/gbnode.c:*			B_reflect_node_propagate():
gnode/gbnode.c:*            reflected (ang 2)    RP->state[1]    incident (ang 1)
gnode/gbnode.c:*	A B_reflect_node is a node at which four curves meet, with two
gnode/gbnode.c:EXPORT	int B_reflect_node_propagate(
gnode/gbnode.c:	O_CURVE		Oldcref, Newcref;	/* reflected curve */
gnode/gbnode.c:	debug_print("B_reflect_node","Entered B_reflect_node_propagate()\n");
gnode/gbnode.c:	if (debugging("B_reflect_node")) 
gnode/gbnode.c:	             B_reflect_node_propagate)
gnode/gbnode.c:	find_curve_with_status(oldn,&Oldcref.curve,&Oldcref.orient,REFLECTED);
gnode/gbnode.c:	             B_reflect_node_propagate)
gnode/gbnode.c:		         B_reflect_node_propagate)
gnode/gbnode.c:		         B_reflect_node_propagate)
gnode/gbnode.c:	if (debugging("B_reflect_node")) 
gnode/gbnode.c:	    (void) printf("\t\tOLD REFLECTED CURVE:\n");
gnode/gbnode.c:	    debug_print("B_reflect_node","oldn propagates into ahead comp\n");
gnode/gbnode.c:	    debug_print("B_reflect_node","oldn propagates out of ahead comp\n");
gnode/gbnode.c:	    debug_print("B_reflect_node","Left B_reflect_node_propagate()\n");
gnode/gbnode.c:	    	debug_print("B_reflect_node","Left B_reflect_node_propagate()\n");
gnode/gbnode.c:	if (is_regular_reflection(node_v,fr,RP)) 
gnode/gbnode.c:	    debug_print("B_reflect_node","\t\tREGULAR REFLECTION:\n");
gnode/gbnode.c:	    modify_regular_reflection_node(&Pc,crossbinc,crossbahead,
gnode/gbnode.c:	    debug_print("B_reflect_node","Left B_reflect_node_propagate()\n");
gnode/gbnode.c:	    g_reflect_node_bifurcation(fr,wave,&Oldcinc,&Newcinc,&Oldcref,
gnode/gbnode.c:	    debug_print("B_reflect_node","Left B_reflect_node_propagate()\n");
gnode/gbnode.c:	if (debugging("B_reflect_node")) 
gnode/gbnode.c:	    (void) printf("\t\tNEW REFLECTED CURVE:\n");
gnode/gbnode.c:	debug_print("B_reflect_node","Left B_reflect_node_propagate()\n");
gnode/gbnode.c:}		/*end B_reflect_node_propagate*/
gnode/gbnode.c:*               modify_regular_reflection_node():
gnode/gbnode.c:*	Given that a regular reflection node remains such a node under
gnode/gbnode.c:*	also propagated with respect to the reflected curve newcref, becoming
gnode/gbnode.c:*	a new point is inserted to maintain the correct reflection angle.
gnode/gbnode.c:LOCAL void modify_regular_reflection_node(
gnode/gbnode.c:	debug_print("B_reflect_node","Entered modify_regular_reflection_node()\n");
gnode/gbnode.c:	debug_print("B_reflect_node","Calling shift_node_past\n");
gnode/gbnode.c:		/* Propagate old reflected curve and insert a new bond */
gnode/gbnode.c:	debug_print("B_reflect_node","Left modify_regular_reflection_node()\n");
gnode/gbnode.c:}		/*end modify_regular_reflection_node*/
gnode/gbnode.c:*			is_regular_reflection():
gnode/gbnode.c:*       This is shock boundary regular reflection.   The input and output
gnode/gbnode.c:*       incident shock, while st_2 is the state behind reflected shock.
gnode/gbnode.c:*       Returns 1 if there is a regular reflection according to the von Neumann
gnode/gbnode.c:*	criteria, and 0 otherwise. If there is a regular reflection, then the
gnode/gbnode.c:*	state st_2 behind the reflected shock is computed and placed in
gnode/gbnode.c:EXPORT int is_regular_reflection(
gnode/gbnode.c:	debug_print("B_reflect_node","Entered is_regular_reflection()\n");
gnode/gbnode.c:	if (debugging("B_reflect_node"))
gnode/gbnode.c:	if (debugging("B_reflect_node"))
gnode/gbnode.c:		debug_print("B_reflect_node","Left is_regular_reflection()\n");
gnode/gbnode.c:		debug_print("B_reflect_node","Left is_regular_reflection()\n");
gnode/gbnode.c:	if (debugging("B_reflect_node"))
gnode/gbnode.c:		print_RP_node_states("States after is_regular_reflection()",
gnode/gbnode.c:			nod_v,RP,B_REFLECT_NODE);
gnode/gbnode.c:	debug_print("B_reflect_node","Left is_regular_reflection()\n");
gnode/gbnode.c:}		/*end is_regular_reflection*/
gnode/gmdnode.c:*	diffraction.  We simply untrack the reflected wave, and change the
gnode/gmdnode.c:		/* untrack any reflected waves */
gnode/gnode.c:*	B_REFLECT_NODE	boundary curve		FIXED	no change of type or
gnode/gnode.c:*			shock			REFLECTED
gnode/gnode.c:*			shock			REFLECTED
gnode/gnode.c:*			shock			REFLECTED
gnode/gnode.c:*			shock			REFLECTED
gnode/gnode.c:*			shock			REFLECTED
gnode/gnode.c:*			rarefaction		REFLECTED
gnode/gnode.c:*			rarefaction		REFLECTED
gnode/gnode.c:*			rarefaction		REFLECTED
gnode/gnode.c:*			rarefaction		REFLECTED
gnode/gnode.c:	case B_REFLECT_NODE:
gnode/gnode.c:	    status = B_reflect_node_propagate(fr,wave,oldn,newn,rp,
gnode/gnodeprotos.h:IMPORT	int	B_reflect_node_propagate(Front*,Wave*,NODE*,NODE*,RPROBLEM**,
gnode/gnodeprotos.h:IMPORT	int	is_regular_reflection(float*,Front*,RP_DATA*);
gnode/gnodeprotos.h:IMPORT	void	ramp_reflection_corner_posn(float*,int,int);
gnode/gnodeprotos.h:IMPORT	int	anomalous_reflection_propagate(Front*,Wave*,NODE*,NODE*,
gnode/gnodesub.c:*		modify_regular_reflection_node()
gnode/gpcnode.c:*	this node is not tracked.  One reflected wave at this cross node
gnode/gpcnode.c:*	internal reflection node with the old reflected rarefaction.  It is
gnode/gpcnode.c:*	RP and newc.  The other reflected wave at the cross node is directed
gnode/gpcnode.c:*	toward the reflected rarefaction, creating an overtake node on the
gnode/gpcnode.c:	                  "precursor with reflected rarefaction.\n");
gnode/gpcnode.c:	if (track_scattered_wave(CROSS_NODE,SHOCK_WAVE,REFLECTED,
gnode/gpcnode.c:*	bifurcation from regular diffraction to precursor with reflected
gnode/gpcnode.c:*	found (if possible).  Then the total internal reflection is set
gnode/gpcnode.c:*	reflection node.  The newn corresponding to the old diffraction node
gnode/gpcnode.c:*	becomes the new total internal reflection node.  After this function
gnode/gpcnode.c:*	the RP structure will reflect this, but not newc.  The old incident
gnode/gpcnode.c:*	total internal reflection doesn't really exist at this point.
gnode/gpcnode.c:		(void) printf("Total internal reflection node position <%g, %g>\n",
gnode/gpcnode.c:*	at the total internal reflection, completing the transfer from
gnode/gpcnode.c:*	diffraction node to total internal reflection.
gnode/gpcnode.c:		end_status(curves[0]) = REFLECTED;
gnode/gpcnode.c:		start_status(curves[1]) = REFLECTED;
gnode/gpcnode.c:*	any reflected waves are deleted (untracked) and components reset if
gnode/gpcnode.c:*	This function installs the tracking for the reflected shock at the
gnode/gpcnode.c:*	cross node which is directed toward the reflected rarefaction.  This
gnode/gpcnode.c:*	creates an overtake node at the point of intersection of the reflected
gnode/gpcnode.c:		start_status(final_shock) = REFLECTED;
gnode/gpcnode.c:		end_status(final_shock) = REFLECTED;
gnode/gpcnode.c:*	at the total internal reflection in the direction of the Mach angle
gnode/gscnode.c:*	reflected wave
gnode/gscnode.c:* (curve 2)  reflected shock or
gnode/gscnode.c:* trailing edges of reflected 
gnode/gscnode.c:*	contact discontinuity producing a transmitted shock and a reflected
gnode/gscnode.c:	bool		is_plus_or, is_reflected_shock;
gnode/gscnode.c:		"REFLECTED RAREFACTION LEADING EDGE",
gnode/gscnode.c:		"REFLECTED SHOCK CURVE",
gnode/gscnode.c:		"REFLECTED RAREFACTION TRAILING EDGE",
gnode/gscnode.c:	                                         &is_reflected_shock,fr,
gnode/gscnode.c:	    case ANOMALOUS_REFLECTION:
gnode/gscnode.c:	            is_reflected_shock = NO;
gnode/gscnode.c:	            status = anomalous_reflection_propagate(fr,wave,oldn,newn,
gnode/gscnode.c:	                                                    &is_reflected_shock,
gnode/gscnode.c:	    case PRECURSOR_WITH_REFLECTED_RAREFACTION:
gnode/gscnode.c:	    case PRECURSOR_WITH_REFLECTED_SHOCK:
gnode/gscnode.c:	               "PRECURSOR_WITH_REFLECTED_SHOCK CODE NEEDED\n");
gnode/gscnode.c:	    if ((oldc[2]->curve != NULL && !is_reflected_shock) ||
gnode/gscnode.c:	        (is_reflected_shock && (oldc[1]->curve || oldc[3]->curve)))
gnode/gscnode.c:		    (void) sprintf(s,"The reflected wave has undergone a "
gnode/gscnode.c:	                             "is_reflected_shock = %s\n",
gnode/gscnode.c:	                             y_or_n(is_reflected_shock));
gnode/gscnode.c:				"REFLECTED RAREFACTION LEADING EDGE",
gnode/gscnode.c:				"REFLECTED SHOCK CURVE",
gnode/gscnode.c:				"REFLECTED RAREFACTION TRAILING EDGE",
gnode/gscnode.c:*	in the rp for a deprecursion from precursor with reflected 
gnode/gscnode.c:*	the deprecursion is being triggered by a total internal reflection
gnode/gscnode.c:	identify_curves_with_status(node,&Ctmp0,&Ctmp1,REFLECTED);
gnode/gscnode.c:	    /* Untracked reflected wave */
gnode/gscnode.c:	    /* Tracked reflected shock */
gnode/gscnode.c:	    /* Tracked reflected rarefaction */
gnode/gscnode.c:		/* disconnect reflected waves from node */
gnode/gscnode.c:		    (void) printf("Unable to propagate reflected "
gnode/gscnode.c:*	contact discontinuity producing a transmitted shock but no reflected
gnode/gscnode.c:*	reflected rarefaction back to a single diffraction node.  We allow
gnode/gscnode.c:		/* Identify cross and total internal reflection nodes */
gnode/gscnode.c:	identify_curves_with_status(new_irnode,&Ctmp0,&Ctmp1,REFLECTED);
gnode/gscnode.c:	identify_curves_with_status(old_irnode,&Ctmp0,&Ctmp1,REFLECTED);
gnode/gscnsts.c:*	either a reflected shock or rarefaction wave. The projection is defined
gnode/gscnsts.c:	bool		*is_reflected_shock,
gnode/gscnsts.c:	    *is_reflected_shock = NO;
gnode/gscnsts.c:	    				       is_reflected_shock,
gnode/gscnsts.c:							is_reflected_shock,
gnode/gscnsts.c:	    *is_reflected_shock = (pstar >= p1) ? YES : NO;
gnode/gscnsts.c:	    	     * value of *is_reflected_shock.  We use the
gnode/gscnsts.c:	    		*is_reflected_shock = NO;
gnode/gscnsts.c:	    	status = (*is_reflected_shock) ?
gnode/gscnsts.c:	    		    PRECURSOR_WITH_REFLECTED_SHOCK :
gnode/gscnsts.c:	    		    PRECURSOR_WITH_REFLECTED_RAREFACTION;
gnode/gscnsts.c:	/* Find state behind reflected wave */
gnode/gscnsts.c:	if (*is_reflected_shock) 
gnode/gscnsts.c:	{	 /* reflected wave is a shock */
gnode/gscnsts.c:	    	                 "s_polar_3() failed for reflected shock","\n",
gnode/gscnsts.c:	    /* Check for sonic transition at reflected shock */
gnode/gscnsts.c:	{	/* Reflected wave is a rarefaction */
gnode/gscnsts.c:				 "for reflected rarefaction","\n",flag);
gnode/gscnsts.c:	    (void) printf("is_reflected_shock = %s\n",
gnode/gscnsts.c:	    			(*is_reflected_shock) ? "YES" : "NO");
gnode/gscnsts.c:	/* check for possible anomalous reflection	*/
gnode/gscnsts.c:	    	/* Previous reflected wave was a shock */
gnode/gscnsts.c:	    	/* Previous reflected wave was a rarefaction */
gnode/gscnsts.c:	    	status = ANOMALOUS_REFLECTION;
gnode/gscnsts.c:	    status = ANOMALOUS_REFLECTION;
gnode/gscnsts.c:	    status = ANOMALOUS_REFLECTION;
gnode/gscnsts.c:	    status = REGULAR_TO_MACH_REFLECTION;
gnode/gscnsts.c:	    status = REGULAR_TO_MACH_REFLECTION;
gnode/gscnsts.c:	bool	  *is_reflected_shock,
gnode/gscnsts.c:	*is_reflected_shock = NO;
gnode/gscnsts.c:	/* Reflected wave is a rarefaction */
gnode/gscnsts.c:	    (void) printf("is_reflected_shock = %s\n",
gnode/gscnsts.c:	    		  (*is_reflected_shock) ? "YES" : "NO");
gnode/gscnsts.c:	bool	                 *is_reflected_shock,
gnode/gscnsts.c:	if (is_small_inc_ang_reg_diff_node(t,RP,is_reflected_shock,flag))
gnode/gscnsts.c:	    *is_reflected_shock = NO;
gnode/gscnsts.c:	    status = PRECURSOR_WITH_REFLECTED_RAREFACTION;
gnode/gscnsts.c:	    *is_reflected_shock = NO;
gnode/gscnsts.c:	    *is_reflected_shock = NO;
gnode/gscnsts.c:	    *is_reflected_shock = YES;
gnode/gscnsts.c:	    status = (p6mta < p1mta) ? PRECURSOR_WITH_REFLECTED_SHOCK :
gnode/gscnsts.c:	    	/* Previous reflected wave was a shock */
gnode/gscnsts.c:	    	status = PRECURSOR_WITH_REFLECTED_SHOCK;
gnode/gscnsts.c:	    	*is_reflected_shock = YES;
gnode/gscnsts.c:	    	/* Previous reflected wave was a rarefaction */
gnode/gscnsts.c:	    	status = PRECURSOR_WITH_REFLECTED_RAREFACTION;
gnode/gscnsts.c:	    	*is_reflected_shock = NO;
gnode/gscnsts.c:	*is_reflected_shock = NO;
gnode/gscnsts.c:	bool	  *is_reflected_shock,
gnode/gscnsts.c:		/* Find state behind reflected wave */
gnode/gscnsts.c:	    *is_reflected_shock = YES;
gnode/gscnsts.c:	    *is_reflected_shock = NO;
gnode/gscnsts.c:	    if (*is_reflected_shock) 
gnode/gscnsts.c:	    	print_angle("Reflected angle = ", RP->ang[2],"\n");
gnode/gscnsts.c:* 	transmitted shock but no reflected shock.
gnode/gsndnode.c:LOCAL	bool	snd_B_reflect_node_propagate(Front*,Front*,
gnode/gsndnode.c:	        case B_REFLECT_NODE:
gnode/gsndnode.c:	            ok = snd_B_reflect_node_propagate(fr,newfr,oldn,
gnode/gsndnode.c:LOCAL bool snd_B_reflect_node_propagate(
gnode/gsndnode.c:}		/*end snd_B_reflect_node_propagate*/
gnode/gsndnode.c:*	handled.  This is the case of no reflected waves (either
gnode/gssnode.c:	O_CURVE		Oldcrefl, Newcrefl;	/* the reflected curve */
gnode/gssnode.c:	    ((start_status(Oldcrefl.curve) != REFLECTED) &&
gnode/gssnode.c:	     (end_status(Oldcrefl.curve) != REFLECTED)))
gnode/gssnode.c:	    /* If the reflected wave is untracked, the adjacent curve
gnode/gssnode.c:*	angles of the reflected shock and slip line. After shifting the node
gnode/gssnode.c:*	turn to the old reflected shock and slipline. Since their node has been
gnode/gssnode.c:*	to maintain the correct reflection angle and assign their states.
gnode/gssnode.c:	    /* Find states on reflected curve */
gnode/gssnode.c:	    /* Propagate old reflected curve and insert a new bond */
gnode/gssnode.c:*			ramp_reflection_corner_posn():
gnode/gssnode.c:*	gets called in a problem other than a ramp reflection.
gnode/gssnode.c:EXPORT	void	ramp_reflection_corner_posn(
gnode/gssnode.c:}		/*end ramp_reflection_corner_posn*/
gnode/gssnode.c:*	caused by an increase in the incident angle, in which case the reflected
gnode/gssnode.c:*	The reflected shock and slip line are untracked, with the incident
gnode/gssnode.c:	    (void) printf("\t\t%s REFLECTED CURVE:\n",string);
gnode/gssnode.c:	    (void) printf("\n\t\t%s REFLECTED wave untracked\n",string);
gnode/gssnode.c:*               reflected wave \  state1   / 
gnode/gssnode.c:*	        reflected shock/ 	   \incident shock
gnode/gssnode.c:*	producing two reflected waves. 
gnode/gssnode.c:	                  "The reflected shock has disappeared, it "
gnode/gssnode.c:*	in the rp for a deprecursion from precursor with reflected 
gnode/gssnode.c:	    (void) printf("OLD REFLECTED CURVE 1:\n");
gnode/gssnode.c:	    (void) printf("OLD REFLECTED CURVE 3:\n");
gnode/gssnode.c:	/* Are any reflected shocks tracked? */
gnode/gssnode.c:	identify_curves_with_status(oldn,oldc[1],oldc[3],REFLECTED);
gnode/gssnode.c:	/* The reflected curves may need to be switched */
gnode/gssnode.c:*	reflected wave
gnode/gssnode.c:* (curve 3)  reflected shock or
gnode/gssnode.c:* trailing edges of reflected 
gnode/gssnode.c:*	reflected wave.
gnode/gssnode.c:	                  "REFLECTED RAREFACTION LEADING EDGE",
gnode/gssnode.c:	                  "REFLECTED SHOCK CURVE",
gnode/gssnode.c:	                  "REFLECTED RAREFACTION TRAILING EDGE",
gnode/gssnode.c:	identify_curves_with_status(oldn,&Ctmp0,&Ctmp1,REFLECTED);
gnode/gssnode.c:	    /* Untracked reflected wave */
gnode/gssnode.c:	    /* Tracked reflected shock */
gnode/gssnode.c:	    /* Tracked reflected rarefaction */
gnode/gssnode.c:	    /* untrack reflected waves and slip */
gnode/gssnsts.c:*	of the reflected shock, the contact and the Mach stem are then 
gnode/gssnsts.c:*	degenerate mach node (one in which the reflected shock is very
gnode/gssnsts.c:*	colliding producing reflected shocks and a contact discontinuity.
gnode/gssnsts.c:	/* Find states behind reflected waves */
gnode/gssnsts.c:	/* Find states behind reflected waves */
gprop/gprop.c:EXPORT	void	reflect_wssten(
gprop/gprop.c:		reflect_state(sl[i],front->interf,rcrds[i],pt,nor);
gprop/gprop.c:		reflect_state(sr[i],front->interf,lcrds[i],pt,nor);
gprop/gprop.c:	    screen("ERROR in reflect_wssten(), "
gprop/gprop.c:}		/*end reflect_wssten*/
gprop/gpropprotos.h:IMPORT	void	reflect_wssten(WSSten*,SIDE,Front*);
gprop/gwspeed.c:	       "of a reflection\n"
gprop/gwspeed.c:	    reflect_wssten(wssten,side,front);
gprt/g2dprint.c:*		show_front_states_for_ramp_reflections():
gprt/g2dprint.c:*	For use in the RAMP_REFLECTION problem.  We assume the physical
gprt/g2dprint.c:*	curves are all oriented positively wrt the reflection point.
gprt/g2dprint.c:EXPORT void show_front_states_for_ramp_reflections(
gprt/g2dprint.c:	ramp_reflection_corner_posn(ss_origin,NO,dim);
gprt/g2dprint.c:				(node_type(*n) == B_REFLECT_NODE))
gprt/g2dprint.c:			/* orient now mach stem wrt reflection point */
gprt/g2dprint.c:		    show_front_states_for_ramp_reflections)
gprt/g2dprint.c:	find_curve_with_status(*n,&bow_shock,&orient,REFLECTED);
gprt/g2dprint.c:			/* orient now bow shock wrt reflection point */
gprt/g2dprint.c:	else	/* (node_type(n) == B_REFLECT_NODE) */
gprt/g2dprint.c:		/* orient is now refl wall wrt reflection point */
gprt/g2dprint.c:			     show_front_states_for_ramp_reflections)
gprt/g2dprint.c:		       "along REFLECTED SHOCK (%llu)",curve_number(bow_shock));
gprt/g2dprint.c:		    (void) sprintf(title,"along REFLECTION WALLS (%llu) ",
gprt/g2dprint.c:		(void) sprintf(title,"along REFLECTION WALLS (%llu) and (%llu)",
gprt/g2dprint.c:		(void) sprintf(title,"along REFLECTION WALLS (%llu) and (%llu)",
gprt/g2dprint.c:			"along REFLECTION WALLS (%llu), (%llu), and (%llu)",
gprt/g2dprint.c:		    show_front_states_for_ramp_reflections)
gprt/g2dprint.c:}		/*end show_front_states_for_ramp_reflection*/
gprt/g2dprint.c:		(void) printf("\t\tREFLECTED CURVE 1:\n");
gprt/g2dprint.c:		(void) printf("\t\tREFLECTED CURVE 1 IS UNTRACKED:\n\n");
gprt/g2dprint.c:		(void) printf("\t\tREFLECTED CURVE 3:\n");
gprt/g2dprint.c:		(void) printf("\t\tREFLECTED CURVE 3 IS UNTRACKED:\n\n");
gprt/g2dprint.c:*	This function is simply a driver for the output for ramp reflection
gprt/g2dprint.c:*	runs.  We just fork here for Mach or regular reflection.
gprt/g2dprint.c:		if (node_type(*n) == B_REFLECT_NODE)
gprt/g2dprint.c:			       &newcrefl->orient,REFLECTED);
gprt/g2dprint.c:	ramp_reflection_corner_posn(corner_posn,NO,dim);
gprt/g2dprint.c:	    /* The timestep test is to give the reflected wave a chance
gprt/g2dprint.c:	     * avoid a kink that appears in the last bond of the reflected
gprt/g2dprint.c:	     * to identify a "valid" kink in the reflected shock.
gprt/g2dprint.c:*	This function provides the diagnostic output for regular reflection
gprt/g2dprint.c:			       &newcrefl->orient,REFLECTED);
gprt/g2dprint.c:	ramp_reflection_corner_posn(corner_posn,NO,dim);
gprt/gprcur.c:*	Currently, this function is only used for RAMP_REFLECTION problems.
gprt/gprcur.c:	    /* mark endpoint (for corner point in reflection problems) */
gprt/gprstate.c:	case B_REFLECT_HSBDRY:
gprt/gprstate.c:	    (void) fprintf(file,"B_REFLECT_%s",suffix);
gprt/gprstate.c:	    hsbdry_type = B_REFLECT_HSBDRY;
gprt/gprstate.c:	case ANOMALOUS_REFLECTION:
gprt/gprstate.c:	    (void) printf("%sANOMALOUS_REFLECTION\n",message);
gprt/gprstate.c:	case PRECURSOR_WITH_REFLECTED_RAREFACTION:
gprt/gprstate.c:	    (void) printf("%sPRECURSOR_WITH_REFLECTED_RAREFACTION\n",message);
gprt/gprstate.c:	case PRECURSOR_WITH_REFLECTED_SHOCK:
gprt/gprstate.c:	    (void) printf("%sPRECURSOR_WITH_REFLECTED_SHOCK\n",message);
gprt/gprstate.c:	case REFLECTED:
gprt/gprstate.c:	    (void) fprintf(file,"%sREFLECTED\n",message);
gprt/gprstate.c:	    status = REFLECTED;
gprt/gprstate.c:	    ang_names[1] = "Reflected rarefaction leading edge angle";
gprt/gprstate.c:	    ang_names[2] = "Reflected shock angle";
gprt/gprstate.c:	    ang_names[3] = "Reflected rarefaction trailing edge angle";
gprt/gprstate.c:	    		"Reflected rarefaction leading edge angle";
gprt/gprstate.c:	    		"Reflected rarefaction trailing edge angle";
gprt/gprstate.c:	    	ang_names[4] = "Reflected shock angle";
gprt/gprstate.c:	    	ang_names[1] = "Reflected shock angle";
gprt/gprstate.c:	    		"Reflected rarefaction leading edge angle";
gprt/gprstate.c:	    		"Reflected rarefaction trailing edge angle";
gprt/gprstate.c:	    	ang_names[1] = "Reflected shock angle";
gprt/gprstate.c:	    	ang_names[3] = "Reflected shock angle";
gprt/gprstate.c:	    		"Reflected rarefaction leading edge angle";
gprt/gprstate.c:	    		"Reflected rarefaction trailing edge angle";
gprt/gprstate.c:	    	ang_names[2] = "Reflected shock angle";
gprt/gprstate.c:	    ang_names[1] = "Reflected shock angle";
gprt/gprstate.c:	case B_REFLECT_NODE:
gprt/gprstate.c:	    ang_names[2] = "Reflected shock angle";
gprt/gprtprotos.h:IMPORT	void	show_front_states_for_ramp_reflections(Grid*,Wave*,Front*,
gprt/grectstat.c:	/* Assumes upper global boundaries are not reflecting or periodic. */
gprt/grectstat.c:		/* Note: assumes top boundary not reflecting or periodic */
gstate/gbstate.c:*	Imposes a reflecting boundary condition at points laying outside of a 
gstate/gbstate.c:*	neumann boundary.  The point (x,y) is reflected through the neumann 
gstate/gbstate.c:*	boundary curve Nbdry, the state at this reflected point (xref, yref) 
gstate/gbstate.c:*	normal component of velocity.  It is assumed that the reflected point 
gstate/gbstate.c:	/* coords_ref is the location of the reflection of the point coords
gstate/gbstate.c:	    || !reflect_pt_about_Nbdry(coords,coords_ref,nor,
gstate/gbstate.c:	       state at coords_tmp, but reflects the velocity of the state at
gstate/gbstate.c:	    /* Reflect the velocity of the state at coords_ref. */
gstate/gbstate.c:	       reflection point and boundary (coords_ref and coords_tmp,
gstate/gbstate.c:	    /* This version reflects the velocity of the state at coords_ref
gstate/gbstate.c:	    /* Reflect the velocity of the state at coords_ref. */
gstate/gbstate.c:	    /* Yet another variant:  Reflect or extrapolate the velocity. */
gstate/gbstate.c:	    if (debugging("nbdry4a"))  /* reflect the velocity */
gstate/gbstate.c:	       this function.  If there is gravity, then reflect_state() calls
gstate/gbstate.c:	       g_reflect_and_stratify_state(), which stratifies the density and
gstate/gbstate.c:	       pressure with respect to the state being reflected. */
gstate/gbstate.c:	    /* coords_tmp is the point on the reflection plane. */
gstate/gbstate.c:	    reflect_state(state,front->interf,coords,coords_tmp,nor);
gstate/gintrp.c:*	may be modified to reflect weights which are better for
gstate/gintrp.c:*	may be modified to reflect weights which are better for cylindrical
gstate/gintrp.c:*	may be modified to reflect weights which are better for
gstate/gintrp.c:*	reflect	weights which are better for cylindrical geometry.  (The 
gstate/goverinterp.c:                if(rect_boundary_type(rt_interf,axis,side) == REFLECTION_BOUNDARY)
gstate/grstate.c:	bool		reflect_velocity[3];
gstate/grstate.c:	    reflect_velocity[i] = NO;
gstate/grstate.c:	        if (rect_boundary_type(intfc,i,0) == REFLECTION_BOUNDARY)
gstate/grstate.c:	            reflect_velocity[i] = YES;
gstate/grstate.c:	        if (rect_boundary_type(intfc,i,1) == REFLECTION_BOUNDARY)
gstate/grstate.c:	            reflect_velocity[i] = YES;
gstate/grstate.c:	    if (reflect_velocity[i] == YES)
gstate/gstate.c:*			g_reflect_state():
gstate/gstate.c:*	Reflects the given state about the plane of symmetry with
gstate/gstate.c:EXPORT void g_reflect_state(
gstate/gstate.c:	Locstate	state,	/* state being reflected */
gstate/gstate.c:	INTERFACE	*intfc, /* interface being reflected */
gstate/gstate.c:	float		*pt,	/* position of state being reflected */
gstate/gstate.c:	float		*p,	/* point on reflection plane */
gstate/gstate.c:	float		*nor)	/* unit normal to reflection plane */
gstate/gstate.c:		(is_bad_state(state,YES,"g_reflect_state")))
gstate/gstate.c:	    screen("ERROR in g_reflect_state(), bad state detected\n");
gstate/gstate.c:}		/*end g_reflect_state*/
gstate/gstate.c:*			g_reflect_and_stratify_state():
gstate/gstate.c:*	Reflects the given state about the plane of symmetry with
gstate/gstate.c:EXPORT void g_reflect_and_stratify_state(
gstate/gstate.c:	Locstate	state,	/* state being reflected */
gstate/gstate.c:	INTERFACE	*intfc, /* interface being reflected */
gstate/gstate.c:	float		*pt,	/* position of state being reflected */
gstate/gstate.c:	float		*p,	/* arbitrary point on reflection plane */
gstate/gstate.c:	float		*nor)	/* unit normal to reflection plane */
gstate/gstate.c:	float	ds;		/* distance from pt to reflection plane */
gstate/gstate.c:	float	pr[MAXD];	/* Coordinates of reflected point */
gstate/gstate.c:	g_reflect_state(state,intfc,pt,p,nor);
gstate/gstate.c:		(is_bad_state(state,YES,"g_reflect_and_stratify_state")))
gstate/gstate.c:	    screen("ERROR in g_reflect_and_stratify_state(), "
gstate/gstate.c:}		/*end g_reflect_and_stratify_state*/
gstate/gstateprotos.h:IMPORT	void	g_reflect_and_stratify_state(Locstate,INTERFACE*,
gstate/gstateprotos.h:IMPORT	void	g_reflect_state(Locstate,INTERFACE*,float*,float*,float*);
gstate/gstglobs.c:	    fuh->_reflect_state = g_reflect_and_stratify_state;
gstate/gstglobs.c:	    f_user_interface(T->interface)._reflect_state =
gstate/gstglobs.c:	    	g_reflect_and_stratify_state;
in/ct_rgb:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/ct_rgb:Wave strength tolerance for tracking reflected shocks at attached boundary
in/ct_rgb:	reflection nodes = always track
in/ct_rgb:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/ct_rgb:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/ct_rgb:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/ct_rgb:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/ct_rgb:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/ct_rgb:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/ct_rgb:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/ct_rgb:Wave strength tolerance for tracking reflected shocks at shock-contact
in/ct_rgb:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/ct_rgb:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/ct_rgb:Don't Turn off tracking of regular reflection node if node propagation fails
in/ct_rgb:	Neumann boundary states are computed by an average of a reflection
in/ct_rgb:Reflect small loop shocks (dflt = no): 
in/ct_rgb:		a ramp reflection problem (RR),
in/ct_rgb:		Obstacle (behind reflecting wall) (O),
in/ct_rgb:                Obstacle (behind reflecting wall) (O),                 
in/ct_rgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/ct_rgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/ct_rgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/ct_rgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/drop:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/drop:Wave strength tolerance for tracking reflected shocks at attached boundary
in/drop:	reflection nodes = always track
in/drop:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/drop:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/drop:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/drop:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/drop:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/drop:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/drop:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/drop:Wave strength tolerance for tracking reflected shocks at shock-contact
in/drop:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/drop:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/drop:Don't Turn off tracking of regular reflection node if node propagation fails
in/drop:	Neumann boundary states are computed by an average of a reflection
in/drop:Reflect small loop shocks (dflt = no): 
in/drop:		a ramp reflection problem (RR),
in/drop:		Obstacle (behind reflecting wall) (O),
in/drop:                Obstacle (behind reflecting wall) (O),                 
in/drop:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/drop:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/drop:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/drop:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/frgb:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/frgb:Wave strength tolerance for tracking reflected shocks at attached boundary
in/frgb:	reflection nodes = always track
in/frgb:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/frgb:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/frgb:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/frgb:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/frgb:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/frgb:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/frgb:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/frgb:Wave strength tolerance for tracking reflected shocks at shock-contact
in/frgb:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/frgb:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/frgb:Don't Turn off tracking of regular reflection node if node propagation fails
in/frgb:	Neumann boundary states are computed by an average of a reflection
in/frgb:Reflect small loop shocks (dflt = no): 
in/frgb:		a ramp reflection problem (RR),
in/frgb:		Obstacle (behind reflecting wall) (O),
in/frgb:                Obstacle (behind reflecting wall) (O),                 
in/frgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/frgb:	for the left boundary in the x direction: Reflecting
in/frgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/frgb:	for the right boundary in the x direction: Reflecting
in/frgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/frgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/hcf_cross:	Neumann boundary states are computed by an average of a reflection
in/hcf_cross:		Obstacle (behind reflecting wall) (O),
in/hcf_cross:		Obstacle (behind reflecting wall) (O),
in/hcf_crossr:	Neumann boundary states are computed by an average of a reflection
in/hcf_jet:	Neumann boundary states are computed by an average of a reflection
in/hcf_jet:		Obstacle (behind reflecting wall) (O),
in/hcf_jet:		Obstacle (behind reflecting wall) (O),
in/hcf_jet:		Obstacle (behind reflecting wall) (O),
in/hcf_jet1:	Neumann boundary states are computed by an average of a reflection
in/hcf_jet1:		Obstacle (behind reflecting wall) (O),
in/hcf_jet1:		Obstacle (behind reflecting wall) (O),
in/hcf_jet1:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2:	Neumann boundary states are computed by an average of a reflection
in/hcf_jet2:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_a:	Neumann boundary states are computed by an average of a reflection
in/hcf_jet2_a:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_a:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_parab:	Neumann boundary states are computed by an average of a reflection
in/hcf_jet2_parab:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_parab:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_parab_160:	Neumann boundary states are computed by an average of a reflection
in/hcf_jet2_parab_160:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_parab_160:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_parab_160_1:	Neumann boundary states are computed by an average of a reflection
in/hcf_jet2_parab_160_1:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_parab_160_1:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_parab_320:	Neumann boundary states are computed by an average of a reflection
in/hcf_jet2_parab_320:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_parab_320:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_parab_bdry:	Neumann boundary states are computed by an average of a reflection
in/hcf_jet2_parab_bdry:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_parab_bdry:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_parab_bdry.r:	Neumann boundary states are computed by an average of a reflection
in/hcf_jet2_parab_new:	Neumann boundary states are computed by an average of a reflection
in/hcf_jet2_parab_new:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_parab_new:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_parab_nosgs_bdry:	Neumann boundary states are computed by an average of a reflection
in/hcf_jet2_parab_nosgs_bdry:		Obstacle (behind reflecting wall) (O),
in/hcf_jet2_parab_nosgs_bdry:		Obstacle (behind reflecting wall) (O),
in/hcf_jetr:	Neumann boundary states are computed by an average of a reflection
in/hcf_mirko:	Neumann boundary states are computed by an average of a reflection
in/hcf_mirko:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_320:	Neumann boundary states are computed by an average of a reflection
in/hcf_mirko_320:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_320:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_nozzle:	Neumann boundary states are computed by an average of a reflection
in/hcf_mirko_nozzle:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_nozzle:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_nozzle_r:	Neumann boundary states are computed by an average of a reflection
in/hcf_mirko_ut:	Neumann boundary states are computed by an average of a reflection
in/hcf_mirko_ut:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_ut:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_ut_full:	Neumann boundary states are computed by an average of a reflection
in/hcf_mirko_ut_full:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_ut_full:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_ut_full_320:	Neumann boundary states are computed by an average of a reflection
in/hcf_mirko_ut_full_320:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_ut_full_320:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_ut_full_sgs_160_new:	Neumann boundary states are computed by an average of a reflection
in/hcf_mirko_ut_full_sgs_160_new:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_ut_full_sgs_160_new:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_ut_full_sgs_160_new~:	Neumann boundary states are computed by an average of a reflection
in/hcf_mirko_ut_full_sgs_160_new~:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_ut_full_sgs_160_new~:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_ut_full_sgs_320_new:	Neumann boundary states are computed by an average of a reflection
in/hcf_mirko_ut_full_sgs_320_new:		Obstacle (behind reflecting wall) (O),
in/hcf_mirko_ut_full_sgs_320_new:		Obstacle (behind reflecting wall) (O),
in/hcf_mungal:	Neumann boundary states are computed by an average of a reflection
in/hcf_mungal:		Obstacle (behind reflecting wall) (O),
in/hcf_mungal:		Obstacle (behind reflecting wall) (O),
in/hcf_mungal_nosgs:	Neumann boundary states are computed by an average of a reflection
in/hcf_mungal_nosgs:		Obstacle (behind reflecting wall) (O),
in/hcf_mungal_nosgs:		Obstacle (behind reflecting wall) (O),
in/hcf_mungal_nosgs_cfl01:	Neumann boundary states are computed by an average of a reflection
in/hcf_mungal_nosgs_cfl01:		Obstacle (behind reflecting wall) (O),
in/hcf_mungal_nosgs_cfl01:		Obstacle (behind reflecting wall) (O),
in/hcf_mungal_nosgs_tvd:	Neumann boundary states are computed by an average of a reflection
in/hcf_mungal_nosgs_tvd:		Obstacle (behind reflecting wall) (O),
in/hcf_mungal_nosgs_tvd:		Obstacle (behind reflecting wall) (O),
in/hcf_mungal_oldmethod:	Neumann boundary states are computed by an average of a reflection
in/hcf_mungal_oldmethod:		Obstacle (behind reflecting wall) (O),
in/hcf_mungal_oldmethod:		Obstacle (behind reflecting wall) (O),
in/hcf_noparab:	Neumann boundary states are computed by an average of a reflection
in/hcf_noparab:		Obstacle (behind reflecting wall) (O),
in/hcf_noparab:		Obstacle (behind reflecting wall) (O),
in/hcf_test:	Neumann boundary states are computed by an average of a reflection
in/hcf_test:		Obstacle (behind reflecting wall) (O),
in/hcf_test:		Obstacle (behind reflecting wall) (O),
in/hcf_test1:	Neumann boundary states are computed by an average of a reflection
in/hcf_test1:		Obstacle (behind reflecting wall) (O),
in/hcf_test1:		Obstacle (behind reflecting wall) (O),
in/hcf_test1:Before construct_reflect_bdry 0, interface is consistent.
in/hcf_test1:After construct_reflect_bdry 0, interface is consistent.
in/hcf_testa:	Neumann boundary states are computed by an average of a reflection
in/hcf_testa:		Obstacle (behind reflecting wall) (O),
in/hcf_testa:		Obstacle (behind reflecting wall) (O),
in/hcf_test_bak:	Neumann boundary states are computed by an average of a reflection
in/hcf_test_bak:		Obstacle (behind reflecting wall) (O),
in/hcf_test_bak:		Obstacle (behind reflecting wall) (O),
in/hcf_testm:	Neumann boundary states are computed by an average of a reflection
in/hcf_testm:		Obstacle (behind reflecting wall) (O),
in/hcf_testm:		Obstacle (behind reflecting wall) (O),
in/hcf_testr:	Neumann boundary states are computed by an average of a reflection
in/imp2d:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/imp2d:Wave strength tolerance for tracking reflected shocks at attached boundary
in/imp2d:	reflection nodes = always track
in/imp2d:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/imp2d:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/imp2d:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/imp2d:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/imp2d:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/imp2d:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/imp2d:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/imp2d:Wave strength tolerance for tracking reflected shocks at shock-contact
in/imp2d:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/imp2d:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/imp2d:Don't Turn off tracking of regular reflection node if node propagation fails
in/imp2d:Wave strength tolerance for tracking reflected shocks at regular reflections,
in/imp2d:Wave strength tolerance for tracking reflected shocks at attached boundary
in/imp2d:	reflection nodes, yes = always track, no = never track, otherwise enter
in/imp2d:Wave strength tolerance for tracking the Mach stem at Mach reflections, yes =
in/imp2d:Wave strength tolerance for tracking reflected shocks at Mach reflections, yes
in/imp2d:Wave strength tolerance for tracking the slip line at Mach reflections, yes =
in/imp2d:Wave strength tolerance for tracking reflected shocks at shock crossings, yes =
in/imp2d:Wave strength tolerance for tracking reflected shocks at shock overtakes, yes =
in/imp2d:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/imp2d:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/imp2d:Wave strength tolerance for tracking reflected shocks at shock-contact
in/imp2d:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/imp2d:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/imp2d:	regular reflection node fails (default = 'n'): 
in/imp2d:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/imp2d:Wave strength tolerance for tracking reflected shocks at attached boundary
in/imp2d:	reflection nodes = always track
in/imp2d:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/imp2d:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/imp2d:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/imp2d:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/imp2d:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/imp2d:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/imp2d:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/imp2d:Wave strength tolerance for tracking reflected shocks at shock-contact
in/imp2d:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/imp2d:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/imp2d:Don't Turn off tracking of regular reflection node if node propagation fails
in/imp2d:	Neumann boundary states are computed by an average of a reflection
in/imp2d:	Neumann boundary states are computed by an average of a reflection
in/imp2d:Reflect small loop shocks (dflt = no): 
in/imp2d:		a ramp reflection problem (RR),
in/imp2d:		Obstacle (behind reflecting wall) (O),
in/imp2d:		Obstacle (behind reflecting wall) (O),
in/imp2d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/imp2d:	for the left boundary in the x direction: Reflecting
in/imp2d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/imp2d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/imp2d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/jet16:	Neumann boundary states are computed by an average of a reflection
in/jet16:          Obstacle (behind reflecting wall) (O),
in/jet16:		Obstacle (behind reflecting wall) (O),
in/jet16:		Obstacle (behind reflecting wall) (O),
in/jet16L:	Neumann boundary states are computed by an average of a reflection
in/jet16L:          Obstacle (behind reflecting wall) (O),
in/jet16L:		Obstacle (behind reflecting wall) (O),
in/jet16L:		Obstacle (behind reflecting wall) (O),
in/jet16LL:	Neumann boundary states are computed by an average of a reflection
in/jet16LL:		Obstacle (behind reflecting wall) (O),
in/jet16LL:		Obstacle (behind reflecting wall) (O),
in/jet16LL:		Obstacle (behind reflecting wall) (O),
in/jet16LLr:	Neumann boundary states are computed by an average of a reflection
in/jet16Lr:	Neumann boundary states are computed by an average of a reflection
in/jet16r:	Neumann boundary states are computed by an average of a reflection
in/jet2d1:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/jet2d1:Wave strength tolerance for tracking reflected shocks at attached boundary
in/jet2d1:	reflection nodes = always track
in/jet2d1:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/jet2d1:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/jet2d1:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/jet2d1:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/jet2d1:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/jet2d1:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/jet2d1:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/jet2d1:Wave strength tolerance for tracking reflected shocks at shock-contact
in/jet2d1:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/jet2d1:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/jet2d1:Don't Turn off tracking of regular reflection node if node propagation fails
in/jet2d1:	Neumann boundary states are computed by an average of a reflection
in/jet2d1:	Neumann boundary states are computed by an average of a reflection
in/jet2d1:Reflect small loop shocks (dflt = no): 
in/jet2d1:		a ramp reflection problem (RR),
in/jet2d1:		Obstacle (behind reflecting wall) (O),
in/jet2d1:		Obstacle (behind reflecting wall) (O),
in/jet2d1:Current choices are                 Obstacle (behind reflecting wall) (O),
in/jet2d1:		Obstacle (behind reflecting wall) (O),
in/jet2d1:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/jet2d1:	for the left boundary in the x direction: Reflecting
in/jet2d1:Enter the boundary type -- Unknown, Periodic, Reflecting,
in/jet2d1:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/jet2d2:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/jet2d2:Wave strength tolerance for tracking reflected shocks at attached boundary
in/jet2d2:	reflection nodes = always track
in/jet2d2:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/jet2d2:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/jet2d2:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/jet2d2:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/jet2d2:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/jet2d2:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/jet2d2:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/jet2d2:Wave strength tolerance for tracking reflected shocks at shock-contact
in/jet2d2:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/jet2d2:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/jet2d2:Don't Turn off tracking of regular reflection node if node propagation fails
in/jet2d2:	Neumann boundary states are computed by an average of a reflection
in/jet2d2:	Neumann boundary states are computed by an average of a reflection
in/jet2d2:Reflect small loop shocks (dflt = no): 
in/jet2d2:		a ramp reflection problem (RR),
in/jet2d2:		Obstacle (behind reflecting wall) (O),
in/jet2d2:		Obstacle (behind reflecting wall) (O),
in/jet2d2:Current choices are                 Obstacle (behind reflecting wall) (O),
in/jet2d2:		Obstacle (behind reflecting wall) (O),
in/jet2d2:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/jet2d2:	for the left boundary in the x direction: Reflecting
in/jet2d2:Enter the boundary type -- Unknown, Periodic, Reflecting,
in/jet2d2:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/jet2d3:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/jet2d3:Wave strength tolerance for tracking reflected shocks at attached boundary
in/jet2d3:	reflection nodes = always track
in/jet2d3:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/jet2d3:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/jet2d3:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/jet2d3:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/jet2d3:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/jet2d3:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/jet2d3:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/jet2d3:Wave strength tolerance for tracking reflected shocks at shock-contact
in/jet2d3:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/jet2d3:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/jet2d3:Don't Turn off tracking of regular reflection node if node propagation fails
in/jet2d3:	Neumann boundary states are computed by an average of a reflection
in/jet2d3:	Neumann boundary states are computed by an average of a reflection
in/jet2d3:Reflect small loop shocks (dflt = no): 
in/jet2d3:		a ramp reflection problem (RR),
in/jet2d3:		Obstacle (behind reflecting wall) (O),
in/jet2d3:		Obstacle (behind reflecting wall) (O),
in/jet2d3:Current choices are                 Obstacle (behind reflecting wall) (O),
in/jet2d3:		Obstacle (behind reflecting wall) (O),
in/jet2d3:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/jet2d3:	for the left boundary in the x direction: Reflecting
in/jet2d3:Enter the boundary type -- Unknown, Periodic, Reflecting,
in/jet2d3:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/jet2d4:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/jet2d4:Wave strength tolerance for tracking reflected shocks at attached boundary
in/jet2d4:	reflection nodes = always track
in/jet2d4:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/jet2d4:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/jet2d4:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/jet2d4:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/jet2d4:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/jet2d4:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/jet2d4:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/jet2d4:Wave strength tolerance for tracking reflected shocks at shock-contact
in/jet2d4:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/jet2d4:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/jet2d4:Don't Turn off tracking of regular reflection node if node propagation fails
in/jet2d4:	Neumann boundary states are computed by an average of a reflection
in/jet2d4:	Neumann boundary states are computed by an average of a reflection
in/jet2d4:Reflect small loop shocks (dflt = no): 
in/jet2d4:		a ramp reflection problem (RR),
in/jet2d4:		Obstacle (behind reflecting wall) (O),
in/jet2d4:		Obstacle (behind reflecting wall) (O),
in/jet2d4:Current choices are                 Obstacle (behind reflecting wall) (O),
in/jet2d4:		Obstacle (behind reflecting wall) (O),
in/jet2d4:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/jet2d4:	for the left boundary in the x direction: Reflecting
in/jet2d4:Enter the boundary type -- Unknown, Periodic, Reflecting,
in/jet2d4:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/jet2d5:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/jet2d5:Wave strength tolerance for tracking reflected shocks at attached boundary
in/jet2d5:	reflection nodes = always track
in/jet2d5:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/jet2d5:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/jet2d5:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/jet2d5:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/jet2d5:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/jet2d5:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/jet2d5:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/jet2d5:Wave strength tolerance for tracking reflected shocks at shock-contact
in/jet2d5:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/jet2d5:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/jet2d5:Don't Turn off tracking of regular reflection node if node propagation fails
in/jet2d5:	Neumann boundary states are computed by an average of a reflection
in/jet2d5:	Neumann boundary states are computed by an average of a reflection
in/jet2d5:Reflect small loop shocks (dflt = no): 
in/jet2d5:		a ramp reflection problem (RR),
in/jet2d5:		Obstacle (behind reflecting wall) (O),
in/jet2d5:		Obstacle (behind reflecting wall) (O),
in/jet2d5:Current choices are                 Obstacle (behind reflecting wall) (O),
in/jet2d5:		Obstacle (behind reflecting wall) (O),
in/jet2d5:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/jet2d5:	for the left boundary in the x direction: Reflecting
in/jet2d5:Enter the boundary type -- Unknown, Periodic, Reflecting,
in/jet2d5:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/jetdf2:	Neumann boundary states are computed by an average of a reflection
in/jetdf2:		Obstacle (behind reflecting wall) (O),
in/jetdf2:		Obstacle (behind reflecting wall) (O),
in/jetdf2:		Obstacle (behind reflecting wall) (O),
in/jetdf3:	Neumann boundary states are computed by an average of a reflection
in/jetdf3:		Obstacle (behind reflecting wall) (O),
in/jetdf3:		Obstacle (behind reflecting wall) (O),
in/jetdf3:		Obstacle (behind reflecting wall) (O),
in/jetdf31:	Neumann boundary states are computed by an average of a reflection
in/jetdf31:		Obstacle (behind reflecting wall) (O),
in/jetdf31:		Obstacle (behind reflecting wall) (O),
in/jetdf31:		Obstacle (behind reflecting wall) (O),
in/jetdf32:	Neumann boundary states are computed by an average of a reflection
in/jetdf32:		Obstacle (behind reflecting wall) (O),
in/jetdf32:		Obstacle (behind reflecting wall) (O),
in/jetdf32:		Obstacle (behind reflecting wall) (O),
in/jetkr:	Neumann boundary states are computed by an average of a reflection
in/jetkr:		Obstacle (behind reflecting wall) (O),
in/jetkr:		Obstacle (behind reflecting wall) (O),
in/jetkr:		Obstacle (behind reflecting wall) (O),
in/jetkr1:	Neumann boundary states are computed by an average of a reflection
in/jetkr1:		Obstacle (behind reflecting wall) (O),
in/jetkr1:		Obstacle (behind reflecting wall) (O),
in/jetkr1:		Obstacle (behind reflecting wall) (O),
in/jetkr2:	Neumann boundary states are computed by an average of a reflection
in/jetkr2:		Obstacle (behind reflecting wall) (O),
in/jetkr2:		Obstacle (behind reflecting wall) (O),
in/jetkr2:		Obstacle (behind reflecting wall) (O),
in/jetkrr:	Neumann boundary states are computed by an average of a reflection
in/jetng:	Neumann boundary states are computed by an average of a reflection
in/jetng:		Obstacle (behind reflecting wall) (O),
in/jetng:		Obstacle (behind reflecting wall) (O),
in/jetng:		Obstacle (behind reflecting wall) (O),
in/jetng10:	Neumann boundary states are computed by an average of a reflection
in/jetng10:		Obstacle (behind reflecting wall) (O),
in/jetng10:		Obstacle (behind reflecting wall) (O),
in/jetng10:		Obstacle (behind reflecting wall) (O),
in/jetng11:	Neumann boundary states are computed by an average of a reflection
in/jetng11:		Obstacle (behind reflecting wall) (O),
in/jetng11:		Obstacle (behind reflecting wall) (O),
in/jetng11:		Obstacle (behind reflecting wall) (O),
in/jetng12:	Neumann boundary states are computed by an average of a reflection
in/jetng12:		Obstacle (behind reflecting wall) (O),
in/jetng12:		Obstacle (behind reflecting wall) (O),
in/jetng12:		Obstacle (behind reflecting wall) (O),
in/jetng12r:	Neumann boundary states are computed by an average of a reflection
in/jetng12ra:	Neumann boundary states are computed by an average of a reflection
in/jetng13:	Neumann boundary states are computed by an average of a reflection
in/jetng13:		Obstacle (behind reflecting wall) (O),
in/jetng13:		Obstacle (behind reflecting wall) (O),
in/jetng13:		Obstacle (behind reflecting wall) (O),
in/jetng131:	Neumann boundary states are computed by an average of a reflection
in/jetng131:		Obstacle (behind reflecting wall) (O),
in/jetng131:		Obstacle (behind reflecting wall) (O),
in/jetng131:		Obstacle (behind reflecting wall) (O),
in/jetng131:Before construct_reflect_bdry 0, interface is consistent.
in/jetng131:no tris, exit construct_reflect_bdry.
in/jetng131:After construct_reflect_bdry 0, interface is consistent.
in/jetng13r:	Neumann boundary states are computed by an average of a reflection
in/jetng14:	Neumann boundary states are computed by an average of a reflection
in/jetng14:		Obstacle (behind reflecting wall) (O),
in/jetng14:		Obstacle (behind reflecting wall) (O),
in/jetng14:		Obstacle (behind reflecting wall) (O),
in/jetng15:	Neumann boundary states are computed by an average of a reflection
in/jetng15:		Obstacle (behind reflecting wall) (O),
in/jetng15:		Obstacle (behind reflecting wall) (O),
in/jetng15:		Obstacle (behind reflecting wall) (O),
in/jetng15r:	Neumann boundary states are computed by an average of a reflection
in/jetng16:	Neumann boundary states are computed by an average of a reflection
in/jetng16:		Obstacle (behind reflecting wall) (O),
in/jetng16:		Obstacle (behind reflecting wall) (O),
in/jetng16:		Obstacle (behind reflecting wall) (O),
in/jetng16r:	Neumann boundary states are computed by an average of a reflection
in/jetng16rt:	Neumann boundary states are computed by an average of a reflection
in/jetng17:	Neumann boundary states are computed by an average of a reflection
in/jetng17:		Obstacle (behind reflecting wall) (O),
in/jetng17:		Obstacle (behind reflecting wall) (O),
in/jetng17:		Obstacle (behind reflecting wall) (O),
in/jetng17r:	Neumann boundary states are computed by an average of a reflection
in/jetng1r:	Neumann boundary states are computed by an average of a reflection
in/jetng2:	Neumann boundary states are computed by an average of a reflection
in/jetng2:		Obstacle (behind reflecting wall) (O),
in/jetng2:		Obstacle (behind reflecting wall) (O),
in/jetng2:		Obstacle (behind reflecting wall) (O),
in/jetng2r:	Neumann boundary states are computed by an average of a reflection
in/jetng3:	Neumann boundary states are computed by an average of a reflection
in/jetng3:		Obstacle (behind reflecting wall) (O),
in/jetng3:		Obstacle (behind reflecting wall) (O),
in/jetng3:		Obstacle (behind reflecting wall) (O),
in/jetng4:	Neumann boundary states are computed by an average of a reflection
in/jetng4:		Obstacle (behind reflecting wall) (O),
in/jetng4:		Obstacle (behind reflecting wall) (O),
in/jetng4:		Obstacle (behind reflecting wall) (O),
in/jetng4r:	Neumann boundary states are computed by an average of a reflection
in/jetng5:	Neumann boundary states are computed by an average of a reflection
in/jetng5:		Obstacle (behind reflecting wall) (O),
in/jetng5:		Obstacle (behind reflecting wall) (O),
in/jetng5:		Obstacle (behind reflecting wall) (O),
in/jetng6:	Neumann boundary states are computed by an average of a reflection
in/jetng6:		Obstacle (behind reflecting wall) (O),
in/jetng6:		Obstacle (behind reflecting wall) (O),
in/jetng6:		Obstacle (behind reflecting wall) (O),
in/jetng7:	Neumann boundary states are computed by an average of a reflection
in/jetng7:		Obstacle (behind reflecting wall) (O),
in/jetng7:		Obstacle (behind reflecting wall) (O),
in/jetng7:		Obstacle (behind reflecting wall) (O),
in/jetng8:	Neumann boundary states are computed by an average of a reflection
in/jetng8:		Obstacle (behind reflecting wall) (O),
in/jetng8:		Obstacle (behind reflecting wall) (O),
in/jetng8:		Obstacle (behind reflecting wall) (O),
in/jetng9:	Neumann boundary states are computed by an average of a reflection
in/jetng9:		Obstacle (behind reflecting wall) (O),
in/jetng9:		Obstacle (behind reflecting wall) (O),
in/jetng9:		Obstacle (behind reflecting wall) (O),
in/jetngr:	Neumann boundary states are computed by an average of a reflection
in/jetnk:	Neumann boundary states are computed by an average of a reflection
in/jetnk:		Obstacle (behind reflecting wall) (O),
in/jetnk:		Obstacle (behind reflecting wall) (O),
in/jetnk:		Obstacle (behind reflecting wall) (O),
in/jetnk1:	Neumann boundary states are computed by an average of a reflection
in/jetnk1:		Obstacle (behind reflecting wall) (O),
in/jetnk1:		Obstacle (behind reflecting wall) (O),
in/jetnk1:		Obstacle (behind reflecting wall) (O),
in/jetnk1r:	Neumann boundary states are computed by an average of a reflection
in/jetnk2:	Neumann boundary states are computed by an average of a reflection
in/jetnk2:		Obstacle (behind reflecting wall) (O),
in/jetnk2:		Obstacle (behind reflecting wall) (O),
in/jetnk2:		Obstacle (behind reflecting wall) (O),
in/jetnk2r:	Neumann boundary states are computed by an average of a reflection
in/jetnk3:	Neumann boundary states are computed by an average of a reflection
in/jetnk3:		Obstacle (behind reflecting wall) (O),
in/jetnk3:		Obstacle (behind reflecting wall) (O),
in/jetnk3:		Obstacle (behind reflecting wall) (O),
in/jetnkr:	Neumann boundary states are computed by an average of a reflection
in/kh3d1:	Neumann boundary states are computed by an average of a reflection
in/kh3d1:		Obstacle (behind reflecting wall) (O),
in/kh3d1:	reflecting boundary state (F, default).
in/kh3d1:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/kh3d1:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/kh3d2:	Neumann boundary states are computed by an average of a reflection
in/kh3d2:		Obstacle (behind reflecting wall) (O),
in/kh3d2:		Obstacle (behind reflecting wall) (O),
in/kh3d2:	reflecting boundary state (F, default).
in/kh3d2:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/kh3d2:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/kh3d32:	Neumann boundary states are computed by an average of a reflection
in/kh3d32:		Obstacle (behind reflecting wall) (O),
in/kh3d32:	reflecting boundary state (F, default).
in/kh3d32:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/kh3d32:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/kh3d34:	Neumann boundary states are computed by an average of a reflection
in/kh3d34:		Obstacle (behind reflecting wall) (O),
in/kh3d34:	reflecting boundary state (F, default).
in/kh3d34:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/kh3d34:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/kh_test:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/kh_test:Wave strength tolerance for tracking reflected shocks at attached boundary
in/kh_test:	reflection nodes = always track
in/kh_test:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/kh_test:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/kh_test:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/kh_test:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/kh_test:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/kh_test:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/kh_test:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/kh_test:Wave strength tolerance for tracking reflected shocks at shock-contact
in/kh_test:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/kh_test:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/kh_test:Don't Turn off tracking of regular reflection node if node propagation fails
in/kh_test:	Neumann boundary states are computed by an average of a reflection
in/kh_test:Reflect small loop shocks (dflt = no): 
in/kh_test:		a ramp reflection problem (RR),
in/kh_test:		Obstacle (behind reflecting wall) (O),
in/kh_test:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/kh_test:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/kh_test:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/kh_test:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/muscl1d:	Neumann boundary states are computed by an average of a reflection
in/muscl1d:	Neumann boundary states are computed by an average of a reflection
in/muscl1d:		Obstacle (behind reflecting wall) (O),
in/muscl1d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/muscl1d:	for the left boundary in the x direction: Reflecting
in/muscl1d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/muscl1d:                           Mixed, Neumann, Reflecting, 
in/muscl1d:	for the right boundary in the x direction: Reflecting
in/new2:	Neumann boundary states are computed by an average of a reflection
in/new2:		Obstacle (behind reflecting wall) (O),
in/new2:		Obstacle (behind reflecting wall) (O),
in/new2~:	Neumann boundary states are computed by an average of a reflection
in/new2~:		Obstacle (behind reflecting wall) (O),
in/new2~:		Obstacle (behind reflecting wall) (O),
in/random2d1:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/random2d1:Wave strength tolerance for tracking reflected shocks at attached boundary
in/random2d1:	reflection nodes = always track
in/random2d1:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/random2d1:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/random2d1:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/random2d1:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/random2d1:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/random2d1:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d1:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d1:Wave strength tolerance for tracking reflected shocks at shock-contact
in/random2d1:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d1:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d1:Don't Turn off tracking of regular reflection node if node propagation fails
in/random2d1:	Neumann boundary states are computed by an average of a reflection
in/random2d1:Reflect small loop shocks (dflt = no): 
in/random2d1:		a ramp reflection problem (RR),
in/random2d1:		Obstacle (behind reflecting wall) (O),
in/random2d1:		Obstacle (behind reflecting wall) (O),
in/random2d1:	reflecting boundary state (F, default).
in/random2d1:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d1:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d10:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/random2d10:Wave strength tolerance for tracking reflected shocks at attached boundary
in/random2d10:	reflection nodes = always track
in/random2d10:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/random2d10:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/random2d10:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/random2d10:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/random2d10:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/random2d10:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d10:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d10:Wave strength tolerance for tracking reflected shocks at shock-contact
in/random2d10:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d10:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d10:Don't Turn off tracking of regular reflection node if node propagation fails
in/random2d10:	Neumann boundary states are computed by an average of a reflection
in/random2d10:Reflect small loop shocks (dflt = no): 
in/random2d10:		a ramp reflection problem (RR),
in/random2d10:		Obstacle (behind reflecting wall) (O),
in/random2d10:		Obstacle (behind reflecting wall) (O),
in/random2d10:	reflecting boundary state (F, default).
in/random2d10:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d10:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d11:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/random2d11:Wave strength tolerance for tracking reflected shocks at attached boundary
in/random2d11:	reflection nodes = always track
in/random2d11:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/random2d11:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/random2d11:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/random2d11:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/random2d11:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/random2d11:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d11:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d11:Wave strength tolerance for tracking reflected shocks at shock-contact
in/random2d11:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d11:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d11:Don't Turn off tracking of regular reflection node if node propagation fails
in/random2d11:	Neumann boundary states are computed by an average of a reflection
in/random2d11:Reflect small loop shocks (dflt = no): 
in/random2d11:		a ramp reflection problem (RR),
in/random2d11:		Obstacle (behind reflecting wall) (O),
in/random2d11:		Obstacle (behind reflecting wall) (O),
in/random2d11:	reflecting boundary state (F, default).
in/random2d11:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d11:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d12:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/random2d12:Wave strength tolerance for tracking reflected shocks at attached boundary
in/random2d12:	reflection nodes = always track
in/random2d12:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/random2d12:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/random2d12:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/random2d12:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/random2d12:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/random2d12:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d12:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d12:Wave strength tolerance for tracking reflected shocks at shock-contact
in/random2d12:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d12:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d12:Don't Turn off tracking of regular reflection node if node propagation fails
in/random2d12:	Neumann boundary states are computed by an average of a reflection
in/random2d12:Reflect small loop shocks (dflt = no): 
in/random2d12:		a ramp reflection problem (RR),
in/random2d12:		Obstacle (behind reflecting wall) (O),
in/random2d12:		Obstacle (behind reflecting wall) (O),
in/random2d12:	reflecting boundary state (F, default).
in/random2d12:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d12:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d1_tvd:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/random2d1_tvd:Wave strength tolerance for tracking reflected shocks at attached boundary
in/random2d1_tvd:	reflection nodes = always track
in/random2d1_tvd:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/random2d1_tvd:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/random2d1_tvd:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/random2d1_tvd:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/random2d1_tvd:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/random2d1_tvd:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d1_tvd:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d1_tvd:Wave strength tolerance for tracking reflected shocks at shock-contact
in/random2d1_tvd:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d1_tvd:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d1_tvd:Don't Turn off tracking of regular reflection node if node propagation fails
in/random2d1_tvd:	Neumann boundary states are computed by an average of a reflection
in/random2d1_tvd:Reflect small loop shocks (dflt = no): 
in/random2d1_tvd:		a ramp reflection problem (RR),
in/random2d1_tvd:		Obstacle (behind reflecting wall) (O),
in/random2d1_tvd:		Obstacle (behind reflecting wall) (O),
in/random2d1_tvd:	reflecting boundary state (F, default).
in/random2d1_tvd:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d1_tvd:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d2:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/random2d2:Wave strength tolerance for tracking reflected shocks at attached boundary
in/random2d2:	reflection nodes = always track
in/random2d2:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/random2d2:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/random2d2:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/random2d2:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/random2d2:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/random2d2:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d2:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d2:Wave strength tolerance for tracking reflected shocks at shock-contact
in/random2d2:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d2:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d2:Don't Turn off tracking of regular reflection node if node propagation fails
in/random2d2:	Neumann boundary states are computed by an average of a reflection
in/random2d2:Reflect small loop shocks (dflt = no): 
in/random2d2:		a ramp reflection problem (RR),
in/random2d2:		Obstacle (behind reflecting wall) (O),
in/random2d2:		Obstacle (behind reflecting wall) (O),
in/random2d2:	reflecting boundary state (F, default).
in/random2d2:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d2:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d3:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/random2d3:Wave strength tolerance for tracking reflected shocks at attached boundary
in/random2d3:	reflection nodes = always track
in/random2d3:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/random2d3:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/random2d3:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/random2d3:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/random2d3:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/random2d3:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d3:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d3:Wave strength tolerance for tracking reflected shocks at shock-contact
in/random2d3:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d3:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d3:Don't Turn off tracking of regular reflection node if node propagation fails
in/random2d3:	Neumann boundary states are computed by an average of a reflection
in/random2d3:Reflect small loop shocks (dflt = no): 
in/random2d3:		a ramp reflection problem (RR),
in/random2d3:		Obstacle (behind reflecting wall) (O),
in/random2d3:		Obstacle (behind reflecting wall) (O),
in/random2d3:	reflecting boundary state (F, default).
in/random2d3:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d3:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d4:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/random2d4:Wave strength tolerance for tracking reflected shocks at attached boundary
in/random2d4:	reflection nodes = always track
in/random2d4:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/random2d4:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/random2d4:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/random2d4:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/random2d4:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/random2d4:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d4:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d4:Wave strength tolerance for tracking reflected shocks at shock-contact
in/random2d4:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d4:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d4:Don't Turn off tracking of regular reflection node if node propagation fails
in/random2d4:	Neumann boundary states are computed by an average of a reflection
in/random2d4:Reflect small loop shocks (dflt = no): 
in/random2d4:		a ramp reflection problem (RR),
in/random2d4:		Obstacle (behind reflecting wall) (O),
in/random2d4:		Obstacle (behind reflecting wall) (O),
in/random2d4:	reflecting boundary state (F, default).
in/random2d4:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d4:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d5:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/random2d5:Wave strength tolerance for tracking reflected shocks at attached boundary
in/random2d5:	reflection nodes = always track
in/random2d5:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/random2d5:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/random2d5:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/random2d5:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/random2d5:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/random2d5:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d5:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d5:Wave strength tolerance for tracking reflected shocks at shock-contact
in/random2d5:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d5:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d5:Don't Turn off tracking of regular reflection node if node propagation fails
in/random2d5:	Neumann boundary states are computed by an average of a reflection
in/random2d5:Reflect small loop shocks (dflt = no): 
in/random2d5:		a ramp reflection problem (RR),
in/random2d5:		Obstacle (behind reflecting wall) (O),
in/random2d5:		Obstacle (behind reflecting wall) (O),
in/random2d5:	reflecting boundary state (F, default).
in/random2d5:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d5:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d6:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/random2d6:Wave strength tolerance for tracking reflected shocks at attached boundary
in/random2d6:	reflection nodes = always track
in/random2d6:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/random2d6:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/random2d6:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/random2d6:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/random2d6:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/random2d6:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d6:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d6:Wave strength tolerance for tracking reflected shocks at shock-contact
in/random2d6:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d6:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d6:Don't Turn off tracking of regular reflection node if node propagation fails
in/random2d6:	Neumann boundary states are computed by an average of a reflection
in/random2d6:Reflect small loop shocks (dflt = no): 
in/random2d6:		a ramp reflection problem (RR),
in/random2d6:		Obstacle (behind reflecting wall) (O),
in/random2d6:		Obstacle (behind reflecting wall) (O),
in/random2d6:	reflecting boundary state (F, default).
in/random2d6:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d6:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d7:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/random2d7:Wave strength tolerance for tracking reflected shocks at attached boundary
in/random2d7:	reflection nodes = always track
in/random2d7:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/random2d7:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/random2d7:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/random2d7:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/random2d7:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/random2d7:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d7:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d7:Wave strength tolerance for tracking reflected shocks at shock-contact
in/random2d7:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d7:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d7:Don't Turn off tracking of regular reflection node if node propagation fails
in/random2d7:	Neumann boundary states are computed by an average of a reflection
in/random2d7:Reflect small loop shocks (dflt = no): 
in/random2d7:		a ramp reflection problem (RR),
in/random2d7:		Obstacle (behind reflecting wall) (O),
in/random2d7:		Obstacle (behind reflecting wall) (O),
in/random2d7:	reflecting boundary state (F, default).
in/random2d7:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d7:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d8:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/random2d8:Wave strength tolerance for tracking reflected shocks at attached boundary
in/random2d8:	reflection nodes = always track
in/random2d8:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/random2d8:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/random2d8:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/random2d8:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/random2d8:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/random2d8:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d8:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d8:Wave strength tolerance for tracking reflected shocks at shock-contact
in/random2d8:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d8:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d8:Don't Turn off tracking of regular reflection node if node propagation fails
in/random2d8:	Neumann boundary states are computed by an average of a reflection
in/random2d8:Reflect small loop shocks (dflt = no): 
in/random2d8:		a ramp reflection problem (RR),
in/random2d8:		Obstacle (behind reflecting wall) (O),
in/random2d8:		Obstacle (behind reflecting wall) (O),
in/random2d8:	reflecting boundary state (F, default).
in/random2d8:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d8:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d9:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/random2d9:Wave strength tolerance for tracking reflected shocks at attached boundary
in/random2d9:	reflection nodes = always track
in/random2d9:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/random2d9:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/random2d9:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/random2d9:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/random2d9:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/random2d9:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d9:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d9:Wave strength tolerance for tracking reflected shocks at shock-contact
in/random2d9:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/random2d9:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/random2d9:Don't Turn off tracking of regular reflection node if node propagation fails
in/random2d9:	Neumann boundary states are computed by an average of a reflection
in/random2d9:Reflect small loop shocks (dflt = no): 
in/random2d9:		a ramp reflection problem (RR),
in/random2d9:		Obstacle (behind reflecting wall) (O),
in/random2d9:		Obstacle (behind reflecting wall) (O),
in/random2d9:	reflecting boundary state (F, default).
in/random2d9:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/random2d9:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/ri3d:	Neumann boundary states are computed by an average of a reflection
in/ri3d:		Obstacle (behind reflecting wall) (O),
in/ri3d:		Obstacle (behind reflecting wall) (O),
in/ri3d:		Obstacle (behind reflecting wall) (O),
in/ri3d3:	Neumann boundary states are computed by an average of a reflection
in/ri3d3:		Obstacle (behind reflecting wall) (O),
in/ri3d3:		Obstacle (behind reflecting wall) (O),
in/ri3d3:		Obstacle (behind reflecting wall) (O),
in/riem1d:	Neumann boundary states are computed by an average of a reflection
in/riem1d:	Neumann boundary states are computed by an average of a reflection
in/riem1d:		Obstacle (behind reflecting wall) (O),
in/riem1d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/riem1d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/rigfm:	Neumann boundary states are computed by an average of a reflection
in/rigfm:		Obstacle (behind reflecting wall) (O),
in/rigfm:		Obstacle (behind reflecting wall) (O),
in/rigfm:		Obstacle (behind reflecting wall) (O),
in/rigfm0:	Neumann boundary states are computed by an average of a reflection
in/rigfm0:		Obstacle (behind reflecting wall) (O),
in/rigfm0:		Obstacle (behind reflecting wall) (O),
in/rigfm0:		Obstacle (behind reflecting wall) (O),
in/rigfm0a:	Neumann boundary states are computed by an average of a reflection
in/rigfm0a:		Obstacle (behind reflecting wall) (O),
in/rigfm0a:		Obstacle (behind reflecting wall) (O),
in/rigfm0a:		Obstacle (behind reflecting wall) (O),
in/rigfma:	Neumann boundary states are computed by an average of a reflection
in/rigfma:		Obstacle (behind reflecting wall) (O),
in/rigfma:		Obstacle (behind reflecting wall) (O),
in/rigfma:		Obstacle (behind reflecting wall) (O),
in/rm2d:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/rm2d:Wave strength tolerance for tracking reflected shocks at attached boundary
in/rm2d:	reflection nodes = always track
in/rm2d:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/rm2d:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/rm2d:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/rm2d:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/rm2d:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/rm2d:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/rm2d:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/rm2d:Wave strength tolerance for tracking reflected shocks at shock-contact
in/rm2d:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/rm2d:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/rm2d:Don't Turn off tracking of regular reflection node if node propagation fails
in/rm2d:Wave strength tolerance for tracking reflected shocks at regular reflections,
in/rm2d:Wave strength tolerance for tracking reflected shocks at attached boundary
in/rm2d:	reflection nodes, yes = always track, no = never track, otherwise enter
in/rm2d:Wave strength tolerance for tracking the Mach stem at Mach reflections, yes =
in/rm2d:Wave strength tolerance for tracking reflected shocks at Mach reflections, yes
in/rm2d:Wave strength tolerance for tracking the slip line at Mach reflections, yes =
in/rm2d:Wave strength tolerance for tracking reflected shocks at shock crossings, yes =
in/rm2d:Wave strength tolerance for tracking reflected shocks at shock overtakes, yes =
in/rm2d:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/rm2d:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/rm2d:Wave strength tolerance for tracking reflected shocks at shock-contact
in/rm2d:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/rm2d:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/rm2d:	regular reflection node fails (default = 'n'): 
in/rm2d:Reflect small loop shocks (dflt = no): 
in/rm2d:		a ramp reflection problem (RR),
in/rm2d:		Obstacle (behind reflecting wall) (O),
in/rm2d:		Obstacle (behind reflecting wall) (O),
in/rm2d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/rm2d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/rm2d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/rm2d1:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/rm2d1:Wave strength tolerance for tracking reflected shocks at attached boundary
in/rm2d1:	reflection nodes = always track
in/rm2d1:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/rm2d1:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/rm2d1:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/rm2d1:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/rm2d1:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/rm2d1:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/rm2d1:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/rm2d1:Wave strength tolerance for tracking reflected shocks at shock-contact
in/rm2d1:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/rm2d1:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/rm2d1:Don't Turn off tracking of regular reflection node if node propagation fails
in/rm2d1:Wave strength tolerance for tracking reflected shocks at regular reflections,
in/rm2d1:Wave strength tolerance for tracking reflected shocks at attached boundary
in/rm2d1:	reflection nodes, yes = always track, no = never track, otherwise enter
in/rm2d1:Wave strength tolerance for tracking the Mach stem at Mach reflections, yes =
in/rm2d1:Wave strength tolerance for tracking reflected shocks at Mach reflections, yes
in/rm2d1:Wave strength tolerance for tracking the slip line at Mach reflections, yes =
in/rm2d1:Wave strength tolerance for tracking reflected shocks at shock crossings, yes =
in/rm2d1:Wave strength tolerance for tracking reflected shocks at shock overtakes, yes =
in/rm2d1:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/rm2d1:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/rm2d1:Wave strength tolerance for tracking reflected shocks at shock-contact
in/rm2d1:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/rm2d1:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/rm2d1:	regular reflection node fails (default = 'n'): 
in/rm2d1:Reflect small loop shocks (dflt = no): 
in/rm2d1:		a ramp reflection problem (RR),
in/rm2d1:		Obstacle (behind reflecting wall) (O),
in/rm2d1:		Obstacle (behind reflecting wall) (O),
in/rm2d1:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/rm2d1:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/rm2d1:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/rt2d1:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/rt2d1:Wave strength tolerance for tracking reflected shocks at attached boundary
in/rt2d1:	reflection nodes = always track
in/rt2d1:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/rt2d1:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/rt2d1:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/rt2d1:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/rt2d1:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/rt2d1:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/rt2d1:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/rt2d1:Wave strength tolerance for tracking reflected shocks at shock-contact
in/rt2d1:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/rt2d1:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/rt2d1:Don't Turn off tracking of regular reflection node if node propagation fails
in/rt2d1:	Neumann boundary states are computed by an average of a reflection
in/rt2d1:Reflect small loop shocks (dflt = no): 
in/rt2d1:		a ramp reflection problem (RR),
in/rt2d1:		Obstacle (behind reflecting wall) (O),
in/rt2d1:		Obstacle (behind reflecting wall) (O),
in/rt2d1:	reflecting boundary state (F, default).
in/rt2d1:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/rt2d1:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/rt3d:	Neumann boundary states are computed by an average of a reflection
in/rt3d:		Obstacle (behind reflecting wall) (O),
in/rt3d:		Obstacle (behind reflecting wall) (O),
in/rt3d:	reflecting boundary state (F, default).
in/rt3d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/rt3d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/rt3d1:		a ramp reflection problem (R),
in/rt3d1:	Obstacle (behind reflecting wall) (O).
in/rt3d1:	Obstacle (behind reflecting wall) (O).
in/rt3d1:	reflecting boundary state (F, default).
in/rt3d1q:		a ramp reflection problem (R),
in/rt3d1q:	Obstacle (behind reflecting wall) (O).
in/rt3d1q:	Obstacle (behind reflecting wall) (O).
in/rt3d1q:	reflecting boundary state (F, default).
in/rt3d2q:		a ramp reflection problem (R),
in/rt3d2q:	Obstacle (behind reflecting wall) (O).
in/rt3d2q:	Obstacle (behind reflecting wall) (O).
in/rt3d2q:	reflecting boundary state (F, default).
in/rt3dm:	Neumann boundary states are computed by an average of a reflection
in/rt3dm:		Obstacle (behind reflecting wall) (O),
in/rt3dm:		Obstacle (behind reflecting wall) (O),
in/rt3dm:	reflecting boundary state (F, default).
in/rt3dm:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/rt3dm:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/rt3d_mp:	Neumann boundary states are computed by an average of a reflection
in/rt3d_mp:		Obstacle (behind reflecting wall) (O),
in/rt3d_mp:	reflecting boundary state (F, default).
in/rt3d_mp:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/rt3d_mp:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/rt3dr:	Neumann boundary states are computed by an average of a reflection
in/shear:	Neumann boundary states are computed by an average of a reflection
in/shear:		Obstacle (behind reflecting wall) (O),
in/shear:		Obstacle (behind reflecting wall) (O),
in/shear:		Obstacle (behind reflecting wall) (O),
in/shear1:	Neumann boundary states are computed by an average of a reflection
in/shear1:		Obstacle (behind reflecting wall) (O),
in/shear1:		Obstacle (behind reflecting wall) (O),
in/shear1:		Obstacle (behind reflecting wall) (O),
in/sh_mrgb:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/sh_mrgb:Wave strength tolerance for tracking reflected shocks at attached boundary
in/sh_mrgb:	reflection nodes = always track
in/sh_mrgb:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/sh_mrgb:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/sh_mrgb:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/sh_mrgb:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/sh_mrgb:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/sh_mrgb:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/sh_mrgb:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/sh_mrgb:Wave strength tolerance for tracking reflected shocks at shock-contact
in/sh_mrgb:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/sh_mrgb:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/sh_mrgb:Don't Turn off tracking of regular reflection node if node propagation fails
in/sh_mrgb:	Neumann boundary states are computed by an average of a reflection
in/sh_mrgb:Reflect small loop shocks (dflt = no): 
in/sh_mrgb:		a ramp reflection problem (RR),
in/sh_mrgb:		Obstacle (behind reflecting wall) (O),
in/sh_mrgb:                Obstacle (behind reflecting wall) (O),                 
in/sh_mrgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/sh_mrgb:	for the left boundary in the x direction: Reflecting
in/sh_mrgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/sh_mrgb:	for the right boundary in the x direction: Reflecting
in/sh_mrgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/sh_mrgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/sh_rgb:Wave strength tolerance for tracking reflected shocks at regular reflections =
in/sh_rgb:Wave strength tolerance for tracking reflected shocks at attached boundary
in/sh_rgb:	reflection nodes = always track
in/sh_rgb:Wave strength tolerance for tracking the slip line at Mach reflections = always
in/sh_rgb:Wave strength tolerance for tracking reflected shocks at Mach reflections =
in/sh_rgb:Wave strength tolerance for tracking the Mach stem at Mach reflections = always
in/sh_rgb:Wave strength tolerance for tracking reflected shocks at shock crossings =
in/sh_rgb:Wave strength tolerance for tracking reflected shocks at shock overtakes =
in/sh_rgb:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/sh_rgb:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/sh_rgb:Wave strength tolerance for tracking reflected shocks at shock-contact
in/sh_rgb:Wave strength tolerance for tracking reflected rarefaction leading edges at
in/sh_rgb:Wave strength tolerance for tracking reflected rarefaction trailing edges at
in/sh_rgb:Don't Turn off tracking of regular reflection node if node propagation fails
in/sh_rgb:	Neumann boundary states are computed by an average of a reflection
in/sh_rgb:Reflect small loop shocks (dflt = no): 
in/sh_rgb:		a ramp reflection problem (RR),
in/sh_rgb:		Obstacle (behind reflecting wall) (O),
in/sh_rgb:                Obstacle (behind reflecting wall) (O),                 
in/sh_rgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/sh_rgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/sh_rgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/sh_rgb:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/test1d:	Neumann boundary states are computed by an average of a reflection
in/test1d:	Neumann boundary states are computed by an average of a reflection
in/test1d:		Obstacle (behind reflecting wall) (O),
in/test1d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/test1d:	for the left boundary in the x direction: Reflecting
in/test1d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/test1d:                           Mixed, Neumann, Reflecting, 
in/test1d:	for the right boundary in the x direction: Reflecting
in/test.r:	Neumann boundary states are computed by an average of a reflection
in/tmp1:	Neumann boundary states are computed by an average of a reflection
in/tmp1:		Obstacle (behind reflecting wall) (O),
in/tmp1:		Obstacle (behind reflecting wall) (O),
in/tmp1~:	Neumann boundary states are computed by an average of a reflection
in/tmp1~:		Obstacle (behind reflecting wall) (O),
in/tmp1~:		Obstacle (behind reflecting wall) (O),
in/tmp1_1:	Neumann boundary states are computed by an average of a reflection
in/tmp1_1:		Obstacle (behind reflecting wall) (O),
in/tmp1_1:		Obstacle (behind reflecting wall) (O),
in/tmp2:	Neumann boundary states are computed by an average of a reflection
in/tmp2:		Obstacle (behind reflecting wall) (O),
in/tmp2:		Obstacle (behind reflecting wall) (O),
in/tmp_2:	Neumann boundary states are computed by an average of a reflection
in/tmp_2:		Obstacle (behind reflecting wall) (O),
in/tmp_2:		Obstacle (behind reflecting wall) (O),
in/tmp2~:	Neumann boundary states are computed by an average of a reflection
in/tmp2~:		Obstacle (behind reflecting wall) (O),
in/tmp2~:		Obstacle (behind reflecting wall) (O),
in/tmp_for_gdb:	Neumann boundary states are computed by an average of a reflection
in/tmp_for_gdb:		Obstacle (behind reflecting wall) (O),
in/tmp_for_gdb:		Obstacle (behind reflecting wall) (O),
in/tmp_for_gdb~:	Neumann boundary states are computed by an average of a reflection
in/tmp_for_gdb~:		Obstacle (behind reflecting wall) (O),
in/tmp_for_gdb~:		Obstacle (behind reflecting wall) (O),
in/tmp_for_gdb2:	Neumann boundary states are computed by an average of a reflection
in/tmp_for_gdb2:		Obstacle (behind reflecting wall) (O),
in/tmp_for_gdb2:		Obstacle (behind reflecting wall) (O),
in/tu_in:	Neumann boundary states are computed by an average of a reflection
in/tu_in:		Obstacle (behind reflecting wall) (O),
in/tu_in:		Obstacle (behind reflecting wall) (O),
in/tvd1d:	Neumann boundary states are computed by an average of a reflection
in/tvd1d:	Neumann boundary states are computed by an average of a reflection
in/tvd1d:		Obstacle (behind reflecting wall) (O),
in/tvd1d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/tvd1d:	for the left boundary in the x direction: Reflecting
in/tvd1d:Enter the boundary type -- Unknown, Periodic, Reflecting, 
in/tvd1d:                           Mixed, Neumann, Reflecting, 
in/tvd1d:	for the right boundary in the x direction: Reflecting
in/vel_test:	Neumann boundary states are computed by an average of a reflection
in/vel_test:		Obstacle (behind reflecting wall) (O),
in/vel_test1:	Neumann boundary states are computed by an average of a reflection
in/vel_test1:		Obstacle (behind reflecting wall) (O),
in/vel_test1r:	Neumann boundary states are computed by an average of a reflection
in/vel_test2:	Neumann boundary states are computed by an average of a reflection
in/vel_test2:		Obstacle (behind reflecting wall) (O),
in/vinay_version:	Neumann boundary states are computed by an average of a reflection
in/vinay_version:		Obstacle (behind reflecting wall) (O),
in/vinay_version:		Obstacle (behind reflecting wall) (O),
in/xu_in:	Neumann boundary states are computed by an average of a reflection
in/xu_in:		Obstacle (behind reflecting wall) (O),
in/xu_in:		Obstacle (behind reflecting wall) (O),
in/xu_in:Before construct_reflect_bdry 0, interface is consistent.
in/xu_in:After construct_reflect_bdry 0, interface is consistent.
in/xu_inc:	Neumann boundary states are computed by an average of a reflection
in/xu_inc:		Obstacle (behind reflecting wall) (O),
in/xu_inc:		Obstacle (behind reflecting wall) (O),
in/xu_inc:Before construct_reflect_bdry 0, interface is consistent.
in/xu_inc:After construct_reflect_bdry 0, interface is consistent.
in/xu_inf:	Neumann boundary states are computed by an average of a reflection
in/xu_inf:		Obstacle (behind reflecting wall) (O),
in/xu_inf:		Obstacle (behind reflecting wall) (O),
in/xu_inf:Before construct_reflect_bdry 0, interface is consistent.
in/xu_inf:After construct_reflect_bdry 0, interface is consistent.
