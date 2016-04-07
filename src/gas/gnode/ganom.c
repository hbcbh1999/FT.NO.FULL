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


#if defined(FULL_PHYSICS) && defined(TWOD)
/*
*
*				ganom.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Transforms a regular shock diffraction into anomalous reflection
*
*/


#include <gdecs/gdecs.h>

typedef struct {
	POINT   *pc;
	BOND    *bc; /* Output */
	O_CURVE *newc6;
	WSSten	 *wssten;
	Locstate rpst0, rpst1;
	double    abs_v[MAXD]; /* Output */

	/* Anomalous node states by mass flux specific data */

	POINT    *p00, *pvertex;
	BOND     *b_ct[2];
	Locstate rpst6;
	double    abs_v00[MAXD];
	double    **t;
	SIDE     inc_side;

	/* Anomalous node states by pressure specific data */

	POINT    *oldp, *newp;
	O_CURVE  *newc0;
	BOND     *binc;
	Locstate osb0, osb;
	Locstate sot;
	double    *oabs_v;
	int      is_plus_or;
} ANOM;

/* Values returned by set_limits_for_anom_find_root() */

enum _SL_AFR {
	UNABLE_TO_SET_LIMITS	= 0,
	LIMITS_SET		= 1,
	SUPERSONIC_BY_MASS_FLUX	= 2
};
typedef enum _SL_AFR SL_AFR;

	/* LOCAL Function Declarations */
LOCAL	ANOM	*allocate_static_ANOM_data_structure(Front*);
LOCAL	SL_AFR	set_limits_for_anom_find_root(NODE*,Front*,double,RP_DATA*,ANOM*,
					      double*,double*,double*,double*);
LOCAL	boolean	f_mf_anom(double,double*,POINTER);
LOCAL	boolean	f_pr_anom(double,double*,POINTER);
LOCAL	int	anom_states_by_mass_flux(NODE*,NODE*,O_CURVE**,O_CURVE**,
					 BOND**,POINT*,Front*,Wave*,
					 RPROBLEM**,double,double*,NODE_FLAG,
					 RP_DATA*,double**);
LOCAL	int	anom_states_by_pressure(NODE*,NODE*,O_CURVE**,O_CURVE**,
					BOND**,POINT*,Front*,Wave*,
					RPROBLEM**,double,double*,NODE_FLAG,
					RP_DATA*,double**);
LOCAL	int	prop_char_up_incident_shock(NODE*,NODE*,O_CURVE**,O_CURVE**,
					    BOND**,POINT*,Front*,Wave*,
					    RPROBLEM**,double,double*,
					    ANGLE_DIRECTION,NODE_FLAG);
LOCAL	void	correct_for_refl_wave_under_shoot(Front*,O_CURVE*,BOND*,POINT*,
						  BOND**,POINT*);
LOCAL	void	init_ANOM_data_structure(ANOM*,NODE*,NODE*,NODE**,O_CURVE**,
					 O_CURVE**,BOND**,BOND**,POINT*,POINT*,
					 POINT*,POINT**,Front*,Wave*,double,
					 NODE_FLAG,RP_DATA*,double*,double**);

#if defined(DEBUG_NODE_PROPAGATE)
LOCAL	void	verbose_print_f_mf_anom_data(double,ANOM*,BOND*,
					     double*,double,double);
#endif /* defined(DEBUG_NODE_PROPAGATE) */



EXPORT	int anomalous_reflection_propagate(
	Front		*fr,
	Wave		*wave,
	NODE		*oldn,
	NODE		*newn,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	boolean		*is_refl_shock,
	RP_DATA		*RP,
	POINT		*pc,
	double		*tcr,
	BOND		**newb)
{
	double		va[MAXD], sinsqr, mf_sqr, qa, qasqr, len, aisqr;
	POINT		*pvertex;
	int		status = ERROR_NODE;
	int		dim = fr->rect_grid->dim;
	int		i;
	static	int	(*anomalous_node_states[2])(NODE*,NODE*,O_CURVE**,
						    O_CURVE**,BOND**,POINT*,
						    Front*,Wave*,RPROBLEM**,
						    double,double*,NODE_FLAG,
						    RP_DATA*,double**);
	static	int	anom_index = 0;
	static	double	**t = NULL;

	debug_print("diffraction","Entered anomalous_reflection_propagate()\n");

	if (oldn == NULL)
	{
	    (void) printf("WARNING in anomalous_reflection_propagate(), "
	                  "Oldn == NULL, double bifurcation in rproblem\n");
	    return status;
	}

	if (t == NULL)
	{
	    anomalous_node_states[0] = anom_states_by_mass_flux;
	    anomalous_node_states[1] = anom_states_by_pressure;
	    bi_array(&t,2,2,FLOAT);
	}

		/*     Propagate leading edge of reflected	*/
		/*     rarefaction wave up the incident shock.	*/

	if (newc[1]->curve != NULL) /* Bifurcation from regular to anomalous */
	{
	    status = refl_curve_overtakes_incident_shock(oldc,newc,newb,pc,fr,
							 wave,rp,dt,dt_frac,1,
							 RP->ang_dir,flag);


	    if (status != GOOD_NODE) 
	    {
#if defined(DEBUG_NODE_PROPAGATE)
	    	if (debugging("diffraction"))
	    	{
	    	    (void) printf("Unable to propagate anomalous reflection\n");
		}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
		goto leave;
	    }
	}
	else
	{
	    status = prop_char_up_incident_shock(oldn,newn,oldc,newc,newb,pc,
						 fr,wave,rp,dt,dt_frac,
						 RP->ang_dir,flag);

	    if (status != GOOD_NODE) 
	    {
#if defined(DEBUG_NODE_PROPAGATE)
	    	if (debugging("diffraction"))
	    	{
	    	    (void) printf("Unable to propagate anomalous reflection\n");
	    	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    	goto leave;
	    }
	}

	/* Check for sonic transition */

	pvertex  = Point_adjacent_to_node(newc[0]->curve,newc[0]->orient);
	for (i = 0; i < dim; i++)
	    t[0][i] = Coords(pvertex)[i] - Coords(pc)[i];
	len = mag_vector(t[0],dim);
	t[0][0] /= len;			t[0][1] /= len;

	for (i = 0; i < dim; i++)
	    va[i] = vel(i,RP->state[0]) - Node_vel(newn)[i];
	qasqr = scalar_product(va,va,dim);	    qa = sqrt(qasqr);
	for (i = 0; i < dim; i++)
	    t[1][i] = -va[i]/qa;

	(void) vector_product(t[1],t[0],&sinsqr,dim); sinsqr = sinsqr*sinsqr;
	mf_sqr = sqr(Dens(RP->state[0])) * qasqr * sinsqr;
	aisqr = acoustic_impedance_squared(RP->state[0]);
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("diffraction"))
	{
	    double ang;
	    (void) printf("Sonic transition test\n pvertex = <%g, %g>\n",
			  Coords(pvertex)[0],Coords(pvertex)[1]);
	    (void) printf("pc = <%g, %g>\n",Coords(pc)[0],Coords(pc)[1]);
	    (void) printf("Scaled separation pvertex to pc = %g\n",
			  scaled_separation(pvertex,pc,fr->rect_grid->h,
			  fr->rect_grid->dim));
	    (void) printf("New incident shock tangent = <%g, %g>,\n",
			  t[0][0],t[0][1]);
	    ang = atan2(t[0][1],t[0][0]);
	    (void) printf("\t\t\tmagnitude = %g, ",mag_vector(t[0],dim));
	    print_angle("angle =",ang,"\n");
	    (void) printf("New incident contact tangent = <%g, %g>,\n",
			  t[1][0],t[1][1]);
	    ang = atan2(t[1][1],t[1][0]);
	    (void) printf("\t\t\tmagnitude = %g, ",mag_vector(t[1],dim));
	    print_angle("angle =",ang,"\n");
	    (void) printf("mf_sqr = %g,  acoustic impedance squared = %g\n",
			  mf_sqr,aisqr);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	if (mf_sqr < aisqr)
	{
	    status = sonic_incident_shock_at_diffraction(newn,oldc,newc,newb,
							 pc,tcr,RP,fr,wave,
							 rp,dt,dt_frac,flag);
	    goto leave;
	}

	for (i = 0; i < 2; i++)
	{
	    status = (*anomalous_node_states[anom_index])(oldn,newn,oldc,newc,
							  newb,pc,fr,
							  (Wave*)wave,rp,
							  dt,dt_frac,flag,RP,t);

	    if (status == GOOD_NODE)
		break;
	    anom_index = (anom_index + 1)%2;
	}

	if (status != GOOD_NODE)
	{
	    (void) printf("WARNING in anomalous_reflection_propagate(), "
	                  "Unable to find state behind incident shock\n");
	    goto leave;
	}

	if (is_regular_diffraction_node(Coords(pc),Node_vel(newn),
					Node_vel(oldn),t,RP,Rp_data(oldn),
					is_refl_shock,fr,DIFFRACTION_NODE,flag)
	    != REGULAR_DIFFRACTION)
	{
	    (void) printf("WARNING in anomalous_reflection_propagate(), "
	                  "is_regular_diffraction_node() fails\n");
	    status = ERROR_NODE;
	    goto leave;
	}

	status = modify_diffraction_node(pc,oldn,newn,oldc,newc,newb,fr,wave,
					 dt,RP,flag);

leave:
	debug_print("diffraction","Left anomalous_reflection_propagate(), ");
	if (debugging("diffraction"))
	    print_node_status("status = ",status,"\n");
	return status;
}		/*end anomalous_reflection_propagate*/


EXPORT	int refl_curve_overtakes_incident_shock(
	O_CURVE		**oldc,
	O_CURVE		**newc,
	BOND		**newb,
	POINT		*pc,
	Front		*fr,
	Wave		*wave,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	int		index,
	ANGLE_DIRECTION	ang_dir,
	NODE_FLAG	flag)
{
	POINT		*pcrs;
	BOND		*newb0, *newbi;
	NODE		*newn, *diffnd;
	CURVE		*newc0 = newc[0]->curve;
	CURVE		**curs;
	INTERFACE	*intfc = newc[0]->curve->interface;
	COMPONENT	left0, right0, left1, right1;
	double		tcrs0, tcrsi;
	double		*nv, *dnv, sav_v[MAXD];
	ORIENTATION	newc0_orient = newc[0]->orient;
	boolean		sav_intrp = interpolate_intfc_states(intfc);
	boolean		sav_copy = copy_intfc_states();
	int		status = ERROR_NODE;
	boolean		c_ext[2];
	int		i, dim = fr->rect_grid->dim;
	boolean		sav_scss = interpolate_states_at_split_curve_node();
	Locstate	lst, rst;

	debug_print("diffraction","Entered refl_curve_overtakes_incident_shock()\n");

	if (newc[index] == NULL || newc[index]->curve == NULL) 
	{
		status = GOOD_NODE;
		debug_print("diffraction",
			"Left refl_curve_overtakes_incident_shock(), ");
		if (debugging("diffraction"))
		    print_node_status("status = ",status,"\n");
		return status;
	}

	pcrs = Point(NULL);
	lst = left_state(pcrs);	rst = right_state(pcrs);

	diffnd = Node_of_o_curve(newc[0]);
	for (i = 0; i < dim; i++) sav_v[i] = Node_vel(diffnd)[i];
	status = cross_or_extend_to_cross_two_propagated_curves(oldc[0],
		    newc[0],oldc[index],newc[index],&pcrs,&newb0,&newbi,
		    &tcrs0,&tcrsi,fr,(POINTER)wave,rp,dt,dt_frac,flag,c_ext);
	
	if (status != GOOD_NODE) 
	{
#if defined(DEBUG_NODE_PROPAGATE)
		if (debugging("diffraction"))
		{
			(void) printf("Unable to propagate ");
			(void) printf("anomalous reflection\n");
			(void) printf("oldc[0]\n");	print_o_curve(oldc[0]);
			(void) printf("oldc[%d]\n",index);
			print_o_curve(oldc[index]);
			(void) printf("newc[0]\n");	print_o_curve(newc[0]);
			(void) printf("newc[%d]\n",index);
			print_o_curve(newc[index]);
		}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
		status = MODIFY_TIME_STEP_NODE;
		interpolate_intfc_states(intfc) = sav_intrp;
		for (i = 0; i < dim; i++) Node_vel(diffnd)[i] = sav_v[i];
		debug_print("diffraction",
			"Left refl_curve_overtakes_incident_shock(), ");
		if (debugging("diffraction"))
		    print_node_status("status = ",status,"\n");
		*dt_frac = min(*dt_frac,Min_time_step_modification_factor(fr));
		return status;
	}

	correct_for_refl_wave_under_shoot(fr,newc[0],newb[0],pc,&newb0,pcrs);

		/*
		*  Disconnect newc[index]->curve from its node and install a
		*  new overtake node at the crossing position between
		*  newc[0]->curve and newc[index]->curve
		*/

	interpolate_intfc_states(intfc) = YES;
	if (insert_point_in_bond(pcrs,newb0,newc0) != FUNCTION_SUCCEEDED)
	{
	    screen("ERROR in refl_curve_overtakes_incident_shock(), "
		   "insert_point_in_bond() failed\n");
	    clean_up(ERROR);
	}
	if (newc0_orient == POSITIVE_ORIENTATION)
	{
		left1  = negative_component(newc0);
		right1 = positive_component(newc0);
		if (ang_dir == COUNTER_CLOCK)
		{
			right0 = positive_component(newc0);
			left0  = (newc[index]->orient==POSITIVE_ORIENTATION) ?
				  negative_component(newc[index]->curve) :
				  positive_component(newc[index]->curve);
		}
		else
		{
			left0  = negative_component(newc0);
			right0 = (newc[index]->orient==POSITIVE_ORIENTATION) ?
				  positive_component(newc[index]->curve) :
				  negative_component(newc[index]->curve);
		}
	}
	else
	{
		left0  = negative_component(newc0);
		right0 = positive_component(newc0);
		if (ang_dir == COUNTER_CLOCK)
		{
			left1  = negative_component(newc[0]->curve);
			right1 = (newc[index]->orient==POSITIVE_ORIENTATION) ?
				  negative_component(newc[index]->curve) :
				  positive_component(newc[index]->curve);
		}
		else
		{
			right1 = positive_component(newc[0]->curve);
			left1  = (newc[index]->orient==POSITIVE_ORIENTATION) ?
				  positive_component(newc[index]->curve) :
				  negative_component(newc[index]->curve);
		}
	}

	set_copy_intfc_states(YES);
	interpolate_intfc_states(intfc) = NO;
	set_interpolate_states_at_split_curve_node(NO);
	curs = split_curve(pcrs,newb0,newc[0]->curve,left0,right0,
				left1,right1);
	set_interpolate_states_at_split_curve_node(sav_scss);
	set_copy_intfc_states(sav_copy);
	interpolate_intfc_states(intfc) = sav_intrp;

	newn = curs[0]->end;
	nv = Node_vel(newn);
	dnv = Node_vel(diffnd);
	for (i = 0; i < dim; i++)
	{
		nv[i] = dnv[i];
		dnv[i] = sav_v[i];
	}
	if (newc0_orient == POSITIVE_ORIENTATION)
	{
		newc0 = newc[0]->curve = curs[0];
		end_status(curs[0]) = TRANSMITTED;
		start_status(curs[1]) = OVERTOOK;
	}
	else
	{
		newc0 = newc[0]->curve = curs[1];
		start_status(curs[1]) = TRANSMITTED;
		end_status(curs[0]) = OVERTOOK;
	}

	delete_interior_points_of_curve(fr,newc0);
	newb[0] = Bond_at_node_of_o_curve(newc[0]);

	left_state_along_bond(tcrsi,newbi,newc[index]->curve,lst);
	right_state_along_bond(tcrsi,newbi,newc[index]->curve,rst);
	change_node_of_curve(newc[index]->curve,newc[index]->orient,newn);
	set_status_at_node(newc[index]->curve,newc[index]->orient,INCIDENT);

	cut_curve(pcrs,newbi,newc[index]->curve,newc[index]->orient,fr,lst,rst);

	clear_state(fr->interf,lst,fr->sizest);
	clear_state(fr->interf,rst,fr->sizest);
	newbi = Bond_at_node_of_o_curve(newc[index]);
	newc[index]->curve = NULL;
	oldc[index]->curve = NULL;
	node_type(newn) = OVERTAKE_NODE;
	propagation_status(newn) = PROPAGATED_NODE;

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("diffraction"))
	{
		CURVE **c;

		(void) printf("NEW OVERTAKE NODE - ");		print_node(newn);
		(void) printf("CURVES AT NEWN\n");

		for (c = newn->in_curves; c && *c; c++)
		{
			(void) printf("In curve %llu - ",curve_number(*c));
			print_curve(*c);
		}
		for (c = newn->out_curves; c && *c; c++)
		{
			(void) printf("Out curve %llu - ",curve_number(*c));
			print_curve(*c);
		}

		(void) printf("END CURVES AT NEWN\n");
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	debug_print("diffraction","Left refl_curve_overtakes_incident_shock(), ");
	if (debugging("diffraction"))
	    print_node_status("status = ",status,"\n");
	return status;
}		/*end refl_curve_overtakes_incident_shock*/


LOCAL	void correct_for_refl_wave_under_shoot(
	Front		*fr,
	O_CURVE		*newc0,
	BOND		*newb00,
	POINT		*pc00,
	BOND		**newb0,
	POINT		*pc0)
{
	POINT		*p_opp;
	double		sc_len, s;
	double		sep;
	double		*h = fr->rect_grid->h;
	ORIENTATION	orient = newc0->orient;
	int		i, dim = fr->interf->dim;
	double		MIN_EXPANSION = 0.6*MIN_SC_SEP(fr->interf);/*TOLERANCE*/

	p_opp = Point_of_bond(newb00,Opposite_orient(orient));
	sep = scaled_separation(pc0,pc00,h,dim);

	if (sep > MIN_EXPANSION)
	{
		if (newb00 == *newb0)
		{
		    if (separation(pc0,p_opp,fr->interf->dim) <
				separation(pc00,p_opp,fr->interf->dim))
		        return;
		}
		else if ((orient == POSITIVE_ORIENTATION &&
			    (bonds_in_strict_order(newb00,*newb0) == YES))
		         ||
		         (orient == NEGATIVE_ORIENTATION &&
			    (bonds_in_strict_order(*newb0,newb00) == YES)))
			return;
	}

	*newb0 = newb00;

	sc_len = scaled_separation(pc00,p_opp,h,dim);

	if (sc_len < 2.0*MIN_EXPANSION)
	{
		for (i = 0; i < dim; i++)
		    Coords(pc0)[i] = 0.5*(Coords(p_opp)[i] + Coords(pc00)[i]);
		return;
	}

	s = MIN_EXPANSION/sc_len;
	for (i = 0; i < dim; i++)
	{
		Coords(pc0)[i] = Coords(pc00)[i] +
					s*(Coords(p_opp)[i] - Coords(pc00)[i]);
	}
}		/*end correct_for_refl_wave_under_shoot*/

LOCAL	int prop_char_up_incident_shock(
	NODE		*oldn,
	NODE		*newn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	BOND		**newb,
	POINT		*pc,
	Front		*fr,
	Wave		*wave,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	ANGLE_DIRECTION	ang_dir,
	NODE_FLAG	flag)
{
	static BOND	Bnew, Bvirtual, Bold;
	static POINT	Pnews, Pnewe, Polds, Polde, Ptmp;
	static POINT	*p0 = NULL, *p0_opp = NULL;
	static POINT	*ptmp = NULL;
	static WSSten	*sten = NULL;
	static const double SMALL_EXT_FAC = 1.5; /*TOLERANCE*/

	POINT		*pcrs;
	BOND		*newb0, *newb1;
	BOND		*oppb0;
	NODE		*oppn0;
	NODE		*interact_nodes[5];
	CURVE		*newc0 = newc[0]->curve;
	INTERFACE	*intfc = newc[0]->curve->interface;
	COMPONENT	compb;
	RECT_GRID	*rgr = fr->rect_grid;
	BOND		Newbdir;
	double		len, v1[MAXD], *nor, t[MAXD], t0[MAXD];
	double		*coords, V[MAXD], dtan;
	double		V0[MAXD];
	double		dt_tmp;
	ORIENTATION	newc0_orient = newc[0]->orient;
	boolean		sav_intrp = interpolate_intfc_states(intfc);
	int		status = ERROR_NODE;
	int		i, dim = intfc->dim;
	Locstate	ostb;
	double		M1, sinA1, cosA1;


	debug_print("diffraction","Entered prop_char_up_incident_shock()\n");

	if (ptmp == NULL)
	{
	    sten = AllocDefaultWSSten(fr);
	    p0 = Static_point(fr->interf);
	    p0_opp = Static_point(fr->interf);
	    ptmp = Static_point(fr->interf);
	    Bnew.start = &Pnews;	Bnew.end = &Pnewe;
	    Bold.start = &Polds;	Bold.end = &Polde;
	}
	else
	    ClearWSStenData(sten);
	sten->front = fr;
	sten->wave = wave;
	sten->w_type = FORWARD_SOUND_WAVE_LE;
	sten->pjump = 0.0;
	sten->hs = NULL;
	sten->dt = dt;

	if (is_forward_wave(wave_type(oldc[0]->curve)))
	{
	    ostb = Left_state_at_node_of_o_curve(oldc[0]);
	    compb = negative_component(oldc[0]->curve);
	}
	else
	{
	    ostb = Right_state_at_node_of_o_curve(oldc[0]);
	    compb = positive_component(oldc[0]->curve);
	}

	M1 = mach_number(ostb,Node_vel(oldn));

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("diffraction"))
	{
	    (void) printf("Old Node Velocity = < %g, %g>\n",
			  Node_vel(oldn)[0],Node_vel(oldn)[1]);
	    (void) printf("M1 = %g, M1 - 1.0 = %g\n",M1,M1-1.0);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	if (M1 < SONIC_MINUS)
	{
	    (void) printf("WARNING in prop_char_up_incident_shock(), "
		          "Subsonic state behind incident shock at oldn\n");
	    M1 = 1.0;
	    /**
	    interpolate_intfc_states(intfc) = sav_intrp;
	    status = ERROR_NODE;
	    debug_print("diffraction","Left prop_char_up_incident_shock(), ");
	    if (debugging("diffraction"))
		print_node_status("status = ",status,"\n");
	    return status;
	    **/
	}
	M1 = max(1.0,M1);

	nor = sten->nor;
	sinA1 = 1.0/M1;		cosA1 = sqrt(1.0 - sqr(sinA1));
	for (i = 0; i < dim; i++)
	    v1[i] = vel(i,ostb) - Node_vel(oldn)[i];
	if (ang_dir == COUNTER_CLOCK)
	{
	    nor[0] = -v1[0]*sinA1 + v1[1]*cosA1;
	    nor[1] = -v1[1]*sinA1 - v1[0]*cosA1;
	    len = mag_vector(nor,dim);
	    nor[0] /= len;			nor[1] /= len;
	    t[0] = -nor[1];			t[1] = nor[0];
	}
	else
	{
	    nor[0] = -v1[0]*sinA1 - v1[1]*cosA1;
	    nor[1] = -v1[1]*sinA1 + v1[0]*cosA1;
	    len = mag_vector(nor,dim);
	    nor[0] /= len;			nor[1] /= len;
	    t[0] = nor[1];			t[1] = -nor[0];
	}
	coords = Coords(Node_of_o_curve(oldc[0])->posn);
	states_near_location(sten,coords,nor,compb,compb,ostb,ostb);
	npt_w_speed(sten,left_state(ptmp),right_state(ptmp),V);
	for (i = 0; i < dim; i++)
	    Coords(ptmp)[i] = coords[i] + V[i]*dt;
	dtan = grid_size_in_direction(t,rgr->h,dim);

	for (i = 0; i < dim; i++)
	{
	    Coords(&Pnews)[i] = Coords(ptmp)[i] + SMALL_EXT_FAC*dtan*t[i];
	    Coords(&Pnewe)[i] = Coords(ptmp)[i] - SMALL_EXT_FAC*dtan*t[i];
	}
	newb1 = &Bnew;
	set_bond_length(newb1,dim);

	init_curve_for_crossing(p0,p0_opp,&Bvirtual,oldc[0],newc[0],
		                &oppn0,&oppb0,fr,wave,dt,V0,flag);
	find_bonds_for_extension_direction(&Bvirtual,oldc[0],newc[0],
		                           &Newbdir,NULL,fr);
	if (newc[0]->orient == POSITIVE_ORIENTATION)
	{
	    for (i = 0; i < dim; i++)
	    {
	        t0[i] = (Coords(Newbdir.start)[i] - Coords(Newbdir.end)[i])/
				bond_length(&Newbdir);
	        Coords(Bvirtual.start)[i] += SMALL_EXT_FAC*dtan*t0[i];
	    }
	}
	else
	{
	    for (i = 0; i < dim; i++)
	    {
	        t0[i] = (Coords(Newbdir.end)[i] - Coords(Newbdir.start)[i])/
				bond_length(&Newbdir);
	        Coords(Bvirtual.end)[i] += SMALL_EXT_FAC*dtan*t0[i];
	    }
	}
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("diffraction"))
	{
	    (void) printf("Before crossing check\n");
	    (void) printf("newb1 - \n");	print_bond(newb1);
	    (void) printf("Bvirtual - \n");	print_bond(&Bvirtual);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	pcrs = Point(NULL);
	if (!intersect_bond_with_curve_segment(newb1,&Bvirtual,
					Bond_at_opp_node_of_o_curve(newc[0]),
					newc[0],&newb0,pcrs,fr->rect_grid))
	{
	    NODE *endn_oldc0 = Opp_node_of_o_curve(oldc[0]);
	    NODE *endn_newc0 = Opp_node_of_o_curve(newc[0]);


#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("diffraction"))
	    {
	    	(void) printf("WARNING in prop_char_up_incident_shock(), "
	    	              "No intersection between propagated "
	    	              "incident shock and "
	    	              "behind characteristic\n");
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	    if ( propagation_status(endn_newc0) == UNPROPAGATED_NODE)
	    {
	    	status = PSEUDOCROSS_NODE_NODE;
	    	goto leave;
	    }

	    if (is_short_curve(oldc[0]->curve,oldc[0]->orient,rgr,1.0))
	    {
	    	for (i = 0; i < dim; i++)
	    	{
	    	    Coords(&Polds)[i] =
				Coords(oldn->posn)[i] + SMALL_EXT_FAC*dtan*t[i];
	    	    Coords(&Polde)[i] =
				Coords(oldn->posn)[i] - SMALL_EXT_FAC*dtan*t[i];
	    	}
	    	if (robust_cross_trace(rgr,endn_oldc0->posn,p0_opp,
				       &Bold,&Bnew,&dt_tmp,&Ptmp))
		{
		    *dt_frac = min(dt_tmp,*dt_frac);
		    interact_nodes[0] = newn;
		    interact_nodes[1] = oldn;
	    	    interact_nodes[2] = endn_newc0;
		    interact_nodes[3] = endn_oldc0;
		    interact_nodes[4] = NULL;
		    augment_rproblem_list(rp,interact_nodes,
					  dt,dt_tmp,oldc[0]->curve->interface,
					  newc[0]->curve->interface,fr,wave);
		    status = CROSS_NODE_NODE;
		    goto leave;
		}
	    }

	    (void) printf("ERROR in prop_char_up_incident_shock(), "
	                  "No intersection between propagated "
	                  "incident shock and "
	                  "behind characteristic, "
	                  "and no detected interaction\n");
	    status = ERROR_NODE;
	    goto leave;
	}
	if (newb0 == &Bvirtual) newb0 = Bond_at_node(newc0,newc0_orient);
	status = GOOD_NODE;
	correct_for_refl_wave_under_shoot(fr,newc[0],newb[0],pc,&newb0,pcrs);
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("diffraction"))
	{
	    (void) printf("Crossing bonds in prop_char_up_incident_shock()\n");
	    (void) printf("newb0 - ");	print_bond(newb0);
	    (void) printf("newb1 - ");	print_bond(newb1);
	    (void) printf("pcrs = <%g, %g>\n",Coords(pcrs)[0],Coords(pcrs)[1]);
	    (void) printf("newc0 before modification - ");
	    print_curve(newc0);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

leave:
	set_point_of_bond(oppn0->posn,oppb0,Opposite_orient(newc0_orient),dim);
	if (status == GOOD_NODE)
	{
	    interpolate_intfc_states(intfc) = YES;
	    if (insert_point_in_bond(pcrs,newb0,newc0) != FUNCTION_SUCCEEDED)
	    {
	        screen("ERROR in prop_char_up_incident_shock(), "
		       "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }
	    while(Point_adjacent_to_node(newc0,newc0_orient) != pcrs)
	    	(void) delete_point_adjacent_to_node(fr,newc0,newc0_orient);
	    newb[0] = Bond_at_node_of_o_curve(newc[0]);
	    if (is_forward_wave(wave_type(newc0)))
	    {
	    	ft_assign(left_state(pcrs),left_state(ptmp),fr->sizest);
	    }
	    else
	    {
	    	ft_assign(right_state(pcrs),left_state(ptmp),fr->sizest);
	    }
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("diffraction"))
	    {
	    	(void) printf("newc0 after modification - ");
	    	print_curve(newc0);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	}
	interpolate_intfc_states(intfc) = sav_intrp;
	debug_print("diffraction","Left prop_char_up_incident_shock(), ");
	if (debugging("diffraction"))
	    print_node_status("status = ",status,"\n");
	return status;
}		/*end prop_char_up_incident_shock*/


LOCAL	ANOM *allocate_static_ANOM_data_structure(
	Front		*fr)
{
	ANOM		*anom;

	scalar(&anom,sizeof(ANOM));

	alloc_state(fr->interf,&anom->osb,fr->sizest);
	scalar(&anom->binc,sizeof(BOND));
	anom->binc->prev = anom->binc->next = NULL;
	scalar(&anom->newp,sizeof(POINT));
	anom->p00 = Static_point(fr->interf);
	anom->wssten = AllocDefaultWSSten(fr);

	return anom;
}		/*end allocate_static_ANOM_data_structure*/

LOCAL	void init_ANOM_data_structure(
	ANOM		*anom,
	NODE		*oldn,
	NODE		*newn,
	NODE		**oppn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	BOND		**newb,
	BOND		**oppb,
	POINT		*pc,
	POINT		*ploc,
	POINT		*ploc_opp,
	POINT		**p_sav,
	Front		*fr,
	Wave		*wave,
	double		dt,
	NODE_FLAG	flag,
	RP_DATA		*RP,
	double		*virtual_len,
	double		**t)
{
	WSSten		*wssten = anom->wssten;
	INTERFACE	*intfc = newn->interface;
	CURVE		*oldc0 = oldc[0]->curve;
	CURVE		*c6;
	POINT		*oldp;
	BOND		Bvirtual, *bvertex;
	double		*nor;		/* components of unit normal to front */
	double		dn;		/* mesh spacing in normal direction */
	double		dir[MAXD], ds;
	double		V[MAXD];
	double		*h = fr->rect_grid->h;
	boolean		sav_intrp = interpolate_intfc_states(intfc);
	ORIENTATION	oldc0_orient = oldc[0]->orient;
	ORIENTATION	c6_orient;
	int		dim = intfc->dim;
	int		i, j;
	size_t		sizest = fr->sizest;

	c6 = newc[6]->curve;	c6_orient = newc[6]->orient;

	anom->pc = pc;	anom->newc6 = newc[6];
	anom->rpst0 = RP->state[0];	anom->rpst1 = RP->state[1];
	wssten->dt = dt;
	wssten->front = fr;
	wssten->wave = wave;
	*p_sav = Point_of_bond(Bond_at_node(c6,c6_orient),c6_orient);
	init_curve_for_crossing(ploc,ploc_opp,&Bvirtual,oldc[6],newc[6],
		                oppn,oppb,fr,wave,dt,V,flag);
	*virtual_len = bond_length(&Bvirtual);
	set_point_of_bond(ploc,Bond_at_node(c6,c6_orient),c6_orient,dim);
	find_tangent_to_curve(ploc,Bond_at_node(c6,c6_orient),c6,c6_orient,
		              dir,fr);
	ds = grid_size_in_direction(dir,h,dim);
	for (i = 0; i < dim; i++)
	    Coords(ploc)[i] -= 2.0*ds*dir[i];

	interpolate_intfc_states(intfc) = YES;
	for (i = 0; i < dim; i++)
	    Coords(anom->p00)[i] = Coords(pc)[i];
	if (insert_point_in_bond(anom->p00,newb[6],newc[6]->curve) !=
	    FUNCTION_SUCCEEDED)
	{
	    screen("ERROR in init_ANOM_data_structure(), "
		   "insert_point_in_bond() failed\n");
	    clean_up(ERROR);
	}
	interpolate_intfc_states(intfc) = sav_intrp;
	anom->b_ct[0] = newb[6];	anom->b_ct[1] = newb[6]->next;
	anom->pvertex = Point_adjacent_to_node(newc[0]->curve,newc[0]->orient);
	bvertex = Bond_at_node_of_o_curve(newc[0]);
	anom->rpst6 = RP->state[6];
	for (i = 0; i < dim; i++)
	    anom->abs_v00[i] = Node_vel(newn)[i];
	anom->t = t;
	anom->inc_side = (curve_ang_oriented_l_to_r(RP->ang_dir,c6_orient)) ?
				POSITIVE_SIDE : NEGATIVE_SIDE;

	anom->oldp = oldp = oldn->posn;
	anom->newc0 = newc[0];
	anom->oabs_v = Node_vel(oldn);
	anom->is_plus_or = (RP->ang_dir == COUNTER_CLOCK) ? YES : NO;
	wssten->w_type = wave_type(oldc0);
	wssten->hs = NULL;
	if (wave_type(newc[0]->curve) == wave_type(oldc[0]->curve))
	{
	    left_state(anom->newp) = Left_state_at_node_of_o_curve(newc[0]);
	    right_state(anom->newp) = Right_state_at_node_of_o_curve(newc[0]);
	}
	else
	{
	    right_state(anom->newp) = Left_state_at_node_of_o_curve(newc[0]);
	    left_state(anom->newp) = Right_state_at_node_of_o_curve(newc[0]);
	}
	if (newc[0]->orient == POSITIVE_ORIENTATION)
	{
	    anom->binc->start = anom->newp;
	    anom->binc->end =
	        Point_adjacent_to_node(newc[0]->curve,newc[0]->orient);
	}
	else
	{
	    anom->binc->end = anom->newp;
	    anom->binc->start =
	        Point_adjacent_to_node(newc[0]->curve,newc[0]->orient);
	}

	/* Calculate normal direction */
	/* and positions of normally displaced states */

	nor = wssten->nor;
	normal(oldp,Hyper_surf_element(Bond_at_node(oldc0,oldc0_orient)),
	       Hyper_surf(oldc0),nor,fr);

	wssten->dn = dn = grid_size_in_direction(nor,h,dim);
	wssten->pcomp = wssten->ncomp = NO_COMP;
	for (i = 0; i < dim; i++)
	{
	    wssten->coords[i] = wssten->lcrds[0][i] = wssten->rcrds[0][i] =
		Coords(oldp)[i];
	}
	for (i = 1; i < wssten->nsts; i++)
	{
	    for (j = 0; j < dim; j++)
	    {
	        wssten->lcrds[i][j] = Coords(oldp)[j] - i*nor[j]*dn;
	        wssten->rcrds[i][j] = Coords(oldp)[j] + i*nor[j]*dn;
	    }
	}


	if (is_forward_wave(wave_type(oldc0)))
	{
	    anom->osb0 = Left_state_at_node(oldc0,oldc0_orient);
	    wssten->sr[0] = Right_state_at_node(oldc0,oldc0_orient);
	    for (i = 1; i < wssten->nsts; i++)
	    {
		wssten->sr[i] = (Locstate) (wssten->sr_store + i*sizest);
	        hyp_solution(wssten->rcrds[i],positive_component(oldc0),
			     Hyper_surf(oldc0),POSITIVE_SIDE,fr,
			     wave,wssten->sr[i],wssten->sr[0]);
	    }
	    /*wssten->sl[0] = Left_state_at_node(oldc0,oldc0_orient);*/
	    for (i = 0; i < wssten->nsts; i++)
	        wssten->sl[i] = anom->osb;
	}
	else
	{
	    anom->osb0 = Right_state_at_node(oldc0,oldc0_orient);
	    wssten->sl[0] = Left_state_at_node(oldc0,oldc0_orient);
	    for (i = 1; i < wssten->nsts; i++)
	    {
	        wssten->sl[i] = (Locstate) (wssten->sl_store + i*sizest);
	        hyp_solution(wssten->lcrds[i],negative_component(oldc0),
			     Hyper_surf(oldc0),NEGATIVE_SIDE,fr,
			     wave,wssten->sl[i],wssten->sl[0]);
	    }
	    /*wssten->sr[0] = Right_state_at_node(oldc0,oldc0_orient);*/
	    for (i = 0; i < wssten->nsts; i++)
	        wssten->sr[i] = anom->osb;
	}

	anom->sot = (is_forward_wave(wave_type(newc[0]->curve))) ?
	    left_state_at_point_on_curve(anom->pvertex,bvertex,newc[0]->curve) :
	    right_state_at_point_on_curve(anom->pvertex,bvertex,newc[0]->curve);
			
}		 /*end init_ANOM_data_structure*/

/*ARGSUSED*/
LOCAL	int anom_states_by_pressure(
	NODE		*oldn,
	NODE		*newn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	BOND		**newb,
	POINT		*pc,
	Front		*fr,
	Wave		*wave,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag,
	RP_DATA		*RP,
	double		**t)
{
	CURVE		*c6;
	POINT		*p_sav, *pt;
	NODE		*oppn;
	BOND		*oppb;
	double		p, pupper, plower, M;
	double		epsilon, delta;
	double		para;
	double		virtual_len;
	double		va[MAXD], qa;
	int		i, dim = fr->interf->dim;
	ORIENTATION	c6_orient;
	static	ANOM	*anom = NULL;
	static	WSSten	*wssten = NULL;
	static	POINT	*ploc = NULL, *ploc_opp = NULL;


	debug_print("diffraction","Entered anom_states_by_pressure()\n");

		/* Allocate storage */

	if (ploc_opp == NULL)
	{
	    anom = allocate_static_ANOM_data_structure(fr);
	    wssten = anom->wssten;
	    ploc = Static_point(fr->interf);
	    ploc_opp = Static_point(fr->interf);
	}
	else
	    ClearWSStenData(anom->wssten);



#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("diffraction"))
	{
	    double pr1, prsonic, M0sq;

	    pr1 = pressure(RP->state[1]);
	    M0sq = mach_number_squared(RP->state[0],
			Node_vel(newn),(double *)NULL);
	    prsonic = pressure_at_sonic_point(M0sq,RP->state[0]);
	    (void) printf("Original behind pressure = %g, ",pr1);
	    (void) printf("sonic point pressure = %g,\n",prsonic);
	    (void) printf("\t\t\t\tdifference = %g\n",pr1 - prsonic);
	    (void) printf("Old node velocity = <%g, %g>\n",
			Node_vel(oldn)[0],Node_vel(oldn)[1]);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	init_ANOM_data_structure(anom,oldn,newn,&oppn,oldc,newc,newb,&oppb,
		                 pc,ploc,ploc_opp,&p_sav,fr,wave,dt,flag,RP,
				 &virtual_len,t);

	pupper = pressure(anom->osb0);

	plower = pressure(Rp_data(oldn)->state[4]);

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("diffraction")) 
		(void) printf("plower = %g, pupper = %g\n",plower,pupper);
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	epsilon = 0.5*SONIC_TOL;	delta = EPS*fabs(pupper - plower);

	(void) f_pr_anom(plower,&M,(POINTER)anom);
	if (M < SONIC_MINUS)
	{
	    if (is_forward_wave(wave_type(oldc[0]->curve)))
	        wssten->sl[0] = anom->osb;
	    else
		wssten->sr[0] = anom->osb;
	}

	if (find_root(f_pr_anom,(POINTER)anom,1.0,&p,plower,
		      pupper,epsilon,delta) == FUNCTION_FAILED)
	{
	    (void) printf("WARNING in anom_states_by_pressure(), "
	                  "Unable to find state behind incident shock\n");
	    debug_print("diffraction","Left anom_states_by_pressure()\n");
	    return ERROR_NODE;
	}

	for (i = 0; i < dim; i++)
	    Node_vel(newn)[i] = anom->abs_v[i];

	if (newc[0]->orient == POSITIVE_ORIENTATION)
	{
	    for (i = 0; i < dim; i++)
	    {
	        t[0][i] = (Coords(anom->binc->end)[i] - 
			   Coords(anom->binc->start)[i]) /
				   bond_length(anom->binc);
	    }
	}
	else
	{
	    for (i = 0; i < dim; i++)
	    {
	        t[0][i] = (Coords(anom->binc->start)[i] -
			   Coords(anom->binc->end)[i]) /
				   bond_length(anom->binc);
	    }
	}

	newb[6] = (anom->bc == anom->b_ct[1]) ? newb[6] : anom->bc;
	c6 = newc[6]->curve;	c6_orient = newc[6]->orient;
	(void) delete_start_of_bond(anom->b_ct[1],c6);
	set_point_of_bond(p_sav,Bond_at_node(c6,c6_orient),c6_orient,dim);
	set_point_of_bond(oppn->posn,oppb,Opposite_orient(c6_orient),dim);

	for (i = 0; i < dim; i++)
	    va[i] = vel(i,RP->state[0]) - Node_vel(newn)[i];
	qa = mag_vector(va,dim);
	for (i = 0; i < dim; i++)
	    t[1][i] = -va[i]/qa;

	pt = Point_of_bond(newb[6],Opposite_orient(newc[6]->orient));
	if (newb[6] != Bond_at_node_of_o_curve(newc[6]))
	    virtual_len = bond_length(newb[6]);
	para = min(separation(pc,pt,dim)/virtual_len,1.0);
	if (newc[6]->orient == POSITIVE_ORIENTATION)
	    para = 1.0 - para;


	if (curve_ang_oriented_l_to_r(RP->ang_dir,c6_orient))
	{
	    left_state_along_bond(para,newb[6],c6,RP->state[6]);
	}
	else
	{
	    right_state_along_bond(para,newb[6],c6,RP->state[6]);
	}

	debug_print("diffraction","Left anom_states_by_pressure()\n");
	return GOOD_NODE;
}		/*end anom_states_by_pressure*/


LOCAL	boolean f_pr_anom(
	double		p,
	double		*M,
	POINTER		panom)
{
	ANOM		*anom = (ANOM *) panom;
	WSSten		*wssten = anom->wssten;
	Front		*fr = wssten->front;
	POINT		*newp = anom->newp;
	POINT		*oldp = anom->oldp;
	POINT		*pc = anom->pc;
	O_CURVE		*newc0 = anom->newc0;
	O_CURVE		*newc6 = anom->newc6;
	BOND		*binc = anom->binc;
	Locstate	osb0 = anom->osb0;
	Locstate	osb = anom->osb;
	Locstate	rpst0 = anom->rpst0, rpst1 = anom->rpst1;
	double		dt = wssten->dt;
	double		*nor = wssten->nor;
	double		*oabs_v = anom->oabs_v;
	int		w_type = wssten->w_type;
	int		is_plus_or = anom->is_plus_or;
	ORIENTATION	c0_or = newc0->orient;
	int		i, dim = fr->interf->dim;
	double		p1;
	double		para;
	double		ang0, ang1, theta;
	double		V[MAXD];
	double		ds, len, dir[MAXD], *s;
	double		vb[MAXD], qbsqr;
	static Locstate stdummy = NULL;

#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("fanom","Entered f_pr_anom(), p = %g\n",p);
#endif /* defined(DEBUG_NODE_PROPAGATE) */


	if (stdummy == NULL)
	{
	    alloc_state(fr->interf,&stdummy,fr->sizest);
	}

	if (!prandtl_meyer_wave(osb0,p,(is_plus_or)?NO:YES,
			           oabs_v,osb,&ang0,&ang1,&theta)) 
	{
#if defined(DEBUG_NODE_PROPAGATE)
	    (void) printf("WARNING in f_pr_anom(), "
	                  "prandtl_meyer_wave() failed\n");
	    debug_print("fanom","Left f_pr_anom()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return FUNCTION_FAILED;
	}

	npt_w_speed(wssten,left_state(newp),right_state(newp),V);

	for (i = 0; i < dim; i++)
	    Coords(newp)[i] = Coords(oldp)[i] + V[i]*dt;

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("fanom"))
	{
	    char	s[80];
	    verbose_print_state("osb0",osb0);
	    (void) printf("State data into npt_w_speed()\n");
	    (void) printf("nor = %g, %g, dt = %g\n",nor[0],nor[1],dt);
	    print_wave_type("w_type = ",w_type,"\n",fr->interf);
	    for (i = wssten->nsts-1; i >= 0; i--)
	    {
		(void)sprintf(s,"sl[%d]",i);
	        verbose_print_state(s,wssten->sl[i]);
	    }
	    for (i = 0; i < wssten->nsts; i++)
	    {
		(void)sprintf(s,"sr[%d]",i);
	        verbose_print_state(s,wssten->sr[i]);
	    }
	    (void) printf("Answer\n");
	    (void) printf("V = %g, %g\n",V[0],V[1]);
	    (void) printf("Coords(newp) = <%g, %g>\n",
	    	          Coords(newp)[0],Coords(newp)[1]);
	    verbose_print_state("left_state(newp)",left_state(newp));
	    verbose_print_state("right_state(newp)",right_state(newp));
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("fanom"))
	{
	    set_bond_length(binc,dim);
	    (void) printf("binc before extension - ");
	    print_bond(binc);
	    print_orientation("c0_or = ",c0_or,"\n");
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	len = separation(binc->start,binc->end,dim);
	s = Coords(binc->start);
	for (i = 0; i < dim; i++)
	    dir[i] = (Coords(binc->end)[i] - s[i])/len;
	ds = grid_size_in_direction(dir,fr->rect_grid->h,dim);

	if (c0_or == POSITIVE_ORIENTATION)
	{
	    for (i = 0; i < dim; i++)
	    	Coords(newp)[i] -= 2.0*ds*dir[i];
	}
	else
	{
	    for (i = 0; i < dim; i++)
	    	Coords(newp)[i] += 2.0*ds*dir[i];
	}
	set_bond_length(binc,dim);

	if (!intersect_bond_with_curve_segment(binc,
	                     Bond_at_node_of_o_curve(newc6),
	                     Bond_at_opp_node_of_o_curve(newc6),
	                     newc6,&anom->bc,pc,fr->rect_grid))
	{
#if defined(DEBUG_NODE_PROPAGATE)
	    (void) printf("WARNING in f_pr_anom(), "
	                  "No intersection of bond with curve\n");
	    (void) printf("binc\n");
	    print_bond(binc);
	    (void) printf("newc6\n");
	    print_o_curve(newc6);
	    debug_print("fanom","Left f_pr_anom()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return FUNCTION_FAILED;
	}
	for (i = 0; i < dim; i++)
	    anom->abs_v[i] = (Coords(pc)[i] - Coords(oldp)[i])/dt;
	para = (scalar_product(Coords(pc),dir,dim) -
	    	scalar_product(s,dir,dim))/len;
	if (para < 0.0)
	    para = 0.0;
	else if (para > 1.0)
	    para = 1.0;

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("fanom"))
	{
	    BOND *b = Bond_at_node_of_o_curve(newc0);
	    CURVE *c = newc0->curve;

	    (void) printf("Cross point <%g, %g>\n",Coords(pc)[0],Coords(pc)[1]);
	    (void) printf("binc after extension - ");
	    print_bond(binc);
	    (void) printf("Crossing bond on contact - ");
	    print_bond(anom->bc);
	    (void) printf("para = %g\n",para);
	    (void) printf("Calculated node velocity = <%g, %g>\n",
			  anom->abs_v[0],anom->abs_v[1]);

	    (void) printf("Bond at node of newc0 - "); print_bond(b);
	    verbose_print_state("Left state b->start",
			        left_state_at_point_on_curve(b->start,b,c));
	    verbose_print_state("Right state b->start",
			        right_state_at_point_on_curve(b->start,b,c));
	    verbose_print_state("Left state b->end",
			        left_state_at_point_on_curve(b->end,b,c));
	    verbose_print_state("Right state b->end",
			        right_state_at_point_on_curve(b->end,b,c));
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	if (is_forward_wave(wave_type(newc0->curve)))
	{
	    right_state_along_bond(para,Bond_at_node_of_o_curve(newc0),
			           newc0->curve,rpst0);
	    left_state_along_bond(para,Bond_at_node_of_o_curve(newc0),
			          newc0->curve,rpst1);
	}
	else
	{
	    left_state_along_bond(para,Bond_at_node_of_o_curve(newc0),
			          newc0->curve,rpst0);
	    right_state_along_bond(para,Bond_at_node_of_o_curve(newc0),
			           newc0->curve,rpst1);
	}
	p1 = pressure(rpst1);

	if (!s_polar_3(rpst0,YES,p1,is_plus_or,NO,anom->abs_v,
		          rpst1,&ang0,&theta))
	{
	    (void) printf("WARNING - in f_pr_anom(), "
	                  "s_polar_3() failed\n");
	    debug_print("fanom","Left f_pr_anom()\n");
	    return FUNCTION_FAILED;
	}
	for (i = 0; i < dim; i++)
	    dir[i] = Coords(pc)[i] - Coords(anom->pvertex)[i];
	len = mag_vector(dir,dim);
	dir[0] /= len;			dir[1] /= len;
	w_speed(Coords(anom->pvertex),anom->sot,rpst1,stdummy,rpst1,V,0.0,dir,
	        FORWARD_SOUND_WAVE_TE,fr);

	for (i = 0; i < dim; i++)
	    vb[i] = vel(i,rpst1) - anom->abs_v[i]; 
	qbsqr = scalar_product(vb,vb,dim);

	*M = sqrt(qbsqr/sound_speed_squared(rpst1));

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("fanom"))
	    (void) printf("p = %g, p1 = %g, M = %g\n",p,p1,*M);
	debug_print("fanom","Left f_pr_anom()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	return FUNCTION_SUCCEEDED;
}		/*end f_pr_anom*/


/*ARGSUSED*/
LOCAL	int anom_states_by_mass_flux(
	NODE		*oldn,
	NODE		*newn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	BOND		**newb,
	POINT		*pc,
	Front		*fr,
	Wave		*wave,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag,
	RP_DATA		*RP,
	double		**t)
{
	NODE		*oppn;
	POINT		*p_sav;
	CURVE		*c6;
	BOND		*oppb;
	double		ds, ds_min, ds_max;
	double		epsilon, delta;
	double		virtual_len;
	int		status;
	ORIENTATION	c6_orient;
	int		i, dim = fr->interf->dim;
	static	ANOM	*anom = NULL;
	static	POINT	*ploc = NULL, *ploc_opp = NULL;

#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("diffraction","Entered anom_states_by_mass_flux()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	if (ploc_opp == NULL)
	{
		anom = allocate_static_ANOM_data_structure(fr);
		ploc = Static_point(fr->interf);
		ploc_opp = Static_point(fr->interf);
	}
	else
	    ClearWSStenData(anom->wssten);

	init_ANOM_data_structure(anom,oldn,newn,&oppn,oldc,newc,newb,&oppb,
		                 pc,ploc,ploc_opp,&p_sav,fr,wave,dt,flag,RP,
				 &virtual_len,t);


	switch(set_limits_for_anom_find_root(oldn,fr,dt,RP,anom,&ds_min,&ds_max,
			                     &epsilon,&delta))
	{
	case SUPERSONIC_BY_MASS_FLUX:
	    status = GOOD_NODE;
	    goto leave;

	case LIMITS_SET:
	    if (find_root(f_mf_anom,(POINTER)anom,1.0,&ds,ds_min,
			  ds_max,epsilon,delta) == FUNCTION_FAILED)
	    {
	        (void) printf("WARNING in anom_states_by_mass_flux(), "
	                      "Unable to find state behind incident shock\n");
	        status = ERROR_NODE;
	        goto leave;
	    }
	    status = GOOD_NODE;
	    break;
	
	case UNABLE_TO_SET_LIMITS:
	default:
	    (void) printf("WARNING in anom_states_by_mass_flux(), "
	                  "Unable to find state behind incident shock\n"
	                  "set_limits_for_anom_find_root() failed\n");
	    status = ERROR_NODE;
	    goto leave;
	}


#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("diffraction"))
	{
		double dang;

		dang = Rp_data(oldn)->ang[3] - Rp_data(oldn)->ang[1];
		dang = (Rp_data(oldn)->ang_dir == CLOCKWISE) ?
			normalized_angle(-dang) : normalized_angle(dang);
		
		(void) printf("dang = %g, dt = %g, ds = %g\n",dang,dt,ds);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	newb[6] = (anom->bc == anom->b_ct[1]) ? newb[6] : anom->bc;
	for (i = 0; i < dim; i++) Node_vel(newn)[i] = anom->abs_v[i];

leave:
	c6 = newc[6]->curve;	c6_orient = newc[6]->orient;
	(void) delete_start_of_bond(anom->b_ct[1],c6);
	set_point_of_bond(p_sav,Bond_at_node(c6,c6_orient),c6_orient,dim);
	set_point_of_bond(oppn->posn,oppb,Opposite_orient(c6_orient),dim);

#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("diffraction","Left anom_states_by_mass_flux(), ");
	if (debugging("diffraction"))
	    print_node_status("status = ",status,"\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	return status;
}		/*end anom_states_by_mass_flux*/


LOCAL	SL_AFR set_limits_for_anom_find_root(
	NODE		*oldn,
	Front		*fr,
	double		dt,
	RP_DATA		*RP,
	ANOM		*anom,
	double		*pds_min,
	double		*pds_max,
	double		*epsilon,
	double		*delta)
{
	static const int	MAX_ITER = 20;
	static const double	EFAC = 0.5;	/* TOLERANCE:  must be < 1 */
	static const double	DFAC = 2.0;	/* TOLERANCE:  should be >= 1 */
	double		Mds0, Mds_mid, Mds_min, Mds_max;
	double		Pmin, Pmid, Pmax, P1ds0, Plower, Pupper;
	double		ds0, ds_mid, ds_min, ds_max, dslen;
	Locstate	rpst1 = anom->rpst1;
	int		Mds_max_set, Mds_min_set;
	int		i, j;
	SL_AFR		status;

#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("fanom","Entered set_limits_for_anom_find_root()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	Plower = pressure(RP->state[0]);
	Pupper = pressure(RP->state[1]);
	ds0 = 0.0;
	if (f_mf_anom(ds0,&Mds0,(POINTER)anom) == FUNCTION_FAILED)
	{
	    (void) printf("WARNING in set_limits_for_anom_find_root(), "
	                  "Unable to evaluate f_mf_anom() for ds = %g\n",ds0);
	    status = UNABLE_TO_SET_LIMITS;
	    goto leave;
	}
	if (Mds0 >= 1.0)
	{
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("fanom"))
	    	(void) printf("M1 > 1 at ds = 0.0\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    status = SUPERSONIC_BY_MASS_FLUX;
	    goto leave;
	}
	P1ds0 = pressure(rpst1);
	if (P1ds0 < Plower)
	{
	    (void) printf("WARNING in set_limits_for_anom_find_root(), "
	                  "Invalid base node velocity\n"
	                  "This should have been previously detected\n");
	    status = UNABLE_TO_SET_LIMITS;
	    goto leave;
	}
	Pupper = max(Pupper,P1ds0);

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("fanom"))
	{
	    (void) printf("\nPlower = %g, Pupper = %g\n",Plower,Pupper);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */


	{
	    double p4, M;
	    double *h = fr->rect_grid->h;
	    int dim = fr->interf->dim;
	    int max_num_iter = 5;/*TOLERANCE*/
	    int i;

	    p4 = pressure(Rp_data(oldn)->state[4]);
	    for (i = 0; i < max_num_iter; i++)
	    {
#if defined(DEBUG_NODE_PROPAGATE)
	    	debug_print("fanom","Calling f_pr_anom() with p4 = %g\n",p4);
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    	if (f_pr_anom(p4,&M,(POINTER)anom)==FUNCTION_SUCCEEDED)
	    		break;
	    	p4 *= 0.95;/*TOLERANCE*/
	    }
	    if (i == max_num_iter)
	    	dslen = fabs(sqr(Mds0) - 1.0) * dt / (1.0 + dt);
	    else
	    	dslen = 0.5*scaled_separation(anom->pc,anom->p00,h,dim);
	}


#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("fanom"))
	{
	    (void) printf("Initial value of dslen = %g\n",dslen);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	ds_min = -dslen;	Mds_min_set = NO;
	ds_max = dslen;		Mds_max_set = NO;
	ds_mid = 0.0;		Mds_mid = Mds0;		Pmid = P1ds0;
	for (i = 0; i < MAX_ITER; i++)
	{
	    for (j = 0; j < MAX_ITER; j++)
	    {
	        if (Mds_max_set) break;
	        if (f_mf_anom(ds_max,&Mds_max,(POINTER)anom) ==
		    FUNCTION_SUCCEEDED)
	        {
	            Pmax = pressure(rpst1);
	            if (Mds_max < Mds_mid)
			break;
	            if (Between(Pmax,Plower,Pupper))
			break;
	        }
	        else
	        {
	            (void) printf("WARNING in set_limits_for_anom_find_root()"
	                          ", Unable to evaluate f_mf_anom() "
	                          "for ds = %g\n",ds_max);
	        }
	        ds_max = 0.5 * (ds_mid + ds_max);
	    }
	    if (j >= MAX_ITER)
	    {
	        (void) printf("WARNING in set_limits_for_anom_find_root(), "
	                      "Unable to find acceptable ds_max after "
	                      "%d iterations\n",j);
	        status = UNABLE_TO_SET_LIMITS;
	        goto leave;
	    }
	    Mds_max_set = YES;
	    for (j = 0; j < MAX_ITER; j++)
	    {
	        if (Mds_min_set) break;
	        if (f_mf_anom(ds_min,&Mds_min,(POINTER)anom) ==
		    FUNCTION_SUCCEEDED)
	        {
	            Pmin = pressure(rpst1);
	            if (Mds_min < Mds_mid) break;
	            if (Between(Pmin,Plower,Pupper)) break;
	        }
	        else
	        {
	            (void) printf("WARNING in "
	                          "set_limits_for_anom_find_root(), "
	                          "Unable to evaluate f_mf_anom() "
	                          "for ds = %g\n",ds_min);
	        }
	        ds_min = 0.5 * (ds_mid + ds_min);
	    }
	    if (j >= MAX_ITER)
	    {
	        (void) printf("WARNING in set_limits_for_anom_find_root(), "
	                      "Unable to find acceptable ds_min after "
	                      "%d iterations\n",j);
	        status = UNABLE_TO_SET_LIMITS;
	        goto leave;
	    }
	    Mds_min_set = YES;

	    if (Mds_mid < 1.0)
	    {
	        if (Mds_min > 1.0)
	        {
	            ds_max = ds_mid;
	            break;
	        }
	        if (Mds_max > 1.0)
	        {
	            ds_min = ds_mid;
	            break;
	        }
	        if (Mds_mid < Mds_min && Mds_mid > Mds_max)
	        {
	            ds_max = ds_mid; Mds_max = Mds_mid; Pmax = Pmid;
	            ds_mid = ds_min; Mds_mid = Mds_min; Pmid = Pmin;
	            ds_min = 2.0*ds_mid - ds_max;
	            Mds_min_set = NO;
	        }
	        else if (Mds_mid > Mds_min && Mds_mid < Mds_max)
	        {
	            ds_min = ds_mid; Mds_min = Mds_mid; Pmin = Pmid;
	            ds_mid = ds_max; Mds_mid = Mds_max; Pmid = Pmax;
	            ds_max = 2.0*ds_mid - ds_min;
	            Mds_max_set = NO;
	        }
	        else
	        {
	            dslen *= 0.5;
	            ds_min = -dslen; Mds_min_set = NO;
	            ds_max = dslen;  Mds_max_set = NO;
	            ds_mid = 0.0;
		    Mds_mid = Mds0;
		    Pmid = P1ds0;
	        }
	    }
	    else if (Mds_min < 1.0)
	    {
	        ds_max = ds_mid;
	        break;
	    }
	    else if (Mds_max < 1.0)
	    {
	        ds_min = ds_mid;
	        break;
	    }
	    else
	    {
	        dslen *= 0.5;
	        ds_min = -dslen; Mds_min_set = NO;
	        ds_max =  dslen; Mds_max_set = NO;
	        ds_mid = 0.0;
		Mds_mid = Mds0;
		Pmid = P1ds0;
	    }
	}
	if (i >= MAX_ITER)
	{
	    (void) printf("WARNING in set_limits_for_anom_find_root(), "
	                  "Uable to find limits after %d iterations\n",i);
	    status = UNABLE_TO_SET_LIMITS;
	    goto leave;
	}

	*pds_min = ds_min;		*pds_max = ds_max;
	*epsilon = EFAC*SONIC_TOL;	*delta = DFAC*EPS*(ds_max - ds_min);
	status = LIMITS_SET;

leave:
#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("fanom","Left set_limits_for_anom_find_root()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	return status;
}		/*end set_limits_for_anom_find_root*/

LOCAL	boolean f_mf_anom(
	double		ds,
	double		*M,
	POINTER		panom)
{
	ANOM		*anom = (ANOM *) panom;
	WSSten		*wssten = anom->wssten;
	BOND		*b_ct = anom->b_ct[0];
	BOND		*curr_b;
	CURVE		*c_ct = anom->newc6->curve;
	CURVE		*curr_c;
	Front		*fr = wssten->front;
	Locstate	rpst0 = anom->rpst0;
	Locstate	rpst6 = anom->rpst6;
	Locstate	rpst1 = anom->rpst1;
	POINT		*p00 = anom->p00;
	POINT		*pvertex = anom->pvertex;
	POINT		*pc = anom->pc;
	POINT		*ptmp;
	double		*abs_v00 = anom->abs_v00;
	double		dt = wssten->dt;
	double		ang0, theta, V[MAXD];
	double		va[MAXD]; 
	double		d[MAXD], len, sinsqr, qasqr, qa, mf_sqr;
	double		abs_v[MAXD];
	double		rhoa;
	double		nor[MAXD];
	double		para;
	ORIENTATION	c_ct_orient = anom->newc6->orient;
	SIDE		inc_side = anom->inc_side;
	int		is_plus_or = anom->is_plus_or;
	int		i, dim = fr->interf->dim;
	static BOND	*bdir = NULL;
	static Locstate stdummy = NULL;
	static double	**t = NULL;

#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("fanom","\nEntered f_mf_anom(), ds = %g\n",ds);
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	if (stdummy == NULL)
	{
		alloc_state(fr->interf,&stdummy,fr->sizest);
		scalar(&bdir,sizeof(BOND));
		bi_array(&t,2,2,FLOAT);
		bdir->start = Static_point(fr->interf);
		bdir->end = Static_point(fr->interf);
	}

	bond_secant_to_curve(p00,b_ct,c_ct,c_ct_orient,bdir,fr,ds);

	curr_c = c_ct;
	if (c_ct_orient == POSITIVE_ORIENTATION)
	{
		ptmp = bdir->end;
		curr_b = bdir->next;
		for (i = 0; i < dim; i++)
			Coords(pc)[i] = Coords(bdir->end)[i];
	}
	else
	{
		ptmp = bdir->start;
		curr_b = bdir->prev;
		for (i = 0; i < dim; i++)
			Coords(pc)[i] = Coords(bdir->start)[i];
	}
	anom->bc = (curr_c == c_ct) ? curr_b :
				Bond_at_node(c_ct,Opposite_orient(c_ct_orient));


	para = (bond_length(curr_b) < 0.00001) ? 0.5 :	/*TOLERANCE*/
		separation(ptmp,curr_b->start,fr->interf->dim) /
			bond_length(curr_b);

	if (inc_side == NEGATIVE_SIDE)
	{
		left_state_along_bond(para,curr_b,curr_c,rpst0);
		right_state_along_bond(para,curr_b,curr_c,rpst6);
	}
	else
	{
		right_state_along_bond(para,curr_b,curr_c,rpst0);
		left_state_along_bond(para,curr_b,curr_c,rpst6);
	}


	for (i = 0; i < dim; i++) t[0][i] = Coords(pvertex)[i] - Coords(pc)[i];
	len = mag_vector(t[0],dim);
	for (i = 0; i < dim; i++)
	{
		t[0][i] /= len;
		anom->t[0][i] = t[0][i];

		d[i] = Coords(pc)[i] - Coords(p00)[i];
		abs_v[i] = abs_v00[i] + d[i]/dt;
		anom->abs_v[i] = abs_v[i];
		va[i] = vel(i,rpst0) - abs_v[i]; 
	}

	qasqr = scalar_product(va,va,dim);    qa = sqrt(qasqr);
	for (i = 0; i < dim; i++)
	{
		t[1][i] = -va[i]/qa;
		anom->t[1][i] = t[1][i];
	}

	(void) vector_product(t[1],t[0],&sinsqr,dim);
	sinsqr = sinsqr*sinsqr;


	rhoa = Dens(rpst0);
	mf_sqr = sqr(rhoa) * qasqr * sinsqr;
	if (!s_polar_5(rpst0,YES,mf_sqr,is_plus_or,NO,
		anom->abs_v,rpst1,&ang0,&theta))
	{
	    (void) printf("WARNING: in f_mf_anom(), s_polar_5() failed\n");
	    debug_print("fanom","Left f_mfanom()\n");
	    return FUNCTION_FAILED;
	}
	for (i = 0; i < dim; i++) nor[i] = -t[0][i];
	w_speed(Coords(pvertex),anom->sot,rpst1,stdummy,rpst1,V,0.0,nor,
		FORWARD_SOUND_WAVE_TE,fr);

	if (!s_polar_3(rpst0,YES,pressure(rpst1),
			is_plus_or,NO,anom->abs_v,rpst1,
			&ang0,&theta))
	{
	    (void) printf("WARNING: in f_mf_anom(), s_polar_3() failed\n");
	    debug_print("fanom","Left f_mfanom()\n");
	    return FUNCTION_FAILED;
	}

	*M = mach_number(rpst1,anom->abs_v);

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("fanom"))
		verbose_print_f_mf_anom_data(mf_sqr,anom,bdir,d,
			qa/sound_speed(rpst0),*M);
	debug_print("fanom","Left f_mf_anom()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	return FUNCTION_SUCCEEDED;
}		/*end f_mf_anom*/

#if defined(DEBUG_NODE_PROPAGATE)
LOCAL	void verbose_print_f_mf_anom_data(
	double		mf_sqr,
	ANOM		*anom,
	BOND		*bdir,
	double		*d,
	double		Ma,
	double		Mb)
{
	WSSten		*wssten = anom->wssten;
	Front		*fr = wssten->front;
	POINT		*p00 = anom->p00;
	POINT		*pvertex = anom->pvertex;
	POINT		*pc = anom->pc;
	Locstate	rpst0 = anom->rpst0;
	Locstate	rpst1 = anom->rpst1;
	double		**t = anom->t;
	double		*abs_v = anom->abs_v;
	double		dt = wssten->dt;
	double		sin_ang, cos_ang, ang;
	double		*h = fr->rect_grid->h;
	BOND		*bc;
	int		dim = fr->interf->dim;

	(void) printf("Tangent uni_arrays\n");
	(void) printf("Incident shock tangent = <%g, %g>,\n",t[0][0],t[0][1]);
	ang = atan2(t[0][1],t[0][0]);
	(void) printf("\t\t\tmagnitude = %g, ",mag_vector(t[0],dim));
	print_angle("angle =",ang,"\n");
	(void) printf("Incident contact tangent = <%g, %g>,\n",t[1][0],t[1][1]);
	ang = atan2(t[1][1],t[1][0]);
	(void) printf("\t\t\tmagnitude = %g, ",mag_vector(t[1],dim));
	print_angle("angle =",ang,"\n");
	
	(void) vector_product(t[0],t[1],&sin_ang,dim);
	cos_ang = scalar_product(t[0],t[1],dim);
	
	ang = atan2(sin_ang,cos_ang);
	(void) printf("Incident angle (shock to contact)\n");
	(void) printf("sin(ang) = %g, cos(ang) = %g, ",sin_ang,cos_ang);
	print_angle("ang =",ang,"\n");
	(void) printf("Node velocity = <%g, %g>,\n",abs_v[0],abs_v[1]);
	ang = atan2(abs_v[1],abs_v[0]);
	(void) printf("\t\t\tmagnitude = %g, ",mag_vector(abs_v,dim));
	print_angle("angle =",ang,"\n");
	(void) printf("mf_sqr = %g\n",mf_sqr);
	(void) printf("Acoustic Impedance Squared for Ahead State = %g\n",
		acoustic_impedance_squared(rpst0));

	(void) printf("p00 = <%g, %g>\n",Coords(p00)[0],Coords(p00)[1]);
	(void) printf("pvertex = <%g, %g>, pc = <%g, %g>,\n",
		Coords(pvertex)[0],Coords(pvertex)[1],
		Coords(pc)[0],Coords(pc)[1]);
	(void) printf("\t\t\t sep %g, scaled sep %g\n",
		separation(pvertex,pc,fr->interf->dim),
		scaled_separation(pvertex,pc,h,dim));

	(void) printf("dt = %g, d = %g, %g\n",dt,d[0],d[1]);

	(void) printf("bdir - ");      print_bond(bdir);
	ang = atan2(Coords(bdir->end)[1] - Coords(bdir->start)[1],
			Coords(bdir->end)[0] - Coords(bdir->start)[0]);
	(void) printf("\tscaled len = %g, ",scaled_bond_length(bdir,h,dim));
	print_angle("ang =",ang,"\n");

	bc = anom->bc;
	(void) printf("bc - ");        print_bond(bc);
	ang = atan2(Coords(bc->end)[1] - Coords(bc->start)[1],
			Coords(bc->end)[0] - Coords(bc->start)[0]);
	(void) printf("\tscaled len = %g, ",scaled_bond_length(bc,h,dim));
	print_angle("ang =",ang,"\n");
	verbose_print_state("Ahead rpst0",rpst0);
	verbose_print_state("Behind rpst1",rpst1);

	(void) printf("Ahead  Mach Number = %g\n",Ma);
	(void) printf("Behind Mach Number = %g\n",Mb);
}		/*end verbose_print_f_mf_anom_data*/
#endif /* defined(DEBUG_NODE_PROPAGATE) */
#endif /* defined(FULL_PHYSICS) && defined(TWOD) */
