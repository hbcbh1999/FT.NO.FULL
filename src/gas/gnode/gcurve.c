/*
*
*				gcurve.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Propagation algorithm for curve hypersurface boundaries in three
*	space dimensions.
*/

#if defined(THREED)

#include <gdecs/gdecs.h>



struct _O_BOND
{
  	BOND 		*old_bond;
 	BOND_TRI 	**_btris;
  	double		*angle;
  	int		*orientation; /*TODO change to type ORIENTATION */
  	int		num_btris;
};
typedef struct _O_BOND O_BOND;

LOCAL  	O_BOND   *make_o_bond(BOND*);
LOCAL  	int      attached_b_curve_propagate(Front*,POINTER,CURVE*,CURVE*,double);
LOCAL   SURFACE  *find_physical_surface_at_curve(CURVE*, ORIENTATION*);
LOCAL   SURFACE  *adjacent_surface(SURFACE*, ORIENTATION, CURVE*, 
                                   ANGLE_DIRECTION, ORIENTATION*);
LOCAL   boolean     matchable_comps(COMPONENT , COMPONENT,INTERFACE*);
LOCAL   void     reset_fixed_btri_states(BOND*, O_BOND*, Front*);
LOCAL   int      fixed_curve_propagate(Front*,POINTER,CURVE*,CURVE*,double);
LOCAL   int      subdomain_curve_propagate(Front*,POINTER,CURVE*,CURVE*,double);
LOCAL	void	 difference3d(double*,double*,double*);
LOCAL	void	 vector_scale3d(double*,double);
LOCAL   int      fixed_curve_debug = NO;
LOCAL   void     set_curve_vel_to_zero(CURVE*);

#define surf_ang_oriented_l_to_r(ang_dir,orient)                             \
        ((((ang_dir) == CLOCKWISE && (orient) == POSITIVE_ORIENTATION) ||     \
        ((ang_dir) == COUNTER_CLOCK && (orient) == NEGATIVE_ORIENTATION)) ?  \
        YES : NO)



EXPORT  void	g_curve_propagate_3d(
	Front		*front,
	POINTER 	wave,
	CURVE		*oldc,
	CURVE		*newc,
	double		dt)
{
	int     status;
		
	DEBUG_ENTER(g_curve_propagate_3d)
        debug_print("curve_propagate","Entered g_curve_propagate_3d()\n");
	
	/* Propagate curve according to its type */

        switch (curve_type(newc))
	{
	    case PASSIVE_CURVE:   /* to do */
		 status = GOOD_CURVE;
		 break;
    	    case DIRICHLET_CURVE:
            case FIXED_CURVE:
                status = fixed_curve_propagate(front,wave,oldc,newc,dt);
                break; 
 	    case ATTACHED_B_CURVE:
	        status = attached_b_curve_propagate(front,wave,oldc,newc,dt);
	        break;
	    case NEUMANN_CURVE_W:
	    case NEUMANN_CURVE_P:
	        status = attached_b_curve_propagate(front,wave,oldc,newc,dt);
	        break;
	    case SUBDOMAIN_HSBDRY:
	        status = subdomain_curve_propagate(front,wave,oldc,newc,dt);
	        break;
	    default:
	        (void) printf("ERROR: g_curve_propagate_3d() oldc = %p ",oldc); 
	        print_hsbdry_type("Unable to process HSBDRY type ",
			      curve_type(oldc),"\n",oldc->interface);
		clean_up(ERROR);
	        break;
	}

        debug_print("curve_propagate","Left g_curve_propagate(), \n");
	DEBUG_LEAVE(g_curve_propagate_3d)
	return;
} 		/*end g_curve_propagate_3d*/

LOCAL   int     attached_b_curve_propagate(
        Front           *fr,
        POINTER         wave,
        CURVE           *oldc,
        CURVE           *newc,
        double           dt)
{
        O_SURFACE       Si, NewSi;
	O_SURFACE       Sa, NewSa;
	O_SURFACE       Sb, NewSb;
	BOND            *oldb, *newb;
        POINT           *oldp, *newp;
        Locstate        ahead, behind;
        Locstate        new_ahead, new_behind;
        double           V[MAXD];
        int             dim = fr->interf->dim;
        int             status;

        debug_print("curve_propagate","Entered attached_b_curve_propagate(), \n");
        
        oldb = oldc->first;
        newb = newc->first;

	Si.surface = find_physical_surface_at_curve(oldc,&Si.orient);
	NewSi.surface = find_physical_surface_at_curve(newc,&NewSi.orient);
	
        if (Si.surface != NULL &&  wave_type(Si.surface) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE)
        {    
             while(oldb)
             {
                 if ((oldb != oldc->last) && (!n_pt_propagated(newb->end)))
                 {
                     point_propagate(fr,wave,oldb->end,newb->end,oldb->next,
                                     oldc,dt,V);
                     n_pt_propagated(newb->end) = YES;
                     if(wave_type(oldc) == FORWARD_SHOCK_WAVE) /*???????? */
                     {
                         ahead = right_state(newb->end);
                         behind = left_state(newb->end);
                     }
                     else
                     {
                         ahead = left_state(newb->end);
                         behind = right_state(newb->end);
                     }
                 }
                 /*if (fr->bond_propagate != NULL) */
                 /*   (*fr->bond_propagate)(fr,wave,oldb,newb,oldc,dt); */
                 else
                     set_bond_length(newb,dim); /* Update new bond length */
                 if (oldb == oldc->last)
                     break;
                 oldb = oldb->next;
                 newb = newb->next;
             }
             if (Si.orient == POSITIVE_ORIENTATION)
             {
                 Sa.surface = adjacent_surface(Si.surface,Si.orient,oldc,
                                               CLOCKWISE,&Sa.orient);
                 Sb.surface = adjacent_surface(Si.surface,Si.orient,oldc,
                                               COUNTER_CLOCK,&Sb.orient);
             }  
             else
             {
                 Sa.surface = adjacent_surface(Si.surface,Si.orient,oldc,
                                           COUNTER_CLOCK,&Sa.orient);
                 Sb.surface = adjacent_surface(Si.surface,Si.orient,oldc,
                                           CLOCKWISE,&Sb.orient);
             }
             
             if (NewSi.orient == POSITIVE_ORIENTATION)
             {
                 NewSa.surface = adjacent_surface(NewSi.surface,NewSi.orient,
                                             newc, CLOCKWISE,&NewSa.orient);
                 NewSb.surface = adjacent_surface(NewSi.surface,NewSi.orient,
                                             newc, COUNTER_CLOCK,&NewSb.orient);
             }
             else
             {
                 NewSa.surface = adjacent_surface(NewSi.surface,NewSi.orient,
                                             newc, COUNTER_CLOCK,&NewSa.orient);
                 NewSb.surface = adjacent_surface(NewSi.surface,NewSi.orient,
                                             newc, CLOCKWISE,&NewSb.orient);
             }
        }
        else
        {
            status = fixed_curve_propagate(fr,wave,oldc,newc,dt);
            return status;
        }
        debug_print("curve_propagate","Left attached_b_curve_propagate(), \n");

}

/*  find the btri according to comp, only right for jet3d */
LOCAL  BOND_TRI * find_correspond_btri(
	BOND_TRI *btri, 
	BOND *bond)
{
	BOND_TRI    **newbtri;
	SURFACE     *newsurf, *surf;

	surf = btri->surface;
	for(newbtri = Btris(bond); newbtri && *newbtri; newbtri++)
	{
	    newsurf = (*newbtri)->surface;
	    if(positive_component(newsurf) == positive_component(surf) &&
	       negative_component(newsurf) == negative_component(surf))
	       return *newbtri;
	}
	return NULL;
}

/*set the intersection states between two btris at the end of a bond  */
LOCAL  void reset_fixed_btri_states(
	BOND 	*bond, 
	O_BOND  *ordered_bond, 
	Front   *fr)
{
int       	i, j, k;
SURFACE   	*surf1, *surf2;
BOND_TRI  	*bt1, *bt2;
Locstate  	state1, state2;
static Locstate tmpst = NULL;
int	  	sizest = fr->sizest;
POINT		*pos;

	if(tmpst == NULL)
	    alloc_state(fr->interf, &tmpst, sizest);

	for (i = 0; i < ordered_bond->num_btris; i++)
	    for(k=0; k<1; k++)
	    {
	        bt1 = find_correspond_btri(ordered_bond->_btris[i], bond);
            	surf1 = bt1->surface;
                
	    	j = (i == ordered_bond->num_btris - 1) ? 0 : i+1;
            	bt2 = find_correspond_btri(ordered_bond->_btris[j], bond);
            	surf2 = bt2->surface;

                if (wave_type(surf1) == SUBDOMAIN_BOUNDARY ||
                    wave_type(surf2) == SUBDOMAIN_BOUNDARY)
                    continue;
	
		if (wave_type(surf1) == PASSIVE_BOUNDARY ||
                    wave_type(surf2) == PASSIVE_BOUNDARY)
                    continue;
		
		/*only consider the end states and ASSUME one pos_comp and one neg_comp */
		/*uniquly determine one surface. */
		if(k == 0)
		{
		    if(positive_component(surf1) == positive_component(surf2))
		    {
		        state1 = right_end_btri_state(bt1);
		        state2 = right_end_btri_state(bt2);
		    }
		    else  if(positive_component(surf1) == negative_component(surf2))
		    {
		        state1 = right_end_btri_state(bt1);
		        state2 = left_end_btri_state(bt2);
		    }
		    else  if(negative_component(surf1) == negative_component(surf2))
		    {
		        state1 = left_end_btri_state(bt1);
		        state2 = left_end_btri_state(bt2);
		    }
		    else  if(negative_component(surf1) == positive_component(surf2))
		    {
		        state1 = left_end_btri_state(bt1);
		        state2 = right_end_btri_state(bt2);
		    }
		    else
		    {
		        printf("ERROR reset_fixed_btri_state, inconsistent two surfaces.\n");
			clean_up(ERROR);
		    }

		    pos = bond->end;
		}
               
		if(wave_type(surf1) == NEUMANN_BOUNDARY &&
                   wave_type(surf2) == DIRICHLET_BOUNDARY)
                {
                    ft_assign(state1,state2,sizest);
                }
                else if (wave_type(surf2) == NEUMANN_BOUNDARY &&
                     wave_type(surf1) == DIRICHLET_BOUNDARY)
                {
                    ft_assign(state2,state1,sizest);
                }
                else if ((wave_type(surf1)  >= FIRST_PHYSICS_WAVE_TYPE) &&
                    (wave_type(surf2) <  FIRST_PHYSICS_WAVE_TYPE))
                {
                    ft_assign(state2,state1,sizest);
                }
                else if ((wave_type(surf2) >= FIRST_PHYSICS_WAVE_TYPE) &&
                         (wave_type(surf1)  <  FIRST_PHYSICS_WAVE_TYPE))
                {
                    ft_assign(state1,state2,sizest);
                }
                else
                {
		    /*printf("#wtype %d  %d\n", wave_type(surf1), wave_type(surf2)); */
		    interpolate_states(fr,0.5,0.5,Coords(pos),state1,
                                   Coords(pos),state2,tmpst);
		    ft_assign(state1,tmpst,sizest);
		    ft_assign(state2,tmpst,sizest);
		}
            }
}


LOCAL   int      fixed_curve_propagate(
        Front     *fr,
        POINTER   wave,
        CURVE     *oldc,
        CURVE     *newc,
        double     dt)
{
BOND        *oldb, *newb;
POINT       *n_pos, *o_pos;
BOND_TRI    *oldbtri, *newbtri, *newbtri1;
O_BOND      *ordered_bond;
Locstate    sl, sr; 
double       V[MAXD];
size_t      sizest = fr->sizest; 
int	    i, sav_mov;

	DEBUG_ENTER(fixed_curve_propagate)

	sav_mov = fr->movingframe;
	fr->movingframe = NO;

	/*printf("#curve type %d  %d\n", curve_type(oldc), curve_type(newc)); */
	if(curve_type(newc) == NEUMANN_CURVE_P)
	{
	    /*printf("#prt curve states\n"); */
	    /*fprint_curve_states(stdout, &oldc); */
	}
	
	for(oldb = oldc->first, newb = newc->first; oldb && newb; 
	    oldb = oldb->next, newb = newb->next)
	{
            ordered_bond = make_o_bond(oldb);
	    for(i=0; i<ordered_bond->num_btris; i++)
	    {
	        oldbtri = ordered_bond->_btris[i];
		newbtri = find_correspond_btri(oldbtri, newb);
		/*this is for moving curve cases, bond_tris are not 1-1  */
		if(newbtri == NULL)
		    continue;
		if(wave_type(newbtri->surface) == PASSIVE_BOUNDARY)
		    continue;

		o_pos = oldb->end;
		n_pos = newb->end;
		normal_at_point(o_pos)[0] = HUGE_VAL;
	    	
	        ft_assign(left_state(o_pos),left_end_btri_state(oldbtri), sizest);
		ft_assign(right_state(o_pos),right_end_btri_state(oldbtri), sizest);
	
	        /*printf("#o_pos %d  %d n_pos %d  %d  %d", Params(left_state(o_pos)), Params(right_state(o_pos)), */
		/*       Params(left_state(n_pos)), Params(right_state(n_pos)), n_pos); */
		
	    	if(curve_type(newc) != NEUMANN_CURVE_P)
		{
		    point_propagate(fr,wave,o_pos,n_pos,oldbtri->tri,oldbtri->surface,dt,V);
		    Coords(n_pos)[0] = Coords(o_pos)[0];
		    Coords(n_pos)[1] = Coords(o_pos)[1];
		    Coords(n_pos)[2] = Coords(o_pos)[2];
		}
		else
		{
		    if(oldb == oldc->last && !is_closed_curve(oldc))
		    {
		        /*printf("#not closed ed\n"); */
			ft_assign(left_state(n_pos),left_state(o_pos), sizest);
		        ft_assign(right_state(n_pos),right_state(o_pos), sizest);
		    }
		    else
		    {
		        point_propagate(fr,wave,o_pos,n_pos,oldbtri->tri,oldbtri->surface,dt,V);
			/*inside the wall, should be moved ouside */
			if(gas_params_for_comp(component(Coords(n_pos), fr->interf), fr->interf) != NULL)
		            point_propagate_along_wall(fr,wave,o_pos,oldb,oldc, 
			         oldbtri->tri,oldbtri->surface,n_pos,dt,V);
		    }
		}

		sl = left_end_btri_state(newbtri);
		sr = right_end_btri_state(newbtri);
		ft_assign(sl, left_state(n_pos), sizest);
		ft_assign(sr, right_state(n_pos), sizest);
	    }
	    if(curve_type(newc) != NEUMANN_CURVE_W && 
	       curve_type(newc) != NEUMANN_CURVE_P)
	        reset_fixed_btri_states(newb, ordered_bond, fr);
	}

	/*for closed curve, the start state of the first bond comes from the end  */
	/*state of the last bond.  */
	/*for none closed curve, ASSUMING the node is always on the subdomain bdry */
	/*so it will be removed after scatter_front, just fill the state with the state */
	/*from last time step. */

        ordered_bond = make_o_bond(oldc->first);
	for(i=0; i<ordered_bond->num_btris; i++)
	{
	    oldbtri = ordered_bond->_btris[i];
	    newbtri = find_correspond_btri(oldbtri, newc->first);
	    
	    /*this is for moving curve cases, bond_tris are not 1-1  */
	    if(newbtri == NULL)
	        continue;
	    if(wave_type(newbtri->surface) == PASSIVE_BOUNDARY)
		continue;

	    newbtri1 = find_correspond_btri(oldbtri, newc->last);

	    /*newbtri: btri of the first bond. */
	    /*newbtri1: btri of the last bond. */

	    if(newbtri1 == NULL)
	    {
	        printf("ERROR: fixed_curve_propagate, interface inconsistent for btri.\n");
		clean_up(ERROR);
	    }
	    if(newbtri1->surface != newbtri->surface)
	    {
	        printf("ERROR: fixed_curve_propagate, interface inconsistent for surface.\n");
		clean_up(ERROR);
	    }

	    if(is_closed_curve(oldc))
	    {
	        sl = left_end_btri_state(newbtri1);
	        sr = right_end_btri_state(newbtri1);
	        ft_assign(left_start_btri_state(newbtri), sl, sizest);
		ft_assign(right_start_btri_state(newbtri), sr, sizest);
	    }
	    else
	    {
	        sl = left_start_btri_state(oldbtri);
	        sr = right_start_btri_state(oldbtri);
	        ft_assign(left_start_btri_state(newbtri), sl, sizest);
		ft_assign(right_start_btri_state(newbtri), sr, sizest);
	    }
	}

	DEBUG_LEAVE(fixed_curve_propagate)

	fr->movingframe = sav_mov;
        
	return GOOD_CURVE; 
}   /*end fixed_curve_propagate */


LOCAL   int      subdomain_curve_propagate(
        Front     *fr,
        POINTER   wave,
        CURVE     *oldc,
        CURVE     *newc,
        double     dt)
{
	BOND        *oldb, *newb;
	BOND_TRI    *oldbtri, *newbtri;
	size_t      sizest = fr->sizest; 
	int	    i, k;
	INTERFACE   *intfc;
	POINT	    *oldp, *newp;
	Locstate    sl, sr;
	double	    V[3], nor[4], bdry_fac = 0.5;

	DEBUG_ENTER(subdomain_curve_propagate)
	
	/*printf("#curve type subdomain %d  %d\n", curve_type(oldc), curve_type(newc)); */
	intfc = fr->interf;

	for(oldb = oldc->first, newb = newc->first; 
	    oldb && newb; 
	    oldb = oldb->next, newb = newb->next)
	{
	    if(Btris(oldb) == NULL || Btris(newb) == NULL)
	    {
	        printf("ERROR: subdomain_curve_propagate. NULL Btris\n");
		clean_up(ERROR);
	    }
	    /*subdomain curves have only one bond tri. */
	    oldbtri = Btris(oldb)[0];
	    newbtri = Btris(newb)[0];
	    if(oldbtri == NULL || newbtri == NULL)
	    {
	        printf("ERROR: subdomain_curve_propagate. NULL bond_tri\n");
		clean_up(ERROR);
	    }
	    /*subdomain curves states will be discarded after scatter_front */
	    ft_assign(left_start_btri_state(newbtri), 
	    	left_start_btri_state(oldbtri), sizest);
	    ft_assign(right_start_btri_state(newbtri), 
	    	right_start_btri_state(oldbtri), sizest);
	    ft_assign(left_end_btri_state(newbtri), 
	    	left_end_btri_state(oldbtri), sizest);
	    ft_assign(right_end_btri_state(newbtri), 
	    	right_end_btri_state(oldbtri), sizest);

	    /*moving the starting point. */
	    oldp = oldb->start;
	    newp = newb->start;
	    if(point_outside_open_bdry(&k, nor, oldp, intfc))
	    {
		point_propagate(fr,wave,oldp,newp,
		    Hyper_surf_element(oldbtri->tri), 
		    Hyper_surf(oldbtri->surface),dt,V);
		slsr(newp, Hyper_surf_element(newbtri->tri), 
			   Hyper_surf(newbtri->surface), &sl, &sr);
		ft_assign(sl,left_state(newp),sizest);
		ft_assign(sr,right_state(newp),sizest);
		
		if(V[k]*nor[k] < 0.0 || nor[3] > bdry_fac)
		    ft_assign(Coords(newp), Coords(oldp), 3*FLOAT);
	    }

	    if(newb->next != NULL)
		continue;
	    if(is_closed_curve(newc))
		continue;
	    
	    /*moving the ending point. */
	    newp = newb->end;
	    oldp = oldb->end;
	    if(point_outside_open_bdry(&k, nor, oldp, intfc))
	    {
		point_propagate(fr,wave,oldp,newp,
		    Hyper_surf_element(oldbtri->tri), 
		    Hyper_surf(oldbtri->surface),dt,V);
		slsr(newp, Hyper_surf_element(newbtri->tri), 
			   Hyper_surf(newbtri->surface), &sl, &sr);
		ft_assign(sl,left_state(newp),sizest);
		ft_assign(sr,right_state(newp),sizest);
		
		if(V[k]*nor[k] < 0.0 || nor[3] > bdry_fac)
		    ft_assign(Coords(newp), Coords(oldp), 3*FLOAT);
	    }
	}

	smooth_curve(newc);

	DEBUG_LEAVE(subdomain_curve_propagate)
		
        return GOOD_CURVE; 
}   /*end fixed_curve_propagate */

LOCAL   SURFACE  *find_physical_surface_at_curve(
	CURVE       *c,
	ORIENTATION *orient)
{
	SURFACE     **s, **surf;
	int         i;
	ORIENTATION s_or;


	for (i = 0, s = c->pos_surfaces, s_or = POSITIVE_ORIENTATION;
	     i < 2;
	     ++i,   s = c->neg_surfaces, s_or = NEGATIVE_ORIENTATION)
	{
    	    for (surf = s; surf && *surf; ++surf)
            {
                 if (wave_type(*surf) >= FIRST_PHYSICS_WAVE_TYPE)
                 {
                      *orient = s_or;
                      return *surf;
                 }
             }
 	
	}
        return NULL; 		
}

LOCAL   SURFACE  *adjacent_surface(
        SURFACE          *surf,
        ORIENTATION      s_orient,
        CURVE            *curve,
        ANGLE_DIRECTION  angle_dir,
        ORIENTATION      *adj_s_orient)
{
        INTERFACE   *intfc = surf->interface;
        BOND        *bond = curve->first;
        CURVE       *c;
        SURFACE     **s, *ans = NULL; 
        BOND_TRI    **btris;
        const TRI   *tri,*tri1,*tri2;
        double       sin12,cos12,oldsin12,oldcos12;
	const double *t1,*t2;
	int         i, dim = curve->interface->dim;
	COMPONENT   test_comp;
        static int nfail = 0;
 
        if (debugging("adjacent_surface"))
        {
            print_surface(surf);
            print_orientation("s_orient = ",s_orient,", ");
            print_angle_direction("angle_dir = ",angle_dir,"\n");
        }
        printf("Entered adjacent_surface()\n");
        test_comp = (surf_ang_oriented_l_to_r(angle_dir,s_orient)) ?
                                positive_component(surf) :
                                negative_component(surf);

	for(btris =  Btris(bond); btris && *btris;++btris)
	{
	   tri = (*btris)->tri;
	   if((*btris)->surface == surf)
	   {
               tri = (*btris)->tri;
	       t1 = Tri_normal(tri);	   
	   }
	}

        for(s == curve->pos_surfaces; s && *s; ++s)
	{
            if (*s == surf && s_orient == NEGATIVE_ORIENTATION)
                continue;

                /* Test for consistent component */
                                                                                
            if (((angle_dir == CLOCKWISE) &&
                    !matchable_comps(positive_component(*s),test_comp,intfc))
                                 ||
                ((angle_dir == COUNTER_CLOCK) &&
                    !matchable_comps(negative_component(*s),test_comp,intfc)))
                     continue;
            for(btris =  Btris(bond); btris && *btris;++btris)
            {
                if((*btris)->surface == *s)
                {
                    tri = (*btris)->tri; 
                    t2 = Tri_normal(tri);
                }
            }
            (void) vector_product(t1,t2,&sin12,dim);
            cos12 = scalar_product(t1,t2,dim);
            if (ans == NULL)
            {
                oldsin12 = sin12;
                oldcos12 = cos12;
                ans = *s;
                *adj_s_orient = NEGATIVE_ORIENTATION;
                printf("in adjacent_surface(), test1\n" );
                print_surface(ans);
                continue;
            }
            if (is_new_angle_smaller(sin12,cos12,oldsin12,oldcos12,angle_dir))
            {
                oldsin12 = sin12;
                oldcos12 = cos12;
                ans = *s;
                *adj_s_orient = NEGATIVE_ORIENTATION;
                printf("in adjacent_surface(), test2\n" );
                print_surface(ans);
            }
        }
        
        for(s == curve->neg_surfaces; s && *s; ++s)
        {
            if (*s == surf && s_orient == POSITIVE_ORIENTATION)
                continue;
                                                                                                                            
                /* Test for consistent component */
                                                                                                                             
            if (((angle_dir == CLOCKWISE) &&
                    !matchable_comps(negative_component(*s),test_comp,intfc))
                                 ||
                ((angle_dir == COUNTER_CLOCK) &&
                    !matchable_comps(positive_component(*s),test_comp,intfc)))
                     continue;
            for(btris =  Btris(bond); btris && *btris;++btris)
            {
                if((*btris)->surface == *s)
                {
                    tri = (*btris)->tri;
                    t2 = Tri_normal(tri);
                }
            }
            (void) vector_product(t1,t2,&sin12,dim);
            cos12 = scalar_product(t1,t2,dim);
            if (ans == NULL)
            {
                oldsin12 = sin12;
                oldcos12 = cos12;
                ans = *s;
                *adj_s_orient = POSITIVE_ORIENTATION;
                continue;
                printf("in adjacent_surface(), test3\n" );
                print_surface(ans);

            }
            if (is_new_angle_smaller(sin12,cos12,oldsin12,oldcos12,angle_dir))
            {
                oldsin12 = sin12;
                oldcos12 = cos12;
                ans = *s;
                *adj_s_orient = POSITIVE_ORIENTATION;
                printf("in adjacent_surface(), test4\n" );
                print_surface(ans);
            }
        }
        
        if (ans == NULL)
        {
            if (nfail++ < 10) /* TOLERANCE */
                (void) printf("WARNING in adjacent_surface(), returning null\n");
            else
            {
                screen("ERROR in adjacent_surface(), "
                       "can't find adjacent surface\n");
                clean_up(ERROR);
            }
        }
        else
            nfail = 0;
        debug_print("adjacent_surface","Leaving adjacent_surface(), ans = %d\n",ans);

        return  ans;           
} 

LOCAL   boolean    matchable_comps(
        COMPONENT comp1,
        COMPONENT comp2,
        INTERFACE *intfc)
{
        if (equivalent_comps(comp1,comp2,intfc))
            return YES;
        else if (is_exterior_comp(comp1,intfc) && is_excluded_comp(comp2,intfc))            return YES;
        else if (is_exterior_comp(comp2,intfc) && is_excluded_comp(comp1,intfc))            return YES;
        else
            return NO;
}               /*end matchable_comps*/

LOCAL  	O_BOND 	*make_o_bond(
	BOND		*b)
{
	double		x_axis[3], y_axis[3], z_axis[3];
	double		vector_in_xz_plane[3], uni_array[3];
	double		*origin = Coords(b->start);
	double		magx, magy;
	double		x_proj, y_proj;
	int		i, j_start, j_end, other;
	int		num_btris = 0;
	int		not_ordered;
	TRI		*tri;
	O_BOND		*ordered_bond;
	static int  	third_point[3][3] = 
	{
	    {-1, 2, 1},
	    { 2,-1, 0},
	    { 1, 0,-1}
	};

	static int 	orient_map[3][3] =
	{
	    { 0, 1,-1},
	    {-1, 0, 1},
	    { 1,-1, 0}
	};

	DEBUG_ENTER(make_o_bond)
	while (Btris(b)[num_btris])
	    num_btris++;

	ordered_bond = (O_BOND*) Store(sizeof(O_BOND));
	ordered_bond->_btris = (BOND_TRI**) Store(num_btris*sizeof(BOND_TRI*));
	ordered_bond->angle = (double*) Store(num_btris*sizeof(double));
	ordered_bond->orientation = (int*) Store(num_btris*sizeof(int));

	for (i = 0; i < num_btris; i++)
	    ordered_bond->_btris[i] = Btris(b)[i];
	ordered_bond->num_btris = num_btris;
	ordered_bond->old_bond = b;

	/* Construct local coordinate system. *
	 * x_axis defined to be zero radians  *
	 *  and in plane of first tri,        *
	 * z_axis defined to be bond b,       *
	 * y_axis = z_axis X x_axis           */
	
	tri = Btris(b)[0]->tri;  
	j_start = j_end = 0;
	
	while (Point_of_tri(tri)[j_start] != b->start && j_start < 3)
	    j_start++;
	while (Point_of_tri(tri)[j_end] != b->end && j_start < 3)
	    j_end++;
	other = third_point[j_start][j_end];

	difference3d(Coords(b->end),origin,z_axis);
	difference3d(Coords(Point_of_tri(tri)[other]),origin,
		     vector_in_xz_plane);
	magy = vector_product(z_axis,vector_in_xz_plane,y_axis,3);
	magx = vector_product(y_axis,z_axis,x_axis,3);
	vector_scale3d(y_axis,(1.0/magy));
	vector_scale3d(x_axis,(1.0/magx));


	for (i = 0; i < num_btris; i++)
	{
	    tri = Btris(b)[i]->tri;
	    j_start = j_end = 0;
	    
	    while (Point_of_tri(tri)[j_start] != b->start && j_start < 3)
	        j_start++;
	    while (Point_of_tri(tri)[j_end] != b->end && j_start < 3)
	        j_end++;
	    ordered_bond->orientation[i] = orient_map[j_start][j_end];
	    other = third_point[j_start][j_end];
	    
	    if (ordered_bond->orientation[i] == 0)
	    {
	        screen("unable to orient tri wrt bond in make_o_bond()\n");
		clean_up(ERROR);
	    }

	    difference3d(Coords(Point_of_tri(tri)[other]),origin,uni_array);

	    x_proj = scalar_product(uni_array,x_axis,3);
	    y_proj = scalar_product(uni_array,y_axis,3);
      
/*	    if (sqr(x_proj) + sqr(y_proj) <= .000001)*/ /*TOLERANCE*/
/*	    {
	        screen("ERROR in make_o_bond(), degenerate TRI\n");
		print_tri(tri,Btris(b)[i]->surface->interface);
		clean_up(ERROR);
	    }
*/
	    ordered_bond->angle[i] = normalized_angle(atan2(y_proj,x_proj));
	}
        
	ordered_bond->angle[0] = 0.0;

	if(fixed_curve_debug)
	{
	    printf("IN make_o_bond(),before order the bond\n");
	    print_bond(b);
	    for(i = 0; i < num_btris; i++)
	    {
		printf("btris[%d] on the surface %llu,angle = %f\n",
			i,surface_number(Btris(b)[i]->surface),ordered_bond->angle[i]);    
		/*print_tri_coords(Btris(b)[i]->tri); */
	    }
	}
	not_ordered = YES;
	while (not_ordered)  /* TODO: this should be rewritten using qsort() */
        {
	    not_ordered = NO;
	    for (i = 0; i < num_btris - 1; i++)
	        if (ordered_bond->angle[i+1] < ordered_bond->angle[i])
		{
		    double ftmp;
		    int   itmp;
		    BOND_TRI *btmp;

		    ftmp = ordered_bond->angle[i];
		    btmp = ordered_bond->_btris[i];
		    itmp = ordered_bond->orientation[i];

		    ordered_bond->angle[i] = ordered_bond->angle[i+1];
		    ordered_bond->_btris[i] = ordered_bond->_btris[i+1];
		    ordered_bond->orientation[i] = 
		        ordered_bond->orientation[i+1];

		    ordered_bond->angle[i+1] = ftmp;
		    ordered_bond->_btris[i+1] = btmp;
		    ordered_bond->orientation[i+1] = itmp;

		    not_ordered = YES;
		} 
	}
        DEBUG_LEAVE(make_o_bond)
	return ordered_bond;
}
		/*end make_o_bond*/

LOCAL	void	difference3d(
	double *a,
	double *b,
	double *c)
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}		/*end difference3d*/

LOCAL	void	vector_scale3d(
	double *a,
	double s)
{
	a[0] *= s;
	a[1] *= s;
	a[2] *= s;
}		/*end vector_scale3d*/

LOCAL   void set_curve_vel_to_zero(CURVE *c)
{
        BOND         *b;
        BOND_TRI     **btris;
        int          i;
        Locstate     start_left,start_right,end_left,end_right;
        INTERFACE    *intfc = c->interface;
                                                                                           
        for(b = c->first; b; b=b->next)
        {
            for(btris = Btris(b); btris && *btris; btris++)
            {
                start_left = left_start_btri_state(*btris);
                start_right = right_start_btri_state(*btris);
                end_left = left_end_btri_state(*btris);
                end_right = right_end_btri_state(*btris);
                for(i = 0; i < 3;i++)
                {
                    Vel(start_left)[i] = 0.0;
                    Vel(start_right)[i] = 0.0;
                    Vel(end_left)[i] = 0.0;
                    Vel(end_right)[i] = 0.0;
                }
            }
        }
}

#endif /* defined(THREED) */
