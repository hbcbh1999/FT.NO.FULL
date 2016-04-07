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
*				gtop.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*			User Supplied Operations
*	for gas dynamical specific interface operations.
*/


#include <gdecs/gdecs.h>

/* LOCAL variables for CC_NODE */
#define MAX_N_NODES  100

/* LOCAL variables for CC_NODE */
LOCAL   int    change_to_param_index = 1;   /* WARNING: THis is a hard wired code for jet */

/* LOCAL functions for CC_NODE */
LOCAL   void    delete_attached_bubbles(Front*);
LOCAL   boolean    intfc_delete_needles_on_node(Front*); 
LOCAL   boolean    bond_is_consecutive(BOND*,BOND*); 
LOCAL   void    open_bubbles(Front*); 
LOCAL   void    merge_slip_with_other(Front*,CURVE*,CURVE*); 
LOCAL   POINT   *pt_at_ipos(int,CURVE*); 
LOCAL   void    open_bubbles_in_parallel(Front*); 
/* Changed to IMPORT
 *  * LOCAL   void    merge_slip_with_other_in_parallel(Front*,CURVE*); 
 *   */
LOCAL   void    join_2curves_del_node(NODE*,const char*); 
LOCAL   void    mark_bubble_in_parallel(Front*); 

EXPORT	void	g_reflect_point(
	POINT		*point,/* point being reflected */
	double		*p,	/* point on reflection plane */
	double		*n,	/* normal to plane */
	INTERFACE	*intfc)	/* interface being reflected */
{
	f_reflect_point(point,p,n,intfc);
#if defined(ONED)
	if (intfc->dim == 1)
	    wave_type(point) = opposite_wave_type(wave_type(point));
#endif /*defined(ONED)*/
}		/*end f_reflect_point*/


#if defined(TWOD)
EXPORT  void    g_reflect_node2d(
	NODE		*node,/* node being reflected */
	double		*p,     /* point on reflection plane */
	double		*n)     /* normal to plane */
{
	double		*pt = Coords(node->posn);
	RP_DATA		*RP;

	if ((RP = Rp_data(node)) != NULL)
	{
	    INTERFACE	*intfc = node->interface;
	    int	i;
	    double	beta;

	    beta = PI + 2.0*angle(n[0],n[1]);
	    RP->ang_dir = Opposite_ang_dir(RP->ang_dir);
	    for (i = 0; i < MAX_N_CURVES; i++)
	    {
	    	reflect_state(RP->state[i],intfc,pt,p,n);
	    	RP->ang[i] = normalized_angle(beta - RP->ang[i]);
	    	RP->theta[i] = normalized_angle(beta - RP->theta[i]);
	    }
	}
	f_reflect_node2d(node,p,n);
}		/*end g_reflect_node2d*/

EXPORT	void g_delete_small_loops(
	Front		*fr)
{
	CURVE		**c;
	INTERFACE	*intfc = fr->interf;
	NODE		*m;
	double		perim_tol;
	double		*h = computational_grid(intfc)->h;

	perim_tol = h[0] + h[1];

label_2:
	for (c = intfc->curves;  c && *c;  c++)
	{
	    if (!is_closed_curve(*c)) continue;
	    if (wave_type(*c) < FIRST_PHYSICS_WAVE_TYPE) continue;

	    if (((*c)->num_points < 5) ||
		    (curve_length(*c) < 4.0*perim_tol)) /*TOLERANCE*/
	    {
	    	if (is_shock_wave(wave_type(*c)))
	    	{
	    	    double area = area_of_closed_curve(*c);
	    	    if ((area > 0.0) && is_backward_wave(wave_type(*c)))
	    	    {
	    	    	wave_type(*c) = FORWARD_SHOCK_WAVE;
	    	    }
	    	    if ((area < 0.0) && is_forward_wave(wave_type(*c)))
	    	    {
	    	    	wave_type(*c) = BACKWARD_SHOCK_WAVE;
	    	    }
	    	    continue;
	    	}
	    	else
	    	{
	    	    m = (*c)->start;
	    	    (void) delete_curve(*c);
	    	    if (node_type(m) == CLOSED_NODE)
	    	    	(void) delete_node(m);
	    	    goto label_2;
	    	}
	    }
	}
}		/*end g_delete_small_loops*/

/**          081004 added
 ** Consider to move to gbifur/guntan2d.c if there is the need.
 **/
EXPORT	boolean g_redist_on_cc_node(
	Front		*fr)
{
        NODE      **nn; 
        int       found_cc_node = NO;
        int       num_inconsistent_nodes = 0; 
        O_NODE    *onode_list;
        CURVE     **cc; 
        INTERFACE *sav_intfc;

        debug_print("gintfc","Entered g_redist_on_cc_node()\n");

        sav_intfc = current_interface();
        set_current_interface(fr->interf);

        for (nn = fr->interf->nodes;  nn && *nn;  ++nn)
        {
            if(node_type(*nn) == CC_NODE)
            {
                found_cc_node = YES;
                break; 
            }
        }
        if(found_cc_node == NO)
        {
            set_current_interface(sav_intfc);  
            return NO;
        }

	/*
        printf("Entered g_redist_on_cc_node, found CC_NODE fr[%d], level[%d]\n",
	       fr->patch_number, fr->patch_level); 
        for (nn = fr->interf->nodes;  nn && *nn;  ++nn)
        {
            if(node_type(*nn) == CC_NODE)
            {
                printf("print cc_node[%d]\n", *nn); 
                print_node(*nn); 
                for(cc = (*nn)->in_curves; cc && *cc; cc++)
                    print_curve(*cc); 
                for(cc = (*nn)->out_curves; cc && *cc; cc++)
                    print_curve(*cc); 
                printf("end of print cc_node[%d]\n", *nn); 
            }
        }
        */

        delete_attached_bubbles(fr);
        open_bubbles(fr);

        found_cc_node = NO;
        for (nn = fr->interf->nodes;  nn && *nn;  ++nn)
        {
            if(node_type(*nn) == CC_NODE)
            {
                found_cc_node = YES;
                break; 
            }
        }
        if(found_cc_node == NO)
        {
            set_current_interface(sav_intfc);
            return NO;
        }

        /* Do a parallel open CC_NODE, which requires communication */
        mark_bubble_in_parallel(fr); 

        /* THis function is not used since we do not 
         * keep track of the CC_NODE
         * intfc_delete_needles_on_node(fr);
         */

        /** If open_bubbles_in_parallel() is called, an 
         ** inconsistency could be generated.
         **/ 
        num_inconsistent_nodes = check_comps_at_nodes(fr->interf,
                                     &onode_list);
        if(num_inconsistent_nodes != 0)
        {
            printf("ERROR g_redist_on_cc_node, node inconsistency found\n"); 
            print_interface(fr->interf); 
            geomview_intfc_plot2d("gview_plot",fr->interf,fr->rect_grid);
            clean_up(ERROR); 
        }

        debug_print("gintfc","Left g_redist_on_cc_node()\n");
        set_current_interface(sav_intfc);
        return YES; 
}

/* Aug 6: Zhiliang: for triple point */
/**  delete_attached_bubbles remove the following figures **/
/**  ------------------------------------- C1
 **              | |  C2
 **               -
 **  Both C1 and C2 are physical waves, C2 attaches C1.
 **/
LOCAL   void delete_attached_bubbles(
        Front           *fr)
{
        CURVE           **c, **c2, **inc, **outc, *newc;
        INTERFACE       *intfc = fr->interf;
        NODE            *sn, *en, **nn;
        int             n_total, n_in, n_out;
        boolean            sav_interp;

        sav_interp = interpolate_intfc_states(intfc);

redo_label:
        for (c = intfc->curves;  c && *c;  ++c)
        {
            if (wave_type(*c) < FIRST_PHYSICS_WAVE_TYPE) continue;

            if((*c)->num_points <= 3)
            {
                for(c2 = intfc->curves;  c2 && *c2;  ++c2)
                {
                    if (wave_type(*c2) < FIRST_PHYSICS_WAVE_TYPE) continue;
                    if ((*c) == (*c2)) continue;

                    /* curve c and c2 form a closed loop */
                    if( (((*c2)->start == (*c)->start &&
                          (*c2)->end == (*c)->end) ||
                         ((*c2)->start == (*c)->end &&
                          (*c2)->end == (*c)->start)) &&
                        (*c2)->num_points <= 5 )
                    {
                        en = (*c)->end;
                        sn = (*c)->start;
                        /* Move curves on c start node to end node */
                        /* This will make a closed curve for c2 */

                        /* for( inc = (*c)->start->in_curves; inc && *inc; ++inc) */
                        for( inc = (*c)->start->in_curves; inc && *inc;)
                        {
                            change_node_of_curve((*inc), NEGATIVE_ORIENTATION,
                               (*c)->end);
                        }
                        /* for( outc = (*c)->start->out_curves; outc && *outc; ++outc) */
                        for( outc = (*c)->start->out_curves; outc && *outc;)
                        {
                            change_node_of_curve((*outc), POSITIVE_ORIENTATION,
                               (*c)->end);
                        }
 
                        /* remove c from interface */
                        delete_curve(*c);

                        /* if sn node become null, remove */
                        if(sn->in_curves == NULL &&
                           sn->out_curves == NULL)
                            delete_node(sn);

                        /* detach c2 from node ((*c)->end) */
                        /**
                        printf("THe bubble curves \n");
                        print_curve(*c2);
                        **/
                        if( (*c2)->num_points <= 3)
                            delete_curve(*c2);
                        else
                        {
                            BOND   *b;
                            NODE   *newn;

                            /**
                            b = (*c2)->first->next;
                            move_closed_loop_node(*c2, b);
                            b = (*c2)->last;
                            delete_start_of_bond(b, *c2);
                            **/
                            /* c2 is assumed to be a closed curve */
                            if((*c2)->start != (*c2)->end)
                            {
                                printf("ERROR delete_attached_bubbles\n");
                                printf("curve c2 should be a closed curve now\n");
                                print_curve(*c2);
                                clean_up(ERROR);
                            }
                            b = (*c2)->first;
                            newn = make_node(b->start);
                            node_type(newn) = CLOSED_NODE;
                            change_node_of_curve((*c2), NEGATIVE_ORIENTATION,newn);
                            change_node_of_curve((*c2), POSITIVE_ORIENTATION,newn);
                            start_status(*c2) = end_status(*c2) = INCIDENT;
                            b = (*c2)->first->next;
                            move_closed_loop_node(*c2, b);
                            b = (*c2)->last;
                            delete_start_of_bond(b, *c2);
                        }

                        n_total = num_curves_at_node(en, &n_in, &n_out);
                        if (n_total == 2 && n_in == 1 && n_out == 1 &&
                            en->in_curves[0] != en->out_curves[0] &&
                            wave_type(en->in_curves[0]) ==
                            wave_type(en->out_curves[0]))
                        {
                            interpolate_intfc_states(intfc) = YES;
                            if(negative_component(en->in_curves[0]) !=
                               negative_component(en->out_curves[0]) ||
                               positive_component(en->in_curves[0]) !=
                               positive_component(en->out_curves[0]))
                            {
                                printf("ERROR, in delete_attached_bubbles\n");
                                printf("curves have different comp\n");
                                clean_up(ERROR);
                            }

                            if(NULL  == (newc = join_curves(en->in_curves[0], 
                                en->out_curves[0],
                                negative_component(en->in_curves[0]),
                                positive_component(en->in_curves[0]),NULL)))
                            {
                                printf("ERROR, in delete_attached_bubbles\n");
                                printf("join_curves failed\n");
                                clean_up(ERROR);
                            }
                            (void) delete_node(en);
                            if(newc->num_points >= 3)
                            {
                                curve_delete_very_short_bonds(newc);
                                equi_curve_redistribute(fr, newc, YES); 
                            }
                        }
                        else if (n_total == 2 && n_in == 2 && n_out == 0 &&
                            en->in_curves[0] != en->in_curves[1] &&
                            wave_type(en->in_curves[0]) ==
                            wave_type(en->in_curves[1]))
                        {
                            interpolate_intfc_states(intfc) = YES;
                            invert_curve(en->in_curves[1]); 

                            if(negative_component(en->in_curves[0]) !=
                               negative_component(en->out_curves[0]) ||
                               positive_component(en->in_curves[0]) !=
                               positive_component(en->out_curves[0]))
                            {
                                printf("ERROR, in delete_attached_bubbles\n");
                                printf("curves have different comp case 2\n");
                                clean_up(ERROR);
                            }
                            if(NULL == (newc = join_curves(en->in_curves[0], 
                                en->out_curves[0],
                                negative_component(en->in_curves[0]),
                                positive_component(en->in_curves[0]),NULL)))
                            {
                                printf("ERROR, in delete_attached_bubbles\n");
                                printf("join_curves failed\n");
                                clean_up(ERROR);
                            }
                            (void) delete_node(en);
                            if(newc->num_points >= 3)
                            {
                                equi_curve_redistribute(fr, newc, YES); 
                                curve_delete_very_short_bonds(newc);
                            }
                        }
                        else if(n_total == 2 && n_in == 0 && n_out == 2 &&
                            en->out_curves[0] != en->out_curves[1] &&
                            wave_type(en->out_curves[0]) ==
                            wave_type(en->out_curves[1]))
                        {
                            interpolate_intfc_states(intfc) = YES;
                            invert_curve(en->out_curves[1]);
                            if(negative_component(en->in_curves[0]) !=
                               negative_component(en->out_curves[0]) ||
                               positive_component(en->in_curves[0]) !=
                               positive_component(en->out_curves[0]))
                            {
                                printf("ERROR, in delete_attached_bubbles\n");
                                printf("curves have different comp case 3\n");
                                clean_up(ERROR);
                            }
                            if(NULL == (newc = join_curves(en->in_curves[0], 
                                en->out_curves[0],
                                negative_component(en->in_curves[0]),
                                positive_component(en->in_curves[0]),NULL)))
                            {
                                printf("ERROR, in delete_attached_bubbles\n");
                                printf("join_curves failed case 3\n");
                                clean_up(ERROR);
                            }
                            (void) delete_node(en);
                            if(newc->num_points >= 3)
                            {
                                curve_delete_very_short_bonds(newc);
                                equi_curve_redistribute(fr, newc, YES); 
                            }
                        }
                        else
                        {
                            printf("WARNING, in delete_attached_bubbles\n");
                            printf("Curve merging is not performed %d %d %d\n",
                               n_total, n_in, n_out);
                            for( inc = en->in_curves; inc && *inc; ++inc)
                                print_curve(*inc);
                            for( inc = en->out_curves; inc && *inc; ++inc)
                                print_curve(*inc);
                            printf("end of print curve remained\n");
                        }
                        goto redo_label;
                    }
                }
            }
        }

redo_label_node: 
        for (nn = intfc->nodes;  nn && *nn;  ++nn)
        {
            if((*nn)->in_curves == NULL &&
               (*nn)->out_curves == NULL)
            {
                delete_node(*nn); 
                goto redo_label_node;
            }
        }

        interpolate_intfc_states(intfc) = sav_interp;
}

LOCAL boolean bond_is_consecutive(
        BOND  *b1,
        BOND  *b2)
{
        
        if(b1->start == b2->end ||
           b1->end == b2->start)
            return YES; 
        return NO; 
}

/* On a CC_NODE, small angles between two bonds can be formed by the 
 * following two cases:
 1)  fold back bonds:    ----------->
                             _______ (node)
                         <---
 2)  bonds both entering (leaving) node :
                         ----------->
                         _______     (node)
                                ---->
 *
 */ 
/* Currently, only CC_NODE generates this configuration
 */
LOCAL boolean intfc_delete_needles_on_node(
        Front        *fr)
{
        INTERFACE    *intfc = fr->interf;
        NODE            **n;
        CURVE           **c;
        int             b_fnd = NO, n_fnd = NO;
        double           min_sc_sep = MIN_SC_SEP(fr->interf);
        double           minus_cos_min_ang, cos_min_ang, min_ang, sp[3], tol[3];
        RECT_GRID       *gr = fr->rect_grid;
        double           *h = gr->h, crds[MAXD];
        int             i, nb = 0, dim = gr->dim;
        int             n_redist_nd = 0, redistributed;
        BOND            *b[3], *shortb, *longb, *thirdb;
        CURVE           *cb[3], *shortcb, *longcb, *thirdcb;
        double           vec1[MAXD], vec2[MAXD]; 
        NODE            *newn, *redist_nd[MAX_N_NODES];
        size_t          sizest = fr->sizest;
        boolean            sv_interp_status; 

        debug_print("gintfc","Entered intfc_delete_needles_on_node()\n");

        min_ang = 0.1 * min(fabs(angle(h[0],h[1])),fabs(angle(h[1],h[0])));
        cos_min_ang = cos(0.1*3.1415926535897932);
        minus_cos_min_ang = -cos_min_ang;

        sv_interp_status = interpolate_intfc_states(intfc); 
        set_size_of_intfc_state(sizest);
        interpolate_intfc_states(intfc) = YES; 

        if(current_interface() != intfc)
            printf("WARNING: intfc_delete_needles_on_node, not current intfc\n"); 

redo_node:
        for(n = intfc->nodes; n && *n; n++)
        {
            redistributed = NO;
            for(i = 0; i < n_redist_nd; i++)
            {
                if((*n) == redist_nd[i])
                {
                    redistributed = YES;
                    break; 
                }
            }
            if(redistributed == YES)
                continue; 

            if(node_type(*n) == CC_NODE)
            {
                nb = 0;
                b_fnd = NO;
                for(c = (*n)->in_curves; c && *c; c++)
                {
                    b[nb] = (*c)->last;
                    cb[nb] = *c;
                    nb++;
                    if(nb > 3 || b[nb-1]->end != (*n)->posn)
                    {
                        printf("ERROR in intfc_delete_needles_on_node\n");
                        printf("nb = %d > 3 ? %d\n", nb, (nb>3) ? 1 : 0 );
                        printf("b[%d]->end = (%"FFMT", %"FFMT")\n", nb-1,
                                Coords(b[nb-1]->end)[0],
                                Coords(b[nb-1]->end)[1]);
                        printf("(*n)->posn = (%"FFMT", %"FFMT")\n",
                                Coords((*n)->posn)[0],
                                Coords((*n)->posn)[1]);
                        print_node_type("node type = ",node_type(*n),"\n",fr->interf);
                        printf("curve c = (*n)->in_curves\n");
                        print_curve(*c);
                        print_interface(intfc);
                        clean_up(ERROR);
                    }
                }
                for(c = (*n)->out_curves; c && *c; c++)
                {
                    b[nb] = (*c)->first;
                    cb[nb] = *c;
                    nb++;
                    if(nb > 3 || b[nb-1]->start != (*n)->posn)
                    {
                        printf("ERROR in intfc_delete_needles_on_node 2\n");
                        clean_up(ERROR);
                    }
                }
                if(cb[0]->num_points <= 3 &&
                   cb[1]->num_points <= 3 &&
                   cb[2]->num_points <= 3 )
                {
                    redist_nd[n_redist_nd] = *n; 
                    n_redist_nd++; 
                    continue; 
                }
           
                if(YES == bond_is_consecutive(b[0],b[1]))
                {
                    for(i = 0; i < dim; i++)
                    {
                        vec1[i] = Coords(b[0]->start)[i] - Coords(b[0]->end)[i]; 
                        vec2[i] = Coords(b[1]->end)[i] - Coords(b[1]->start)[i]; 
                    }
                    sp[0] = scalar_product(vec1, vec2, dim)/
                                   (bond_length(b[0]) * bond_length(b[1]));
                    tol[0] = bond_length(b[0]) * bond_length(b[1]) * cos_min_ang;
                }
                else
                {
                    sp[0] = scalar_product_on_bonds(b[0],b[1],dim)/
                                (bond_length(b[0]) * bond_length(b[1]));
                    tol[0] = bond_length(b[0]) * bond_length(b[1]) * minus_cos_min_ang;
                }

                if(YES == bond_is_consecutive(b[1],b[2]))
                {
                    for(i = 0; i < dim; i++)
                    {
                        vec1[i] = Coords(b[1]->start)[i] - Coords(b[1]->end)[i]; 
                        vec2[i] = Coords(b[2]->end)[i] - Coords(b[2]->start)[i]; 
                    }
                    sp[1] = scalar_product(vec1, vec2, dim)/
                                (bond_length(b[1]) * bond_length(b[2]));
                    tol[1] = bond_length(b[1]) * bond_length(b[2]) * cos_min_ang;
                }
                else
                {
                    sp[1] = scalar_product_on_bonds(b[1],b[2],dim)/
                                (bond_length(b[1]) * bond_length(b[2]));
                    tol[1] = bond_length(b[1]) * bond_length(b[2]) * minus_cos_min_ang;
                }

                if(YES == bond_is_consecutive(b[0],b[2]))
                {
                    for(i = 0; i < dim; i++)
                    {
                        vec1[i] = Coords(b[0]->start)[i] - Coords(b[0]->end)[i]; 
                        vec2[i] = Coords(b[2]->end)[i] - Coords(b[2]->start)[i]; 
                    }
                    sp[2] = scalar_product(vec1, vec2, dim)/
                                (bond_length(b[0]) * bond_length(b[2]));
                    tol[2] = bond_length(b[0]) * bond_length(b[2]) * cos_min_ang;
                }
                else
                {
                    sp[2] = scalar_product_on_bonds(b[0],b[2],dim)/
                                (bond_length(b[0]) * bond_length(b[2]));
                    tol[2] = bond_length(b[0]) * bond_length(b[2]) * minus_cos_min_ang;
                }

                if(sp[0] >= sp[1] && sp[0] >= sp[2])
                {
                    if(sp[0] > cos_min_ang)
                        b_fnd = YES;
                    
                    if(debugging("delete_needles_on_node"))
                    {
                        printf("Found on node[%d][%g %g]\n",
                           *n, Coords((*n)->posn)[0], Coords((*n)->posn)[1]);
                        printf("The min angle[%g] is between b0 -- b1 angle %g\n",
                            acos(cos_min_ang)*180/3.14159, acos(sp[0])*180/3.14159); 
                        print_bond(b[0]);
                        print_bond(b[1]);
                        print_bond(b[2]);
                        printf("\n"); 
                    }
                    if(bond_length(b[0]) < bond_length(b[1]))
                    {
                        shortb = b[0];
                        shortcb = cb[0]; 
                        longb = b[1];
                        longcb = cb[1]; 
                    }
                    else
                    {
                        shortb = b[1];
                        shortcb = cb[1]; 
                        longb = b[0];
                        longcb = cb[0]; 
                    }
                    thirdb = b[2]; 
                    thirdcb = cb[2]; 
                } 
                else if(sp[1] >= sp[0] && sp[1] >= sp[2])
                {
                    if(sp[1] > cos_min_ang)
                        b_fnd = YES;
                    if(debugging("delete_needles_on_node"))
                    {
                        printf("Found on node[%d][%g %g]\n",
                           *n, Coords((*n)->posn)[0], Coords((*n)->posn)[1]);
                        printf("The min angle[%g] is between b1 -- b2 angle %g\n",
                            acos(cos_min_ang)*180/3.14159,acos(sp[1])*180/3.14159); 
                        print_bond(b[0]);
                        print_bond(b[1]);
                        print_bond(b[2]);
                        printf("\n"); 
                    }
                    if(bond_length(b[1]) < bond_length(b[2]))
                    {
                        shortb = b[1];
                        shortcb = cb[1]; 
                        longb = b[2];
                        longcb = cb[2]; 
                    }
                    else
                    {
                        shortb = b[2];
                        shortcb = cb[2]; 
                        longb = b[1];
                        longcb = cb[1]; 
                    }
                    thirdb = b[0]; 
                    thirdcb = cb[0]; 
                }
                else if(sp[2] >= sp[0] && sp[2] >= sp[1])
                {
                    if(sp[2] > cos_min_ang)
                        b_fnd = YES;
                    if(debugging("delete_needles_on_node"))
                    {
                        printf("Found on node[%d][%g %g]\n",
                           *n, Coords((*n)->posn)[0], Coords((*n)->posn)[1]);
                        printf("The min angle[%g] is between b0 -- b2 angle %g\n",
                            acos(cos_min_ang)*180/3.14159,acos(sp[2])*180/3.14159); 
                        print_bond(b[0]);
                        print_bond(b[1]);
                        print_bond(b[2]);
                        printf("\n"); 
                    }
                    if(bond_length(b[0]) < bond_length(b[2]))
                    {
                        shortb = b[0];
                        shortcb = cb[0]; 
                        longb = b[2];
                        longcb = cb[2]; 
                    }
                    else
                    {
                        shortb = b[2];
                        shortcb = cb[2]; 
                        longb = b[0];
                        longcb = cb[0]; 
                    }
                    thirdb = b[1]; 
                    thirdcb = cb[1]; 
                }
                
                if(b_fnd == YES)
                { 
                    if(shortb->start == (*n)->posn)
                    {
                        newn = make_node(Point(Coords(shortb->end)));

                        if(debugging("delete_needles_on_node"))
                        {                        
                            /**
                            printf("intfc_delete_needles_on_node:\n");
                            printf("Replace curve start state with bond st\n");
                            printf("Curve left, right start st:"); 
                            g_verbose_print_state(left_start_state(shortcb)); 
                            g_verbose_print_state(right_start_state(shortcb)); 

                            printf("Bond end left, right start st:"); 
                            g_verbose_print_state(left_state(shortb->end)); 
                            g_verbose_print_state(right_state(shortb->end)); 
                            **/
                        }

                        if(! is_obstacle_state(left_state(shortb->end)))
                        {
                            ft_assign(left_start_state(shortcb),
                               left_state(shortb->end),sizest);
                            ft_assign(right_start_state(shortcb),
                               right_state(shortb->end),sizest);
                        }
                        delete_end_of_bond(shortb,shortcb);
                    }
                    else
                    {
                        newn = make_node(Point(Coords(shortb->start)));

                        if(debugging("delete_needles_on_node"))
                        {
                            /**
                            printf("intfc_delete_needles_on_node:\n");
                            printf("Replace curve end state with bond st\n");
                            printf("Curve left, right end st:");
                            g_verbose_print_state(left_end_state(shortcb));
                            g_verbose_print_state(right_end_state(shortcb));

                            printf("Bond start left, right start st:");
                            g_verbose_print_state(left_state(shortb->start));
                            g_verbose_print_state(right_state(shortb->start));
                            **/ 
                        }

                        if(! is_obstacle_state(left_state(shortb->start)))
                        {
                            ft_assign(left_end_state(shortcb),
                               left_state(shortb->start),sizest);
                            ft_assign(right_end_state(shortcb),
                               right_state(shortb->start),sizest);
                        }
                        delete_start_of_bond(shortb,shortcb);
                    }

                    node_type(newn) = node_type(*n);
                    redist_nd[n_redist_nd] = newn; 
                    n_redist_nd++; 

                    for( c = (*n)->in_curves; c && *c;)
                    {
                        change_node_of_curve((*c), NEGATIVE_ORIENTATION, newn);
                    }
                    for( c = (*n)->out_curves; c && *c;)
                    {
                        change_node_of_curve((*c), POSITIVE_ORIENTATION, newn);
                    }
                    if(longcb->end == newn)
                    {
                        if(debugging("delete_needles_on_node"))
                        {
                            /**
                            printf("intfc_delete_needles_on_node:\n");
                            printf("Replace curve end state with bond st\n");
                            printf("Curve left, right end st:");
                            g_verbose_print_state(left_end_state(longcb));
                            g_verbose_print_state(right_end_state(longcb));

                            printf("Bond start left, right start st:");
                            g_verbose_print_state(left_state(longb->start));
                            g_verbose_print_state(right_state(longb->start));
                            **/
                        }

                        if(! is_obstacle_state(left_state(longb->start)))
                        {
                            ft_assign(left_end_state(longcb),
                               left_state(longb->start),sizest);
                            ft_assign(right_end_state(longcb),
                               right_state(longb->start),sizest);
                        }
                        delete_start_of_bond(longcb->last,longcb);
                    }
                    else
                    {
                        if(debugging("delete_needles_on_node"))
                        {
                            /**
                            printf("intfc_delete_needles_on_node:\n");
                            printf("Replace long curve start state with bond st\n");
                            printf("Curve left, right start st:");
                            g_verbose_print_state(left_start_state(longcb));
                            g_verbose_print_state(right_start_state(longcb));

                            printf("Bond end left, right start st:");
                            g_verbose_print_state(left_state(longb->end));
                            g_verbose_print_state(right_state(longb->end));
                            **/
                        }

                        if(! is_obstacle_state(left_state(longb->end)))
                        {
                            ft_assign(left_start_state(longcb),
                               left_state(longb->end),sizest);
                            ft_assign(right_start_state(longcb),
                               right_state(longb->end),sizest);
                        }
                        delete_end_of_bond(longcb->first,longcb);
                    }

                    for(i = 0; i < dim; i++)
                        crds[i] = 0.5*(Coords(thirdb->start)[i] + Coords(thirdb->end)[i]);
                   
                    /* if (insert_point_in_bond(Point(Coords((*n)->posn)),thirdb,thirdcb) != */
                    /*     FUNCTION_SUCCEEDED) */
                    if (insert_point_in_bond(Point(crds),thirdb,thirdcb) !=
                        FUNCTION_SUCCEEDED)
                    {
                        screen("ERROR in intfc_delete_needles_on_node(), "
                               "insert_point_in_bond() failed\n");
                        clean_up(ERROR);
                    }

                    delete_node(*n);
                    goto redo_node;
                }
            }
        }

        interpolate_intfc_states(intfc) = sv_interp_status;
        debug_print("gintfc","Left intfc_delete_needles_on_node()\n");
        return YES;
}

/** WARING: hardwired code inside: curve-component **/
LOCAL   void open_bubbles(
        Front           *fr)
{
        CURVE           **c, **c2, **inc, **outc, *newc;
        INTERFACE       *intfc = fr->interf;
        NODE            *sn, *en, **nn;
        int             n_total, n_in, n_out;
        boolean            sav_interp;
        double           area = 0.0, min_area = 0.0; 
        RECT_GRID       *gr = fr->rect_grid;

        min_area = 16.0*gr->h[0]*gr->h[1];  /* 4 by 4 size */

remove_min_area: 
        for (c = intfc->curves;  c && *c;  ++c)
        {
            if (wave_type(*c) < FIRST_PHYSICS_WAVE_TYPE) continue;

                for(c2 = intfc->curves;  c2 && *c2;  ++c2)
                {
                    if (wave_type(*c2) < FIRST_PHYSICS_WAVE_TYPE) continue;
                    if ((*c) == (*c2)) continue;

                    /* curve c and c2 form a closed loop */
                    if( (((*c2)->start == (*c)->start &&
                          (*c2)->end == (*c)->end) ||
                         ((*c2)->start == (*c)->end &&
                          (*c2)->end == (*c)->start)))
                    {
                        en = (*c)->end;
                        sn = (*c)->start;
                        /***
                        area = f_area_of_loop(*c,POSITIVE_ORIENTATION,*c2); 
                        if(area == -1.0)
                        {
                            printf("WARNING: open_bubbles\n");
                            printf("area is -1.0\n"); 
                            print_curve(*c); 
                            print_curve(*c2); 
                        }
                        {
                            print_curve(*c); 
                            print_curve(*c2); 
                            printf("closed loop area %g, min_area %g\n", area, min_area); 
                            printf("area/min_area = %g\n", area/min_area); 
                        }
                        ***/

                        /* TEST : open the bubbles no matter the size.  */
                        /* Starting from step = 66750, sp.20.4.c.lf.r3 */
                        /* if(area < 2.0*min_area) */
                        {
                            if(start_status(*c) == end_status(*c) &&
                               start_status(*c) == SLIP &&
                               start_status(*c2) != SLIP &&
                               end_status(*c2) != SLIP &&
                               ((positive_component(*c2) == 2 && 
                                 negative_component(*c2) == 4) ||
                                (positive_component(*c2) == 4 &&
                                 negative_component(*c2) == 2)
                               )
                              )
                            {
                                merge_slip_with_other(fr, *c, *c2); 
                                goto remove_min_area; 
                            }
                            else if(start_status(*c2) == end_status(*c2) &&
                               start_status(*c2) == SLIP &&
                               start_status(*c) != SLIP &&
                               end_status(*c) != SLIP &&
                               ((positive_component(*c) == 2 &&
                                 negative_component(*c) == 4) ||
                                (positive_component(*c) == 4 &&
                                 negative_component(*c) == 2)
                               )
                                   )
                            {
                                merge_slip_with_other(fr, *c2, *c); 
                                goto remove_min_area; 
                            }
                            /***
                            else
                            {
                                printf("WARNING: open_bubbles\n");
                                printf(" curves have other status\n"); 
                                print_curve(*c); 
                                print_curve(*c2); 
                                clean_up(ERROR); 
                            }
                            ***/
                        }
                    }
                }
        }

}

/** Do the following: 
 ** ---------------->-------  slipc
 ** |                      |
 ** |--------------->------|   c 
 **/
LOCAL void merge_slip_with_other(
       Front   *fr,
       CURVE   *slipc,
       CURVE   *c)
{
       BOND *cb, *slipb;
       int nump_c, nump_slipc, i, j;
       POINT *cpt, *slippt;
       size_t      sizest = fr->sizest;
       boolean        sav_interp;
       NODE        *sn, *en; 
       int         n_total, n_in, n_out;  
       CURVE       *newc; 

       sav_interp = interpolate_intfc_states(fr->interf);
       interpolate_intfc_states(fr->interf) = YES;
       nump_c = c->num_points; 
       nump_slipc = slipc->num_points; 
 
       if( c->start == slipc->end &&
           c->end == slipc->start )
       {
           invert_curve(slipc); 
       }

       if(negative_component(c) == 
          positive_component(slipc))
       {
           NULL;
       }
       else if(positive_component(c) ==
          negative_component(slipc))
       {
           NULL;
       }
       else
       {
           printf("ERROR merge_slip_with_other\n");
           printf("closed loop not consistent\n");
           clean_up(ERROR); 
       }

       /* printf("merge_slip_with_other called to open bubble\n"); */
       /* printf("Show slip curve status before change\n");  */
       /* show_curve_states(slipc);  */
      
       /***
       ft_assign(left_state(slipc->first->start),
           left_start_state(slipc),sizest);
       ft_assign(right_state(slipc->first->start),
           right_start_state(slipc),sizest);
       ft_assign(left_state(slipc->last->end),
           left_end_state(slipc),sizest);
       ft_assign(right_state(slipc->last->end),
           right_end_state(slipc),sizest);
       ***/

       cb = c->first; 
       i = 2; 
       while(cb)
       {
           cpt = cb->end; 
           j = (int)(1.0*i/nump_c*nump_slipc);  
           slippt = pt_at_ipos(j, slipc);

           if(slippt == NULL)
           {
               printf("ERROR merge_slip_with_other\n"); 
               printf("NULL point found\n");
               clean_up(ERROR); 
           }

           if(negative_component(c) == 
              positive_component(slipc))
           {
               if(nump_slipc == 2)
               {
                   /**
                   if(slippt == slipc->start->posn)
                       ft_assign(left_state(cpt),
                          left_start_state(slipc),sizest);
                   else
                       ft_assign(left_state(cpt),
                          left_end_state(slipc),sizest);
                   **/
                   /* Do not modify states */
                   Set_params(left_state(cpt), left_start_state(slipc)); 
                   /* Hard wired code */
                   if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                   {
                       if(Params(left_start_state(slipc))->n_comps != 1)
                       {
                           pdens(left_state(cpt))[0] = 0.0;
                           pdens(left_state(cpt))[1] = Dens(left_state(cpt));
                       }
                   }
               }
               else
               {
                   /**
                   ft_assign(left_state(cpt),
                      left_state(slippt),sizest);
                   **/
                   /* Do not modify states */
                   if(slippt == slipc->start->posn)
                   {
                       Set_params(left_state(cpt), left_start_state(slipc)); 
                       /* Hard wired code */
                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       {
                           if(Params(left_start_state(slipc))->n_comps != 1)
                           {
                               pdens(left_state(cpt))[0] = 0.0;
                               pdens(left_state(cpt))[1] = Dens(left_state(cpt));
                           }
                       }
                   }
                   else if(slippt == slipc->end->posn)
                   {
                       Set_params(left_state(cpt), left_end_state(slipc)); 
                       /* Hard wired code */
                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       {
                           if(Params(left_end_state(slipc))->n_comps != 1)
                           {
                               pdens(left_state(cpt))[0] = 0.0;
                               pdens(left_state(cpt))[1] = Dens(left_state(cpt));
                           }
                       }
                   }
                   else
                   {
                       Set_params(left_state(cpt), left_state(slippt)); 
                       /* Hard wired code */
                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       {
                           if(Params(left_state(slippt))->n_comps != 1)
                           {
                               pdens(left_state(cpt))[0] = 0.0;
                               pdens(left_state(cpt))[1] = Dens(left_state(cpt));
                           }
                       }
                   }
               }
           }
           else
           {
               if(nump_slipc == 2)
               {
                   /**
                   if(slippt == slipc->start->posn)
                       ft_assign(right_state(cpt),
                          right_start_state(slipc),sizest);
                   else
                       ft_assign(right_state(cpt),
                          right_end_state(slipc),sizest);
                   **/
                   /* Do not modify states */
                   Set_params(right_state(cpt), right_end_state(slipc)); 
                   /* Hard wired code */
                   if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                   {
                       if(Params(right_end_state(slipc))->n_comps != 1)
                       {
                           pdens(right_state(cpt))[0] = 0.0;
                           pdens(right_state(cpt))[1] = Dens(right_state(cpt));
                       }
                   }
               }
               else
               {
                   /**
                   ft_assign(right_state(cpt),
                      right_state(slippt),sizest);
                   **/
                   /* Do not modify states */
                   if(slippt == slipc->start->posn)
                   {
                       Set_params(right_state(cpt), right_start_state(slipc)); 
                       /* Hard wired code */
                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       {
                           if(Params(right_start_state(slipc))->n_comps != 1)
                           {
                               pdens(right_state(cpt))[0] = 0.0;
                               pdens(right_state(cpt))[1] = Dens(right_state(cpt));
                           }
                       }
                   }
                   else if(slippt == slipc->end->posn)
                   {
                       Set_params(right_state(cpt), right_end_state(slipc)); 
                       /* Hard wired code */
                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       {
                           if(Params(right_end_state(slipc))->n_comps != 1)
                           {
                               pdens(right_state(cpt))[0] = 0.0;
                               pdens(right_state(cpt))[1] = Dens(right_state(cpt));
                           }
                       }
                   }
                   else
                   {
                       Set_params(right_state(cpt), right_state(slippt)); 
                       /* Hard wired code */
                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       {
                           if(Params(right_state(slippt))->n_comps != 1)
                           {
                               pdens(right_state(cpt))[0] = 0.0;
                               pdens(right_state(cpt))[1] = Dens(right_state(cpt));
                           }
                       }
                   }
               }
           }

           cb = cb->next; 
           i++; 
       } 

       if(negative_component(c) == 
          positive_component(slipc))
       {
            ft_assign(left_start_state(c),
               left_start_state(slipc),sizest);
            ft_assign(left_end_state(c),
               left_end_state(slipc),sizest);
            ft_assign(left_state(c->first->start),
               left_state(slipc->first->start),sizest);
            negative_component(c) = negative_component(slipc); 
       }
       else
       {
            ft_assign(right_start_state(c),
               right_start_state(slipc),sizest);
            ft_assign(right_end_state(c),
               right_end_state(slipc),sizest);
            ft_assign(right_state(c->first->start),
               right_state(slipc->first->start),sizest);
            positive_component(c) = positive_component(slipc); 
       }

       /* TMP */
       /* printf("Show curve in merge_slip_with_other()\n");  */
       /* print_curve(c); */
       /* printf("Show curve status after change\n");  */
       /* show_curve_states(c);  */

       delete_curve(slipc); 

       en = c->end;
       sn = c->start;

       n_total = num_curves_at_node(en, &n_in, &n_out);

       if (n_total == 2 && n_in == 2 && n_out == 0 &&
           en->in_curves[0] != en->in_curves[1])
       {
           if(wave_type(en->in_curves[0]) !=
              wave_type(en->in_curves[1]))
           {
               printf("ERROR merge_slip_with_other\n");
               printf("Different curves 1\n"); 
               clean_up(ERROR); 
           }
           if(c == en->in_curves[1])
               invert_curve(en->in_curves[1]);
           else
               invert_curve(en->in_curves[0]);
       }
       else if (n_total == 2 && n_in == 0 && n_out == 2 &&
           en->out_curves[0] != en->out_curves[1])
       {
           if(wave_type(en->out_curves[0]) !=
              wave_type(en->out_curves[1]))
           {
               printf("ERROR merge_slip_with_other\n");
               printf("Different curves 2\n"); 
               clean_up(ERROR); 
           }
           if(c == en->out_curves[1])
               invert_curve(en->out_curves[1]);
           else
               invert_curve(en->out_curves[0]);
       }

       n_total = num_curves_at_node(en, &n_in, &n_out);

       join_2curves_del_node(en,"en_merge_slip_with_other");
       /**
       if (n_total == 2 && n_in == 1 && n_out == 1 &&
            en->in_curves[0] != en->out_curves[0] &&
            wave_type(en->in_curves[0]) ==
            wave_type(en->out_curves[0]))
       {
            if(negative_component(en->in_curves[0]) !=
               negative_component(en->out_curves[0]) ||
               positive_component(en->in_curves[0]) !=
               positive_component(en->out_curves[0]))
            {
                printf("ERROR, in merge_slip_with_other\n");
                printf("curves have different comp 1\n");
                clean_up(ERROR);
            }
                                                                                                      
            if(NULL  == (newc = join_curves(en->in_curves[0],
               en->out_curves[0],
               negative_component(en->in_curves[0]),
               positive_component(en->in_curves[0]),NULL)))
           {
               printf("ERROR, in merge_slip_with_other\n");
               printf("join_curves failed 1\n");
               clean_up(ERROR);
           }
           (void) delete_node(en);
       }
       else
       {
           printf("ERROR merge_slip_with_other\n");
               printf("Different curves 3\n");
               clean_up(ERROR);
       }
       **/

       /*** On the start node side ***/
       n_total = num_curves_at_node(sn, &n_in, &n_out);

       if (n_total == 2 && n_in == 2 && n_out == 0 &&
           sn->in_curves[0] != sn->in_curves[1])
       {
           if(wave_type(sn->in_curves[0]) !=
              wave_type(sn->in_curves[1]))
           {
               printf("ERROR merge_slip_with_other\n");
               printf("Different curves s1\n");
               clean_up(ERROR);
           }
           invert_curve(sn->in_curves[1]);
       }
       else if (n_total == 2 && n_in == 0 && n_out == 2 &&
           sn->out_curves[0] != sn->out_curves[1])
       {
           if(wave_type(sn->out_curves[0]) !=
              wave_type(sn->out_curves[1]))
           {
               printf("ERROR merge_slip_with_other\n");
               printf("Different curves s2\n");
               clean_up(ERROR);
           }
           invert_curve(sn->out_curves[1]);
       }

       n_total = num_curves_at_node(sn, &n_in, &n_out);

       join_2curves_del_node(sn,"sn_merge_slip_with_other");
       /**
       if (n_total == 2 && n_in == 1 && n_out == 1 &&
            sn->in_curves[0] != sn->out_curves[0] &&
            wave_type(sn->in_curves[0]) ==
            wave_type(sn->out_curves[0]))
       {
            if(negative_component(sn->in_curves[0]) !=
               negative_component(sn->out_curves[0]) ||
               positive_component(sn->in_curves[0]) !=
               positive_component(sn->out_curves[0]))
            {
                printf("ERROR, in merge_slip_with_other\n");
                printf("curves have different comp s1\n");
                clean_up(ERROR);
            }
                                                                                                      
            if(NULL  == (newc = join_curves(sn->in_curves[0],
               sn->out_curves[0],
               negative_component(sn->in_curves[0]),
               positive_component(sn->in_curves[0]),NULL)))
           {
               printf("ERROR, in merge_slip_with_other\n");
               printf("join_curves failed s1\n");
               clean_up(ERROR);
           }
           (void) delete_node(sn);
       }
       else
       {
           printf("ERROR merge_slip_with_other\n");
               printf("Different curves s3\n");
               clean_up(ERROR);
       }
       **/

       interpolate_intfc_states(fr->interf) = sav_interp;
}

LOCAL POINT *pt_at_ipos(
       int ipos,
       CURVE *c)
{
       int nump_c, i;
       POINT *cpt;
       BOND  *b;

       nump_c = c->num_points; 

       if(ipos == 1)
           return c->first->end;
       if(ipos == 0)
           return c->first->start;
 
       b = c->first;
       i = 2;
       while(b)
       {
           if(i == ipos)
               return b->end;
           i++; 
           b = b->next; 
       }

       return c->end->posn; 
}

LOCAL   void open_bubbles_in_parallel(
        Front           *fr)
{
        CURVE       **c, **c2, *newc;
        INTERFACE   *intfc = fr->interf;
        boolean        sav_interp;
        NODE        *sn, *en, **nn;
        int         n_total, n_in, n_out;

        sav_interp = interpolate_intfc_states(fr->interf);
        interpolate_intfc_states(fr->interf) = YES;

change_curve_type: 
        for (c = intfc->curves;  c && *c;  ++c)
        {
            if (wave_type(*c) < FIRST_USER_BOUNDARY_TYPE) continue;
            if(change_curve_comp(*c) == YES)
                merge_slip_with_other_in_parallel(fr, *c); 
        }

        /* NEED to delete SLIP curves here */
        for (c = intfc->curves;  c && *c;  ++c)
        {
            if (wave_type(*c) < FIRST_PHYSICS_WAVE_TYPE) continue;
            if(start_status(*c) == end_status(*c) &&
               start_status(*c) == SLIP)
            {
                sn = (*c)->start; 
                en = (*c)->end; 
               
                /* TMP */
                /**
                printf("\nWARNING: in open_bubbles_in_parallel\n");
                printf("SLIP curve remains\n"); 
                print_curve(*c); 
                **/

                delete_curve(*c);

                n_total = num_curves_at_node(en, &n_in, &n_out);
             
                /**
                print_node(en); 
                printf("n_total of end %d, n_in %d, n_out %d\n", 
                     n_total, n_in, n_out); 
                **/

                /* DO merge curve case by case  */
                join_2curves_del_node(en, "open_bubbles_in_parallel_0"); 

                n_total = num_curves_at_node(sn, &n_in, &n_out);

                /**
                print_node(sn); 
                printf("n_total of start %d, n_in %d, n_out %d\n", 
                     n_total, n_in, n_out); 
                **/

                join_2curves_del_node(sn, "open_bubbles_in_parallel_1"); 
            }
        }

        interpolate_intfc_states(fr->interf) = sav_interp;
 
        /* TMP */
        /* printf("print interface in open_bubbles_in_parallel\n"); */
        /* print_interface(intfc); */
 
        /**
        printf("ERROR: clean_up in open_bubbles_in_parallel\n"); 
        clean_up(ERROR);
        **/
}

/* WARNING: change_to_param_index is hardwired */
EXPORT void merge_slip_with_other_in_parallel(
       Front   *fr,
       CURVE   *c)
{
       BOND        *cb;
       POINT       *cpt;
       Gas_param   *change_to_param;

       /* printf("merge_slip_with_other_in_parallel called to open bubble in parallel\n"); */
       /* printf("Show slip curve status before change\n");  */
      
       change_to_param = Params_of_index(change_to_param_index);  

       cb = c->first; 
       while(cb)
       {
           cpt = cb->end; 
           if(change_curve_comp_on_side(c) == NEGATIVE_SIDE)
           {
               /* Set_params(left_state(cpt), left_start_state(slipc));  */
               Params(left_state(cpt)) = change_to_param; 
               /* Hard wired code */
               if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
               {
                   if((change_to_param)->n_comps != 1)
                   {
                       pdens(left_state(cpt))[0] = 0.0;
                       pdens(left_state(cpt))[1] = Dens(left_state(cpt));
                   }
               }
           }
           else
           {
               /* Set_params(right_state(cpt), right_start_state(slipc));  */
               Params(right_state(cpt)) = change_to_param; 
               /* Hard wired code */
               if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
               {
                   if((change_to_param)->n_comps != 1)
                   {
                       pdens(right_state(cpt))[0] = 0.0;
                       pdens(right_state(cpt))[1] = Dens(right_state(cpt));
                   }
               }
           }
           cb = cb->next; 
       } 

       if(change_curve_comp_on_side(c) == NEGATIVE_SIDE)
       {
            /**
            Set_params(left_start_state(c), left_start_state(slipc)); 
            Set_params(left_end_state(c), left_end_state(slipc)); 
            negative_component(c) = negative_component(slipc); 
            **/
            Params(left_start_state(c)) = change_to_param; 
            Params(left_end_state(c)) = change_to_param; 
               /* Hard wired code */
               if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
               {
                   if((change_to_param)->n_comps != 1)
                   {
                       pdens(left_start_state(c))[0] = 0.0;
                       pdens(left_start_state(c))[1] = Dens(left_start_state(c));
                       pdens(left_end_state(c))[0] = 0.0;
                       pdens(left_end_state(c))[1] = Dens(left_end_state(c));
                   }
               }
            negative_component(c) = change_curve_comp_by(c); 
       }
       else
       {
            /**
            Set_params(right_start_state(c), right_start_state(slipc)); 
            Set_params(right_end_state(c), right_end_state(slipc)); 
            positive_component(c) = positive_component(slipc); 
            **/
            Params(right_start_state(c)) = change_to_param; 
            Params(right_end_state(c)) =  change_to_param; 
               /* Hard wired code */
               if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
               {
                   if((change_to_param)->n_comps != 1)
                   {
                       pdens(right_start_state(c))[0] = 0.0;
                       pdens(right_start_state(c))[1] = Dens(right_start_state(c));
                       pdens(right_end_state(c))[0] = 0.0;
                       pdens(right_end_state(c))[1] = Dens(right_end_state(c));
                   }
               }
            positive_component(c) = change_curve_comp_by(c); 
       }
}

LOCAL void join_2curves_del_node(
       NODE        *n,
       const char  *msg)
{
       int   n_total, n_in, n_out;
       CURVE  *newc;

       n_total = num_curves_at_node(n, &n_in, &n_out);

       /* DO merge curve case by case  */
       if (n_total == 2 && n_in == 1 && n_out == 1 &&
           n->in_curves[0] != n->out_curves[0] &&
           wave_type(n->in_curves[0]) ==
           wave_type(n->out_curves[0]))
       {
           if(negative_component(n->in_curves[0]) !=
              negative_component(n->out_curves[0]) ||
              positive_component(n->in_curves[0]) !=
              positive_component(n->out_curves[0]))
           {
               printf("ERROR, in join_2curves_del_node called by %s\n", msg);
               printf("curves have different comp 1\n");
               clean_up(ERROR);
           }

           if(NULL  == (newc = join_curves(n->in_curves[0],
              n->out_curves[0],
              negative_component(n->in_curves[0]),
              positive_component(n->in_curves[0]),NULL)))
           {
               printf("ERROR, in open_bubbles_in_parallel %s\n", msg);
               printf("join_curves failed 1\n");
               print_curve(n->in_curves[0]); 
               print_curve(n->out_curves[0]); 
               clean_up(ERROR);
           }
           (void) delete_node(n);
       }
       else if(n_total == 1 && n_in == 1 && n_out == 1 &&
               n->in_curves[0] == n->out_curves[0])
       {
           /* closed curve detected, do nothing */
           node_type(n) = CLOSED_NODE;
           start_status(n->in_curves[0]) = 
           end_status(n->in_curves[0]) = INCIDENT;
       }
       else
       {
           CURVE **cc;
           printf("ERROR join_2curves_del_node called by %s\n", msg);
           printf("Different curve case\n");
           printf("n_total %d, n_in %d, n_out %d\n", n_total, n_in, n_out);
           printf("print in curves\n");
           for(cc = n->in_curves; cc && *cc; cc++)
               print_curve(*cc);
           printf("print out curves\n");
           for(cc = n->out_curves; cc && *cc; cc++)
               print_curve(*cc);
           clean_up(ERROR);
       }
}

LOCAL   void mark_bubble_in_parallel(
        Front           *fr)
{
        CURVE       **c, **c2;
        INTERFACE   *intfc = fr->interf;
        CURVE       *changed_c[100];
        int         change_n = 0;  

        for (c = intfc->curves;  c && *c;  ++c)
        {
            if (wave_type(*c) < FIRST_PHYSICS_WAVE_TYPE) continue;
            for(c2 = intfc->curves;  c2 && *c2;  ++c2)
            {
                if (wave_type(*c2) < FIRST_PHYSICS_WAVE_TYPE) continue;
                if ((*c) == (*c2)) continue;

                /* find curves need to change type, */
                /* hardwired code */
                if(start_status(*c) == end_status(*c) &&
                   start_status(*c) == SLIP)
                {
                    if((*c2)->start == (*c)->start ||
                       (*c2)->end == (*c)->end ||
                       (*c2)->start == (*c)->end ||
                       (*c2)->end == (*c)->start)
                    {
                        if((negative_component(*c2) == 4 &&
                            positive_component(*c2) == 2) ||
                           (negative_component(*c2) == 2 &&
                            positive_component(*c2) == 4))
                        {
                            /***
                            printf("WARNING: in mark_bubble_in_parallel\n"); 
                            printf("WARNING: change curve %llu comp\n", curve_number(*c2)); 
                            print_curve(*c); 
                            print_curve(*c2); 
                            ***/

                            change_curve_comp(*c2) = YES;
                            if(negative_component(*c2) == 4) 
                                change_curve_comp_on_side(*c2) = NEGATIVE_SIDE;
                            else
                                change_curve_comp_on_side(*c2) = POSITIVE_SIDE; 

                            if(negative_component(*c2) ==
                               positive_component(*c))
                                change_curve_comp_by(*c2) = negative_component(*c); 
                            else if(positive_component(*c2) ==
                                    negative_component(*c))
                                change_curve_comp_by(*c2) = positive_component(*c); 
                            else if(positive_component(*c2) ==
                                    positive_component(*c)) 
                                change_curve_comp_by(*c2) = negative_component(*c); 
                            else if(negative_component(*c2) ==
                                    negative_component(*c))
                                change_curve_comp_by(*c2) = positive_component(*c); 
                            
                            changed_c[change_n] = *c2;
                            change_n++; 
                        } 
                    }
                }
            }
        }
}

/**          11082004 added
 ** Consider to move to gbifur/guntan2d.c if there is the need.
 **/
EXPORT	boolean g_parallel_redist_on_cc_node(
	Front		*fr)
{
        INTERFACE  *sav_intfc;

        sav_intfc = current_interface();
        set_current_interface(fr->interf); 

        open_bubbles_in_parallel(fr); 

        set_current_interface(sav_intfc); 
        return YES; 
}
#endif /* defined(TWOD) */
