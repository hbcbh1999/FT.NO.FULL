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
*                               giboone.c:
*
*       Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*       Initialized a neutrino booster.
*/

#if defined(TWOD) || defined(THREED)

#define	MAX_NUM_PMT 10
#include <ginit/ginit.h>

enum { INTACT = 0, SWEEP };
enum { START_CRACK = 2 };

struct pmt_params
{
    double   rt;     /* transverse radius of the PMT bulb */
    double   rl;     /* longitudinal radius of the PMT bulb */
    double   rb;     /* radius of the bottom part of the PMT */
    double   rs;     /* radius of the shell of the PMT, 0 if no shell */
    double   th;	    /* thickness of the metal cover of the lower part */
    double   in;     /* indentation at the beginning of cracking */
    double   h;      /* height of PMT from bottom to the center of bulb */
    double   hcyl;   /* height of the cylindrical part of the PMT */
    double   d;	    /* distance between the adjacent PMTs */
    double   surf;   /* surface tension between ambient fluid and broken PMT */
    double   hr;	    /* height of compressed fluid in cracked PMT */
    double   pr;	    /* pressure of compressed fluid in cracked PMT */
    double   cc;     /* center of the cracked PMT, assuming it's unique */
    NODE    *upnd;  /* upper node of PMT */
    NODE    *lownd; /* lower node of PMT on the wall */
    NODE    *cracktopnd; /* the node at the top of the cracked PMT */
};

LOCAL struct pmt_params pmt;

LOCAL	void	make_pmt(
	double			center,
	int			is_crack,
        Front 			*front);

LOCAL	void	round_corner(
	double			center[2],
	NODE			*ns,
	NODE			*ne,
	COMPONENT               cmpl,
	COMPONENT               cmpr,
	int                     cndtype,
	int                     wtype,
	int                     wstatus,
	double			arm,
	int			connect_type,
	Front                   *front);
	
LOCAL	double	arg(	/* argument of point (x,y) in (-PI,PI] */
	double			x,
	double			y);

EXPORT  void init_neutrino_booster_detector(
        INIT_DATA       *init,
        INIT_PHYSICS    *ip)
{
        Front           *front = ip->root->front;
        INTERFACE       *intfc = front->interf;
        COMPONENT       ecomp = exterior_component(intfc);
        Gas_param       *paramobst, *paramsa, *paramsb;
        double           *L = front->rect_grid->GL;
        double           *U = front->rect_grid->GU;
        double		coords[MAXD];
	char            s[1024];
        int             i,num_pmt;      /* number of PMT's */
        double		firstpmt;
        int		is_crack;	/* flag */
        NODE		*tlnd,*trnd,*blnd,*brnd,*nd;
	CURVE		*lcv,*cv; /* left boundary */
		
	coords[0] = L[0], coords[1] = L[1];
        nd = make_node(Point(coords));
        set_is_bdry(nd);
        node_type(nd) = FIXED_NODE;
        blnd = nd;
        coords[0] = U[0], coords[1] = L[1];
        nd = make_node(Point(coords));
        set_is_bdry(nd);
        node_type(nd) = FIXED_NODE;
        brnd = nd;
        coords[0] = L[0], coords[1] = U[1];
        nd = make_node(Point(coords));
        set_is_bdry(nd);
        node_type(nd) = FIXED_NODE;
        tlnd = nd;
        coords[0] = U[0], coords[1] = U[1];
        nd = make_node(Point(coords));
        set_is_bdry(nd);
        node_type(nd) = FIXED_NODE;
        trnd = nd;

	debug_print("boone","Entered init_neutrino_booster_detector()\n");

	screen("\nChoose the experiment setup among BooNE(B),\n");
	screen("\tSuperK(S) and other(o): ");
	Gets(s);
	if (s[0] == 'B' || s[0] == 'b') { /* BooNE with frame radius 15.5 */
	    pmt.rt = 10,   pmt.rl = 7,   pmt.rb = 4.2;
	    pmt.rs = 0,    pmt.th = 0.1, pmt.in = 5, pmt.h = 22;
	    pmt.d = 55,    pmt.surf = 0, pmt.hr = 0.5;
	    pmt.pr = 100;
	    pmt.upnd = NULL;
	}
	else if (s[0] == 'S' || s[0] == 's'){ /* SuperK */
	    pmt.rt = 25.4, pmt.rl = 19,  pmt.rb = 12.7;
	    pmt.rs = 0,    pmt.th = 0.1, pmt.in = 2, pmt.h = 42;
	    pmt.d = 75,    pmt.surf = 0, pmt.hr = 0.5;
	    pmt.pr = 100;
	    coords[0] = L[0], coords[1] = L[1] + pmt.h;
            nd = make_node(Point(coords));
            set_is_bdry(nd);
            node_type(nd) = FIXED_NODE;
	    pmt.upnd = nd;
	}
	else
	{ /* other setup */
	    screen("\nDefine the shape of the photomultiplier tube.\n");
            screen("Enter the transverse radius of the bulb of the PMT: ");
            Scanf("%f\n",&pmt.rt);
            screen("Enter the longitudinal radius of the bulb of the PMT: ");
            Scanf("%f\n",&pmt.rl);
            screen("Enter the radius of the bottom part of the PMT: ");
            Scanf("%f\n",&pmt.rb);
	    screen("Enter the radius of the shell of PMT (0 for no shell): ");
	    Scanf("%f\n",&pmt.rs);
	    screen("Enter the thickness of metal covering the lower part: ");
	    Scanf("%f\n",&pmt.th);
	    screen("Enter the indent of the interface to start cracking: ");
	    Scanf("%f\n",&pmt.in);
            screen("Enter the height of PMT from bottom to center of bulb: ");
            Scanf("%f\n",&pmt.h);
	    screen("Enter the distance between the adjacent PMTs: ");
	    Scanf("%f\n",&pmt.d);
            screen("Enter the surface tension on the cracked PMT: ");
            Scanf("%f\n",&pmt.surf);
    	    screen("Enter the height of compressed fluid in cracked PMT: ");
	    Scanf("%f\n",&pmt.hr);
	    screen("Enter the pressure of compressed fluid in cracked PMT: ");
	    Scanf("%f\n",&pmt.pr);
	    screen("Enter if there is a shield connecting the PMTs(dflt=no): ");
    	    Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y'){
		coords[0] = L[0], coords[1] = L[1] + pmt.h;
                nd = make_node(Point(coords));
                set_is_bdry(nd);
                node_type(nd) = FIXED_NODE;
                pmt.upnd = nd;
	    }
	    else
		pmt.upnd = NULL;    
	}

        if (pmt.h + pmt.rl > U[1] - L[1])
        {
            screen("ERROR: PMT height exceeds domain.\n");
            clean_up(ERROR);
        }
	if ( pmt.d <= 2 * max(pmt.rt,pmt.rs) )
	{
	    screen("ERROR: PMTs overlap.\n");
	    clean_up(ERROR);
	}
	
        pmt.hcyl = pmt.h - pmt.rl * sqrt( 1 - sqr(pmt.rb/pmt.rt) );
        pmt.lownd = blnd;

	pmt.cc = 0; /* assume the cracked PMT is at left side */
	
        /* generate the curves along the side where PMT's are installed */
	screen("\nEnter the number of PMT's on the wall at bottom: ");
        Scanf("%d\n",&num_pmt);
	if ( max(pmt.rt,pmt.rs) + (num_pmt-1)*pmt.d >= U[0] - L[0] )
	{
	    screen("ERROR: Total width of PMTs exceeds domain.\n");
	    clean_up(ERROR);
	}
	else
	    firstpmt = L[0]; 
	
        for (i=0; i<num_pmt ;i++)
        {
            screen("Is PMT %d initially cracked (no, yes or start): ",i+1);
            Gets(s);
            if (s[0] == 's' || s[0] == 'S') is_crack = START_CRACK;
	    else if (s[0] == 'y' || s[0] == 'Y' ) is_crack = YES;
	    else is_crack = NO;
	    make_pmt(firstpmt+i*pmt.d,is_crack,front);
	    /* generate left boundary */
	    if (i ==0 )
	    {
		if (is_crack)
		{
		    cv = make_curve(ecomp,COMPA,blnd,tlnd);
		}
		else{
		    cv = make_curve(ecomp,(is_crack)?COMPB:COMPOBST,
		    		blnd,pmt.cracktopnd);
		    set_is_bdry(cv);
		    wave_type(cv) = (is_crack)?NEUMANN_BOUNDARY:
		    		PASSIVE_BOUNDARY;
		    start_status(cv) = end_status(cv) = FIXED;
		    cv = make_curve(ecomp,COMPA,pmt.cracktopnd,tlnd);
		}   
		set_is_bdry(cv);
		wave_type(cv) = NEUMANN_BOUNDARY;
		start_status(cv) = end_status(cv) = FIXED;
		if (is_crack)
		    rect_boundary_type(intfc,0,0) = NEUMANN_BOUNDARY;
		else
		    rect_boundary_type(intfc,0,0) = MIXED_TYPE_BOUNDARY;
	    }
	}
	cv = make_curve(ecomp,(pmt.upnd)?COMPOBST:COMPA,brnd,pmt.lownd);
	set_is_bdry(cv);
	wave_type(cv) = (pmt.upnd)?PASSIVE_BOUNDARY:NEUMANN_BOUNDARY;
	start_status(cv) = end_status(cv) = FIXED;
        rect_boundary_type(intfc,1,0) = MIXED_TYPE_BOUNDARY;
        
	/* generate right boundary */
	if(pmt.upnd){ /* with shield */
	    coords[0] = U[0], coords[1] = L[1] + pmt.h;
            nd = make_node(Point(coords));
            set_is_bdry(nd);
            node_type(nd) = FIXED_NODE;
	    cv = make_curve(COMPOBST,COMPA,nd,pmt.upnd);
	    wave_type(cv) = NEUMANN_BOUNDARY;
	    start_status(cv) = end_status(cv) = FIXED;
	    cv = make_curve(ecomp,COMPOBST,nd,brnd);
	    set_is_bdry(cv);
	    wave_type(cv) = PASSIVE_BOUNDARY;
	    start_status(cv) = end_status(cv) = FIXED;
	    cv = make_curve(ecomp,COMPA,trnd,nd);
	}
	else
	    cv = make_curve(ecomp,COMPA,trnd,brnd);

	set_is_bdry(cv);
	wave_type(cv) = DIRICHLET_BOUNDARY;
	start_status(cv) = end_status(cv) = FIXED;
        bstate_index(cv) = prompt_for_boundary_state(
            wave_type(cv),"right",NULL,COMPA,-1,Hyper_surf(cv),init,ip);
	rect_boundary_type(intfc,0,1) = MIXED_TYPE_BOUNDARY;

        /* generate top boundary */
        cv = make_curve(ecomp,COMPA,tlnd,trnd);
        set_is_bdry(cv);
        wave_type(cv) = DIRICHLET_BOUNDARY;
        start_status(cv) = end_status(cv) = FIXED;
        bstate_index(cv) = prompt_for_boundary_state(
            wave_type(cv),"top",NULL,COMPA,-1,Hyper_surf(cv),init,ip);
        rect_boundary_type(intfc,1,1) = MIXED_TYPE_BOUNDARY;

	/* prompt for EOS parameters */
        prompt_for_eos_params(init,ip,YES,"");
        paramobst = prompt_for_eos_params(init,ip,NO,
			"\n\tfor uncracked detector wall");
        paramsa = prompt_for_eos_params(init,ip,YES,
			"\n\tfor neutrino booster medium fluid");
        paramsb = prompt_for_eos_params(init,ip,YES,
			"\n\tfor gas inside cracking detector");
	prompt_for_ambient_state(comp_type(COMPA),paramsa,
			" for neutrino booster medium fluid",front,init);
	prompt_for_ambient_state(comp_type(COMPB),paramsb,
			" for gas inside cracking detector",front,init);
        set_obstacle_comp_type(comp_type(COMPOBST),front);

	if (debugging("boone"))
        {
            printf("\nInterface in init_neutrino_booster_detector():\n");
            print_interface(intfc);
        }
        
        debug_print("boone","Left init_neutrino_booster_detector()\n");
}       /* end init_neutrino_booster_detector */


LOCAL	void	make_pmt(
	double			center,
	int			is_crack,
        Front 			*front)
{
        COMPONENT       ecomp = exterior_component(front->interf);
        COMPONENT	ambcomp = COMPA;
        COMPONENT	pmtcomp = \
		is_crack ? ( (is_crack == YES) ? COMPA : COMPB ) : COMPOBST;
        COMPONENT	obscomp = COMPOBST;
	double           *L = front->rect_grid->GL;
        double           *U = front->rect_grid->GU;
        double		left=L[0],right=U[0],bottom=L[1];
        NODE		*sup_l,*sup_r,*bulb_l,*bulb_r,*cyl_l,*cyl_r;
        ELLIPSOID       Ellip;
        double           coords[MAXD];
        NODE            *nd,*incorner;
        CURVE           *curve,*cv;
	double		arm = 1; /* roundness of the interior corners of PMT */
	double		elev = 5; /* twisted end of photocathode */
	double		theta,r;

        /* Curve of photocathode and nodes of two ends */
	if (is_crack==YES){
	    if(center>left){
		coords[0] = center - pmt.rt, coords[1] = bottom + pmt.h;
		bulb_l = make_node(Point(coords));
		node_type(bulb_l) = FIXED_NODE;
	    }
            coords[0] = center + pmt.rt, coords[1] = bottom + pmt.h;
            bulb_r = make_node(Point(coords));
            node_type(bulb_r) = FIXED_NODE;
	}
	else
	{
	    coords[0] = center, coords[1] = bottom + pmt.h + elev;
	    pmt.cracktopnd = make_node(Point(coords));
	    set_is_bdry(pmt.cracktopnd);
	    coords[0] = center + pmt.rt, coords[1] = bottom + pmt.h + elev;
	    nd = make_node(Point(coords));
            curve = make_curve(pmtcomp,ambcomp,nd,pmt.cracktopnd);
	    if (is_crack)
	    	wave_type(curve) = CONTACT;
	    else
	    	wave_type(curve) = NEUMANN_BOUNDARY;
	    start_status(curve) = end_status(curve) = 
	    		is_crack ? INCIDENT : FIXED;
      	    node_type(curve->start) = node_type(curve->end) = 
			is_crack ? NEUMANN_NODE : FIXED_NODE;
            surface_tension(curve) = pmt.surf;
            bulb_l = curve->end, bulb_r = curve->start;

	    theta = atan( (pmt.h - pmt.hcyl) / (pmt.rt - pmt.rb) );
	    r = ( center + pmt.rb - bulb_r->posn->_coords[0] ) + 
		( bulb_r->posn->_coords[1] - bottom - pmt.hcyl ) * 
		( pmt.rt - pmt.rb ) / ( pmt.h - pmt.hcyl );
	    if(center>left){
		zero_scalar(&Ellip,sizeof(ELLIPSOID));
		Ellip.cen[0] = bulb_l->posn->_coords[0] - r/2;
		Ellip.cen[1] = bulb_l->posn->_coords[1];
		Ellip.rad[0] = r/2, Ellip.rad[1] = r/2;
		Ellip.ThetaS[0] = 0, Ellip.ThetaE[0] = - theta;
		Ellip.closed = NO;
		Ellip.scale = 1;
		Ellip.nor_orient = NEGATIVE_ORIENTATION;
		Ellip.dim = 2;
		Ellip.untracked = NO;
		Ellip.compout = pmtcomp;
		Ellip.compin = ambcomp;
		if (is_crack)
		    Ellip.wv_type = CONTACT;
		else
		    Ellip.wv_type = NEUMANN_BOUNDARY;
		cv = Curve_of_hs(g_make_ellipse(&Ellip,ambcomp,pmtcomp,front));
		nd = cv->start;
		change_node_of_curve(curve,NEGATIVE_ORIENTATION,nd);	
		curve = join_curves(curve,cv,pmtcomp,ambcomp,NULL);
		delete_node(nd);

		coords[0] = bulb_l->posn->_coords[0] - 
			r/2 * (1 - cos(theta) ) * ( 2 + cos(theta) );
		coords[1] = bulb_l->posn->_coords[1] -
		       	r/2 * sin(theta) * ( 1 + cos(theta) );
		nd = make_node(Point(coords));
		node_type(nd) = is_crack ? NEUMANN_NODE : FIXED_NODE;                                                        
		delete_node(bulb_l), bulb_l = curve->end;
		cv = make_curve(pmtcomp,ambcomp,bulb_l,nd);                                                                  
		if (is_crack)
		    wave_type(cv) = CONTACT;
		else
		    wave_type(cv) = NEUMANN_BOUNDARY;
		
		curve = join_curves(curve,cv,pmtcomp,ambcomp,NULL);
		delete_node(bulb_l), bulb_l = curve->end;
	    }
	    zero_scalar(&Ellip,sizeof(ELLIPSOID));
	    Ellip.cen[0] = bulb_r->posn->_coords[0] + r/2;  
	    Ellip.cen[1] = bulb_r->posn->_coords[1];
	    Ellip.rad[0] = r/2, Ellip.rad[1] = r/2;
	    Ellip.ThetaS[0] = PI + theta, Ellip.ThetaE[0] = PI;
	    Ellip.closed = NO;
	    Ellip.scale = 1;
	    Ellip.nor_orient = NEGATIVE_ORIENTATION;
	    Ellip.dim = 2;
	    Ellip.untracked = NO;
	    Ellip.compout = pmtcomp;
	    Ellip.compin = ambcomp;
	    if (is_crack)
	    	Ellip.wv_type = CONTACT;
	    else
	    	Ellip.wv_type = NEUMANN_BOUNDARY;
	    cv = Curve_of_hs(g_make_ellipse(&Ellip,ambcomp,pmtcomp,front));
	    nd = cv->end;
	    change_node_of_curve(curve,POSITIVE_ORIENTATION,nd);
	    curve = join_curves(cv,curve,pmtcomp,ambcomp,NULL);
	    delete_node(nd);
	    
	    coords[0] = bulb_r->posn->_coords[0] + 
		    r/2 * (1 - cos(theta) ) * ( 2 + cos(theta) );
	    coords[1] = bulb_r->posn->_coords[1] -
		    r/2 * sin(theta) * ( 1 + cos(theta) );
	    nd = make_node(Point(coords));
	    node_type(nd) = is_crack ? NEUMANN_NODE : FIXED_NODE;
	    delete_node(bulb_r), bulb_r = curve->start;
	    cv = make_curve(pmtcomp,ambcomp,nd,bulb_r);
	    if (is_crack)
	    	wave_type(cv) = CONTACT;
	    else
	    	wave_type(cv) = NEUMANN_BOUNDARY;
	    
	    curve = join_curves(cv,curve,pmtcomp,ambcomp,NULL);
	    if (is_crack)
	    	wave_type(curve) = CONTACT;
	    else
	    	wave_type(curve) = NEUMANN_BOUNDARY;
	    start_status(curve) = end_status(curve) = is_crack ? 
	    			INCIDENT : FIXED;
	    surface_tension(curve) = pmt.surf;
	    delete_node(bulb_r), bulb_r = curve->start;
 	}

	/* Create the four nodes at the bottom and the support of the PMT */
	if (pmt.upnd==NULL){ /* without shield */
	    if(center>left){
	    	coords[0] = center - ( pmt.rs ? pmt.rs : (pmt.rb + pmt.th) );
		coords[1] = bottom;
 	    	nd = make_node(Point(coords));
 	    	node_type(nd) = FIXED_NODE;
	    	set_is_bdry(nd);
	    	sup_l = nd;
	    }
	    else
		sup_l = pmt.lownd;
	    coords[0] = center + ( pmt.rs ? pmt.rs : (pmt.rb + pmt.th) );
	    coords[1] = bottom;
	    nd = make_node(Point(coords));
	    node_type(nd) = FIXED_NODE;
	    set_is_bdry(nd);
	    sup_r = nd;
	}
	if( is_crack ){
	    if(center>left){
	    	coords[0] = center - pmt.rb, coords[1] = bottom;
	    	nd = make_node(Point(coords));
	    	node_type(nd) = FIXED_NODE;
	    	set_is_bdry(nd);
	    	cyl_l = nd;
	    }
	    else
		cyl_l = pmt.lownd;
	    coords[0] = center + pmt.rb, coords[1] = bottom;
	    nd = make_node(Point(coords));
	    node_type(nd) = FIXED_NODE;
	    set_is_bdry(nd);
	    cyl_r = nd;
	}

 	/* lower boundary and shield connecting previous PMT and current one */
	if (pmt.upnd){ /* with shield */
	    if(center>left){
	    	curve = make_curve(obscomp,ambcomp,bulb_l,pmt.upnd);
	    	wave_type(curve) = NEUMANN_BOUNDARY;
	    	start_status(curve) = end_status(curve) = FIXED;
	    }
	    else
		delete_node(pmt.upnd);
	    pmt.upnd = bulb_r;
	    if (is_crack){
		if(center>left){
		    curve = make_curve(ecomp,obscomp,cyl_l,pmt.lownd);
		    wave_type(curve) = PASSIVE_BOUNDARY;
		    start_status(curve) = end_status(curve) = FIXED;
		    set_is_bdry(curve);
		}
		pmt.lownd = cyl_r;
	    }
	}	
	else{ /* without shield */
            if(center>left){
		curve = make_curve(ecomp,ambcomp,sup_l,pmt.lownd);
                wave_type(curve) = NEUMANN_BOUNDARY;
            	start_status(curve) = end_status(curve) = FIXED;
            	set_is_bdry(curve);
	    }
	    pmt.lownd = sup_r;
	}
        
        /* shells of the PMT */
        if (pmt.upnd==NULL){ /* without shield */
	    if(pmt.rs){ /* with shells */
		if(center>left){
		    curve = make_curve(obscomp,ambcomp,bulb_l,sup_l);
		    wave_type(curve) = NEUMANN_BOUNDARY;
		    start_status(curve) = end_status(curve) = FIXED;
		}
		curve = make_curve(obscomp,ambcomp,sup_r,bulb_r);
		wave_type(curve) = NEUMANN_BOUNDARY;
		start_status(curve) = end_status(curve) = FIXED;
	    }
	    else{ /* without shells */
		double theta = atan( (pmt.rt - pmt.rb) / (pmt.h - pmt.hcyl) );
		
		if(center>left){
		    coords[0] = bulb_l->posn->_coords[0] - pmt.in * sin(theta);
		    coords[1] = bulb_l->posn->_coords[1] + pmt.in * cos(theta);
		    incorner = make_node(Point(coords));
		    node_type(incorner) = FIXED_NODE;
		    curve = make_curve(obscomp,ambcomp,bulb_l,incorner);
		    wave_type(curve) = NEUMANN_BOUNDARY;
		    start_status(curve) = end_status(curve) = FIXED;
		    coords[0] = bulb_l->posn->_coords[0] -
			    pmt.in * sin(theta) - pmt.th * cos(theta);
	       	    coords[1] = bulb_l->posn->_coords[1] +
			    pmt.in * cos(theta) - pmt.th * sin(theta);
		    nd = make_node(Point(coords));
		    node_type(nd) = FIXED_NODE;
		    curve = make_curve(obscomp,ambcomp,incorner,nd);
		    wave_type(curve) = NEUMANN_BOUNDARY;
		    start_status(curve) = end_status(curve) = FIXED;
		    coords[0] = center - pmt.rb - pmt.th;
		    coords[1] = bottom + pmt.hcyl - pmt.th * tan(theta/2);
		    round_corner(coords,nd,sup_l,obscomp,ambcomp,
		    	FIXED_NODE,NEUMANN_BOUNDARY,FIXED,arm,INTACT,front);
		}
		coords[0] = bulb_r->posn->_coords[0] + pmt.in * sin(theta);
		coords[1] = bulb_r->posn->_coords[1] + pmt.in * cos(theta);
		incorner = make_node(Point(coords));
		node_type(incorner) = FIXED_NODE;
		curve = make_curve(obscomp,ambcomp,incorner,bulb_r);
		wave_type(curve) = NEUMANN_BOUNDARY;
		start_status(curve) = end_status(curve) = FIXED;
		coords[0] = bulb_r->posn->_coords[0] + pmt.in * 
				sin(theta) + pmt.th * cos(theta);
		coords[1] = bulb_r->posn->_coords[1] + pmt.in * 
				cos(theta) - pmt.th * sin(theta);
		nd = make_node(Point(coords));
		node_type(nd) = FIXED_NODE;
		curve = make_curve(obscomp,ambcomp,nd,incorner);
		wave_type(curve) = NEUMANN_BOUNDARY;
		start_status(curve) = end_status(curve) = FIXED;
		coords[0] = center + pmt.rb + pmt.th;
		coords[1] = bottom + pmt.hcyl - pmt.th * tan(theta/2);
		round_corner(coords,sup_r,nd,obscomp,ambcomp,
		    FIXED_NODE,NEUMANN_BOUNDARY,FIXED,arm,INTACT,front);
	    }
	    if(is_crack)
	    {
		if(center>left)
		{
        	    curve = make_curve(ecomp,obscomp,cyl_l,sup_l);
        	    wave_type(curve) = PASSIVE_BOUNDARY;
    		    start_status(curve) = end_status(curve) = FIXED;
        	    set_is_bdry(curve);
		}
		curve = make_curve(ecomp,obscomp,sup_r,cyl_r);
        	wave_type(curve) = PASSIVE_BOUNDARY;
       		start_status(curve) = end_status(curve) = FIXED;
        	set_is_bdry(curve);
	    }
	    else{
		curve = make_curve(ecomp,obscomp,sup_r,sup_l);
		wave_type(curve) = PASSIVE_BOUNDARY;
		start_status(curve) = end_status(curve) = FIXED;
		set_is_bdry(curve);
	    }
	}
        
        /* Lower part of the cracked or cracking PMT */
	if( is_crack ){
	    curve = make_curve(ecomp,pmtcomp,cyl_r,cyl_l);
	    wave_type(curve) = NEUMANN_BOUNDARY;
	    start_status(curve) = end_status(curve) = FIXED;
	    set_is_bdry(curve);
				
	    if(center>left){
            	coords[0] = center - pmt.rb, coords[1] = bottom + pmt.hcyl;
            	round_corner(coords,cyl_l,bulb_l,obscomp,pmtcomp,
		    FIXED_NODE,NEUMANN_BOUNDARY,FIXED,arm,SWEEP,front);
	    }
	    coords[0] = center + pmt.rb, coords[1] = bottom + pmt.hcyl;
	    round_corner(coords,bulb_r,cyl_r,obscomp,pmtcomp,
		FIXED_NODE,NEUMANN_BOUNDARY,FIXED,arm,SWEEP,front);
	}
}	/* end make_pmt */

LOCAL   void    round_corner(
        double                   center[2],
	NODE			*ns,
	NODE			*ne,
	COMPONENT		cmpl,
	COMPONENT		cmpr,
	int			cndtype,
	int			wtype,
	int			wstatus,
        double                   arm,
        int                     connect_type,
	Front			*front)
{
	double 		ts,te,tcenter,halfangle;
	double 		*ps = ns->posn->_coords, *pe = ne->posn->_coords;
	NODE 		*nd, *ncs, *nce; /* corner */
	CURVE		*curve, *cvs, *cve, *cvc; /* corner */
	ELLIPSOID       Ellip;
        double 		d; 
	d = min( sqrt( sqr(center[0] - ps[0]) + sqr(center[1] - ps[1]) ) , \
		 sqrt( sqr(center[0] - pe[0]) + sqr(center[1] - pe[1]) ) );
	d = min(arm,d);
	ts = arg(ps[0] - center[0], ps[1] - center[1]);
	te = arg(pe[0] - center[0], pe[1] - center[1]);

	if (fabs(te-ts) == PI){
		if (connect_type == INTACT){ /* 1 piece */
			curve = make_curve(cmpl,cmpr,ns,ne);
			wave_type(curve) = wtype;
			start_status(curve) = end_status(curve) = wstatus;
			return;
		}
		else{ /* SWEEP: 2 segment */
			nd = make_node(Point(center));
			node_type(nd) = cndtype;
			curve = make_curve(cmpl,cmpr,ns,nd);
			wave_type(curve) = wtype;
			start_status(curve) = end_status(curve) = wstatus;
			curve = make_curve(cmpl,cmpr,nd,ne);
			wave_type(curve) = wtype;
			start_status(curve) = end_status(curve) = wstatus;
			return;
		}
	}

	if (fabs(te-ts)>PI) te+=2*PI*((te>ts)?-1:1); /* so that |te-ts| < PI */
	tcenter = (ts+te)/2;
	halfangle = fabs(te-ts)/2;

	zero_scalar(&Ellip,sizeof(ELLIPSOID));
	Ellip.cen[0] = center[0] + arm/cos(halfangle)*cos(tcenter);
        Ellip.cen[1] = center[1] + arm/cos(halfangle)*sin(tcenter);
	Ellip.rad[0] = Ellip.rad[1] = arm*tan(halfangle);
	Ellip.ThetaS[0] = ts + ( (te>ts) ? 3*PI/2 :  PI/2 );
	Ellip.ThetaE[0] = te + ( (te>ts) ?   PI/2 :3*PI/2 );
	Ellip.closed = NO;
	Ellip.scale = 1;
	Ellip.nor_orient = (te>ts)?NEGATIVE_ORIENTATION:POSITIVE_ORIENTATION;
	Ellip.dim = 2;
	Ellip.compin = (te>ts)?cmpr:cmpl;
	Ellip.compout = (te>ts)?cmpl:cmpr;
	Ellip.wv_type = wtype;

        /* Not boundary nodes on both ends */
	cvc = Curve_of_hs(g_make_ellipse(&Ellip,Ellip.compin,
			Ellip.compout,front));
	nce = cvc->end, ncs = cvc->start;	
	cvs = make_curve(cmpl,cmpr,ns,ncs);
	cve = make_curve(cmpl,cmpr,nce,ne);
	wave_type(cvs) = wave_type(cvc) = wave_type(cve) = wtype;

	if (connect_type == INTACT){ /* 1 piece */
		curve = join_curves(cvs,cvc,cmpl,cmpr,NULL);
		curve = join_curves(curve,cve,cmpl,cmpr,NULL);
		wave_type(curve) = wtype;
		start_status(curve) = end_status(curve) = wstatus;
		delete_node(ncs), delete_node(nce);
	}
	else{ /* SWEEP: 3 segments */
		start_status(cvs) = end_status(cvs) = wstatus;
                start_status(cvc) = end_status(cvc) = wstatus;
                start_status(cve) = end_status(cve) = wstatus;
		node_type(ncs) = node_type(nce) = cndtype;
	}
} /* End round_corner */

LOCAL	double	arg(
	double			x,
	double			y)
{
	return (x==0) ? ( (y>=0)?PI/2:-PI/2 ) : \
		( (x>0) ? atan(y/x) : \
		  ( (y>=0) ? atan(y/x)+PI : atan(y/x)-PI ) );
} /* End arg */

EXPORT	void	set_state_cracked_pmt_bottom(
		double		*coords,
		COMPONENT	comp,
		Locstate	state)
{
    debug_print("init_states","Entered set_state_cracked_pmt_bottom()\n");
    if (comp==COMPA && coords[1] < pmt.hr &&  /* assume L[1] = 0 */
	coords[0] > pmt.cc - pmt.rb && coords[0] < pmt.cc + pmt.rb)
    {
	Gas_param       *params = Params(Ambient(comp_type(comp)));
	Locstate        tstate;
	
	(*params->_alloc_state)(&tstate,params->sizest);
	Init_params(tstate,params);
	set_type_of_state(tstate,TGAS_STATE);
	Dens(tstate) = Dens(state); /* assume density increment negligible */
	Press(tstate) = pmt.pr;
	Vel(tstate)[0] = Vel(tstate)[1] = 0;
	set_state(state,GAS_STATE,tstate);
	free(tstate);
    }
    debug_print("init_states","Left set_state_cracked_pmt_bottom()\n");
} /* End set_state_cracked_pmt_bottom */

#endif /* defined(TWOD) */
