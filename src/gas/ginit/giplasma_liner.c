
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
*                               giplasma.c:
*
*       Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*       Initialized a neutrino booster.
*/


#include <ginit/ginit.h>

struct PLASMA_LINER
{
    int		n_jet;	 
    double	r_jet;    
    double	d_jet;
};


EXPORT  void init_plasma_jet(
        INIT_DATA       *init,
        INIT_PHYSICS    *ip)
{
        Front           *front = ip->root->front;
        INTERFACE       *intfc = front->interf;
	NODE            *nd, *n1, *n2;
	CURVE           **cjet;
        Gas_param       *paramsa, *paramsb, *paramsc;
	char		s[Gets_BUF_SIZE];
        double           *L = front->rect_grid->GL;
        double           *U = front->rect_grid->GU;
	double           coords[MAXD], dphi,theta,dtheta;
	int		i, num_jets;
	struct          PLASMA_LINER pl;
	double		diag_angle;

	/*TMP*/
	printf("Entered init_plasma_jet()\n");

	/*This should be odd number, otherwise the formular is not correct */
	screen("Enter the number of the jet: ");
	Scanf("%d\n",&pl.n_jet);

	screen("Enter the radius of the jet: ");
	Scanf("%f\n",&pl.r_jet);

        screen("Enter the distance of the jet to the center: ");
        Scanf("%f\n",&pl.d_jet);

	uni_array(&cjet,pl.n_jet,sizeof(CURVE*));

        dphi = 0.5*PI/(pl.n_jet-1);
	dtheta = asin(pl.r_jet/pl.d_jet);
	diag_angle = atan((U[1] - L[1])/(U[0] - L[0]));
	printf("diag_angle = %f\n",diag_angle);
	printf("dphi = %f  dtheta = %f\n",dphi,dtheta);
	num_jets = pl.n_jet;

	coords[0] = L[0] + pl.d_jet; 
	coords[1] = L[1];
	n1 = make_node(Point(coords));
	set_is_bdry(n1); 
	node_type(n1) = DIRICHLET_NODE;
	coords[0] = U[0]; 
	coords[1] = L[1] + pl.r_jet;
	n2 = make_node(Point(coords));
	set_is_bdry(n2);
	node_type(n2) = SUBDOMAIN_NODE;

	cjet[0] = make_curve(COMPA,COMPB,n1,n2);
	wave_type(cjet[0]) = CONTACT;
	start_status(cjet[0]) = end_status(cjet[0]) = INCIDENT;

	coords[0] = L[0] + pl.d_jet; 
	coords[1] = L[1] + pl.r_jet; 
        insert_point_in_bond(Point(coords),cjet[0]->last,cjet[0]);


        for (i = 1; i < num_jets-1; i++)
        {
	    theta = dphi*i;
	    if (theta - dtheta < diag_angle)
	    {
                coords[0] = U[0];
                coords[1] = L[1] + (U[0] - L[0] - 
				pl.r_jet/sin(theta))*tan(theta);
		n1 = make_node(Point(coords));
		node_type(n1) = DIRICHLET_NODE;
	    }
	    else
	    {
	    	coords[0] = L[0] + (U[1] - L[1] + 
				pl.r_jet/cos(theta))*tan(PI/2 - theta);
		coords[1] = U[1];
		n1 = make_node(Point(coords));
		node_type(n1) = DIRICHLET_NODE;
	    }
	    if (theta + dtheta < diag_angle)
	    {
                coords[0] = U[0];
                coords[1] = L[1] + (U[0] - L[0] + 
				pl.r_jet/sin(theta))*tan(theta);
		n2 = make_node(Point(coords));
		node_type(n2) = DIRICHLET_NODE;
	    }
	    else
	    {
	    	coords[0] = L[0] + (U[1] - L[1] - 
				pl.r_jet/cos(theta))*tan(PI/2 - theta);
		coords[1] = U[1];
		n2 = make_node(Point(coords));
		node_type(n2) = DIRICHLET_NODE;
	    }
	    cjet[i] = make_curve(COMPA,COMPB,n1,n2);
	    wave_type(cjet[i]) = CONTACT;
	    start_status(cjet[i]) = end_status(cjet[i]) = INCIDENT;

	    coords[0] = pl.d_jet*cos(theta - dtheta) + L[0];
	    coords[1] = pl.d_jet*sin(theta - dtheta) + L[1];
            insert_point_in_bond(Point(coords),cjet[i]->last,cjet[i]);
	    coords[0] = pl.d_jet*cos(theta + dtheta) + L[0];
	    coords[1] = pl.d_jet*sin(theta + dtheta) + L[1];
            insert_point_in_bond(Point(coords),cjet[i]->last,cjet[i]);

        }
	    
	coords[0] = L[0] + pl.r_jet;
	coords[1] = U[1]; 
	n1 = make_node(Point(coords));
	set_is_bdry(n1); 
	node_type(n1) = DIRICHLET_NODE;
	coords[0] = L[0]; 
	coords[1] = L[1] + pl.d_jet;
	n2 = make_node(Point(coords));
	set_is_bdry(n2);
	node_type(n2) = SUBDOMAIN_NODE;

	cjet[num_jets-1] = make_curve(COMPA,COMPB,n1,n2);
	wave_type(cjet[num_jets-1]) = CONTACT;
	start_status(cjet[num_jets-1]) = end_status(cjet[num_jets-1]) 
			= INCIDENT;
	
	coords[0] = L[0] + pl.r_jet;
	coords[1] = L[1] + pl.d_jet;
        insert_point_in_bond(Point(coords),cjet[num_jets-1]->last,
				cjet[num_jets-1]);

	/*printf("Test before xgraph_2d_intfc_within_range()\n"); */
	/*xgraph_2d_intfc_within_range("test_llwu",intfc,L,10.0); */

	/* prompt for EOS parameters */
        prompt_for_eos_params(init,ip,YES,"");
        paramsa = prompt_for_eos_params(init,ip,YES,
			"\n\tfor the target + vacuum");
        paramsb = prompt_for_eos_params(init,ip,YES,
			"\n\tfor the plasma liner");
	/*paramsc = prompt_for_eos_params(init,ip,YES, */
	/*			"\n\tfor the outside vacuum"); */

	prompt_for_ambient_state(comp_type(COMPA),paramsa,
			" for the target + vacuum",front,init);
	prompt_for_ambient_state(comp_type(COMPB),paramsb,
			" for the plasma liner",front,init);
	/*prompt_for_ambient_state(comp_type(COMPC),paramsc, */
	/*		 " for the outside vacuum",front,init); */

	/*g_init_ext_deposition(ip->root->wave,ip->root->front); */
}      


