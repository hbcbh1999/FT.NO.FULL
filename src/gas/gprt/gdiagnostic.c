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
*				gdiagnositic.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains printing routines for gas dynamics.
*
*/


#include <gdecs/gdecs.h>

EXPORT 	void 	summary_of_curve_states3d(
	CURVE 		*c)
{
  	BOND   	  *b;
	int    	  b_cnt;
	BOND_TRI  **btris, *bt;
	INTERFACE *intfc = c->interface;
	Locstate  ls, rs;

	(void) printf("\n");
	(void) printf(" start of summary_of_curve_states3d() curve = %p\n",c);
	(void) printf("\n");

	for (b = c->first,b_cnt = 0; b != c->last; b = b->next,b_cnt++)
	{
	    if (21*(b_cnt/21) == b_cnt)
	    {
	        (void) printf("\n");
		(void) printf("%4s %9s %2s | %23s | %7s %7s %7s %7s %7s | "
			      "%7s %7s %7s %7s %7s\n",
			      "cnt","bond ","s","coords end      ",
			      "dens","engy","Mx ","My ","Mz ",
			      "dens","engy","Mx ","My ","Mz ");
		(void) printf("%4s %9s %2s | %23s | %7s %7s %7s %7s %7s | "
			      "%7s %7s %7s %7s %7s\n",
			      "====","=========","==",
			      "=======================","=======","=======",
			      "=======","=======","=======","=======",
			      "=======","=======","=======","=======");
	    }

	    for (btris = Btris(b); btris && *btris; btris++)
	    {
	        bt = *btris;
		ls = left_end_btri_state(bt);
		rs = right_end_btri_state(bt);
		(void) printf("%4d %llu %2d | %g %g %g |",b_cnt,
			      bond_tri_number(bt,intfc),
			      index_of_pointer(
				  (POINTER*)c->interface->surfaces,
				  (POINTER)Surface_of_tri(bt->tri)),
			      Coords(b->end)[0],Coords(b->end)[1],
			      Coords(b->end)[2]);
		if (is_obstacle_state(ls))
		    (void) printf(" %40s|",
				  "is_obstacle_state( left state )      ");
		else
		    (void) printf(" %g %g %g %g %g |",
				  Dens(ls),Energy(ls),
				  Mom(ls)[0],Mom(ls)[1],Mom(ls)[2]);

		if (is_obstacle_state(rs))
		    (void) printf(" %40s\n",
				  "is_obstacle_state( right state )     ");
		else
		    (void) printf(" %g %g %g %g %g\n",
				  Dens(rs),Energy(rs),
				  Mom(rs)[0],Mom(rs)[1],Mom(rs)[2]);
	    }
	    (void) printf("\n");
	}

	(void) printf("\n");
	(void) printf(" end of summary_of_curve_states3d() curve = %p\n",c);
	(void) printf("\n");
} 		/*end summary_of_curve_states3d*/

EXPORT  void 	summarize_front_states3d(
	Front 		*front)
{
	INTERFACE		*intfc = front->interf;
	double			*coords;
	int			cnt = 0;
	POINT			*p;
	HYPER_SURF_ELEMENT 	*hse;
	HYPER_SURF		*hs;
	Locstate		sl, sr;

	(void) printf("\n");
	(void) printf("summarize_front_states3d()  front = %p\n",front);

	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (30*(cnt/30) == cnt)
	    {
		(void) printf("\n");
		(void) printf("cnt = %d\n",cnt);
		(void) printf("%5s %5s %5s %8s %8s %8s | "
			      "%8s %4s %4s %4s %4s %4s %8s | "
			      "%8s %4s %4s %4s %4s %4s %8s\n",
			      "x ","y ","z ","p ","hse ","hs ","l_state","dens",
			      "engy","mx ","my ","mz ","Params","r_state",
			      "dens","engy","mx ","my ","mz ","Params");
		(void) printf("%5s %5s %5s %8s %8s %8s | "
			      "%8s %4s %4s %4s %4s %4s %8s | "
			      "%8s %4s %4s %4s %4s %4s %8s\n",
			      "====","====","====","========","========",
			      "========","========","====","====","====","====",
			      "====","========","========","====","====","====",
			      "====","====","========");
	    }

	    slsr(p,hse,hs,&sl,&sr);  
	    coords = Coords(p);

	    (void) printf("%g %g %g %llu %llu %llu | ",
			  coords[0],coords[1],coords[2],
			  point_number(p),
			  hypersurface_element_number(hse,intfc),
			  hypersurface_number(hs));
	    cnt++;
	    (void) printf("%llu ",ptr2ull(sl));
	    if (is_obstacle_state(sl))
	        (void) printf("%32s | ","---------OBSTACLE STATE--------");
	    else
	        (void) printf("%g %g %g %g %g %llu | ",
			      Dens(sl),
			      Energy(sl),
			      Mom(sl)[0],Mom(sl)[1],Mom(sl)[2],
			      gas_param_number(Params(sl)));

	    (void) printf("%llu ",ptr2ull(sr));
	    if (is_obstacle_state(sr))
	        (void) printf("%32s\n","---------OBSTACLE STATE--------");
	    else
	        (void) printf("%g %g %g %g %g %llu\n",
			      Dens(sr),
			      Energy(sr),
			      Mom(sr)[0],Mom(sr)[1],Mom(sr)[2],
			      gas_param_number(Params(sr)));
	}
	(void) printf("END of summarize_front_states3d() front = %p  "
		      "num_points processed = %d\n",front,cnt);
}		/*end summarize_front_states3d*/
