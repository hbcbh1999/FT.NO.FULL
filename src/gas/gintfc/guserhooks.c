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
*                               guserhooks.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Gas Dynamic Extensions to Interface User Supplied Operations
*/

#include <gdecs/gdecs.h>

EXPORT	void w_speed(
	double		*pt,
	Locstate	sl,
	Locstate	sr,
	Locstate	ansl,
	Locstate	ansr,
	double		*W,
	double		pjump,
	double		*nor,
	int		w_type,
	Front		*front)
{
        INTERFACE       *infc;

	if (front == NULL)
	{
            screen("\nERROR in w_speed(), NULL front !!!\n");
            clean_up(ERROR);
	}
        infc = front->interf;
	if (g_user_interface(infc)._w_speed == NULL)
	{
                screen("\nERROR in w_speed(), ");
                screen("Function pointer is empty!!!\n");
                clean_up(ERROR);
	}
	if (ansl == NULL)
	{
	    static Locstate temp_ansl = NULL;
	    if (temp_ansl == NULL)
	        alloc_state(infc,&temp_ansl,front->sizest);
	    ansl = temp_ansl;
	}
	if (ansr == NULL)
	{
	    static Locstate temp_ansr = NULL;
	    if (temp_ansr == NULL)
	        alloc_state(infc,&temp_ansr,front->sizest);
	    ansr = temp_ansr;
	}
	if (W == NULL)
	{
	    static double temp_W[3];
	    W = temp_W;
	}

	(*g_user_interface(infc)._w_speed)(pt,sl,sr,ansl,ansr,
					   W,pjump,nor,w_type,front);
}		/*end w_speed*/


EXPORT  void npt_w_speed(
	WSSten		*sten,
        Locstate        ansl,
        Locstate        ansr,
        double           *W)
{
        INTERFACE       *infc;

	if (sten == NULL)
	{
            screen("\nERROR in npt_w_speed(), NULL WSSten!!!\n");
            clean_up(ERROR);
	}
	if (sten->front == NULL)
	{
            screen("\nERROR in npt_w_speed(), NULL front in WSSten!!!\n");
            clean_up(ERROR);
	}
        infc = sten->front->interf;
	if (g_user_interface(infc)._npt_w_speed == NULL)
	{
            screen("\nERROR in npt_w_speed(), Function pointer is empty!!!\n");
            clean_up(ERROR);
	}

	if (ansl == NULL)
	{
	    static Locstate temp_ansl = NULL;
	    if (temp_ansl == NULL)
	        alloc_state(infc,&temp_ansl,sten->front->sizest);
	    ansl = temp_ansl;
	}
	if (ansr == NULL)
	{
	    static Locstate temp_ansr = NULL;
	    if (temp_ansr == NULL)
	        alloc_state(infc,&temp_ansr,sten->front->sizest);
	    ansr = temp_ansr;
	}
	if (W == NULL)
	{
	    static double temp_W[3];
	    W = temp_W;
	}

	(*g_user_interface(infc)._npt_w_speed)(sten,ansl,ansr,W);
}	/* end npt_w_speed() */

EXPORT  void unsplit_w_speed2d(
	USWSSten2d	*uswssten,
	Locstate	ansl,
	Locstate	ansr,
        double           *W)
{
        INTERFACE       *infc = uswssten->fr->interf;

/*
*
*	6---------4---------8	---> tan_fr(uswssten)[-1]
*       |         |         | 
*       |         |         | 
*       |         |         |	   
*       |         |         |	
*       |         |         |
*	|         |         |
*	1---->----0---->----2   ---> tan_fr(uswssten)[0]
*	|         |         |
*       |         |         |      W is the front speed of point 0.
*       |         |         | 
*	|	  |         |    uswssten->nor_ans->leftst[0],
*       |         |         |           uswssten->nor_ans->rightst[0]
*       |         |         |           state in point 0 at next time level.
*	5---------3---------7	---> tan_fr(uswssten)[1]
*
*	
*/

        if (g_user_interface(infc)._unsplit_w_speed2d == NULL)
        {
                screen("\nERROR in unsplit_w_speed2d(), ");
                screen("Function pointer is empty!!!\n");
                clean_up(ERROR);
        }
 
        (*g_user_interface(infc)._unsplit_w_speed2d)(uswssten,ansl,ansr,W);
}       /* end unsplit_w_speed2d */
