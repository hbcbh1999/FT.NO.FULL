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
*				guserrp.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*			User Supplied Operations
*/

#if !defined(_USERRP_H)
#define _USERRP_H

#include <gdecs/gdecs.h>

struct _G_RPROBLEM {
	F_RPROBLEM frproblem;
	int _num_incident;
        int _num_refl;
	int _num_transm;
	int _num_inc_shock;
	int _num_contact;
};
typedef struct _G_RPROBLEM G_RPROBLEM;

	/* G_RPROBLEM access macros */
#define g_rproblem(rp)  ((G_RPROBLEM *) (rp))
#define	rp_num_transm(rp)		(g_rproblem(rp)->_num_transm)
#define	rp_num_inc_shock(rp)		(g_rproblem(rp)->_num_inc_shock)
#define	rp_num_contact(rp)		(g_rproblem(rp)->_num_contact)


struct _G_RP_NODE {
	RP_NODE rpn;
	O_CURVE_FAMILY *_inc_shock1;
	O_CURVE_FAMILY *_reflected1;
	O_CURVE_FAMILY *_transmitted;
	O_CURVE_FAMILY *_inc_shock2;
	O_CURVE_FAMILY *_reflected2;
	O_CURVE_FAMILY *_contact1;
	O_CURVE_FAMILY *_contact2;
};
typedef struct _G_RP_NODE G_RP_NODE;

	/* G_RPROBLEM access macros */
#define g_rp_node(rpn)  ((G_RP_NODE *) (rpn))
#define	rpn_inc_shock1(rpn)	(g_rp_node(rpn)->_inc_shock1)
#define	rpn_transmitted(rpn)	(g_rp_node(rpn)->_transmitted)
#define	rpn_inc_shock2(rpn)	(g_rp_node(rpn)->_inc_shock2)
#define	rpn_contact1(rpn)	(g_rp_node(rpn)->_contact1)
#define	rpn_contact2(rpn)	(g_rp_node(rpn)->_contact2)



#endif /* !defined(_USERRP_H) */
