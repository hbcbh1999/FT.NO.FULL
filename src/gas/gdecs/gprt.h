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
*				gprt.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if !defined(_GPRT_H)
#define _GPRT_H

struct _G_Printplot {
	Printplot Prt;

	int	_problem;
};
typedef struct _G_Printplot G_Printplot;

#define g_printplot(prt)	((G_Printplot *)(prt))
#define prt_problem_type(prt)	g_printplot(prt)->_problem

	/* Data structures for additional GAS diagnostics */

typedef struct {
	OUTPUT_DATA	odata;
	RECT_GRID 	rect_grid;
	double		min_var[5], max_var[5], spread_var[5];
	int		clip_var[5], use_logs[5];
	int		frame_number;
	int		compress;
} Woodward_Porter_data;

typedef struct {
	OUTPUT_DATA	odata;
	RECT_GRID	rect_grid;
} Error_Grid_data;

typedef struct {
	OUTPUT_DATA	odata;
#if defined(TWOD) || defined(THREED)
	SINE_PERT	*pert;
#endif /* defined(TWOD) || defined(THREED) */
	int             p;
	int		ratio;
	int		stype;
	void		(*_alternate_state)(Locstate,double*,Front*,Wave*,int);
} Lp_Diff_data;

struct _PROBE {
	OUTPUT_DATA odata;
	int         dim;
	double       c[3];
	double       r[3];
	double       L;
	int         n[3];
	int         N;
	int         line;
	Locstate    st[3];
	Locstate    dst, Ndst;
	double       time[3];
	COMPONENT   *pcomp, *pcompstore;
	Locstate    *pstate, *pstatestore;
	double       *pcoords, *pcoordsstore;
	double       ptol;
	double       wave_tag; /* 1 = inside compression,
	                        -1 = inside rarefaction,
				 0 = otherwise */
};
typedef struct _PROBE PROBE;
#define	probe_index(i,j,k,p)						\
  ((p)->n[2]+(k)+(2*(p)->n[2]+1)*(j+(p)->n[1]+(2*(p)->n[1]+1)*((i)+(p)->n[0])))
#define probe_comp(i,j,k,p)						\
	(p)->pcomp[probe_index(i,j,k,p)]
#define probe_state(i,j,k,p)						\
	(p)->pstate[probe_index(i,j,k,p)]
#define probe_coords(i,j,k,p)						\
	((p)->pcoords+3*probe_index(i,j,k,p))
	

#define EG_data(data)		((Error_Grid_data *)data)
#define WP_data(data)		((Woodward_Porter_data *)data)
#define Lp_data(data)		((Lp_Diff_data *)data)

#if defined(TWOD)
#define print_RP_DATA(RP,v)		fprint_RP_DATA(stdout,RP,v)
#define	print_curve_status(message,status)				\
		fprint_curve_status(stdout,message,status)
#endif /* defined(TWOD) */
#define print_gas_data(state,state_type)	fprint_gas_data(stdout,state)

#endif /* !defined(_GPRT_H) */
