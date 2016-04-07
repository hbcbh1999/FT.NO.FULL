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
*				damr.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Structures for the automatic mesh refinement code.
*/

#if !defined(_DAMR_H)
#define _DAMR_H

#include <driver/ddecs.h>

#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif



/*
*	A CHART is a structure for automatic mesh refinement by means of
*	overlapping coordinate systems, each defined through an affine-linear
*	coordinate transformation and a rectangle with sides parallel to the
*	coordinate system in the new coordinates.
*/

struct _CHART {
	struct _CHART *parent;
	struct _CHART *prev_chart;	/* Prev chart at level */
	struct _CHART *next_chart;	/* Next chart at level */

	struct _AFLIN *to_root;
	struct _AFLIN *from_root;

	struct _LEVEL *level;

	int dynamic;			/* YES if chart can be remeshed */
	int is_old_chart;               /* YES if might be removed after
					   regridding  */

	Grid	*grid;
	Front	*front;
	Wave	*wave;
	Front	*newfront;
	Wave	*newwave;
	struct _Printplot *prt;
	POINTER	el_map;

#if defined(USE_OVERTURE)
        Overparam                 *overparam;
        POINTER                   old_cg_over;  /* pointer of CompositeGrid */
        POINTER                   cg_over_function; /* pointer of doubleCompositeGridFunction */
        int                       totalNumberOfPatches;
        int                       use_overture_state;
#endif /* defined(USE_OVERTURE) */

	/* basic one chart integration step */
	int (*hyp_solver)(double,double*,Wave*,Front*);
	int (*parab)(double,double*,Wave*,Front*);
#if defined(SUBGRID)
        void (*subgrid)(double,Front*,Wave*);
#endif /* defined SUBGRID */
#if defined (USE_OVERTURE)
        void (*parab_npt)(double,Front*,Wave*,Wave*);
#endif /* if defined (USE_OVERTURE) */
	void (*ellip)(INIT_DATA*,struct _CHART*,const IO_TYPE*);
	void (*bc_propagate)(Grid*);
	int  (*make_bubbles)(Wave*, Front*);
} ;

typedef struct _CHART CHART;

struct _LEVEL {
	int level_num;			/* 0 is root (coarsest) level */
	int *max_num_levels;		/* num ref. levels currently in use */

		/* Neighboring levels and charts at level */

	struct _LEVEL *next_coarser_level;
	struct _LEVEL *next_finer_level;
	struct _CHART *first,*last;	/* first and last on level */

		/* Common information for all charts at level */

	int steps_before_synchronize;
	int steps_before_regrid;
	int steps_before_level_change;
	int steps_between_synchronize;
	int steps_between_regrid;
	int refinement_ratio;

	int step;
	double t;
	double dt;
} ;

typedef struct _LEVEL LEVEL;


#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

#endif /* !defined(_DAMR_H) */


