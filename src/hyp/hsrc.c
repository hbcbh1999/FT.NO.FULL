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
*				hsrc.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains routines for the inclusion of source terms in the
*	hyperbolic equations.
*/


#include <hyp/hdecs.h>


	/* LOCAL Function Declarations */
LOCAL	void	source_intfc_node(double*,int,INTERFACE*,RECT_GRID*,NODE**);


/*
*		set_pt_source_interior_vectors():
*
*	This routine sets up the uni_arrays iscrv, isrc_min, isrc_max
*	in the wave structure.  These are used to speed up the 
*	calculation in is_source_block().
*
*	A point source/sink influences a row (column) of the mesh if
*		1) the source/sink lies within the row (column), or
*		2) the midpoint of the row (column) is within the 
*		   rectangle of influence around the source/sink.
*	When a source/sink lies on the boundary of a row (column),
*	the source/sink is associated with the row (column) with
*	the larger index (unless this is outside the mesh).
*/

#define SET_VECTORS(_zpoint,_ZL,_hz,_rz,_zmax,_izv,_iz_min,_iz_max)	\
{									\
	int	_i, i_min, i_max, int_i;				\
	double	real_i;							\
									\
		/* find mesh index of source/sink */			\
	real_i = (_zpoint - _ZL) / _hz;					\
	int_i = (int)(real_i);						\
									\
		/* may happen for source/sink on or near boundary */	\
	if (int_i < 0     ) int_i = 0;					\
	if (int_i >= _zmax) int_i = _zmax-1;				\
									\
		/* find range of influence of source/sink */		\
	i_min =  (int)(ceil(real_i - (_rz + 0.5)));			\
	if (i_min < 0)	i_min = 0;					\
	i_max = (int)(floor(real_i + (_rz - 0.5)));			\
	if (i_max >= _zmax) i_max = _zmax-1;				\
									\
		/* source rect does not overlap any mesh midpoints, */	\
	if (i_min > i_max) i_max = i_min = int_i;			\
									\
	for (_i = i_min; _i <= i_max; ++_i) ++_izv[_i];			\
	_iz_min = i_min;		_iz_max = i_max;		\
}

EXPORT	void set_pt_source_interior_vectors(
	Wave		*wave)
{
	int **isrcv, **isrc_min, **isrc_max;
	int nsrc, num, i;
	int j, dim = wave->rect_grid->dim;

	nsrc = wave->num_point_sources;
	if (nsrc <= 0)
	{
	    wave->isrcv = NULL;
	    wave->isrc_min = wave->isrc_max = NULL;
	    return;
	}

	uni_array(&isrcv,dim,sizeof(int *));
	for (j = 0; j < dim; ++j)
	    uni_array(isrcv + j,wave->rect_grid->gmax[j],INT);

	bi_array(&isrc_min,dim,nsrc,INT);
	bi_array(&isrc_max,dim,nsrc,INT);

	wave->isrcv    = isrcv;
	wave->isrc_min = isrc_min;
	wave->isrc_max = isrc_max;

	for (j = 0; j < dim; ++j)
	    for (i = 0;  i < wave->rect_grid->gmax[j];  ++i)
	    	isrcv[j][i] = 0;

	for (num = 0;  num < nsrc;  ++num)
	{
	    double	X, h, r;
	    int	imax;
	    for (j = 0; j < dim; ++j)
	    {
	    	X = wave->rect_grid->L[j];
	    	h = wave->rect_grid->h[j];
	    	/* ratio (half) source rect side to mesh side */
	    	r = 0.5*wave->pt_src_diam[j][num] / h;
	    	imax = wave->rect_grid->gmax[j];
	    	SET_VECTORS(wave->srcpt[j][num],X,h,r,imax,isrcv[j],
			    isrc_min[j][num],isrc_max[j][num])
	    }
	}
	if (debugging("sources"))
	{
	    int	i_stop = 0;

	    for (j = 0; j < dim; ++j)
	    	i_stop = max(i_stop,wave->rect_grid->gmax[j]);
	    (void) printf("\nset_pt_source_interior_vectors()\n"
	                  "column, row influence flags: isrcv[j][]\n");
	    for (i = 0;  i < i_stop;  ++i)
	    {
	    	(void) printf("%4d  ",i);
	    	for (j = 0; j < dim; ++j)
	    	{
	    	    if (i < wave->rect_grid->gmax[j])
	    	    	(void) printf("%4d  ",isrcv[j][i]);
	    	    else
	    	    	(void) printf("      ");
	    	}
	    }
	    (void) printf("source influence indices: \n");
	    (void) printf("source isrc_min isrc_max\n");
	    for (num = 0;  num < nsrc;  ++num)
	    {
	    	(void) printf("%d ",num);
	    	for (j = 0; j < dim; ++j)
	    	    (void) printf("%d->%d  ",isrc_min[j][num],isrc_max[j][num]);
	    	(void) printf("\n");
	    }
	}
}		/*end set_pt_source_interior_vectors*/

/*
*			free_pt_source_interior_vectors():
*
*	Frees the uni_arrays allocated by set_pt_source_interior_vectors().
*/

EXPORT	void free_pt_source_interior_vectors(
	Wave		*wave)
{
	int		j, dim = wave->rect_grid->dim;

	if (wave->isrcv != NULL)
	{
	    for (j = 0; j < dim; ++j)
	    	if (wave->isrcv[j] != NULL)
		    free(wave->isrcv[j]);
	    free(wave->isrcv);
	}
	if (wave->isrc_min != NULL)
	    free(wave->isrc_min);
	if (wave->isrc_max != NULL)
	    free(wave->isrc_max);
}		/*end free_pt_source_interior_vectors*/


/*
*			is_source_block():
*
*	Given a mesh block icoords, this routine determines whether a
*	point source/sink injecting/producing component comp influences
*	the mesh block and returns the index of the point source (as
*	recorded in the wave structure).  An index of -1 indicates no 
*	source/sink.
*
*	It is currently a fatal error if more than one source/sink influences
*	the well.
*/

EXPORT	int is_source_block(
	Wave		*wave,
	INTERFACE	*intfc,
	COMPONENT	comp,
	int		*icoords)
{
	NODE		  *n;
	double		  coords[MAXD];
	int		  num, stype;
	int		  i, j, dim = intfc->dim;
	int		  comps_match;
	int		  deb_sources = debugging("sources");
	int		  index;

	if (deb_sources) 
	{
	    (void) printf("is_source_block(): ");
	    for (j = 0; j < dim; ++j)
	    	(void) printf(" %d ",icoords[j]);
	    (void) printf("comp %d\n",comp);
	}

	index = -1;
	if (wave->num_point_sources <= 0)
	    return index;
	for (j = 0; j < dim; ++j)
	{
	    if ((wave->isrcv[j][icoords[j]] == 0)) 
	    {
	        if (deb_sources) 
	        {
	    	    (void) printf("block is not in column or row list, "
	    	                  "return: index %d\n",index);
	        }
	        return index;
	    }
	}

	for (num = 0;  num < wave->num_point_sources;  ++num)
	{
			/* check if source influences block icoords */

	    for (j = 0; j < dim; ++j)
	    {
	        coords[j] = wave->srcpt[j][num];
	        if (deb_sources) 
	        {
	    	    (void) printf("testing source %d  ",num);
	    	    (void) printf("coord dir %d  rows %d->%d\n",
	    		          j,wave->isrc_min[j][num],
	    		          wave->isrc_max[j][num]);
	    	    if ((wave->isrc_min[j][num] > icoords[j]) ||
	    	        (wave->isrc_max[j][num] < icoords[j]))
	    	    {
	    	        (void) printf(" - does not influence block\n");
	            }
	        }
		if ((wave->isrc_min[j][num] > icoords[j]) ||
		    (wave->isrc_max[j][num] < icoords[j]))
		    break;
	    }
	    if (j < dim)
		continue;

			/* find the interface node at the source */

	    stype = (wave->source_type[num]==SOURCE) ? SOURCE_NODE : SINK_NODE;
	    source_intfc_node(coords,stype,intfc,wave->rect_grid,&n);
	    if (deb_sources) 
	    {
	    	for (j = 0; j < dim; ++j)
	    	    (void) printf("s[%d] %g",j,coords[j]);
	    	(void) printf(" stype %s\n",(stype == SOURCE_NODE) ?
	    		      "SOURCE_NODE" : "SINK_NODE");
	    	(void) printf("source_intfc_node\n");
	    	print_node(n);
	    }
	    if (n == NULL)
	    {
		screen("ERROR in is_source_block(), No intfc node "
		       "corresponding to %s node found at ",
		       (stype == SOURCE_NODE) ? "source":"sink");
	    	for (j = 0; j < dim; ++j)
	    	    screen("%g ",coords[j]);
	    	screen("\n");
	    	print_pt_sources(wave);
	    	print_interface(intfc);
	    	clean_up(ERROR);
	    }

		/* check for component match */

	    comps_match = NO;
	     if ((n->in_curves == NULL) && (n->out_curves == NULL))
	    {
	        if (deb_sources)
	    	    (void) printf("checking component at node position\n");
	        if (component(coords,intfc) == comp)
		    comps_match = YES;
	    }
#if defined(TWOD)
	    else
	    {
	        int i;
	        CURVE **c, **c_beg;

	        if (deb_sources)
	            (void) printf("checking components of curves at node\n");
		for (i = 0,  c_beg = n->in_curves; i < 2; 
		     ++i,    c_beg = n->out_curves)
		{
		    for (c = c_beg;  c && *c;  ++c)
		    {
			if ((negative_component(*c) == comp) ||
			    (positive_component(*c) == comp))
			{
			    comps_match = YES;
			    break;
			}
		    }
		    if (comps_match)
			break;
		}
	    }
#endif /* defined(TWOD) */
	    if (deb_sources)
	        (void) printf("comps_match %d\n",comps_match);
	    if (!comps_match)
		continue;
		
			/* found a valid source/sink influencing block */

	    if (index == -1)
	        index = num;
	    else
	    {
	        screen("ERROR in is_source_block(), ");
	        for (i = 0; i < dim; ++i)
	            screen("i[%d] %d ",i,icoords[i]);
	        screen(" comp %d",comp);
	        screen(" - more than one source/sink influences block\n");
	        print_pt_sources(wave);
	        clean_up(ERROR);
	    }

	    	/* may have to check for more than 1 source/sink */
	    for (j = 0; j < dim; ++j)
	    {
	        if (wave->isrcv[j][icoords[j]] <= 1)
	    	    break;
	    }
	    if (j == dim)
	        continue;
	    else
	        break;

	}

	/* Reach here if					*/
	/*  1)	valid source in block found,			*/
	/*  2)	no sources in block had correct component, or 	*/
	/*  3)	there was not really a source in the block	*/
	/* 	due to "cast_shadow" nature of using ixv[] iyv[]*/
	/* 	instead of ixyv[][] to flag source/sink blocks.	*/
	/* 	Mesh square coords was lying in the "shadow" of	*/
	/*      two other sources.				*/

	if (deb_sources)
	    (void) printf("return: index %d\n",index);
	return index;
}		/*end is_source_block*/

/*
*			source_intfc_node():
*
*	Returns the interface node corresponding to a given source/sink
*/

LOCAL void source_intfc_node(
	double		*coords,
	int		stype,
	INTERFACE	*intfc,
	RECT_GRID	*gr,
	NODE		**nwell)
{
	NODE		**n;
	double		tol_dist;
	double		old_dist, dist;
	int		i, dim = gr->dim;

	*nwell = NULL;
	tol_dist = gr->h[0];
	for (i = 1; i < dim; ++i)
	    if (tol_dist > gr->h[i])
	        tol_dist = gr->h[i];
	old_dist = tol_dist;/*TOLERANCE*/
	for (n = intfc->nodes;  n && *n;  ++n)
	{
		if (node_type(*n) != stype) continue;
		dist = distance_between_positions(coords,
						  Coords((*n)->posn),dim);
		if (dist < old_dist)
		{
			old_dist = dist;
			*nwell = *n;
		}
	}
		/* check if any node close enough */
	if (old_dist > 0.01 * tol_dist) *nwell = NULL;
}		/*end source_intfc_node*/
