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
*				ghypvec.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains the drivers for the uni_arrayized hyperbolic scheme:
*/


#include <ghyp/ghyp.h> 
#include <gdecs/vecdecs.h>


/*
*			  vector_FD()
*
*	The subroutine first loads the state variables into vst, the source
*	variables to src, then passes to the oned_interior_scheme().
*/

/*ARGSUSED*/
EXPORT	void vector_FD(
	int		swp_num,
	int		*iperm,
	double		*dir,
	Wave		*wv,
	Wave		*newwv,
	Front		*fr,
	Front		*newfr,
	int		*icoords,
	int		imin,
	int		imax,
	double		dt,
	double		dh)
{
	Vec_Gas		*vst = g_wave_vgas(wv);
	Vec_Src		*src = g_wave_vsrc(wv);
	int		idir = iperm[swp_num];
	int		dim = wv->rect_grid->dim;
	int		nrad = vsten_radius(wv);
	int		vsize;
	static	double	alpha;
	static	boolean	first = YES;

#if defined(TIME_HYPVEC)
	start_clock("vector_FD");
#endif /* defined(TIME_HYPVEC) */

	debug_print("vector_FD","Entered vector_FD(), dir = %d\n",idir);
	if (is_rotational_symmetry() && first == YES)
	{
	    first = NO;
	    alpha = rotational_symmetry();
	}

	vsize = imax - imin;
	if (vsize < wv->npts_vsten) return;
	if (load_state_vectors(swp_num,iperm,vst,0,vsize,wv,newwv,
			       icoords,imin) == CONSTANT_IN_TIME)
	{
#if defined(TIME_HYPVEC)
	    stop_clock("vector_FD");
#endif /* defined(TIME_HYPVEC) */
	    return;
	}

	if (is_rotational_symmetry() &&  (alpha > 0.0) && (iperm[swp_num] == 0))
 	{
 	    RECT_GRID *rgr = wv->rect_grid;
	    double     *radii = src->radii;
	    int       j;
 
	    src->rmin = pos_radius(0.0,rgr);
 	    for (j = 0; j < vsize; j++)
 	    	radii[j] = pos_radius(cell_center(j+imin,0,rgr),rgr);
 	}

	oned_interior_scheme(swp_num,iperm,icoords,wv,newwv,fr,newfr,NULL,
		             0,vsize,vst,src,dt,dh,dim);

	assign_wave_state_vectors(swp_num,iperm,wv,newwv,vst,nrad,
		                  vsize-nrad,icoords,imin);
	check_and_correct_bad_state(swp_num,iperm,icoords,fr,newfr,wv,newwv,
					dh,dt,dir,nrad,vsize-nrad,imin);
	if (debugging("my_debug"))
	{
	    remove_from_debug("my_debug");
	}

	debug_print("vector_FD","Left vector_FD(), dir = %d\n",idir);
#if defined(TIME_HYPVEC)
	stop_clock("vector_FD");
#endif /* defined(TIME_HYPVEC) */
}		/*end vector_FD*/


EXPORT	void copy_states_to_new_time_level(
	int		*idirs,
	Wave		*wv,
	Wave		*newwv,
	int		imin,
	int		imax,
	int		*icoords,
	int		pbuf)
{
	int		i;
	size_t		sizest = wv->sizest;

	for (i = imin; i < imax; ++i)
	{
	    icoords[idirs[0]] = i + pbuf;
	    ft_assign(Rect_state(icoords,newwv),Rect_state(icoords,wv),sizest);
	}
}		/*end copy_states_to_new_time_level*/


/*	The following function uses Lax-Friedrich scheme to correct
*       states not properly solved by high order numerical schemes.
*       If debugging flag is turned on, it will report the bad state
*       and call clean_up().
*/

EXPORT	void check_and_correct_bad_state(
	int swp_num,
	int *iperm,
	int *icoords,
	Front *fr,
	Front *newfr,
	Wave *wv,
	Wave *newwv,
	double dh,
	double dt,
	const double *dir,
	int imin,
	int imax,
	int pbuf)
{
	int i,j,k,idirs[MAXD];
	int flag;
	static Stencil *sten;
	int icrds_sten[MAXD];
	double *coords;
	int dim = fr->rect_grid->dim;
	Locstate state;
	const char *message = "In check_and correct_bad_state(), LF failed\n";

	if (sten == NULL) 
	    sten = alloc_stencil(3,fr); /*allocate stencil for LF*/
	sten->fr = fr; 		sten->newfr = newfr;
	sten->wave = wv; 	sten->newwave = newwv;

	if (debugging("bad_state")) flag = YES;
	else flag = NO;
        for (i = 0; i < dim; ++i)
            idirs[i] = iperm[(i+swp_num)%dim];
        for (i = imin; i < imax; ++i)
        {
	    icoords[idirs[0]] = i + pbuf;
	    if (is_obstacle_state(Rect_state(icoords,newwv)))
		continue;
	    state = Rect_state(icoords,newwv);
	    if (is_bad_state(state,flag,"check_and correct_bad_state"))
	    {
		for (k = 0; k < dim; k++)
		    icrds_sten[k] = icoords[k];
		for (j = -1; j <= 1; ++j)
		{
		    icrds_sten[idirs[0]] = icoords[idirs[0]] + j;
		    sten->st[j] = Rect_state(icrds_sten,wv);
		    coords = Rect_coords(icrds_sten,wv);
		    for (k = 0; k < dim; ++k)
		    	Coords(sten->p[j])[k] = coords[k];
		}
		LF(dh,dt,state,dir,swp_num,iperm,NULL,sten);
		if (is_bad_state(state,YES,message))
		    clean_up(ERROR);
	    }
	}
	    
}	/* end check_and correct_bad_state */
