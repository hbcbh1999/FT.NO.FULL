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
*				hnpt.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains n-point stencil code for finite difference solution
*	of the hyperbolic equations within a single component of a tracked
*	interface problem.
*	Point sources in the region are found and passed appropriately.
*/


#include <hyp/hdecs.h>

enum _CRX_SIDE {
	PREV = -1,
	NEXT =  1
};
typedef enum _CRX_SIDE CRX_SIDE;


LOCAL	int	gmax[MAXD], vsten_rad, is_vec_method;
LOCAL	int	dim;
LOCAL   int	lbuf[MAXD], ubuf[MAXD];
LOCAL	double	h[MAXD];
LOCAL	double	L[MAXD], VL[MAXD], VU[MAXD];

	/* Redefine rect grid macros for local efficiency */
	/* these will have to be modified to support variable mesh grids */
#undef cell_center
#define cell_center(indx,i,gr)	(L[i] + ((indx) + 0.5)*h[i])

	/* LOCAL Function Declarations */
LOCAL	COMPONENT	set_stencil_comps_and_crosses(int,Stencil*);
LOCAL	boolean	interior_sten_reg(int,int,Stencil*);
LOCAL	void	copy_states(int,CRX_SIDE,int,int,Stencil*,CRXING*);
LOCAL	void	find_outer_states(int,int,Stencil*,CRX_SIDE);
LOCAL	void	set_dirichlet_bdry_sten_states(int,HYPER_SURF*,Stencil*,
                                               CRX_SIDE);
LOCAL	void	set_neumann_bdry_sten_states(int,int,int,COMPONENT,CRXING*,
					     Stencil*,CRX_SIDE);
LOCAL	void	update_reg_grid_state(double,double,int,int*,int,
				      Stencil*,POINT*,COMPONENT);
LOCAL   boolean Debug;

#if defined(TWOD) || defined(THREED)
LOCAL	void	set_static_coords(int**,POINT*,int,int,RECT_GRID*);
#endif /* defined(TWOD) || defined(THREED) */

EXPORT	void set_hyp_npt_globals(
	Wave		*wave)
{
	RECT_GRID	*r_grid;
	int		i;

	r_grid = wave->rect_grid;
	dim = r_grid->dim;
	for (i = 0; i < dim; ++i)
	{
	    gmax[i] = r_grid->gmax[i];
	    h[i] = r_grid->h[i];
	    L[i] = r_grid->L[i];
	    VL[i] = r_grid->VL[i];
	    VU[i] = r_grid->VU[i];
	    lbuf[i] = r_grid->lbuf[i];
	    ubuf[i] = r_grid->ubuf[i];
	}
	vsten_rad = vsten_radius(wave);
	is_vec_method = is_vector_method(wave);
}		/*end set_hyp_npt_globals*/

/*
*	The following globals are shared between hyp_npt() and
*	the subroutine update_reg_grid_state().  They are used
*	to pass quantities that do not change inside the inner loop
*	of the update step.  The are passed as globals for efficiency.
*/

LOCAL	COMPONENT	ext_comp;
LOCAL	double		hdir, dir[MAXD];
LOCAL	GRID_DIRECTION	prev_side, next_side;
LOCAL	INTERFACE	*intfc = NULL;
LOCAL	int		endpt;


/*
*			hyp_npt():
*
*	Single direction sweep for a n-point scheme.
*/

EXPORT void hyp_npt(
	int		swp_num,
	int		*iperm,	    /* sweep based permutation of coord dir */
	double		ds,
	double		dt,
	Wave		*wave,
	Wave		*newwave,
	Front		*fr,
	Front		*newfr,	    /* newfr needed if hlw is to support */
		    		    /*	changing topologies	     */
	COMPONENT	max_comp)
{
	POINT		Basep;
	int 		i, idirs[MAXD];
	int		imin[3], imax[3];
	static	Stencil	*sten = NULL;
	static	int	num_pts;
#if defined(TWOD) || defined(THREED)
	static	int	**icoords = NULL;
	RECT_GRID	*gr = wave->rect_grid;
#endif /* defined(TWOD) || defined(THREED) */

	for (i = 0; i < dim; ++i)
	{
	    idirs[i] = iperm[(i+swp_num)%dim];
	    dir[i] = 0.0;
	}
	dir[idirs[0]] = 1.0;
	debug_print("hyp_npt","Entered hyp_npt(), swp_num = %d, dir = %d\n",
			swp_num,idirs[0]);

	intfc = fr->interf;
	ext_comp = exterior_component(intfc);

	if (sten == NULL)
	{
	    /* This assumes that the stencil size is constant. */
	    num_pts = wave->npts_sten;
	    endpt = stencil_radius(wave);

	    sten = alloc_stencil(num_pts,fr);
#if defined(TWOD) || defined(THREED)
	    icoords = sten->icoords;
#endif /* defined(TWOD) || defined(THREED) */
	}

	hdir = ds;		/* grid spacing in sweep direction */
	switch (idirs[0])
	{
	case 0:	/*X SWEEP */
	    prev_side = WEST;
	    next_side = EAST;
	    if (debugging("hyp_npt"))
	    	(void) printf("prev_side = WEST, next_side = EAST\n");
	    break;
	case 1: /*Y SWEEP */
	    prev_side = SOUTH;
	    next_side = NORTH;
	    if (debugging("hyp_npt"))
	    	(void) printf("prev_side = SOUTH, next_side = NORTH\n");
	    break;
	case 2: /* Z SWEEP */
	    prev_side = LOWER;
	    next_side = UPPER;
	    if (debugging("hyp_npt"))
	    	(void) printf("prev_side = LOWER, next_side = UPPER\n");
	    break;
	}

	sten->fr   = fr;	sten->newfr   = newfr;
	sten->wave = wave;	sten->newwave = newwave;

	set_sweep_limits(wave,swp_num,idirs,imin,imax);
	/*#bjet2 */
	set_limits_for_open_bdry(wave,fr,swp_num,idirs,imin,imax);

	switch (dim)
	{
#if defined(ONED)
	case 1:
	{
	    int	i0;
	    int	i0max;
	    int	i0min;

	    i0min = imin[0];	i0max = imax[0]; 
	    sten->reg_stencil = NO;
	    for (i0 = i0min; i0 < i0max; ++i0) 
	    {
	    	update_reg_grid_state(ds,dt,swp_num,iperm,
				      i0,sten,&Basep,max_comp);
	    }
	    break;
	}
#endif /* defined(ONED) */
#if defined(TWOD)
	case 2:
	{
	    int	i0, i1;
	    int	i0max, i1max;
	    int	i0min, i1min;

	    i0min = imin[0];	i0max = imax[0]; 
	    i1min = imin[1];	i1max = imax[1]; 
	    for (i1 = i1min; i1 < i1max; ++i1)
	    {
	    	set_static_coords(icoords,&Basep,idirs[1],i1,gr);
	    	sten->reg_stencil = NO;
	    	for (i0 = i0min; i0 < i0max; ++i0) 
	    	{
	    	    update_reg_grid_state(ds,dt,swp_num,iperm,
	    	    		          i0,sten,&Basep,max_comp);
	    	}
	    }
	    break;
	}
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:
	{
	    int	i0, i1, i2;
	    int	i0max, i1max, i2max;
	    int	i0min, i1min, i2min;

	    i0min = imin[0];	i0max = imax[0]; 
	    i1min = imin[1];	i1max = imax[1]; 
	    i2min = imin[2];	i2max = imax[2];

	    for (i2 = i2min; i2 < i2max; ++i2)
	    {
	    	set_static_coords(icoords,&Basep,idirs[2],i2,gr);
	    	for (i1 = i1min; i1 < i1max; ++i1)
	    	{
	    	    set_static_coords(icoords,&Basep,idirs[1],i1,gr);
		    sten->reg_stencil = NO;
		    for (i0 = i0min; i0 < i0max; ++i0) 
		    {
			/*if(i0 == 36 && i1 == 23 && i2 == 5) */
			/*{ */
			/*    add_to_debug("npt_3comp"); */
			/*    printf("#posn b %d  %d  %d \n", i0, i1, i2); */
			/*} */
			update_reg_grid_state(ds,dt,swp_num,iperm,i0,
					      sten,&Basep,max_comp);
			if(debugging("npt_3comp"))
			{
			    printf("#posn e\n");
			    remove_from_debug("npt_3comp");
			}
		    }
		}
	    }
	    break;
	}
#endif /* defined(THREED) */
	}
	debug_print("hyp_npt","Leaving hyp_npt(), dir = %d\n",idirs[0]);
}		/*end hyp_npt*/


#if defined(TWOD) || defined(THREED)
/*
*			set_static_coords():
*/

/*ARGSUSED*/
LOCAL	void set_static_coords(
	int		**icoords,
	POINT		*basep,
	int		idir,
	int		is,
	RECT_GRID	*gr)
{
	int		i;

	for (i = -endpt; i <= endpt; ++i)
	    icoords[i][idir] = is;
	Coords(basep)[idir] = cell_center(is,idir,gr);
}		/*end set_static_coords*/
#endif /* defined(TWOD) || defined(THREED) */


/*
*			interior_sten_reg():
*
*	If called after a vector solver, decides whether the stencil at a given
*	point was regular for the vector solver.  In other words, decides
*	whether solution needs to be recomputed at a given point.
*/

LOCAL	boolean interior_sten_reg(
	int		is,
	int		idir,
	Stencil		*sten)
{
	Wave		*wave = sten->wave;
	Front		*newfr = sten->newfr;
	int 		i;
	int 		**icoords = sten->icoords, icrds[MAXD];
	COMPONENT 	cmp;
	COMPONENT	new_comp = sten->newcomp, *comp = sten->comp;
	CRXING		*cross;

	if (ComponentIsFlowSpecified(new_comp,newfr))
	    return NO;
	for (i = -endpt; i <= endpt; ++i)
	    if ((equivalent_comps(comp[i],new_comp,newfr->interf) == NO) ||
		(sten->nc[i] != 0))
	    return NO;
	if (endpt >= vsten_rad)
	    return YES;

	for (i = 0; i < dim; ++i)
	    icrds[i] = icoords[0][i];
	for (i = endpt+1; i <= vsten_rad; ++i)
	{
	    icrds[idir] = is - i;
	    cmp = Find_rect_comp(idir,icrds,wave);
	    cross = Rect_crossing(icrds,next_side,wave);
	    if ((equivalent_comps(cmp,new_comp,newfr->interf) == NO) ||
		(cross != NULL))
	    	return NO;
		
	    icrds[idir] = is + i;
	    cmp = Find_rect_comp(idir,icrds,wave);
	    cross = Rect_crossing(icrds,prev_side,wave);
	    if ((equivalent_comps(cmp,new_comp,newfr->interf) == NO) ||
		(cross != NULL))
		return NO;
	}
	return YES;
}		/*end interior_sten_reg*/

/*
*			update_reg_grid_state():
*
*	Computes the solution at a given point based on at stencil of
*	surrounding points.
*/

LOCAL	void update_reg_grid_state(
	double		ds,
	double		dt,
	int		swp_num,
	int		*iperm,
	int		is,		/*icoords index in sweep direction */
	Stencil		*sten,
	POINT		*basep,
	COMPONENT	max_comp)
{
	Wave		*wave = sten->wave;
	Wave		*newwave = sten->newwave;
	Front		*fr = sten->fr;
	Front		*newfr = sten->newfr;
	COMPONENT	*comp;
	COMPONENT	new_comp;		/* mid comp wrt newfr */
	HYPER_SURF	*hs;
	int		i, j, index, idir = iperm[swp_num];
	int		**icoords = sten->icoords;

	/*TMP */
/*#define DEBUG_HYP_NPT	 */
#if defined(DEBUG_HYP_NPT)
	int		deb_hyp_npt = NO;

static	char fname[3][11] = {"xhyp_npt()","yhyp_npt()","zhyp_npt()"};
	if (debugging("hyp_npt"))			deb_hyp_npt = YES;
	else if (debugging("xhyp_npt") && idir == 0)	deb_hyp_npt = YES;
	else if (debugging("yhyp_npt") && idir == 1)	deb_hyp_npt = YES;
	else if (debugging("zhyp_npt") && idir == 2)	deb_hyp_npt = YES;
	else						deb_hyp_npt = NO;

#endif /* defined(DEBUG_HYP_NPT) */

	for (i = -endpt; i <= endpt; ++i)
	    icoords[i][idir] = is + i;
	sten->prev_reg_stencil = sten->reg_stencil;

		/* Find components and crosses*/

	new_comp = set_stencil_comps_and_crosses(idir,sten);
	comp = sten->comp;

#if defined(DEBUG_HYP_NPT)
	if (icoords[0][0] == 6 && 
	   	    icoords[0][1] == 35)
	{
	    deb_hyp_npt = YES;
	    add_to_debug("mymuscl");
	}
	if (deb_hyp_npt)
	{
	    (void) printf("\n%s - ",fname[idir]);
	    print_int_vector("icoords[0] = ",icoords[0],dim,"\n");
	    (void) printf("\n comps: ");
	    for (i = -endpt; i <= endpt; ++i)
	    	(void) printf("comp[%d] %d ",i,comp[i]);
	    (void) printf("new_comp %d\n",new_comp);
	    /*TMP */
            sten->prev_reg_stencil = NO;
	    goto jump_is_vec_method;
	}
#endif /* defined(DEBUG_HYP_NPT) */

	if (is_vec_method)
	{
	    if (interior_sten_reg(is,idir,sten))
	    	return;		/*Don't overwrite vector solution*/
	    else
	    	sten->prev_reg_stencil = NO;
	}

/*TMP */
jump_is_vec_method:

	for (j = 0; j < dim; ++j)
	{
	    if (j == idir)
	    {
	    	Coords(sten->p[0])[idir] = cell_center(is,idir,wave->rect_grid);
		for (i = 1; i <= endpt; ++i)
		{
		    Coords(sten->p[-i])[idir] = Coords(sten->p[0])[idir]-i*hdir;
		    Coords(sten->p[i])[idir]  = Coords(sten->p[0])[idir]+i*hdir;
		}
	    }
	    else
	    {
	    	for (i = -endpt; i <= endpt; ++i)
		    Coords(sten->p[i])[j] = Coords(basep)[j];
	    }
	}
	sten->reg_stencil = YES;
	if (RegionIsFlowSpecified(Rect_state(icoords[0],newwave),
			          Rect_state(icoords[0],wave),
			          Coords(sten->p[0]),new_comp,comp[0],fr))
	{
	    sten->reg_stencil = NO;
	    return;
	}

	if ((new_comp > max_comp) && (swp_num > 0))
	{
	    sten->reg_stencil = NO;
	    ft_assign(Rect_state(icoords[0],newwave),
	    	   Rect_state(icoords[0],wave),fr->sizest);
	    return;
	}

	if (new_comp > max_comp)
	{
	    sten->reg_stencil = NO;
	    nearest_intfc_state_and_pt(Coords(sten->p[0]),new_comp,newfr,fr,
				       sten->worksp[0],Coords(sten->p[0]),&hs);
	    ft_assign(Rect_state(icoords[0],newwave),sten->worksp[0],fr->sizest);
	    return;
	}

		/* Find mid state */
	if (equivalent_comps(comp[0],new_comp,newfr->interf) == YES)
	    sten->st[0] = Rect_state(icoords[0],wave);
	else 
	{
	    if (!nearest_crossing_state_with_comp(icoords[0],
	    		Coords(sten->p[0]),new_comp,
			wave_tri_soln(newwave)->tri_grid,&sten->st[0]))
	    {
	    	nearest_intfc_state_and_pt(Coords(sten->p[0]),new_comp,fr,
				       newfr,sten->worksp[0],
				       Coords(sten->p[0]),&hs);
	    	sten->st[0] = sten->worksp[0];
	    }
	    detect_and_load_mix_state(iperm[swp_num],sten,0);
	    sten->reg_stencil = NO;
	}
	ft_assign(left_state(sten->p[0]),sten->st[0],fr->sizest);
	ft_assign(right_state(sten->p[0]),sten->st[0],fr->sizest);

	if(debugging("npt_3comp"))
	    printf("#outer bf\n");
	
	find_outer_states(is,idir,sten,PREV);
	find_outer_states(is,idir,sten,NEXT);
	
	if(debugging("npt_3comp"))
	    printf("#outer af\n");

#if defined(DEBUG_HYP_NPT)
	if (deb_hyp_npt)
	{
	    (void) printf("STENCIL: "); print_Stencil(sten);
	    for (i = -endpt; i <= endpt; ++i)
	    {
	    	(void) printf("state[%d]:	",i);
	    	(*fr->print_state)(sten->st[i]);
	    }
	}
#endif /* defined(DEBUG_HYP_NPT) */

 		/* Check for point source/sink in icoords */

	index = (wave->num_point_sources) ?
	    is_source_block(wave,newfr->interf,new_comp,icoords[0]) : -1;

 		/* Update state */

	/*point_FD */
	npt_solver(ds,dt,Rect_state(icoords[0],newwave),dir,swp_num,iperm,
		   &index,sten,wave);

#if defined(DEBUG_HYP_NPT)
	if (deb_hyp_npt)
	{
	    (void) printf("ANSWER: ");
	    (*fr->print_state)(Rect_state(icoords[0],newwave));
	    remove_from_debug("mymuscl");
	}
#endif /* defined(DEBUG_HYP_NPT) */
}		/*end update_reg_grid_state*/

/*
*			set_stencil_comps_and_crosses():
*
*	This function loads the components of each point and any crosses
*	that occur in the interior of the stencil.  For the crosses,
*	each stencil entry contains the list of crosses lying on its interior
*	side, and these lists are ordered towards the outer points.
*	Thus sten->crx[0] is meaningless, and sten->crx[1] contains the list
*	of crosses (if any) in order FROM point 0 TO point 1.  Periodic
*	boundaries are essentially invisible (cf Find_rect_comp).
*/

LOCAL	COMPONENT set_stencil_comps_and_crosses(
	int		idir,
	Stencil		*sten)
{
	Wave		*wave = sten->wave;
	Wave		*newwave = sten->newwave;
	int		**icoords = sten->icoords;
	TRI_GRID	*tg = wave_tri_soln(wave)->tri_grid;
	int		i, *nc = sten->nc;
	CRXING		***crx = sten->crx;
	COMPONENT	*comp = sten->comp;

	sten->newcomp = Rect_comp(icoords[0],newwave);
	comp[0] = Rect_comp(icoords[0],wave);
	crx[0][0] = NULL;
	nc[0] = 0;
	sten->hs[0] = NULL;/*TODO REMOVE*/

	for (i = 1; i <= endpt; ++i)
	{
	    sten->hs[-i] = sten->hs[i] = NULL;/*TODO REMOVE*/
	    comp[-i] = Find_rect_comp(idir,icoords[-i],wave);
	    nc[-i] = crossings_in_direction(crx[-i],icoords[-i+1],prev_side,tg);
	    comp[i] = Find_rect_comp(idir,icoords[i],wave);
	    nc[i] = crossings_in_direction(crx[i],icoords[i-1],next_side,tg);
	}
	return sten->newcomp;
}		/*set_stencil_comps_and_crosses*/

/*
*			Find_rect_comp():
*
*	Finds the component of a given point using the macro Rect_comp.
*/

EXPORT	COMPONENT	Find_rect_comp(
	int		idir,
	int		*icoords,
	Wave		*wave)
{
	int		is, imax;
	
	is = icoords[idir];
	imax = gmax[idir] + ubuf[idir];

	if ((is < -lbuf[idir]) || (is >= imax))
	    return ext_comp;
	else
	    return Rect_comp(icoords,wave);
		
}		/*Find_rect_comp*/

/*
*			find_outer_states():
*
*	Loads appropriate states for the outer indices of a stencil on either
*	the prev or next sides.  If an interface is crossed, an appropriate
*	state is picked off the interface, and copied out to the end of the
*	stencil (cf copy_states).  Neumann boundaries are also handled
*	appropriately (cf set_neumann_bdry_sten_states).
*/

LOCAL	void find_outer_states(
	int	 is,
	int	 idir,
	Stencil	 *sten,
	CRX_SIDE side)
{
	Wave	   *wave = sten->wave;
	Front	   *fr = sten->fr;
	Front	   *newfr = sten->newfr;
	COMPONENT  *comp = sten->comp, new_comp = sten->newcomp;
	int	   **icoords = sten->icoords;
	int	   i, indx;
	CRXING	   *crx;
	HYPER_SURF *hs;
	CRXING	   *new_crx[MAX_NUM_CRX];

	crx = NULL;
	for  (i = 1; i <= endpt; ++i)
	{
	    indx = i * side;
	    if (equivalent_comps(comp[indx],new_comp,newfr->interf) == YES)
	    {
	        sten->st[indx] = Rect_state(icoords[indx],wave);
	        continue;
	    }
	    if (detect_and_load_mix_state(idir,sten,indx))
	        continue;
	    sten->reg_stencil = NO;
	    crx = sten->crx[indx][0];
	    if (crx != NULL)
	    {
	        /*there is a cross between indx and in_one*/

	        (void) crossings_in_direction(new_crx,icoords[(i-1)*side],
	                               (side == PREV) ? prev_side : next_side,
	                               wave_tri_soln(sten->newwave)->tri_grid);

	        hs = crx->hs;
		if ((new_crx[0] == NULL) ||
		    (wave_type(hs) != wave_type(new_crx[0]->hs)))
		{
		   /* crx has moved out of stencil OR
		    * another curve has moved between indx and the
		    * old crx over the course of the time step. */

	            copy_states(i,side,endpt,new_comp,sten,crx);
	            return;
		}
	        if ((wave_type(hs) == NEUMANN_BOUNDARY) &&
	            (fr->neumann_bdry_state))
	        {
	            set_neumann_bdry_sten_states(i,is,idir,new_comp,
	                                         crx,sten,side);
	            return;
	        }
	        else if (wave_type(hs) == DIRICHLET_BOUNDARY)
	        {
	            set_dirichlet_bdry_sten_states(i,hs,sten,side);
	            return;
	        }
	        else
	        {  
	            /* crx is interior */

	            copy_states(i,side,endpt,new_comp,sten,crx);
	            return;
	        }
	    }
	    else
	    {
		int		   gr_side = ERROR;
	        double		   *coords;
	        static const double OFFSET = 0.01; /* TOLERANCE */

	        coords = Coords(sten->p[indx]);
	        if (is+indx < -lbuf[idir])
		{
		    gr_side = 0;
	            coords[idir] = VL[idir] + hdir*OFFSET;
	        }
		else if (is+indx >= gmax[idir]+ubuf[idir])
		{
		    gr_side = 1;
	            coords[idir] = VU[idir] - hdir*OFFSET;
	        }
	
		if(gr_side != ERROR &&
		   rect_boundary_type(fr->interf, idir, gr_side) == 
			OPEN_BOUNDARY)
		{
		    sten->st[indx] = sten->st[indx-side];
		    copy_states(i,side,endpt,new_comp,sten,crx);
		    return;
		}
		
		sten->st[indx] = sten->worksp[indx];
		if (!nearest_crossing_state_with_comp(icoords[indx],
	    		Coords(sten->p[indx]),new_comp,
			wave_tri_soln(sten->newwave)->tri_grid,
			&sten->st[indx]))
		{
	            nearest_intfc_state_and_pt(coords,new_comp,
	                                   sten->fr,sten->newfr,
					   sten->st[indx],
					   Coords(sten->p[indx]),
	                                   &sten->hs[indx]);
		}
	        
		copy_states(i,side,endpt,new_comp,sten,crx);
	        return;
	    }
	}
}		/*end find_outer_states*/


/*
*			set_neumann_bdry_sten_states():
*
*	Sets stencil states using reflection through a Neumann boundary.
*/

LOCAL	void	set_neumann_bdry_sten_states(
	int	  i,
	int	  is,
	int	  idir,
	COMPONENT new_comp,
	CRXING	  *nbdry,
	Stencil	  *sten,
	CRX_SIDE  side)
{
    	HYPER_SURF *hs = nbdry->hs;
	Wave	   *wave = sten->wave;
	Front	   *fr = sten->fr;
	CRXING	   *rcrx;
	double 	   coords[MAXD], coordsref[MAXD], n[MAXD];
	int	   j, indx, ri;
	size_t	   sizest = fr->sizest;

	for (; i <= endpt; ++i)
	{
	    indx = i * side;
	    if ((*fr->neumann_bdry_state)(Coords(sten->p[indx]),new_comp,
					  nbdry->pt,hs,fr,
					  (POINTER)wave,sten->worksp[indx]))
	    {
	        sten->st[indx] = sten->worksp[indx];
	        for (j = 0; j < dim; ++j)
	    	    Coords(sten->p[indx])[j] = Coords(nbdry->pt)[j];
	        ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
	        ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);
	    }
	    else 
	    {
	        if (sten->crx[indx][0] != nbdry)
	    	    copy_states(i,side,endpt,new_comp,sten,NULL);
	        else
	    	    copy_states(i,side,endpt,new_comp,sten,nbdry);
	        return;
	    }

	    if (is_bdry_hs(hs))
	    {
	        /* TODO:
		 * Currently this section of set_neumann_bdry_sten_states() is
	         * only implemented for rectangular Neumann boundaries
	         * as indicated by the above test.  Interior
	         * or oblique boundaries do not use this boundary
	         * state function. */

	        /* Have set the first state, now need to check for another
	         * cross on the other side of the Neumann boundary. Some
	         * technical issues arise due to the way neumann_bdry_state()
	         * reflects coords.
	         */

	        if (i == endpt)
		    break;
	        ri = (side == PREV) ? -indx - 2*is - 1 :
			              -indx + 2*(gmax[idir] - is) - 1;
	        if (side == PREV)
	        {
	    	    rcrx = (ri >= 0) ? sten->crx[ri+1][0] :
				       sten->crx[ri][sten->nc[ri]-1];
	        }
	        else if (side == NEXT)
	        {
	    	    rcrx = (ri <= 0) ? sten->crx[ri-1][0] :
	    			       sten->crx[ri][sten->nc[ri]-1];
	        }
	        if (rcrx != NULL) 
	        {
	    	    indx = (++i) * side;
	    	    for (j = 0; j < dim; ++j) 
	    	        coords[j] = Coords(rcrx->pt)[j];
	    	    if ((reflect_pt_about_Nbdry(coords,coordsref,n,new_comp,
					        nbdry->hs,fr))
		        && 
		        (*fr->neumann_bdry_state)(coordsref,new_comp,nbdry->pt,
					          hs,fr,(POINTER)wave,
					          sten->worksp[indx]))
		    {
		        sten->st[indx] = sten->worksp[indx];
		        for (j = 0; j < dim; ++j)
		    	    Coords(sten->p[indx])[j] = coordsref[i];
		        ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
		        ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);
		    }
		    copy_states(i,side,endpt,new_comp,sten,NULL);
		    return;
	        }
	    }
	}
}		/*end set_neumann_bdry_sten_states*/

/*
*			set_dirichlet_bdry_sten_states():
*
*	Sets stencil states through a Dirichlet boundary.
*/

/*#bjet2 */
LOCAL	void	set_dirichlet_bdry_sten_states(
	int	   i,
	HYPER_SURF *hs,
	Stencil	   *sten,
	CRX_SIDE   side)
{
	Wave		*wave = sten->wave;
	Front		*fr = sten->fr;
	int		indx; 
	size_t		sizest = fr->sizest;

	for (; i <= endpt; ++i)
	{
	    indx = i * side;
	    sten->st[indx] = sten->worksp[indx];
	    evaluate_dirichlet_boundary_state(Coords(sten->p[indx]),hs,fr,wave,
	    				     sten->st[indx]);
	    ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
	    ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);
	}
}		/*end set_dirichlet_bdry_sten_states*/

/*
*			copy_states():
*
*	We have crossed a non-reflecting curve and are ready to simply copy an
*       appropriate state out to the end of the stencil.  If the cross is not
*	NULL, then it is assumed that the start state has not been set, and
*	an appropriate state is loaded and then copied.
*/

LOCAL	void copy_states(
	int	 start,
	CRX_SIDE side,
	int	 endpt,
	int	 new_comp,
	Stencil	 *sten,
	CRXING	 *cross)
{
	int		i, j, indx, in_one;	
	size_t		sizest = sten->fr->sizest;

	if (cross != NULL)
	{
	    indx = start*side;
	    sten->st[indx] = state_with_comp(cross->pt,cross->hs,new_comp);
	    if (sten->st[indx] != NULL)
	    {
	        for (j = 0; j < dim; ++j)
	    	    Coords(sten->p[indx])[j] = Coords(cross->pt)[j];
	        ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
	        ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);
	    }
	    else
	    {
	        sten->st[indx] = sten->worksp[indx];
	        nearest_intfc_state_and_pt(Coords(sten->p[indx]),new_comp,
					   sten->fr,sten->newfr,sten->st[indx],
					   Coords(sten->p[indx]),
					   &sten->hs[indx]);
	    }
	}

	for (i = start+1; i <= endpt; ++i)
	{
	    indx = i*side;
	    in_one = indx - side;
	    sten->st[indx] = sten->st[in_one];
	    for (j = 0; j < dim; ++j)
	    	Coords(sten->p[indx])[j] = Coords(sten->p[in_one])[j];
	    ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
	    ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);
	}
	return;
}		/*end copy_states*/



#if defined(TWOD)

/*                  h_symmetry_states_printout
*
*       This function print out the front states of two symmetric
*       points w.r.t. the center line on the interface. It's for
*       the symmetric testing of the single mode interface.
*       First we find the two bonds which cross these two symmetric
*       vertical lines, then locate the coordinates of the two points
*       of the intersections.  The left and right states are printed
*       on these two points, in order to check the symmetry on the front.
*/

LOCAL   int	 h_symmetry_states_printout(
        Front           *front,
        Wave		*wave)
{
        int             i, k, dim = wave->rect_grid->dim;
        int             wave_type;
        int             imin[MAXD], imax[MAXD];
        HYPER_SURF      *hs;
        CURVE           *c;
        BOND            *b;
        double           lp[2], rp[2], s0, s1, s2, s3;
        double           dh[MAXD], sgn_x_l, sgn_x_r;
        Locstate        sl_l, sl_r, sr_l, sr_r;
        INTERFACE       *interf = front->interf;
        COMPONENT       pcomp = positive_component(hs);
        COMPONENT       ncomp = negative_component(hs);

        alloc_state(intfc,&sl_l,front->sizest);
        alloc_state(intfc,&sl_r,front->sizest);
        alloc_state(intfc,&sr_l,front->sizest);
        alloc_state(intfc,&sr_r,front->sizest);

        for (i = 0; i < dim; i++)
        {
            dh[i] = wave->rect_grid->h[i];
            imin[i] = 0;        imax[i] = wave->rect_grid->gmax[i];
            if (wave->rect_grid->lbuf[i] == 0) imin[i] += 1;
            if (wave->rect_grid->ubuf[i] == 0) imax[i] -= 1;
        }

        printf("dh[0] = %lf, dh[1] = %lf\n",dh[0],dh[1]);

        if (make_bond_comp_lists(interf) == FUNCTION_FAILED)
        {
                (void) printf("WARNING in make_bond_lists(), ");
                (void) printf("make_bond_listsd() failed\n");
                return FUNCTION_FAILED;
        }

        k = 0;      

        lp[0] = ((imin[0] + imax[0]) / 4) * (dh[0]);

        rp[0] = (((imin[0] + imax[0]) * 3 / 4) - 1) * (dh[0]);


        (void) next_bond(interf,NULL,NULL);

        while (next_bond(interf,&b,&c))
        {
                if(wave_type(c) < FIRST_SCALAR_PHYSICS_WAVE_TYPE) continue;

	        hs = Hyper_surf(c);
	       	pcomp = positive_component(hs);
         	ncomp = negative_component(hs);

                s0 = (Coords(b->end)[0] - lp[0]);
                s1 = (Coords(b->start)[0] - lp[0]);
                s2 = (Coords(b->end)[0] - rp[0]);
                s3 = (Coords(b->start)[0] - rp[0]);

                sgn_x_l = s0 * s1;
                sgn_x_r = s2 * s3;

                if(sgn_x_l > 0.0 && sgn_x_r > 0.0) continue;

                if(s0 == 0.0)
                {
                  lp[1] = Coords(b->end)[1];

                  k = k + 1;
                }

                if(sgn_x_l < 0.0)
                {
                  lp[1] = ((Coords(b->end)[1] - Coords(b->start)[1]) *
                  (Coords(b->start)[0] - lp[0]) /
                  (Coords(b->end)[0] - Coords(b->start)[0])) +
                   Coords(b->start)[1];

                  k = k + 1;
                }

                if(s2 == 0.0)
                {
                  rp[1] = Coords(b->end)[1];

                  k = k + 1;
                }

                if(sgn_x_r < 0.0)
                {
                  rp[1] = ((Coords(b->end)[1] - Coords(b->start)[1]) *
                  (Coords(b->start)[0] - rp[0]) /
                  (Coords(b->end)[0] - Coords(b->start)[0])) +
                   Coords(b->start)[1];

                  k = k + 1;
                }

        }

        printf("left point = (%lf,%lf)\n",lp[0],lp[1]);

        printf("right point = (%lf,%lf)\n",rp[0],rp[1]);

        printf("k = %d\n",k);

        hyp_solution(lp,ncomp,hs,NEGATIVE_SIDE,front,wave,sl_l,NULL);

        hyp_solution(lp,pcomp,hs,POSITIVE_SIDE,front,wave,sl_r,NULL);

        hyp_solution(rp,ncomp,hs,NEGATIVE_SIDE,front,wave,sr_l,NULL);

        hyp_solution(rp,pcomp,hs,POSITIVE_SIDE,front,wave,sr_r,NULL);

        (*front->print_state)(sl_l);
        (*front->print_state)(sl_r);
        (*front->print_state)(sr_l);
        (*front->print_state)(sr_r);

}

#endif /* defined(TWOD) */



