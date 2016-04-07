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
*				hscatter.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains function for copying states across distinct processors
*/

#define DEBUG_STRING	"hscatter"
#include <hyp/hdecs.h>

#define state_id(i)     (STATE_ID + (i+1))
LOCAL size_t BLOCK_SIZE = 0; /*TOLERANCE - TODO: what is a good value*/

	/* LOCAL Function Declarations */
LOCAL	POINTER	buffer_storage(size_t,int*,int*,int,POINTER,size_t*,size_t*);
LOCAL	void	copy_states_across_domain(int*,int,int,Front*,Wave*);
LOCAL	void	pp_receive_interior_states(int*,int,int*,int,Front*,Wave*);
LOCAL	void	pp_send_interior_states(int*,int,int*,int,Front*,Wave*);
LOCAL	void	set_receive_domain(int*,int*,int*,int,int,RECT_GRID*);
LOCAL	void	set_send_domain(int*,int*,int*,int,int,RECT_GRID*);

EXPORT	void	SetHypPPBlockSize(
	size_t	new_BLOCK_SIZE)
{
	BLOCK_SIZE = new_BLOCK_SIZE;
	EnsureSufficientMessageBufferSize(BLOCK_SIZE);
}		/*end SetHypPPBlockSize*/

EXPORT	void	SetDefaultHypPPBlockSize(
	int	dim)
{
	if (BLOCK_SIZE != 0) /*BLOCK_SIZE explicitly set on command line*/
	    return;

	switch (dim)
	{
	case 1:
	case 2:
	    BLOCK_SIZE = 9992;
	    break;
	case 3:
	    BLOCK_SIZE = 1000000;
	    break;
	default:
	    screen("ERROR in SetDefaultHypPPBlockSize(), "
		   "invalid dimension %d\n",dim);
	    clean_up(ERROR);
	}
	EnsureSufficientMessageBufferSize(BLOCK_SIZE);
}		/*end SetDefaultHypPPBlockSize*/

/*
*			h_scatter_states():
*
*	Main control function for parallel communication of rectangular
*	lattice states.  The states in the buffer zone regions are shipped
*	to adjacent subdomains in this function.  This is conducted in
*	dim*2 passes.  As state information in each direction is transmitted,
*	the domain of validity is increased.  This removes the need to
*	transmit corner information.  The basic scheme for two dimensions
*	and iperm = [0 1],  is outlined below.
*
*	gr = front->rect_grid
*	Vl = -gr->lbuf	Vu = gr->gmax + gr->ubuf
*	l = 0,  u = gr->gmax
*
*	i = 0
*
*	                 l[0]+lbuf[0]     u[0]-ubuf[0]
*                            |                |
*			    \|/              \|/
*	  l[1] |------|------------------------------|------|
*	       |      |      |                |      |      |
*	       |  R   |      |                |      |  R   |
*	       |  E   |      |                |      |  E   |
*	       |  C   |  S   |                |  S   |  C   |
*	       |  E   |  E   |                |  E   |  E   |
*	       |  I   |  N   |                |  N   |  I   |
*	       |  V   |  D   |                |  D   |  V   |
*	       |  E   |      |                |      |  E   |
*	       |      |      |                |      |      |
*	       |      |      |                |      |      |
*	       |      |      |                |      |      |
*	   l[1]|------|------------------------------|------|
*         Vl[0]      l[0]                           u[0]  Vu[0]
*
*	i = 1
*
*         Vu[1]----------------------------------------------
*	       |                                            |
*	       |             RECEIVE                        |
*	  l[1] |--------------------------------------------|
*	       |                                            |
*	       |               SEND                         |
*	       |--------------------------------------------|
*	       |                                            |
*	       |                                            |
*	       |                                            |
*	       |                                            |
*	       |                                            |
*	       |                                            |
*	       |--------------------------------------------|
*	       |                                            |
*	       |               SEND                         |
*	   l[1]|--------------------------------------------|
*	       |                                            |
*	       |             RECEIVE                        |
*         Vl[1]----------------------------------------------
*         Vl[0]                                           Vu[0]
*/

EXPORT 	boolean	h_scatter_states(
	Wave		*wave,
	Front		*front,
	int             *iperm,
	int		swp)
{
	PP_GRID		*pp_grid = wave->pp_grid;
	int		myid;
	int		me[MAXD];
	int		i,side, dim;

	DEBUG_ENTER(h_scatter_states)

	if (wave->sizest == 0)
	{
	    DEBUG_LEAVE(h_scatter_states)
	    return FUNCTION_SUCCEEDED;
	}

	start_clock("h_scatter_states");

	myid = pp_mynode();
	dim = wave->rect_grid->dim;

	find_Cartesian_coordinates(myid,pp_grid,me);

	for (side = 0; side < 2; ++side)
	{
	    pp_gsync();
	    pp_send_interior_states(me,swp,iperm,side,front,wave);
	    pp_receive_interior_states(me,swp,iperm,(side+1)%2,front,wave);
	}

	stop_clock("h_scatter_states");
	DEBUG_LEAVE(h_scatter_states)
	return FUNCTION_SUCCEEDED;
}		/*end h_scatter_states*/

/*
*			pp_send_interior_states():
*
*	Sends state information in a single buffer region.
*
*/

LOCAL	void	pp_send_interior_states(
	int		*me,
	int		swp,
	int		*iperm,
	int		side,
	Front		*front,
	Wave		*wave)
{
	INTERFACE     *intfc = front->interf;
	PP_GRID	      *pp_grid = front->pp_grid;
        RECT_GRID     *gr = front->rect_grid;
        byte	      *ps;
	int	      L[MAXD], U[MAXD];
	int	      i;
	int	      him[MAXD];
	int	      myid, dst_id;
        int	      dim = gr->dim;
	size_t        len;
	static byte   *storage = NULL;
	static size_t alloc_len = 0;

	DEBUG_ENTER(pp_send_interior_states)

	if (rect_boundary_type(intfc,iperm[swp],side) == REFLECTION_BOUNDARY)
	{
	    reflect_states_across_domain(iperm,side,swp,front,wave);
	    DEBUG_LEAVE(pp_send_interior_states)
	    return;
	}
	if (rect_boundary_type(intfc,iperm[swp],side) != SUBDOMAIN_BOUNDARY)
	{
	    DEBUG_LEAVE(pp_send_interior_states)
	    return;
	}

	myid = pp_mynode();
	dst_id = neighbor_id(him,me,iperm[swp],side,pp_grid);
	if (myid == dst_id)
	{
	    copy_states_across_domain(iperm,side,swp,front,wave);
	    DEBUG_LEAVE(pp_send_interior_states)
	    return;
	}

	set_send_domain(L,U,iperm,side,swp,gr);
	storage = (byte*) buffer_storage(front->sizest,L,U,dim,
					 (POINTER)storage,&len,&alloc_len);
	(*wave->bundle_states)(L,U,wave,storage);
	start_clock("pp_send_states");
	for (ps = storage, i = 0; len >= BLOCK_SIZE; 
				len -= BLOCK_SIZE, ps += BLOCK_SIZE, ++i)
	{
	    pp_send(state_id(i),(POINTER)ps,BLOCK_SIZE,dst_id);
	}
	if (len != 0)
	    pp_send(state_id(i),(POINTER)ps,len,dst_id);
	stop_clock("pp_send_states");
	DEBUG_LEAVE(pp_send_interior_states)
}		/*end pp_send_interior_states*/

/*
*			pp_receive_interior_states():
*
*	Receives state information in a single buffer region.
*
*/

LOCAL	void	pp_receive_interior_states(
	int		*me,
	int		swp,
	int		*iperm,
	int		side,
	Front		*front,
	Wave		*wave)
{
	INTERFACE     *intfc = front->interf;
	PP_GRID	      *pp_grid = front->pp_grid;
        RECT_GRID     *gr = front->rect_grid;
        byte	      *ps;
	int	      L[MAXD], U[MAXD];
	int	      him[MAXD];
	int	      myid, src_id;
        int	      dim = gr->dim;
	int	      i;
	size_t        len;
	static byte   *storage = NULL;
	static size_t alloc_len = 0;

	DEBUG_ENTER(pp_receive_interior_states)

	if (rect_boundary_type(intfc,iperm[swp],side) != SUBDOMAIN_BOUNDARY)
	{
	    DEBUG_LEAVE(pp_receive_interior_states)
	    return;
	}

	myid = pp_mynode();
	src_id = neighbor_id(him,me,iperm[swp],side,pp_grid);
	if (myid == src_id)
	{
	    DEBUG_LEAVE(pp_receive_interior_states)
	    return; /* Already done */
	}

	set_receive_domain(L,U,iperm,side,swp,gr);
	storage = (byte*) buffer_storage(front->sizest,L,U,dim,
					 (POINTER)storage,&len,&alloc_len);
	start_clock("pp_receive_states");
	for (ps = storage, i = 0; len >= BLOCK_SIZE;
				len -= BLOCK_SIZE, ps += BLOCK_SIZE,++i)
	{
	    pp_recv(state_id(i),src_id,(POINTER)ps,BLOCK_SIZE);
	}
	if (len != 0)
	    pp_recv(state_id(i),src_id,(POINTER)ps,len);
	stop_clock("pp_receive_states");
	(*wave->unbundle_states)(L,U,wave,storage);
	DEBUG_LEAVE(pp_receive_interior_states)
}		/*end pp_receive_interior_states*/

/*
*		set_send_domain():
*
*	Sets the domain limits for a buffer zone to be transmitted.
*
*/

LOCAL	void	set_send_domain(
	int		*L,
	int		*U,
	int		*iperm,
	int		side,
	int		swp,
	RECT_GRID	*gr)
{
        int		dim = gr->dim;
        int		*lbuf = gr->lbuf;
        int		*ubuf = gr->ubuf;
	int		*gmax = gr->gmax;
	int		j;

	DEBUG_ENTER(set_send_domain)
	for (j = 0; j < swp; ++j)
	{
	    L[iperm[j]] = -lbuf[iperm[j]];
	    U[iperm[j]] = gmax[iperm[j]] + ubuf[iperm[j]];
	}
	if (side == 0)
	{
	    L[iperm[swp]] =  0;
	    U[iperm[swp]] = lbuf[iperm[swp]];
	}
	else
	{
	    L[iperm[swp]] = gmax[iperm[swp]] - ubuf[iperm[swp]];
	    U[iperm[swp]] = gmax[iperm[swp]];
	}
	for (j = swp+1; j < dim; ++j)
	{
	    L[iperm[j]] = -lbuf[iperm[j]];
	    U[iperm[j]] = gmax[iperm[j]] + ubuf[iperm[j]];
	}
	DEBUG_LEAVE(set_send_domain)
}		/*end set_send_domain*/

/*
*		set_receive_domain():
*
*	Sets the domain limits for a buffer zone to be received.
*
*/

LOCAL	void	set_receive_domain(
	int		*L,
	int		*U,
	int		*iperm,
	int		side,
	int		swp,
	RECT_GRID	*gr)
{
        int		dim = gr->dim;
        int		*lbuf = gr->lbuf;
        int		*ubuf = gr->ubuf;
	int		*gmax = gr->gmax;
	int		j;

	DEBUG_ENTER(set_receive_domain)
	for (j = 0; j < swp; ++j)
	{
	    L[iperm[j]] = -lbuf[iperm[j]];
	    U[iperm[j]] = gmax[iperm[j]] + ubuf[iperm[j]];
	}
	if (side == 0)
	{
	    L[iperm[swp]] = -lbuf[iperm[swp]];
	    U[iperm[swp]] = 0;
	}
	else
	{
	    L[iperm[swp]] = gmax[iperm[swp]];
	    U[iperm[swp]] = gmax[iperm[swp]] + ubuf[iperm[swp]];
	}
	for (j = swp+1; j < dim; ++j)
	{
	    L[iperm[j]] = -lbuf[iperm[j]];
	    U[iperm[j]] = gmax[iperm[j]] + ubuf[iperm[j]];
	}
	if (DEBUG)
	{
	    (void) printf("swp = %d, side = %d, ",swp,side);
	    print_int_vector("iperm = ",iperm,dim,"\n");
	    print_int_vector("L = ",L,dim,", ");
	    print_int_vector("U = ",U,dim,"\n");
	}
	DEBUG_LEAVE(set_receive_domain)
}		/*end set_receive_domain*/

/*
*			buffer_storage():
*
*	Allocates storage buffer for buffer zone communication.
*/

LOCAL	POINTER buffer_storage(
	size_t  sizest,
	int     *L,
	int     *U,
	int     dim,
	POINTER storage,
	size_t  *plen,
	size_t  *palloc_len)
{
	int    i;
	size_t len;

	DEBUG_ENTER(buffer_storage)

	for (i = 0, len = 1; i < dim; ++i)
	    len *= U[i] - L[i];
	len *= sizest;
	*plen = len;

	if (len > *palloc_len)
	{
	    if (storage != NULL)
		free(storage);
	    scalar(&storage,len);
	    *palloc_len = len;
	}
	DEBUG_LEAVE(buffer_storage)
	return storage;
}		/*end buffer_storage*/

/*
*			copy_states_across_domain():
*
*	This special purpose function is provided to improve performance
*	in the case where the sending and receiving processors are the
*	same.  It copies the state in one region of the rectangular
*	lattice into another region.
*/

LOCAL	void	copy_states_across_domain(
	int		*iperm,
	int		side,
	int		swp,
	Front		*front,
	Wave		*wave)
{
        RECT_GRID	*gr = front->rect_grid;
	int		src_ic[MAXD], dst_ic[MAXD];
	int		src_L[MAXD], src_U[MAXD];
	int		dst_L[MAXD], dst_U[MAXD];
	size_t		sizest = front->sizest;

	DEBUG_ENTER(copy_states_across_domain)
	set_send_domain(src_L,src_U,iperm,side,swp,gr);
	set_receive_domain(dst_L,dst_U,iperm,(side+1)%2,swp,gr);

	switch (gr->dim)
	{
#if defined(ONED)
	case 1:
		{
		    int src_ix, dst_ix;
		    for (src_ix = src_L[0], dst_ix = dst_L[0];
					src_ix < src_U[0]; ++src_ix, ++dst_ix)
		    {
		        src_ic[0] = src_ix;
		        dst_ic[0] = dst_ix;
		        ft_assign(Rect_state(dst_ic,wave),
				Rect_state(src_ic,wave),sizest);
		    }
		}
		break;
#endif /* defined(ONED) */
#if defined(TWOD)
	case 2:
		{
		    int src_ix, dst_ix, src_iy, dst_iy;
		    for (src_iy = src_L[1], dst_iy = dst_L[1];
					src_iy < src_U[1]; ++src_iy, ++dst_iy)
		    {
		        src_ic[1] = src_iy;
		        dst_ic[1] = dst_iy;
		        for (src_ix = src_L[0], dst_ix = dst_L[0];
					src_ix < src_U[0]; ++src_ix, ++dst_ix)
		        {
		            src_ic[0] = src_ix;
		            dst_ic[0] = dst_ix;
		            ft_assign(Rect_state(dst_ic,wave),
				Rect_state(src_ic,wave),sizest);
		        }
		    }
		}
		break;
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:
		{
		    int src_ix, dst_ix, src_iy, dst_iy, src_iz, dst_iz;
		    for (src_iz = src_L[2], dst_iz = dst_L[2];
					src_iz < src_U[2]; ++src_iz, ++dst_iz)
		    {
		        src_ic[2] = src_iz;
		        dst_ic[2] = dst_iz;
		        for (src_iy = src_L[1], dst_iy = dst_L[1];
					src_iy < src_U[1]; ++src_iy, ++dst_iy)
		        {
		            src_ic[1] = src_iy;
		            dst_ic[1] = dst_iy;
		            for (src_ix = src_L[0], dst_ix = dst_L[0];
					src_ix < src_U[0]; ++src_ix, ++dst_ix)
		            {
		                src_ic[0] = src_ix;
		                dst_ic[0] = dst_ix;
		                ft_assign(Rect_state(dst_ic,wave),
				    Rect_state(src_ic,wave),sizest);
		            }
		        }
		    }	
		}
		break;
#endif /* defined(THREED) */
	}
	DEBUG_LEAVE(copy_states_across_domain)
}		/*end copy_states_across_domain*/

/*
*		reflect_states_across_domain():
*
*	Copies the rect state storage across a reflection boundary.
*	States locations are reflected with respect to the reflection
*	boundary and the corresponding states are reflected.
*/

EXPORT	void	reflect_states_across_domain(
	int		*iperm,
	int		side,
	int		swp,
	Front		*front,
	Wave		*wave)
{
        RECT_GRID	*gr = front->rect_grid;
	Locstate	dst_st, src_st;
	double		*n;
	double		*dst_crds, *src_crds, p[MAXD];
	int		i, dim = gr->dim;
	int		dir = iperm[swp];
	int		src_ic[MAXD], dst_ic[MAXD];
	int		src_L[MAXD], src_U[MAXD];
	int		dst_L[MAXD], dst_U[MAXD];
	size_t		sizest = front->sizest;
	static double	nors[] = {  1.0,  0.0,  0.0,
				    0.0,  1.0,  0.0,
				    0.0,  0.0,  1.0,
				   -1.0,  0.0,  0.0,
				    0.0, -1.0,  0.0,
				    0.0,  0.0, -1.0};

	DEBUG_ENTER(reflect_states_across_domain)
	set_send_domain(src_L,src_U,iperm,side,swp,gr);
	set_receive_domain(dst_L,dst_U,iperm,side,swp,gr);

#define refl_ic(ic,dir,sL,dU)	((dU)[dir] - 1 + ((sL)[dir] - (ic)[dir]))

	n = nors+3*dir+9*side;
	switch (dim)
	{
#if defined(ONED)
	case 1:
	    {
	        int src_ix;
	        for (src_ix = src_L[0]; src_ix < src_U[0]; ++src_ix)
	        {
	            src_ic[0] = src_ix;
	    	    dst_ic[0] = refl_ic(src_ic,0,src_L,dst_U);
	    	    dst_st = Rect_state(dst_ic,wave);
	    	    src_st = Rect_state(src_ic,wave);
	    	    dst_crds = Rect_coords(dst_ic,wave);
	    	    src_crds = Rect_coords(src_ic,wave);
	    	    for (i = 0; i < dim; ++i)
	    		p[i] = 0.5*(dst_crds[i]+src_crds[i]);
	            ft_assign(dst_st,src_st,sizest);
	    	    reflect_state(dst_st,front->interf,src_crds,p,n);
	        }
	    }
	    break;
#endif /* defined(ONED) */
#if defined(TWOD)
	case 2:
	    {
	        int src_ix, dst_ix, src_iy, dst_iy;
	        for (src_iy = src_L[1], dst_iy = dst_L[1];
	    			src_iy < src_U[1]; ++src_iy, ++dst_iy)
	        {
	            src_ic[1] = src_iy;
	            dst_ic[1] = dst_iy;
	            for (src_ix = src_L[0], dst_ix = dst_L[0];
	    			src_ix < src_U[0]; ++src_ix, ++dst_ix)
	            {
	                src_ic[0] = src_ix;
	                dst_ic[0] = dst_ix;
	    	        dst_ic[dir] = refl_ic(src_ic,dir,src_L,dst_U);
	    	        dst_st = Rect_state(dst_ic,wave);
	    	        src_st = Rect_state(src_ic,wave);
	    	        dst_crds = Rect_coords(dst_ic,wave);
	    	        src_crds = Rect_coords(src_ic,wave);
	    	        for (i = 0; i < dim; ++i)
	    		    p[i] = 0.5*(dst_crds[i]+src_crds[i]);
	                ft_assign(dst_st,src_st,sizest);
	    	        reflect_state(dst_st,front->interf,src_crds,p,n);
	            }
	        }
	    }
	    break;
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:
	    {
	        int src_ix, dst_ix, src_iy, dst_iy, src_iz, dst_iz;
	        for (src_iz = src_L[2], dst_iz = dst_L[2];
	    			src_iz < src_U[2]; ++src_iz, ++dst_iz)
	        {
	            src_ic[2] = src_iz;
	            dst_ic[2] = dst_iz;
	            for (src_iy = src_L[1], dst_iy = dst_L[1];
	    			src_iy < src_U[1]; ++src_iy, ++dst_iy)
	            {
	                src_ic[1] = src_iy;
	                dst_ic[1] = dst_iy;
	                for (src_ix = src_L[0], dst_ix = dst_L[0];
	    			src_ix < src_U[0]; ++src_ix, ++dst_ix)
	                {
	                    src_ic[0] = src_ix;
	                    dst_ic[0] = dst_ix;
	    	            dst_ic[dir] = refl_ic(src_ic,dir,src_L,dst_U);
	    	            dst_st = Rect_state(dst_ic,wave);
	    	            src_st = Rect_state(src_ic,wave);
	    	            dst_crds = Rect_coords(dst_ic,wave);
	    	            src_crds = Rect_coords(src_ic,wave);
	                    for (i = 0; i < dim; ++i)
	        	        p[i] = 0.5*(dst_crds[i]+src_crds[i]);
	                    ft_assign(dst_st,src_st,sizest);
	    	            reflect_state(dst_st,front->interf,src_crds,p,n);
	                }
	            }
	        }	
	    }
	    break;
#endif /* defined(THREED) */
	}
#undef refl_ic
	DEBUG_LEAVE(reflect_states_across_domain)
}		/*end reflect_states_across_domain*/

#if defined(__MPI__)
EXPORT 	boolean 	h_iscatter_states(
	Wave		*wave,
	Front		*front,
	int             *iperm,
        int             swp)
{
	INTERFACE     *intfc = front->interf;
	MPI_Request    ireq[2], oreq[2];
	PP_GRID	       *pp_grid = wave->pp_grid;
        RECT_GRID      *gr = front->rect_grid;
	boolean           received[2], sent[2];
	int	       iL[2][3], iU[2][3];
	int	       oL[2][3], oU[2][3];
	int            itag[2], otag[2];
	int	       myid;
	int	       me[3], him[3];
	int	       i, j, dim;
	int            src_id, dst_id;
	size_t         ilen[2], olen[2];
        static POINTER ibuf[2], obuf[2];
	static size_t  alloc_ilen[2], alloc_olen[2];

	DEBUG_ENTER(h_iscatter_states)

	if (wave->sizest == 0)
	{
	    DEBUG_LEAVE(h_iscatter_states)
	    return YES;
	}

	start_clock("h_iscatter_states");

	myid = pp_mynode();
	dim = wave->rect_grid->dim;

	find_Cartesian_coordinates(myid,pp_grid,me);

	if (DEBUG)
	{
	    (void) printf("step = %d, ",front->step);
	    print_int_vector("iperm = ",iperm,dim,", ");
	    print_int_vector("me = ",me,dim,"\n");
	}


	/* Set domain sizes */
	i = iperm[swp];
	for (j = 0; j < 2; ++j)
	{
	    itag[j] = state_id(2*i+j);
	    otag[j] = state_id(2*i+(j+1)%2);
	    if (rect_boundary_type(intfc,i,j) == REFLECTION_BOUNDARY)
	    {
	        reflect_states_across_domain(iperm,j,swp,front,wave);
		received[j] = YES;
		sent[j] = YES;
	    }
	    else if (rect_boundary_type(intfc,i,j) != SUBDOMAIN_BOUNDARY)
	    {
		received[j] = YES;
		sent[j] = YES;
	    }
	    else
	    {
		received[j] = NO;
	        set_receive_domain(iL[j],iU[j],iperm,j,swp,gr);
		sent[j] = NO;
	        set_send_domain(oL[j],oU[j],iperm,j,swp,gr);
	        ibuf[j] = buffer_storage(front->sizest,iL[j],iU[j],dim,
					     ibuf[j],ilen+j,alloc_ilen+j);
	        obuf[j] = buffer_storage(front->sizest,oL[j],oU[j],dim,
					     obuf[j],olen+j,alloc_olen+j);
	        (*wave->bundle_states)(oL[j],oU[j],wave,(byte*)obuf[j]);
	        src_id = neighbor_id(him,me,i,j,pp_grid);
		if (DEBUG)
		{
			(void) printf("Receiving data in direction %d side %d "
				      "from processor %d with tag %d (%d)\n",
				      i,j,src_id,itag[j],itag[j]-state_id(0));
			(void) printf("iL[%d] = ",j);
			print_int_vector("",iL[j],dim,", ");
			(void) printf("iU[%d] = ",j);
			print_int_vector("",iU[j],dim,"\n");
			(void) printf("ilen[%d] = %d, %d states\n",j,
				      (int)ilen[j],
				      (int)(ilen[j]/front->sizest));
		}
	        pp_irecv(itag[j],src_id,ibuf[j],ilen[j],ireq+j);
	    }
	}

	pp_gsync();

	/* Submit nonblocking sends */
	for (j = 0; j < 2; ++j)
	{
	    if (sent[j] == NO)
	    {
	        dst_id = neighbor_id(him,me,i,j,pp_grid);
		if (DEBUG)
		{
		    (void) printf("Sending data in direction %d side %d "
				      "to processor %d with tag %d (%d)\n",
				      i,j,dst_id,otag[j],otag[j]-state_id(0));
		    (void) printf("oL[%d] = ",j);
		    print_int_vector("",oL[j],dim,", ");
		    (void) printf("oU[%d] = ",j);
		    print_int_vector("",oU[j],dim,"\n");
		    (void) printf("olen[%d] = %d, %d states\n",
				      j,(int)olen[j],
				      (int)(olen[j]/front->sizest));
		}
	        pp_isend(otag[j],obuf[j],olen[j],dst_id,oreq+j);
	    }
	}

	/* Process received state data */
	for (j = 0; j < 2; ++j)
	{
	    if (received[j] == NO)
	    {
		pp_wait(ireq+j);
		received[j] = YES;
	        (*wave->unbundle_states)(iL[j],iU[j],wave,(byte*)ibuf[j]);
	    }
	}

	/* Wait for sends to complete */
	for (j = 0; j < 2; ++j)
	{
	    if (sent[j] == NO)
	    {
		pp_wait(oreq+j);
		sent[j] = YES;
	    }
	}
	pp_gsync();

	stop_clock("h_iscatter_states");
	DEBUG_LEAVE(h_iscatter_states)
	return YES;
}		/*end h_iscatter_states*/
#endif /* defined(__MPI__) */

EXPORT void pp_send_large_data(
        byte          *storage,
        size_t        len,
        int           dst_id)
{
        byte          *ps;
        int           i;

        DEBUG_ENTER(pp_send_large_data)

        for (ps = storage, i = 0; len >= BLOCK_SIZE;
                   len -= BLOCK_SIZE, ps += BLOCK_SIZE, ++i)
        {
            pp_send(state_id(i),(POINTER)ps,BLOCK_SIZE,dst_id);
        }
        if (len != 0)
            pp_send(state_id(i),(POINTER)ps,len,dst_id);

        DEBUG_LEAVE(pp_send_large_data)
}


EXPORT void pp_receive_large_data(
        byte          *storage,
        size_t        len,
        int           src_id)
{
        byte          *ps;
        int           i;

        DEBUG_ENTER(pp_send_large_data)

        for (ps = storage, i = 0; len >= BLOCK_SIZE;
                   len -= BLOCK_SIZE, ps += BLOCK_SIZE, ++i)
        {
            pp_recv(state_id(i),src_id,(POINTER)ps,BLOCK_SIZE);
        }
        if (len != 0)
            pp_recv(state_id(i),src_id,(POINTER)ps,len);

        DEBUG_LEAVE(pp_send_large_data)
}





















