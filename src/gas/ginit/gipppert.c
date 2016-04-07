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
*				gipppert.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Routines parallel computer support for initializing and analyzing
*	normal modes for perturbation problems.
*
*	Multiple mode interfaces (random surfaces) are also treated.
*/

#if defined(TWOD) || defined(THREED)

#include <ginit/ginit.h>

	/* LOCAL Function Declarations */
LOCAL	void	recv_nm(NORMAL_MODE*,LIN_PERT*,int);
LOCAL	void	recv_pts(LIN_PERT*,int);
LOCAL	void	send_nm(NORMAL_MODE*,LIN_PERT*,int,int);
LOCAL	void	send_pts(LIN_PERT*,int);

#if !defined(MASTER)
enum { MASTER = 0 };
#endif /* !defined(MASTER) */

enum {
	ND_ID1	= USER_MIN_MESG_ID + 1000,
	ND_ID2	= USER_MIN_MESG_ID + 2000,
	M_NUM	= USER_MIN_MESG_ID + 3000,
	UN_M	= USER_MIN_MESG_ID + 4000,
	FLAG	= USER_MIN_MESG_ID + 5000,
	N_M	= USER_MIN_MESG_ID + 6000,
	PTS	= USER_MIN_MESG_ID + 7000
};



EXPORT	int pp_list_modes(
	NORMAL_MODE	**normal_mode,
	LIN_PERT	*lin_pert,
	int		n_modes,
	int		flag,
	int		unstable_modes)
{
	static	int	tot_modes = 0;
	int		node_id;
	int		n_nodes;
	NORMAL_MODE	*n_m = lin_pert->normal_mode;

	debug_print("lin_pert","Entered pp_list_modes()\n");

	n_nodes = pp_numnodes();

	if (n_nodes <= 2 || n_modes == 1)
	{
		if (flag == NO)
		{
			free(n_m);
			lin_pert->normal_mode = NULL;
			debug_print("lin_pert","Left pp_list_modes()\n");
			return	unstable_modes;
		}
		normal_mode[unstable_modes] = n_m; 
		debug_print("lin_pert","Left pp_list_modes()\n");
		return	++unstable_modes;
	}

	node_id = pp_mynode();
	if (node_id == MASTER)
	{
	    if (tot_modes >= n_modes || !pp_iprobe_any(ND_ID2))
	    {
	    	debug_print("lin_pert","Left pp_list_modes()\n");
	    	return	unstable_modes;
	    }
	    ++tot_modes;
	    pp_recv_any(ND_ID2,(POINTER)&node_id,INT);
	    pp_send(UN_M,(POINTER)&unstable_modes,INT,node_id);
	    pp_recv(FLAG,node_id,(POINTER)&flag,INT);
	    if (flag == NO)
	    {
	    	debug_print("lin_pert","Left pp_list_modes()\n");
	    	return	unstable_modes;
	    }
	    if (unstable_modes == 0)
	    {
	        if (debugging("lin_pert"))
	            (void) printf("\nMaster is ready to receive pts. ");
	        recv_pts(lin_pert,node_id);
	        if (debugging("lin_pert"))
	            (void) printf("\nMaster has already received pts.\n");
	    }
	    if (debugging("lin_pert"))
	        (void) printf("\nMaster is ready to receive normal_mode. ");
	    recv_nm(normal_mode[unstable_modes],lin_pert,node_id);
	    if (debugging("lin_pert"))
	        (void) printf("\nMaster has already received normal_mode.\n");
	}
	else
	{
	    pp_send(ND_ID2,(POINTER)&node_id,INT,MASTER);
	    pp_recv(UN_M,MASTER,(POINTER)&unstable_modes,INT);
	    pp_send(FLAG,(POINTER)&flag,INT,MASTER);
	    if (flag == NO)
	    {
	    	debug_print("lin_pert","Left pp_list_modes()\n");
	    	return	unstable_modes;
	    }
	    if (unstable_modes == 0)
	    {
	    	if (debugging("lin_pert"))
	    	    (void) printf("\nNode %d is ready to send pts. ",node_id);
	    	send_pts(lin_pert,MASTER);
	    	if (debugging("lin_pert"))
	    	    (void) printf("\nNode %d has already sent pts.\n",node_id);
	    }
	    if (debugging("lin_pert"))
	    	(void) printf("\nNode %d is ready to send normal_mode.\n",
			      node_id);
	    send_nm(n_m,lin_pert,MASTER,YES);
	    if (debugging("lin_pert"))
	    	(void) printf("\nNode %d has already sent normal_mode.\n",
	    		      node_id);
	}
	debug_print("lin_pert","Left pp_list_modes()\n");
	return	++unstable_modes;
}		/*end pp_list_modes*/


EXPORT	int pp_get_next_mode(
	int		*ith_mode,
	int		n_modes,
	LIN_PERT	*lin_pert)
{
	int		node_id, n_nodes;

	debug_print("lin_pert","Entered pp_get_next_mode()\n");
	if ((lin_pert->layer_sys->front->pp_grid->nn == 1) 
		 || (n_nodes = pp_numnodes()) <= 2 || n_modes == 1)
	{
		debug_print("lin_pert","Left pp_get_next_mode()\n");
		return 	++(*ith_mode);
	}
	
	node_id = pp_mynode();
	if (node_id == MASTER)
	{
		if (pp_iprobe_any(ND_ID1))
		{
			++(*ith_mode);
			pp_recv_any(ND_ID1,(POINTER)&node_id,INT);
			pp_send(M_NUM,(POINTER)ith_mode,INT,node_id);
		}
		debug_print("lin_pert","Left pp_get_next_mode()\n");
		return	*ith_mode < n_modes+n_nodes-2 ? -1 : n_modes;
	}
	else
	{
		pp_send(ND_ID1,(POINTER)&node_id,INT,MASTER);
		pp_recv(M_NUM,MASTER,(POINTER)ith_mode,INT);
		debug_print("lin_pert","Left pp_get_next_mode()\n");
		return	*ith_mode;
	}
}		/*end pp_get_next_mode*/


EXPORT	void	brdcst_info(
	LIN_PERT	*lin_pert,
	NORMAL_MODE	**normal_mode,
	int		num_modes,
	int		*pum)
{
	int	i, j, n_nodes, node_id, unstable_modes;

	if (lin_pert->layer_sys->front->pp_grid->nn == 1) return;

	debug_print("lin_pert","Entered brdcst_info()\n");
	if (debugging("lin_pert"))
		(void) printf("\nNode %d is ready for the broadcasting step. ",
			pp_mynode());
	n_nodes = pp_numnodes();
	if (n_nodes <= 2)
	{
		debug_print("lin_pert","Left brdcst_info()\n");
		return;
	}

	node_id = pp_mynode();
	if (node_id == MASTER)
		pp_send(UN_M,(POINTER)pum,INT,-1);
	else
		pp_recv(UN_M,MASTER,(POINTER)pum,INT);
	unstable_modes = *pum;
	if (unstable_modes == 0)
	{
		debug_print("lin_pert","Left brdcst_info()\n");
		return;
	}

	if (node_id == MASTER)
	{
		for (j = 0; j < n_nodes; j++)
			send_pts(lin_pert,j);
		for (i = 0; i < unstable_modes; i++)
		{
			for (j = 0; j < n_nodes; j++)
				send_nm(normal_mode[i],lin_pert,j,NO);
		}
	}
	else
	{
		recv_pts(lin_pert,MASTER);
		for (i = 0; i < unstable_modes; i++)
			recv_nm(normal_mode[i],lin_pert,MASTER);
	}
	for (i=unstable_modes; i<num_modes; i++)
		free(normal_mode[i]);

	if (debugging("lin_pert"))
		(void) printf("\nNode %d has finished the broadcasting step.\n\n",
			pp_mynode());
	debug_print("lin_pert","Left brdcst_info()\n");
}		/*end brdcst_info*/


LOCAL	void	send_nm(
	NORMAL_MODE	*n_m,
	LIN_PERT	*lin_pert,
	int		to,
	int		flag_free)
{
	LAYER_SYS	*layer_sys = lin_pert->layer_sys;
	int		nl = layer_sys->num_layers;
	int		i, j, k, nstep;
	int 		dim = layer_sys->front->rect_grid->dim;
	int		m = lin_pert->tot_pts;
	int		bufsize;
	double		*sbuf;

	debug_print("lin_pert","Entered send_nm()\n");
	bufsize = (dim-1)+5+2*m;
	uni_array(&sbuf,bufsize,FLOAT);

	for (i=0; i<dim-1; i++)
		sbuf[i] = n_m->wv_num[i];
	sbuf[i] = n_m->phase;
	sbuf[++i] = n_m->amp_max;
	sbuf[++i] = n_m->sigma_r;
	sbuf[++i] = n_m->sigma_i;
	sbuf[++i] = n_m->ksq;
	for (j=1; j<=nl; j++)
	{
		nstep = Rt_perturbed(comp_type(layer_sys->layer[j]->comp))
				->lin_pert_intvl;
		for (k=0; k<=nstep; k++)
		{
			sbuf[++i] = n_m->a[j][k];
			sbuf[++i] = n_m->b[j][k];
		}
		if (flag_free == YES)
		{
			free(n_m->a[j]);
			free(n_m->b[j]);
		}
	}
	pp_send(N_M,(POINTER)sbuf,FLOAT*bufsize,to);
	if (flag_free == YES)
	{
		free(n_m->a);
		free(n_m->b);
	}
	free(sbuf);
	debug_print("lin_pert","Left send_nm()\n");
}		/*end send_nm*/


LOCAL	void	recv_nm(
	NORMAL_MODE	*n_m,
	LIN_PERT	*lin_pert,
	int		from)
{
	LAYER_SYS	*layer_sys = lin_pert->layer_sys;
	int		nl = layer_sys->num_layers;
	int		i, j, k, nstep;
	int 		dim = layer_sys->front->rect_grid->dim;
	int		m = lin_pert->tot_pts;
	double		*rbuf;
	int		bufsize;

	debug_print("lin_pert","Entered recv_nm()\n");
	bufsize = (dim-1)+5+2*m;
	uni_array(&rbuf,bufsize,FLOAT);

	pp_recv(N_M,from,(POINTER)rbuf,FLOAT*bufsize);
	for (i=0; i<dim-1; i++)
		n_m->wv_num[i] = rbuf[i];
	n_m->phase = rbuf[i];
	n_m->amp_max = rbuf[++i];
	n_m->sigma_r = rbuf[++i];
	n_m->sigma_i = rbuf[++i];
	n_m->ksq = rbuf[++i];
	uni_array(&(n_m->a),nl+1,sizeof(double *));
	uni_array(&(n_m->b),nl+1,sizeof(double *));
	for (j=1; j<=nl; j++)
	{
		nstep = Rt_perturbed(comp_type(layer_sys->layer[j]->comp))
				->lin_pert_intvl;
		uni_array(&(n_m->a[j]),nstep+1,FLOAT);
		uni_array(&(n_m->b[j]),nstep+1,FLOAT);
		for (k=0; k<=nstep; k++)
		{
			n_m->a[j][k] = rbuf[++i];
			n_m->b[j][k] = rbuf[++i];
		}
	}
	free(rbuf);
	debug_print("lin_pert","Left recv_nm()\n");
}		/*end recv_nm*/


LOCAL	void	send_pts(
	LIN_PERT	*lin_pert,
	int		to)
{
	LAYER_SYS	*layer_sys = lin_pert->layer_sys;
	int		nl = layer_sys->num_layers;
	int		*sbuf, i;

	debug_print("lin_pert","Entered send_pts()\n");
	uni_array(&sbuf,nl+1,INT);
	sbuf[0] = lin_pert->tot_pts;
	for (i=1; i<=nl; i++)
	{
		sbuf[i] = Rt_perturbed(comp_type(layer_sys->layer[i]->comp))
				->lin_pert_intvl;
	}

	pp_send(PTS,(POINTER)sbuf,INT*(nl+1),to);
	free(sbuf);
	debug_print("lin_pert","Left send_pts()\n");
}		/*end send_pts*/


LOCAL	void	recv_pts(
	LIN_PERT	*lin_pert,
	int		from)
{
	LAYER_SYS	*layer_sys = lin_pert->layer_sys;
	int		nl = layer_sys->num_layers;
	int		*rbuf, i;

	debug_print("lin_pert","Entered recv_pts()\n");
	uni_array(&rbuf,nl+1,INT);
	pp_recv(PTS,from,(POINTER)rbuf,INT*(nl+1));
	lin_pert->tot_pts = rbuf[0];
	for (i=1; i<=nl; i++)
	{
		Rt_perturbed(comp_type(layer_sys->layer[i]->comp))
			->lin_pert_intvl = rbuf[i];
	}
	free(rbuf);
	debug_print("lin_pert","Left recv_pts()\n");
}		/*end recv_pts*/

#undef	FLAG
#undef	N_M
#undef	PTS
#endif /* defined(TWOD) || defined(THREED) */
