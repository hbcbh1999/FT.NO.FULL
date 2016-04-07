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
*				gbifur.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This file contains the structures for gas dynamics states.
*/

#if !defined(_GBIFUR_H)
#define _GBIFUR_H

enum {
	UNTRACK_RECURSE_AT_START_NODE = 0x1,
	UNTRACK_RECURSE_AT_END_NODE   = 0x2,
	UNTRACK_RECURSE_MONO_COMP     = 0x4
};

#define untrack_set_or_recurse(flag,orient,status)			  \
	switch (orient)							  \
	{								  \
	case POSITIVE_ORIENTATION:					  \
	    if (status == YES)						  \
	    	(flag).user_untrack_flag |= UNTRACK_RECURSE_AT_START_NODE;\
	    else							  \
	    	(flag).user_untrack_flag &= ~UNTRACK_RECURSE_AT_START_NODE;\
	    break;							  \
	case NEGATIVE_ORIENTATION:					  \
	    if (status == YES)						  \
	    	(flag).user_untrack_flag |=  UNTRACK_RECURSE_AT_END_NODE; \
	    else							  \
	    	(flag).user_untrack_flag &= ~UNTRACK_RECURSE_AT_END_NODE; \
	    break;							  \
	case ORIENTATION_NOT_SET:					  \
	    screen("ERROR in untrack_set_or_recurse() macro, "            \
		   "invalid orientation\n");				  \
	    clean_up(ERROR);						  \
	    break;							  \
	}

#define untrack_set_mono_comp_recurse(flag,status)			\
	if (status == YES)						\
	    (flag).user_untrack_flag |= UNTRACK_RECURSE_MONO_COMP;	\
	else								\
	    (flag).user_untrack_flag &= ~UNTRACK_RECURSE_MONO_COMP;

#define untrack_recurse(flag,orient)					\
	(((orient) == POSITIVE_ORIENTATION) ?				\
	    (flag).user_untrack_flag & UNTRACK_RECURSE_AT_START_NODE :	\
	((orient) == NEGATIVE_ORIENTATION) ?				\
	    (flag).user_untrack_flag & UNTRACK_RECURSE_AT_END_NODE : NO)

#define untrack_mono_comp_recurse(flag)					\
	(flag).user_untrack_flag & UNTRACK_RECURSE_MONO_COMP

#define set_untrack_flag(flag,orient,s1,s2,s3,s4,s5)			\
	{								\
	    ORIENTATION opp_or = Opposite_orient(orient);		\
	    set_states_set_at_node_flag((flag),(orient),(s1));		\
	    set_states_set_at_node_flag((flag),(opp_or),(s2));		\
	    (flag).user_untrack_flag = 0;				\
	    untrack_set_or_recurse((flag),(orient),(s3));		\
	    untrack_set_or_recurse((flag),(opp_or),(s4));		\
	    untrack_set_mono_comp_recurse((flag),(s5));			\
	}

struct	_G_WAVE_CAPTURE {
	F_WAVE_CAPTURE	F_wave_capture;
	double	_init_threshold;
	double	_expand_threshold;
	double	_contract_threshold;
};
typedef struct _G_WAVE_CAPTURE G_WAVE_CAPTURE;

#define	g_wave_capture(fr)	((G_WAVE_CAPTURE*)(f_wave_capture(fr)))
#define	init_threshold(fr)	g_wave_capture(fr)->_init_threshold
#define	expand_threshold(fr)	g_wave_capture(fr)->_expand_threshold
#define	contract_threshold(fr)	g_wave_capture(fr)->_contract_threshold

#endif /* !defined(_GBIFUR_H) */
