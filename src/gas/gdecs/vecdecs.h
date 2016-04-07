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
*				vecdecs.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains definitions for uni_arrayization in C.
*/

#if !defined(_VECDECS_H)
#define _VECDECS_H

		/* DOUBLE PRECISION VERSION */

#define inner_product dinner_product

#define norm_vector dnorm_vector

#define zero_vector dzero_vector

#define subtract_vector_from_vector dsubtract_vector_from_vector

#define add_scalar_times_vector_to_vector dadd_scalar_times_vector_to_vector

#define add_vector_times_vector_to_vector dadd_vector_times_vector_to_vector

#define copy_vector_to_vector dcopy_vector_to_vector

#define multiply_vector_by_scalar dmultiply_vector_by_scalar

#define multiply_vector_by_vector dmultiply_vector_by_vector

#define divide_vector_by_vector ddivide_vector_by_vector

#define vector_equals_scalar dvector_equals_scalar

#define add_scalar_to_vector dadd_scalar_to_vector

#define min_max_vector dmin_max_vector

#define sum_elements_of_vector dsum_elements_of_vector

#define vector_equals_scalar_plus_vector dvector_equals_scalar_plus_vector

#define vector_equals_scalar_times_vector dvector_equals_scalar_times_vector

#define vector_equals_vector_plus_scalar_times_vector			\
		dvector_equals_vector_plus_scalar_times_vector

#define vector_equals_vector_minus_vector dvector_equals_vector_minus_vector


	/* Functions not in toblas */

#define add_vector_to_vector(n,v1,v2)					\
		dadd_scalar_times_vector_to_vector(n,1.0,v1,v2)

#define vector_equals_vector_plus_vector(n,v1,v2,v3)			\
{									\
	dcopy_vector_to_vector(n,v2,v1);				\
	dadd_scalar_times_vector_to_vector(n,1.0,v3,v1);		\
}

#define vector_equals_vector_times_vector(n,v1,v2,v3)			\
{									\
	dcopy_vector_to_vector(n,v2,v1);				\
	dmultiply_vector_by_vector(n,v1,v3);				\
}

#define vector_equals_vector_divided_by_vector(n,v1,v2,v3)		\
{									\
	dcopy_vector_to_vector(n,v2,v1);				\
	ddivide_vector_by_vector(n,v1,v3);				\
}

#define vector_equals_scalar_times_sum_of_vectors(n,v1,a,v2,v3)		\
{									\
	dvector_equals_vector_plus_vector(n,v1,v2,v3);			\
	dmultiply_vector_by_scalar(n,v1,a);				\
}

#define absolute_value_of_vector(n,v)					\
{									\
	int i;								\
									\
	for (i = 0; i < n; i++) v[i] = fabs(v[i]);			\
}


		/* MIXED MODE FUNCTIONS */
	/* Functions not in toblas */

#define dadd_vector_to_vector(n,v1,v2)					\
		dadd_scalar_times_vector_to_vector(n,1.0,v1,v2)

#define dvector_equals_vector_plus_vector(n,v1,v2,v3)			\
{									\
	dcopy_vector_to_vector(n,v2,v1);				\
	dadd_scalar_times_vector_to_vector(n,1.0,v3,v1);		\
}

#define dvector_equals_vector_times_vector(n,v1,v2,v3)			\
{									\
	dcopy_vector_to_vector(n,v2,v1);				\
	dmultiply_vector_by_vector(n,v1,v3);				\
}

#define dvector_equals_vector_divided_by_vector(n,v1,v2,v3)		\
{									\
	dcopy_vector_to_vector(n,v2,v1);				\
	ddivide_vector_by_vector(n,v1,v3);				\
}

#define dvector_equals_scalar_times_sum_of_vectors(n,v1,a,v2,v3)	\
{									\
	dvector_equals_vector_plus_vector(n,v1,v2,v3);			\
	dmultiply_vector_by_scalar(n,v1,a);				\
}

#define dabsolute_value_of_vector(n,v)					\
{									\
	int i;								\
									\
	for (i = 0; i < n; i++) v[i] = fabs(v[i]);			\
}


#endif /* !defined(_VECDECS_H) */
