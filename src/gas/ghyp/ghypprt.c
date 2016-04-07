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
*				ghypprt.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Printing and debugging routines for the ghyp solvers.
*/

#include <ghyp/ghyp.h>
#include <gdecs/vecdecs.h>

EXPORT	 void g_print_internal_energy(
	const char	*mesg,
	double		**ucon,
	Vec_Muscl	*vmuscl,
	int		start,
	int		end)
{
	int		i, j, dim = vmuscl->dim;
	double		int_e, ke;
	double		E, m[MAXD], rho;
	Locstate	*state = vmuscl->vst->state + vmuscl->offset;
	static Locstate st = NULL;

	if (st == NULL)
	{
	    g_alloc_state(&st,vmuscl->sizest);
	}

	(void) printf("%s\n",mesg);
	(void) printf("%-4s %-14s %-14s\n","n","Int_energy","Pressure");
	for (j = start; j < end; ++j)
	{
	    Dens(st) = rho = ucon[0][j];
	    Energy(st) = E = ucon[1][j];
	    ke = 0.0;
	    for (i = 0; i < dim; ++i)
	    {
	    	Mom(st)[i] = m[i] = ucon[2+i][j];
	    	ke += 0.5*sqr(m[i])/rho;
	    }
	    Set_params(st,state[j]);
	    set_type_of_state(st,GAS_STATE);
	    reset_gamma(st);
	    int_e = (E - ke)/rho;
	    (void) printf("%-4d %-14g %-14g",j,int_e,pressure(st));
	    if (int_e < 0.0)
	    	(void) printf("\tNEGATIVE INTERNAL ENERGY\n");
	    else
	    	(void) printf("\n");
	}
}		/*end g_print_internal_energy*/


EXPORT	void	g_printout_vec_data(
	const char	*mesg,
	double		*u1,
	double		*u2,
	double		**u3,
	int		dim,
	int		start,
	int		end,
	const char	*label)
{
	const char *name[5];

	if (strcmp(label,"src") == 0) /*Source uni_arrays*/
	{
	    name[0] = "mass";
	    name[1] = "energy";
	    name[2] = "mom[0]";
	    name[3] = "mom[1]";
	    name[4] = "mom[2]";
	}
	else if (strcmp(label,"muncons") == 0) 
	{				/*Unconserved (MUSCL) state vars*/
	    name[0] = "density";
	    name[1] = "int_en";
	    name[2] = "vel[0]";
	    name[3] = "vel[1]";
	    name[4] = "vel[2]";
	}
	else if (strcmp(label,"cons") == 0) /*Conserved state vars*/
	{
	    name[0] = "density";
	    name[1] = "en_den";
	    name[2] = "mom[0]";
	    name[3] = "mom[1]";
	    name[4] = "mom[2]";
	}
	else
	{
	    name[0] = "uni_array1";
	    name[1] = "uni_array2";
	    name[2] = "uni_array3";
	    name[3] = "uni_array4";
	    name[4] = "uni_array5";
	}

	(void) printf("%s\n",mesg);
	g_print_state_vectors(start,end,name,u1,u2,u3,dim);
}		/*end g_printout_vec_data*/

EXPORT void g_print_state_vectors(
	int		n1,
	int		n2,
	const char	**name,
	double		*v1,
	double		*v2,
	double		**v3,
	int		dim)
{
	int		i,j;

	(void) printf("%-4s %-14s %-14s","n",name[0],name[1]);
	for (j = 0; j < dim; ++j)
	    (void) printf(" %-14s",name[2+j]);
	(void) printf("\n");
	for (i = n1; i < n2; ++i) 
	{
	    (void) printf("%-4d %-22.16g %-22.16g",i,v1[i],v2[i]);
	    for (j = 0; j < dim; ++j)
	    	(void) printf(" %-22.16g",v3[j][i]);
	    (void) printf("\n");
	}
	(void) printf("\n");
}		/*end g_print_state_vectors*/


