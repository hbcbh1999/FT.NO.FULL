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
*				gcheck3d.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Checks interface consistency.
*/

#if defined(THREED)

#include <gdecs/gdecs.h>

EXPORT	boolean	g_consistent_interface(
	INTERFACE *intfc)
{
	COMPONENT          mincomp, maxcomp;
	COMPONENT          comp, ncomp, pcomp;
	Gas_param          *nparams, *pparams;
	HYPER_SURF         *hs;
	HYPER_SURF_ELEMENT *hse;
	Locstate           sl, sr;
	POINT              *p;
	SURFACE            **s;
	TRI                *t;
	boolean               status;
	boolean               test_for_bad;
	int                i, ncomps;
	static Gas_param   **params = NULL;
	static int         maxlen;
	const char         *warn = "WARNING in g_consistent_interface(), ";

	if (intfc->dim != 3)
	    return NO;

	test_for_bad = debugging("bad_state");
	status = f_consistent_interface(intfc);
	mincomp = min_component(intfc);
	maxcomp = max_component(intfc);
	ncomps = maxcomp - mincomp + 1;
	
	if (ncomps > maxlen)
	{
	    if (params != NULL)
		free(params);
	    maxlen = ncomps;
	    uni_array(&params,maxlen,sizeof(Gas_param*));
	}
	
	for (comp = mincomp; comp <= maxcomp; comp++)
	{
	    params[comp-mincomp] = gas_params_for_comp(comp,intfc);
	    /*printf("#comp= %d  %d  %d    %d  %d\n",  */
	    /*    mincomp, maxcomp, comp, params[comp-mincomp],  */
	    /*gas_param_number(params[comp-mincomp])); */
	}

	/*#bjet2 */

	for (s = intfc->surfaces; s && *s; s++)
	{
	    hs = Hyper_surf(*s);
	    ncomp = negative_component(hs);
	    nparams = params[ncomp-mincomp];
	    pcomp = positive_component(hs);
	    pparams = params[pcomp-mincomp];
	    
	    /*printf("#comp= %d  %d  %d  %d   %d  %d\n",  */
	    /*    mincomp, maxcomp, ncomp, pcomp, nparams, pparams); */

	    for (t = first_tri(*s); !at_end_of_tri_list(t,*s); t = t->next)
	    {
		hse = Hyper_surf_element(t);
		for (i = 0; i < 3; i++)
		{
		    p = Point_of_tri(t)[i];
	            slsr(p,hse,hs,&sl,&sr);
		   
	            if (sl == NULL)
	            {
	                (void) printf("%s left state at point is NULL\n",warn);
		        (void) printf("component = %d\n",ncomp);
		        (void) printf("p %llu %g %g %g\n",point_number(p),
			              Coords(p)[0],Coords(p)[1],Coords(p)[2]);
		        (void) printf("surface = %llu\n",surface_number(*s));
		        print_tri(t,intfc);
		        status = NO;
	            }
	            else if (Params(sl) != nparams)
	            {
	                (void) printf("%s left params are inconsistent\n",warn);
		        (void) printf("p %llu %g %g %g\n",point_number(p),
			              Coords(p)[0],Coords(p)[1],Coords(p)[2]);
		        (void) printf("surface = %llu\n",surface_number(*s));
		        print_tri(t,intfc);

			printf("sl=%d  sr=%d\n", sl, sr);
			(void) printf("Params(sl) = %llu, "
				      "params of comp(%d) = %llu\n",
			              gas_param_number(Params(sl)),ncomp,
			              gas_param_number(nparams));
			(void) printf("Params(sr) = %llu, "
				      "params of comp(%d) = %llu\n",
			              gas_param_number(Params(sr)),pcomp,
			              gas_param_number(pparams));
		        status = NO;
	            }
	            if (sr == NULL)
	            {
	                (void) printf("%s right state at point is NULL\n",warn);
		        (void) printf("component = %d\n",pcomp);
		        (void) printf("p %llu %g %g %g\n",point_number(p),
			              Coords(p)[0],Coords(p)[1],Coords(p)[2]);
		        (void) printf("surface = %llu\n",surface_number(*s));
		        print_tri(t,intfc);
		        status = NO;
	            }
	            else if (Params(sr) != pparams)
	            {
	                (void) printf("%s right params are "
				      "inconsistent\n",warn);
		        (void) printf("p %llu %g %g %g\n",point_number(p),
			              Coords(p)[0],Coords(p)[1],Coords(p)[2]);
		        (void) printf("surface = %llu\n",surface_number(*s));
		        print_tri(t,intfc);
			
			printf("sl=%d  sr=%d\n", sl, sr);
			(void) printf("Params(sr) = %llu, "
				      "params of comp(%d) = %llu\n",
			              gas_param_number(Params(sr)),pcomp,
			              gas_param_number(pparams));
			(void) printf("Params(sl) = %llu, "
				      "params of comp(%d) = %llu\n",
			              gas_param_number(Params(sl)),ncomp,
			              gas_param_number(nparams));

		        status = NO;
	            }
	            if (test_for_bad)
	            {
	                if (is_bad_state(sl,YES,"g_consistent_interface"))
	                {
	                    (void) printf("%s left state is bad\n",warn);
		            (void) printf("Params(sl) = %llu, "
				          "params of comp(%d) = %llu\n",
			                  gas_param_number(Params(sl)),ncomp,
			                  gas_param_number(nparams));
		            fprint_raw_gas_data(stdout,sl,3);
		            (void) printf("p %llu %g %g %g\n",point_number(p),
			                  Coords(p)[0],Coords(p)[1],
					  Coords(p)[2]);
		            (void) printf("surface = %llu\n",surface_number(*s));
		            print_tri(t,intfc);
		            status = NO;
	                }
	                if (is_bad_state(sr,YES,"g_consistent_interface"))
	                {
	                    (void) printf("%s right state is bad\n",warn);
		            (void) printf("Params(sr) = %llu, "
				          "params of comp(%d) = %llu\n",
			                  gas_param_number(Params(sr)),pcomp,
			                  gas_param_number(pparams));
		            fprint_raw_gas_data(stdout,sr,3);
		            (void) printf("p %llu %g %g %g\n",point_number(p),
			                  Coords(p)[0],Coords(p)[1],
					  Coords(p)[2]);
		            (void) printf("surface = %llu\n",surface_number(*s));
		            print_tri(t,intfc);
		            status = NO;
	                }
	            }
	        }
	    }
	}
	if (status == NO)
	{
	    (void) printf("Inconsistent interface\n");
	    print_interface(intfc);
	}
	return status;
}		/*end g_consistent_interface*/

#endif /* defined(THREED) */
