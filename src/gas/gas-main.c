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
*
*				gas-main.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <ginit/ginit.h>

LOCAL	const char *__FronTier_VERSION__ = "20020808";

EXPORT	int C_MAIN_PROGRAM(
	int			argc,
	char			**argv)
{
	static G_INIT_DATA	Ginit;
	static G_INIT_PHYSICS	GIphys;
	static Grid		grid;
	static CHART		root;
	static G_Front		gfront;
	static G_Wave		gwave;
	static G_Printplot	gprt;
	INIT_DATA		*init = init_data(&Ginit);
	const char		*executable;
	char                    title[1024];
	int			i;

	/*Set executable name*/
	executable = NULL;
	(void) sprintf(title,"\n\t\tWELCOME TO WIND TUNNEL CALCULATIONS\n"
	                     "\t\tFronTier Gas version %s\n",
			     __FronTier_VERSION__);

	title(init) = title;
	for (i = 0; i < argc; ++i)
	{
	    if (strncmp(argv[i],"-spolars",8)==0)
	    {
	    	executable = "spolars";
	    	break;
	    }
	    if (strncmp(argv[i],"-gsolvers",9)==0)
	    {
	    	executable = "gsolvers";
	    	break;
	    }
	    if (strncmp(argv[i],"-eosplot",8)==0)
	    {
	    	executable = "eosplot";
	    	break;
	    }
	    if (strncmp(argv[i],"-rage90",7)==0)
	    {
	    	executable = "rage90";
	    	break;
	    }
	}
	if (executable == NULL)
	    executable = strdup(basename(argv[0]));
	if (i < argc)
	{
	    /*remove executable option from argv*/
	    for (argc--; i < argc; ++i)
	    	argv[i] = argv[i+1];
	}

	root.front = (Front*)&gfront;
	root.wave = (Wave*)&gwave;
	wave_of_front(root.front) = root.wave;

	chart_of_wave(&gwave) = &root;
	chart_of_front(&gfront) = &root;

	root.grid = &grid;
	GIphys.Iphys.root = &root;
	GIphys.Iphys.prt = (Printplot*)&gprt;
	root.prt = GIphys.Iphys.prt;

	ex_name(init) = executable;

	set_gas_hooks(init,&GIphys.Iphys);

#if defined(TWOD)
	if (strcmp(executable,"spolars") == 0)
	    return sp_main(argc,argv,init,&GIphys.Iphys);
	if (strcmp(executable,"gsolvers") == 0)
	    return gsolv_main(argc,argv,init,&GIphys.Iphys);
	if (strcmp(executable,"eosplot") == 0)
	    return eosp_main(argc,argv,init,&GIphys.Iphys);
#endif /* defined(TWOD) */
#if defined(RAGE90)
	if (strcmp(executable,"rage90") == 0)
	    return rage90_main(argc,argv,init,&GIphys.Iphys);
#endif /* defined(RAGE90) */

	return dmain(argc,argv,init,&GIphys.Iphys);
}		/*end main*/

#if defined(mips)
EXPORT	void cvd_bug_workaround(
	Gas         *gas,
	HYPER_SURF  *hs,
	G_INIT_DATA *ginit)
{
	F_HYPER_SURF Fhs;
	F_POINT Fp;

	gas->rho = ERROR_FLOAT;
	bstate_index(hs) = -1;
	ginit->_prompt_for_ref_state = NULL;
	Fhs._wave_type = CONTACT;
	wave_type(hs) = Fhs._wave_type;
	Fp.point._coords[0] = 0.0;
	Fp._left_state = Fp._right_state = NULL;
	print_general_vector("Fp",Fp.point._coords,1,"\n");
}		/*end cvd_bug_workaround*/
#endif /* defined(mips) */
