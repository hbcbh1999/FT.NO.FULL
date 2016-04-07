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
*			spolars.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#if defined(TWOD) || defined(THREED)
#include <gdecs/gdecs.h>
#include <ginit/ginit.h>
#include <sys/types.h>
#include <time.h>
#include <sys/stat.h>

#if defined(__cplusplus)
LOCAL const int NUM_POINTS = 100;
LOCAL const int MAX_NUM_POINTS = 400;

LOCAL const unsigned int DONT_PLOT	  = 0x00;
LOCAL const unsigned int PLOT_DENSITY     = 0x01;
LOCAL const unsigned int PLOT_PRESSURE    = 0x02;
LOCAL const unsigned int PLOT_VELOCITY    = 0x04;
LOCAL const unsigned int PLOT_MOMENTUM    = 0x08;
LOCAL const unsigned int PLOT_ENERGY      = 0x10;

#else /* defined(__cplusplus) */
enum {
	NUM_POINTS = 100,
	MAX_NUM_POINTS = 400};

enum {
	DONT_PLOT	 = 0x00,
	PLOT_DENSITY     = 0x01,
	PLOT_PRESSURE    = 0x02,
	PLOT_VELOCITY  	 = 0x04,
	PLOT_MOMENTUM    = 0x08,
	PLOT_ENERGY      = 0x10};
#endif /* defined(__cplusplus) */

#define plot_density(plt_ctrl)		((plt_ctrl) & PLOT_DENSITY)
#define plot_pressure(plt_ctrl)	        ((plt_ctrl) & PLOT_PRESSURE)
#define plot_velocity(plt_ctrl)		((plt_ctrl) & PLOT_VELOCITY)
#define plot_momemtum(plt_ctrl)		((plt_ctrl) & PLOT_MOMEMTUM)
#define plot_energy(plt_ctrl)		((plt_ctrl) & PLOT_ENERGY)

LOCAL	size_t sizest;

LOCAL	void	testsolver_clean_up(void);
LOCAL	void	start_up(int,char**,INIT_DATA*);

LOCAL char *temporary_input_file = NULL;

/*ARGSUSED*/
int gsolv_main(
	int		argc,
	char**		argv,
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	IMPORT	boolean	suppress_prompts;
	static		RECT_GRID	Rgr;
	char		s[Gets_BUF_SIZE];
	Locstate	state[3];
	int		nfloats;
	int		dim = 2;

	suppress_prompts = NO;
	start_up(argc,argv,init);

	ip->root->grid->rect_grid = &Rgr;
	ip->root->front->rect_grid = &Rgr;
	ip->root->wave->rect_grid = &Rgr;
	Rgr.dim = dim;

	g_compute_sizest(PURE_NON_REACTIVE,&sizest,&nfloats,dim);
	g_set_sizeof_state(ip->root,sizest,nfloats);

	set_interface_hooks(dim,init);
	i_intfc(init) = make_interface(dim);

	set_binary_output(NO);
	clean_up(0);
	return 0;
}		 /*end main*/

/*
*				start_up():
*
*	Calls Routines to handle system error messages, and prints
*	Current messages at top of output or on the screen.
*
*/

#define SHIFT (--argc,++argv)

LOCAL	void start_up(
	int	     argc,
	char	     **argv,
	INIT_DATA    *init)
{
	IMPORT	boolean	suppress_prompts;

	init_clean_up(testsolver_clean_up,NULL);

	SHIFT;
	while (argc && argv[0][0]=='-')
	{
	    if (strncmp(argv[0],"-i",2) == 0)
	    {
	    	static char infile[1024];
	    	SHIFT;
		temporary_input_file = infile;
	    	stripcomm(infile,argv[0]);
	    	if (freopen(infile,"r",stdin) == NULL)
	    	{
	    	    screen("ERROR in start_up(), "
	    	           "can't reopen %s to %s\n","stdin",infile);
	    	    clean_up(ERROR);
	    	}
	    	suppress_prompts = YES;
	    }
	    else if (strncmp(argv[0],"-o",2) == 0)
	    {
	    	SHIFT;
	    	if (freopen(argv[0],"w",stdout) == NULL)
	    	{
	    	    screen("ERROR in start_up(), "
	    	           "can't reopen %s to %s\n","stdout",argv[0]);
	    	    clean_up(ERROR);
	    	}
	    }
	    else if (strncmp(argv[0],"-e",2) == 0)
	    {
	    	SHIFT;
	    	if (freopen(argv[0],"w",stderr) == NULL)
	    	{
	    	    screen("ERROR in start_up(), ");
	    	    screen("can't reopen %s to %s\n","stderr",argv[0]);
	    	    clean_up(ERROR);
	    	}
	    }
	    SHIFT;
	}

	set_error_immediate(stdout);
	record_print_version(stdout);
	print_title(stdout,title(init));
	init_prompting_and_debugging(init);
}		/*end start_up*/

LOCAL	void testsolver_clean_up(void)
{
	if ((temporary_input_file != NULL) && (is_io_node(pp_mynode())))
	    (void) unlink(temporary_input_file);
}		/*end spolars_clean_up*/

#endif /* defined(TWOD) || defined(THREED) */
