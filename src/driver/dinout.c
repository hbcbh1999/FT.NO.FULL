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
*				dinout.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Provides input and output functions for state variables
*	associated with an interface and an underlying rectangular
*	grid.
*
*
*	Usage for Output:
*	----------------
*
*	The very first line in the output file must be generated
*	by the call
*
*		record_print_version().
*
*	The format for printing defaults to the latest version, but it
*	may be modified by a call to
*
*		set_print_version(print_version).
*
*	BEFORE the call to record_print_version().  (At present the print
*	format applies to the whole printout file, so modification after
*	the call to record_print_version() will ruin the output file.)
*
*	It is desirable for some programs that read the output file
*	that the underlying rectangular grid for the interfaces be
*	printed before any interface or state information is printed.
*	This can be accomplished with a call to print_rect_grid().
*	(Actually, this is unnecessary for print_version >= 1.)
*
*	The interface should be printed before the printouts for the
*	various state variables by calling print_interface().
*
*	The state information is printed by calling
*
*		print_states(file,prt).
*
*	print_states uses the values n_vars = prt->n_rect_state_vars,  and
*	os = prt->output_soln. Here "n_vars" is the number of variables to
*	be printed, and "output_soln" is a vector of pointers to 
*	OUTPUT_SOLN structures with length "n_vars".  Each member
*	os = output_soln[i] must have been initialized as follows:
*
*	os->name	-- a char * pointing to the string to be printed out
*			   to identify the state variable being printed
*	os->fit		-- an int with value either CONSTANT or LINEAR
*	os->smoothness	-- an int with value either SMOOTH or SINGULAR
*
*	os->grid	-- a RECT_GRID * that defines the rectangular grid;
*			   used ONLY IF os->intfc == NULL 
*	os->intfc	-- the interface overlayed on the rectangular grid;
*			   used ONLY IF smoothness== SINGULAR.
*
*	os->solution	-- a pointer to a function that returns the solution to
*			   be printed at the centers of grid blocks of the
*			   rectangular grid;
*	os->mv_solution	-- a pointer to a function that returns the multivalued
*			   solution for the given mesh block and component wrt
*			   os->intfc;
*			   used ONLY IF fit==CONSTANT and smoothness==SINGULAR.
*	os->intfc_solution	-- a pointer to a function that sets the
*			   states on the two sides of os->intfc at the given
*			   point; used ONLY IF fit==LINEAR and
*			   smoothness==SINGULAR.
*
*	os->var		-- an integer that is passed to the solution functions.
*	os->extra	-- a POINTER that is passed to the solution functions.
*
*	The information to be printed is specified by the fit and smoothness
*	flags in the following way:
*
*	CONSTANT/SMOOTH:
*		(*os->solution)(os,coords,icoords) is called for
*	coords at the centers of the grid blocks of
*	gr = computational_grid(os->intfc) (or gr = os->grid if
*	os->intfc == NULL or doesn't have a fixed grid). These centers are
*	located at
*
*		cell_center(i,j,gr)
*
*	for i = 0, ..., gr->gmax[j] - 1, j = 0, ..., dim-1
*	The double value returned by (*os->solution)() is printed.
*
*	CONSTANT/SINGULAR:
*		In addition to the CONSTANT/SMOOTH printout,
*	(*os->mv_solution)(os,ix,iy,comp) is called for
*	ix = 0, ..., gr->gmax[0] - 1, iy = 0, ..., gr->gmax[1] - 1, and comp
*	running through the components in the ix,iy grid block wrt os->intfc.
*	(If os->intfc doesn't have a fixed_grid, os->grid is used; if
*	necessary, the bond and comp lists are constructed.)  The double
*	value returned by (*os->mv_solution)() is printed.
*
*	LINEAR/SMOOTH:
*		(*os->solution)(os,coords,icoords) is called for
*	coords at the crossings of grid lines of computational_grid(os->intfc)
*	(or os->grid if os->intfc == NULL or doesn't have a fixed grid).
*	These crossings are located at
*
*		cell_edge(i,j,gr)
*
*	for i = 0, ..., gr->gmax[j], j = 0, ..., dim-1
*	The double value returned by (*os->solution)() is printed.
*
*	LINEAR/SINGULAR:
*		In addition to the LINEAR/SMOOTH printout,
*	(*os->intfc_solution)(os,point,hse,hs,&left,&right)
*	is called for every point on os->intfc.  The double values loaded
*	into left and right are printed as the states on the left and right
*	side of the interface.
*
*
*	Usage for Input:
*	----------------
*
*	...
*
*
*	BUGS:
*
*	1. Presently LINEAR/SINGULAR uses the CONSTANT/SMOOTH printout
*	instead of the LINEAR/SMOOTH printout unless STATES_AT_NODES is
*	#define'd.  This means that input_solution() is only piecewise constant.
*	This is bad for restarts using 2d Lax-Wendroff.
*
*	TODO:
*
*	1. Write more documentation.
*
*	2. Write a function
*
*		input_solution_function(var,x,y,comp,input_soln,answer)
*
*	for returning the answer at an arbitrary point.  Perhaps trisoln.c
*	should be used to write a more accurate solution function.
*
*	3. Once 2 is done, rewrite input_solution(), assuming that it is
*	being called as though from the printout routines, thus gaining
*	efficiency.  Change code that uses these functions to use
*	input_solution_function() instead.
*
*		Additional (helpful (hopefully)) notes:
*	1. Data structure Printplot defined in ddecs.h
*		Let	os = prt->output_soln[i], 0 <= i < prt->n_vars
*			rs = prt->restart_soln
*	2. os->smoothness	= SMOOTH	variable continuous across intfc
*				= SINGULAR 	variable (state) discontinuous
*						across interface.
*						uses mv_solution if CONSTANT
*						uses intfc_solution if LINEAR
*	3. os->fit	= CONSTANT	variable at cell center
*					(xmax,ymax)
*			= LINEAR	variable at vertices of grid
*					(xmax+1,ymax+1) if SMOOTH
*					(xmax,ymax) if SINGULAR
*	4. driver is dsub.c/print_front_and_wave()
*		sets following variables
*		POINTER		os->extra[2]
*		extra[0] = front
*		extra[1] = wave
*		os->intfc = front->interf
*	5. restart
*		prt->n_restart_vars <= prt->n_vars
*		rs->name same as os->name
*		rs->fit  same as os->fit
*		rs->smoothness same as os->smoothness
*		variable names for rs must be in same ascending order as os
*/


#include <driver/ddecs.h>

	/* LOCAL Function Declarations */
LOCAL	boolean	allocate_input_soln(INTERFACE*,INPUT_SOLN*,int);
LOCAL	boolean	read_intfc_states(const IO_TYPE*,int,INPUT_SOLN*);
LOCAL	boolean	read_printout_of_states(const IO_TYPE*,int,INPUT_SOLN*);
LOCAL	boolean	read_states(const IO_TYPE*,INTERFACE*,int,INPUT_SOLN*,int);
LOCAL	boolean    read_state_variables1d(const IO_TYPE*,int,INTERFACE*,
                                       INPUT_SOLN**);
LOCAL	uint64_t	hs_point_identifier(HYPER_SURF*);
LOCAL	uint64_t	hs_curve_identifier(HYPER_SURF*);
LOCAL	uint64_t	hs_surface_identifier(HYPER_SURF*);
LOCAL	void	bin_prt_soln_value(FILE*,OUTPUT_VALUE*);
LOCAL	void	bin_prt_soln_values(FILE*,OUTPUT_VALUE*,OUTPUT_VALUE*);
LOCAL	void	prt_soln_value(FILE*,OUTPUT_VALUE*);
LOCAL	void	prt_soln_values(FILE*,OUTPUT_VALUE*,OUTPUT_VALUE*,const char*);
LOCAL	void	print_array_of_states1d(FILE*,RECT_GRID*,OUTPUT_SOLN*);
LOCAL	void	print_array_of_states2d(FILE*,RECT_GRID*,OUTPUT_SOLN*);
LOCAL	void	print_array_of_states3d(FILE*,RECT_GRID*,OUTPUT_SOLN*);
LOCAL	void	print_intfc_states(FILE*,OUTPUT_SOLN*);
LOCAL	void	print_states0(FILE*,OUTPUT_SOLN*);
LOCAL	void	print_states1(FILE*,OUTPUT_SOLN*);
LOCAL	void	print_states_format(FILE*,RECT_GRID*);
LOCAL	void	print_states_format_with_bdry(FILE*,INTERFACE*,RECT_GRID*);
LOCAL	void	read_array_of_states1d(const IO_TYPE*,int*,INPUT_SOLN*);
LOCAL	void	read_array_of_states2d(const IO_TYPE*,int*,INPUT_SOLN*);
LOCAL	void	read_array_of_states3d(const IO_TYPE*,int*,INPUT_SOLN*);
LOCAL	void	read_grid_fit_and_smoothness(FILE*,REMAP*,
                                             RECT_GRID*,int*,int*,int);
LOCAL	void	set_grid(const char*,INTERFACE*,int,int,RECT_GRID*,int);


LOCAL int print_version = 1;

LOCAL void (*print_array_of_states[4])(FILE*,RECT_GRID*,OUTPUT_SOLN*) =
					{NULL,
					print_array_of_states1d,
					print_array_of_states2d,
					print_array_of_states3d};

LOCAL void (*read_array_of_states[4])(const IO_TYPE*,int*,INPUT_SOLN*) =
					{NULL,
					read_array_of_states1d,
					read_array_of_states2d,
					read_array_of_states3d};

LOCAL const char *ffmt	= "%- "FFMT"%s";
LOCAL const char *ifmt	= "%-10d%s";
LOCAL const char *pfmt	= "%-10p%s";
LOCAL const char *lfmt	= "%-10ld%s";
LOCAL const char *ufmt	= "%-10u%s";
LOCAL const char *ulfmt	= "%-10llu%s";

LOCAL	void	bin_prt_soln_value(
	FILE		*file,
	OUTPUT_VALUE	*value)
{
	switch (value->utype)
	{
	case Float:
	    (void) fwrite((const void *)&value->uval.fval,sizeof(double),1,file);
	    break;
	case Int:
	    (void) fwrite((const void *)&value->uval.ival,sizeof(int),1,file);
	    break;
	case Pointer:
	    (void) fwrite((const void *)&value->uval.pval,
			  sizeof(POINTER),1,file);
	    break;
	case Long:
	    (void) fwrite((const void *)&value->uval.lval,
			  sizeof(long int),1,file);
	    break;
	case Unsigned:
	    (void) fwrite((const void *)&value->uval.uval,
			  sizeof(unsigned int),1,file);
	    break;
	case ULong:
	    (void) fwrite((const void *)&value->uval.ulval,
			  sizeof(value->uval.ulval),1,file);
	    break;
	default:
	    screen("ERROR in bin_prt_soln_value(), unknown data type %d\n",
		   value->utype);
	    clean_up(ERROR);
	}
}		/*end bin_prt_soln_value*/

LOCAL	void	prt_soln_value(
	FILE		*file,
	OUTPUT_VALUE	*value)
{
	switch (value->utype)
	{
	case Float:
	    (void) fprintf(file,ffmt,value->uval.fval," ");
	    break;
	case Int:
	    (void) fprintf(file,ifmt,value->uval.ival," ");
	    break;
	case Pointer:
	    (void) fprintf(file,pfmt,value->uval.pval," ");
	    break;
	case Long:
	    (void) fprintf(file,lfmt,value->uval.lval," ");
	    break;
	case Unsigned:
	    (void) fprintf(file,ufmt,value->uval.uval," ");
	    break;
	case ULong:
	    (void) fprintf(file,ulfmt,value->uval.ulval," ");
	    break;
	default:
	    screen("ERROR in prt_soln_value(), unknown data type %d\n",
		   value->utype);
	    clean_up(ERROR);
	}
}		/*end prt_soln_value*/


LOCAL	void	bin_prt_soln_values(
	FILE		*file,
	OUTPUT_VALUE	*left,
	OUTPUT_VALUE	*right)
{
	(void) fprintf(file,"\f%c",2);
	switch (left->utype)
	{
	case Float:
	    (void) fwrite((const void *)&left->uval.fval,sizeof(double),1,file);
	    (void) fwrite((const void *)&right->uval.fval,sizeof(double),1,file);
	    break;
	case Int:
	    (void) fwrite((const void *)&left->uval.ival,sizeof(int),1,file);
	    (void) fwrite((const void *)&right->uval.ival,sizeof(int),1,file);
	    break;
	case Pointer:
	    (void) fwrite((const void *) &left->uval.pval,
			  sizeof(POINTER),1,file);
	    (void) fwrite((const void *) &right->uval.pval,
			  sizeof(POINTER),1,file);
	    break;
	case Long:
	    (void) fwrite((const void *) &left->uval.lval,
			  sizeof(long int),1,file);
	    (void) fwrite((const void *) &right->uval.lval,
			  sizeof(long int),1,file);
	    break;
	case Unsigned:
	    (void) fwrite((const void *) &left->uval.uval,
			  sizeof(unsigned int),1,file);
	    (void) fwrite((const void *) &right->uval.uval,
			  sizeof(unsigned int),1,file);
	case ULong:
	    (void) fwrite((const void *) &left->uval.ulval,
			  sizeof(left->uval.ulval),1,file);
	    (void) fwrite((const void *) &right->uval.ulval,
			  sizeof(right->uval.ulval),1,file);
	    break;
	default:
	    screen("ERROR in bin_prt_soln_values(), unknown data type %d\n",
		   left->utype);
	    clean_up(ERROR);
	}
}		/*end bin_prt_soln_values*/

LOCAL	void	prt_soln_values(
	FILE		*file,
	OUTPUT_VALUE	*left,
	OUTPUT_VALUE	*right,
	const char	*end)
{
	switch (left->utype)
	{
	case Float:
	    (void) fprintf(file,ffmt,left->uval.fval," ");
	    (void) fprintf(file,ffmt,right->uval.fval,end);
	    break;
	case Int:
	    (void) fprintf(file,ifmt,left->uval.ival," ");
	    (void) fprintf(file,ifmt,right->uval.ival,end);
	    break;
	case Pointer:
	    (void) fprintf(file,pfmt,left->uval.pval," ");
	    (void) fprintf(file,pfmt,right->uval.pval,end);
	    break;
	case Long:
	    (void) fprintf(file,lfmt,left->uval.lval," ");
	    (void) fprintf(file,lfmt,right->uval.lval,end);
	    break;
	case Unsigned:
	    (void) fprintf(file,ufmt,left->uval.uval," ");
	    (void) fprintf(file,ufmt,right->uval.uval,end);
	    break;
	case ULong:
	    (void) fprintf(file,ulfmt,left->uval.ulval," ");
	    (void) fprintf(file,ulfmt,right->uval.ulval,end);
	    break;
	default:
	    screen("ERROR in prt_soln_values(), unknown data type %d\n",
		   left->utype);
	    clean_up(ERROR);
	}
}		/*end prt_soln_values*/


EXPORT void set_print_version(
	int		version)
{
	print_version = version;
}		/*end set_print_version*/



EXPORT void record_print_version(
	FILE		*file)
{
	(void) foutput(file);
	(void) fprintf(file,"print version %d\n",print_version);
}		/*end record_print_version*/





EXPORT void d_print_states(
	FILE		*file,
	Printplot	*prt)
{
	int		n_vars = prt->n_rect_state_vars;
	OUTPUT_SOLN	**output_soln = prt->output_soln;
	int		var;

	switch (print_version)
	{
	default:
	case 0:
	    for (var = 0; var < n_vars; ++var)
	    	print_states0(file,output_soln[var]);
	    break;
	case 1:
	    for (var = 0; var < n_vars; ++var)
	    	print_states1(file,output_soln[var]);
	    break;
	}
}		/*end print_states*/


/*ARGSUSED*/
EXPORT void d_print_states1d(
	FILE		*file,
	Printplot	*prt)
{
	OUTPUT_SOLN	**os = prt->output_soln;
	OUTPUT_VALUE	*value;
	OUTPUT_VALUE	left, right;
	POINT		*p;
	RECT_GRID	*c_gr, Gr, *gr = &Gr;
	double		crds[MAXD];
	int		n_vars = prt->n_rect_state_vars;
	int		xstart, xmax;
	int		ix, icrds[MAXD];
	int		i, dim, num_points;
	int		var;

	if ((os[0]->fit | os[0]->smoothness) != (LINEAR_FIT | SINGULAR))
	{
	    d_print_states(file,prt);
	    return;
	}

	if (os[0]->intfc != NULL)
	    c_gr = computational_grid(os[0]->intfc);
	else
	    c_gr = os[0]->grid;

	copy_rect_grid(gr,c_gr);
	dim = gr->dim;

	(void) fprintf(file,"\n\n");
	(void) foutput(file);
	(void) fprintf(file,"Printout of one dimensional states\n");
	for (i = 0; i < dim; ++i)
	    gr->L[i] = cell_center(0,i,gr);
	print_states_format(file,c_gr);
	(void) fprintf(file,"linear singular\n\n");
	(void) foutput(file);
	(void) fprintf(file,"%-15s","POSITION");
	for (var = 0; var < n_vars; ++var)
	    (void) fprintf(file," %-16s",os[var]->name);
	(void) fprintf(file,"\n");

	xstart = 0;
	xmax = gr->gmax[0];
	if (debugging("buffer"))
	{
	    xstart -= c_gr->lbuf[0];
	    xmax += c_gr->ubuf[0];
	}
	if (is_binary_output() == YES)
	{
	    for (ix = xstart; ix < xmax; ++ix)
	    {
	        (void) fprintf(file,"\f%c",n_vars+1);
	        icrds[0] = ix;
	        crds[0] = cell_edge(ix,0,gr);
	        (void) fwrite((const void *)crds,sizeof(double),1,file);
	        for (var = 0; var < n_vars; ++var)
	        {
	            value = (*os[var]->solution)(os[var],crds,icrds);
		    bin_prt_soln_value(file,value);
	        }
	    }
	    (void) fprintf(file,"\n");
	}
	else
	{
	    for (ix = xstart; ix < xmax; ++ix)
	    {
	        icrds[0] = ix;
	        crds[0] = cell_edge(ix,0,gr);
	        (void) fprintf(file,ffmt,crds[0],"");
	        for (var = 0; var < n_vars; ++var)
	        {
	            (void) fprintf(file," ");
	            value = (*os[var]->solution)(os[var],crds,icrds);
		    prt_soln_value(file,value);
	        }
	        (void) fprintf(file,"\n");
	    }
	}
	(void) fprintf(file,"\n");

	(void) fprintf(file,"\nStart of intfc states\n");
	(void) fprintf(file,"%-12s","POSITION");
	for (var = 0; var < n_vars; ++var)
	    (void) fprintf(file," L_%-12s R_%-12s",
	    	os[var]->name,os[var]->name);
	(void) fprintf(file,"\n");
	num_points = os[0]->intfc->num_points;
	if (is_binary_output() == YES)
	{
	    for (i = 0; i < num_points; ++i)
	    {
	        p = os[0]->intfc->points[i];
	        (void) fprintf(file,"\f%c",2*n_vars+1);
	        (void) fwrite((const void *)Coords(p),sizeof(double),1,file);
	        for (var = 0; var < n_vars; ++var)
	        {
		    p = os[var]->intfc->points[i];
		    (*os[var]->intfc_solution)(os[var],p,NULL,Hyper_surf(p),
					       &left,&right);
		    bin_prt_soln_value(file,&left);
		    bin_prt_soln_value(file,&right);
		}
	    }
	    (void) fprintf(file,"\n");
	}
	else
	{
	    for (i = 0; i < num_points; ++i)
	    {
	        p = os[0]->intfc->points[i];
	        (void) fprintf(file,ffmt,Coords(p)[0],"");
	        for (var = 0; var < n_vars; ++var)
	        {
		    p = os[var]->intfc->points[i];
		    (*os[var]->intfc_solution)(os[var],p,NULL,Hyper_surf(p),
					       &left,&right);
		    (void) fprintf(file," ");
		    prt_soln_values(file,&left,&right," ");
		}
	        (void) fprintf(file,"\n");
	    }
	}

	(void) fprintf(file,"\nEnd of intfc states\n\n");
	(void) fprintf(file,"\n\n");
	(void) foutput(file);
	(void) fprintf(file,"End of printout of one dimensional states\n");
}		/*end d_print_states1d*/

LOCAL void print_states0(
	FILE		*file,
	OUTPUT_SOLN	*os)
{
	RECT_GRID	*c_gr, Gr, *gr = &Gr;
	int		i, dim;
	
	if (os->intfc != NULL)
	    c_gr = computational_grid(os->intfc);
	else
	    c_gr = os->grid;

	copy_rect_grid(gr,c_gr);
	dim = gr->dim;

	(void) fprintf(file,"\n\n");
	(void) foutput(file);

	switch (os->fit | os->smoothness)
	{
	case LINEAR_FIT | SMOOTH:
	    (void) fprintf(file,"bi_array ");
	    for (i = 0; i < dim; ++i)
	    {
	    	gr->gmax[i] = gr->gmax[i] + 1;
	    	(void) fprintf(file,ifmt,gr->gmax[i]," ");
	    }
	    (void) fprintf(file,"%s\n\n",os->name);
	    (*print_array_of_states[dim])(file,gr,os);
	    (void) fprintf(file,"\n");
	    break;

	case LINEAR_FIT | SINGULAR:
#if defined(STATES_AT_NODES)
	    (void) fprintf(file,"bi_array ");
	    for (i = 0; i < dim; ++i)
	    {
	    	gr->gmax[i] = gr->gmax[i] + 1;
	    	(void) fprintf(file,ifmt,gr->gmax[i]," ");
	    }
	    (void) fprintf(file,"%s\n\n",os->name);
	    (*print_array_of_states[dim])(file,gr,os);

#else /* defined(STATES_AT_NODES) */

	    (void) fprintf(file,"bi_array ");
	    for (i = 0; i < dim; ++i)
	    {
	    	(void) fprintf(file,ifmt,gr->gmax[i]," ");
	    	gr->L[i] = cell_center(0,i,gr);
	    }
	    (void) fprintf(file,"%s\n\n",os->name);
	    (*print_array_of_states[dim])(file,gr,os);
#endif /* defined(STATES_AT_NODES) */
	    print_intfc_states(file,os);
	    (void) fprintf(file,"\n\n");
	    (void) foutput(file);
	    (void) fprintf(file,"End of %s\n",os->name);
	    break;

	default:
	case CONSTANT_FIT | SMOOTH:
	    (void) fprintf(file,"bi_array ");
	    for (i = 0; i < dim; ++i)
	    {
	    	(void) fprintf(file,ifmt,gr->gmax[i]," ");
	    	gr->L[i] = cell_center(0,i,gr);
	    }
	    (void) fprintf(file,"%s\n\n",os->name);
	    (*print_array_of_states[dim])(file,gr,os);
	    (void) fprintf(file,"\n");
	    break;

	}
}		/*end print_states0*/



LOCAL void print_states1(
	FILE		*file,
	OUTPUT_SOLN	*os)
{
	RECT_GRID	*c_gr, Gr, *gr = &Gr;
	int		i, dim;

	if (os->intfc != NULL)	c_gr = computational_grid(os->intfc);
	else			c_gr = os->grid;

	copy_rect_grid(gr,c_gr);
	dim = gr->dim;

	(void) fprintf(file,"\n\n");
	(void) foutput(file);
	(void) fprintf(file,"\t\t\t%s\n",os->name);

	switch (os->fit | os->smoothness)
	{
	case LINEAR_FIT | SMOOTH:
	    for (i = 0; i < dim; ++i)
	        gr->gmax[i] = gr->gmax[i] + 1;
	    (void) fprintf(file,"linear smooth\n\n");
	    (*print_array_of_states[dim])(file,gr,os);
	    (void) fprintf(file,"\n");
	    break;
	case LINEAR_FIT | SINGULAR:
#if defined(STATES_AT_NODES)
	    for (i = 0; i < dim; ++i)
	        gr->gmax[i] = gr->gmax[i] + 1;
	
	    (void) fprintf(file,"linear singular\n\n");
	    (*print_array_of_states[dim])(file,gr,os);

#else /* defined(STATES_AT_NODES) */
	    for (i = 0; i < dim; ++i)
	        gr->L[i] = cell_center(0,i,gr);
	    
	    print_states_format_with_bdry(file,os->intfc,c_gr);
	    
	    (void) fprintf(file,"linear singular\n\n");
	    (*print_array_of_states[dim])(file,gr,os);
#endif /* defined(STATES_AT_NODES) */
	    print_intfc_states(file,os);
	    (void) fprintf(file,"\n\n");
	    (void) foutput(file);
	    (void) fprintf(file,"End of %s\n",os->name);
	    break;
	case CONSTANT_FIT | SMOOTH:
	default:
	    for (i = 0; i < dim; ++i)
	        gr->L[i] = cell_center(0,i,gr);
	    print_states_format(file,c_gr);
	    (void) fprintf(file,"constant smooth\n\n");
	    (*print_array_of_states[dim])(file,gr,os);
	    (void) fprintf(file,"\n");
	    break;

	}
}		/*end print_states1*/

LOCAL	void print_states_format_with_bdry(
	FILE		*file,
	INTERFACE	*intfc,
	RECT_GRID	*gr)
{
	RECT_GRID	rgr;
	int		i;

	copy_rect_grid(&rgr, gr);
	for(i=0; i<rgr.dim; i++)
	{
	    if(rect_boundary_type(intfc,i,0) == OPEN_BOUNDARY)
	    {
		rgr.L[i] = rgr.VL[i];
		rgr.gmax[i] += rgr.lbuf[i];
		rgr.lbuf[i] = 0;
	    }
	    if(rect_boundary_type(intfc,i,1) == OPEN_BOUNDARY)
	    {
		rgr.U[i] = rgr.VU[i];
		rgr.gmax[i] += rgr.ubuf[i];
		rgr.ubuf[i] = 0;
	    }
	}
	
	print_states_format(file,&rgr);
}

LOCAL	void print_states_format(
	FILE		*file,
	RECT_GRID	*gr)
{
	double	   *L = gr->L, *U = gr->U;
	int	   *gmax = gr->gmax, dim = gr->dim;
	int	   j;
	const char **Dnm = gr->Remap.Dnm;
	const char **dnm = gr->Remap.dnm;

	for (j = 0; j < dim; ++j)
	    (void) fprintf(file,"%sL = %-"FFMT" ",Dnm[j],L[j]);
	for (j = 0; j < dim; ++j)
	    (void) fprintf(file,"%sU = %-"FFMT" ",Dnm[j],U[j]);
	for (j = 0; j < dim; ++j)
	    (void) fprintf(file,"%smax = %-6d%s",dnm[j],gmax[j],
			(j==dim-1) ? "\n" : " ");
}		/*end print_states_format*/


/*ARGSUSED*/
LOCAL void print_array_of_states1d(
	FILE		*file,
	RECT_GRID	*gr,
	OUTPUT_SOLN	*os)
{

	OUTPUT_VALUE	*value;
	OUTPUT_VALUE	*(*solution)(OUTPUT_SOLN*,double*,int*) = os->solution;
	double		coords[MAXD];
	int		xmin, xmax;
	int		ix, icoords[MAXD];

	(void) foutput(file);
	(void) fprintf(file,"%-9s %-10s\n","POSITION",os->name);
	xmin = 0;	xmax = gr->gmax[0];
	if (debugging("buffer"))
	{
	    xmin -= gr->lbuf[0];	xmax += gr->ubuf[0];
	}
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",2*(xmax-xmin));
	    for (ix = xmin; ix < xmax; ++ix)
	    {
	    	icoords[0] = ix;
	    	coords[0] = cell_edge(ix,0,gr);
	    	value = (*solution)(os,coords,icoords);
	    	(void) fwrite((const void *)coords,sizeof(double),1,file);
		bin_prt_soln_value(file,value);
	    }
	}
	else
	{
	    for (ix = xmin; ix < xmax; ++ix)
	    {
	        icoords[0] = ix;
	        coords[0] = cell_edge(ix,0,gr);
	        value = (*solution)(os,coords,icoords);
	        (void) fprintf(file,ffmt,coords[0]," ");
		prt_soln_value(file,value);
	        (void) fprintf(file,"\n");
	    }
	}
}		/*end print_array_of_states1d*/


/*ARGSUSED*/
LOCAL void print_array_of_states2d(
	FILE		*file,
	RECT_GRID	*gr,
	OUTPUT_SOLN	*os)
{

	OUTPUT_VALUE	*(*solution)(OUTPUT_SOLN*,double*,int*) = os->solution;
	OUTPUT_VALUE	*value;
	double		coords[MAXD];
	int		xmax, ymax;
	int		xmin, ymin;
	int		ix,iy,icoords[MAXD];

	xmin = 0;	xmax = gr->gmax[0];
	ymin = 0;	ymax = gr->gmax[1];
	if (debugging("buffer"))
	{
	    xmin -= gr->lbuf[0];	xmax += gr->ubuf[0];
	    ymin -= gr->lbuf[1];	ymax += gr->ubuf[1];
	}
	if (is_binary_output() == YES)
	{
	    for (iy = ymax - 1; iy >= ymin; --iy)
	    {
	    	icoords[1] = iy;
	    	coords[1] = cell_edge(iy,1,gr);
	    	(void) fprintf(file,"\f%c",xmax-xmin);
	    	for (ix = xmin; ix < xmax; ++ix)
	    	{
	    	    icoords[0] = ix;
	    	    coords[0] = cell_edge(ix,0,gr);
	    	    value = (*solution)(os,coords,icoords);
		    bin_prt_soln_value(file,value);
	    	}
	    }
	}
	else
	{
	    for (iy = ymax - 1; iy >= ymin; --iy)
	    {
	        icoords[1] = iy;
	        coords[1] = cell_edge(iy,1,gr);
	        for (ix = xmin; ix < xmax; ++ix)
	        {
	            icoords[0] = ix;
	            coords[0] = cell_edge(ix,0,gr);
	    	    value = (*solution)(os,coords,icoords);
		    prt_soln_value(file,value);
	        }
	        (void) fprintf(file,"\n");
	    }
	}
}		/*end print_array_of_states2d*/


/*ARGSUSED*/
LOCAL void print_array_of_states3d(
	FILE		*file,
	RECT_GRID	*gr,
	OUTPUT_SOLN	*os)
{

	OUTPUT_VALUE	*(*solution)(OUTPUT_SOLN*,double*,int*) = os->solution;
	OUTPUT_VALUE	*value;
	double		coords[MAXD];
	int		xmin, xmax, ymin, ymax, zmin, zmax;
	int		ix,iy,iz,icoords[MAXD];
	int		*icrds_new;

	if (os->repart_at_end == YES)
	{
	    icrds_new = os->icrds_new;
	    xmin = gr->gmax[0]*icrds_new[0];
	    xmax = gr->gmax[0]*(icrds_new[0]+1);
	    ymin = gr->gmax[1]*icrds_new[1];
	    ymax = gr->gmax[1]*(icrds_new[1]+1);
	    zmin = gr->gmax[2]*icrds_new[2];
	    zmax = gr->gmax[2]*(icrds_new[2]+1);
	}
	else
	{
	    xmin = 0;	xmax = gr->gmax[0];
	    ymin = 0;	ymax = gr->gmax[1];
	    zmin = 0;	zmax = gr->gmax[2];
	}
	if (debugging("buffer"))
	{
	    xmin -= gr->lbuf[0];	xmax += gr->ubuf[0];
	    ymin -= gr->lbuf[1];	ymax += gr->ubuf[1];
	    zmin -= gr->lbuf[2];	zmax += gr->ubuf[2];
	}
	
	if(rect_boundary_type(os->intfc, 0, 0) == OPEN_BOUNDARY)
	    xmin = -gr->lbuf[0];
	if(rect_boundary_type(os->intfc, 1, 0) == OPEN_BOUNDARY)
	    ymin = -gr->lbuf[1];
	if(rect_boundary_type(os->intfc, 2, 0) == OPEN_BOUNDARY)
	    zmin = -gr->lbuf[2];
	if(rect_boundary_type(os->intfc, 0, 1) == OPEN_BOUNDARY)
	    xmax = gr->gmax[0] + gr->ubuf[0];
	if(rect_boundary_type(os->intfc, 1, 1) == OPEN_BOUNDARY)
	    ymax = gr->gmax[1] + gr->ubuf[1];
	if(rect_boundary_type(os->intfc, 2, 1) == OPEN_BOUNDARY)
	    zmax = gr->gmax[2] + gr->ubuf[2];

	if (is_binary_output() == YES)
	{
	    for (iz = zmin; iz < zmax; ++iz)
	    {
	    	icoords[2] = iz;
	    	coords[2] = cell_edge(iz,2,gr);
	    	for (iy = ymax - 1; iy >= ymin; --iy)
	    	{
	    	    icoords[1] = iy;
	    	    coords[1] = cell_edge(iy,1,gr);
	    	    (void) fprintf(file,"\f%c",xmax-xmin);
	    	    for (ix = xmin; ix < xmax; ++ix)
	    	    {
	    		icoords[0] = ix;
	    		coords[0] = cell_edge(ix,0,gr);
	    		value = (*solution)(os,coords,icoords);
		        bin_prt_soln_value(file,value);
	    	    }
	    	}
	    }
	}
	else
	{
	    for (iz = zmin; iz < zmax; ++iz)
	    {
	        icoords[2] = iz;
	        coords[2] = cell_edge(iz,2,gr);
	        for (iy = ymax - 1; iy >= ymin; --iy)
	        {
	            icoords[1] = iy;
	            coords[1] = cell_edge(iy,1,gr);
	            for (ix = xmin; ix < xmax; ++ix)
	            {
	                icoords[0] = ix;
	                coords[0] = cell_edge(ix,0,gr);
	                value = (*solution)(os,coords,icoords);
			prt_soln_value(file,value);
	            }
	            (void) fprintf(file,"\n");
	        }
	        (void) fprintf(file,"\n\n");
	    }
	}
}		/*end print_array_of_states3d*/


/*ARGSUSED*/
LOCAL void print_intfc_states(
	FILE		*file,
	OUTPUT_SOLN	*os)
{
	HYPER_SURF_ELEMENT   *hse;
	HYPER_SURF	     *hs, *hslast;
	INTERFACE	     *intfc;
	OUTPUT_VALUE	     left,right;
	POINT		     *p;
	int		     i, dim;
	const char	     *hsname;
	static const char           *dhsname[] = {"point","curve","surface"};
	static uint64_t (*hs_identifier)(HYPER_SURF*);
	static uint64_t (*hs_identifier_choices[3])(HYPER_SURF*) = {
					hs_point_identifier,
					hs_curve_identifier,
					hs_surface_identifier};

	if ((intfc=os->intfc) == NULL)
	    return;
	dim = intfc->dim;
	hsname = dhsname[dim-1];
	hs_identifier = hs_identifier_choices[dim-1];
	(void) fprintf(file,"\nStart of intfc states");
	if (is_binary_output() == YES)
	{
	    hslast = NULL;
	    /*reset point sort status */
	    (void) next_point(intfc,NULL,NULL,NULL);
	    while (next_point(intfc,&p,&hse,&hs))
	    {
	    	if (hs != hslast)
	    	{
	    	    hslast = hs;
		    (void) fprintf(file,"\n%s %llu\n",hsname,
				   (*hs_identifier)(hs));
		}
		(*os->intfc_solution)(os,p,hse,hs,&left,&right);
		bin_prt_soln_values(file,&left,&right);
	    }
	}
	else
	{
	    (void) fprintf(file,"\n");
	    hslast = NULL;
	    (void) next_point(intfc,NULL,NULL,NULL);
	    while (next_point(intfc,&p,&hse,&hs))
	    {
	        if (hs != hslast)
	        {
	            hslast = hs;
	            (void) fprintf(file,"%s %llu\n",hsname,
		                   (*hs_identifier)(hs));
	        }
	        (*os->intfc_solution)(os,p,hse,hs,&left,&right);
	        for (i = 0; i < dim; ++i)
	            (void) fprintf(file,ffmt,Coords(p)[i]," ");
	        prt_soln_values(file,&left,&right,"\n");
	    }
	}
	(void) fprintf(file,"End of intfc states\n\n");
}		/*end print_intfc_states*/


/*ARGSUSED*/
LOCAL	uint64_t	hs_point_identifier(
	HYPER_SURF	*hs)
{
	POINT	*p = NULL;
	if (hs != NULL)
	    p =  Point_of_hs(hs);
	return point_number(p);
}		/*end hs_point_identifier*/

/*ARGSUSED*/
LOCAL	uint64_t	hs_curve_identifier(
	HYPER_SURF	*hs)
{
	CURVE	*c = NULL;
	if (hs != NULL)
	    c = Curve_of_hs(hs);
	return curve_number(c);
	return INT_MAX;
}		/*end hs_curve_identifier*/

/*ARGSUSED*/
LOCAL	uint64_t	hs_surface_identifier(
	HYPER_SURF	*hs)
{
	SURFACE	*s = (hs != NULL) ? Surface_of_hs(hs) : NULL;
	return surface_number(s);
	return INT_MAX;
}		/*end hs_surface_identifier*/


LOCAL int read_version;


EXPORT void determine_read_version(
	FILE		*file)
{
	static OUTPUT *oput = NULL;
	const char    *line;

	oput = save_read_file_variables(file,oput);
	rewind_read_file(file,NULL);
	if ((line = next_output_line_containing_string(file,"print version")))
	{
	    (void) sscanf(line,"%*s%*s%d",&read_version);
	    if (read_version < 0 || read_version > 1)
	    	read_version = 0;
	}
	else
	    read_version = 0;
	reset_read_file_variables(oput);

	return;
}		/*end determine_read_version*/

EXPORT	double is_state(
	INPUT_SOLN	*is,
	int		*icrds)
{
	switch (is->grid.dim)
	{
	case 1:
	    return is->states1d[icrds[0]];
	case 2:
	    return is->states2d[icrds[1]][icrds[0]];
	case 3:
	    return is->states3d[icrds[2]][icrds[1]][icrds[0]];
	}
	return ERROR_FLOAT;
}		/*end is_state*/
	

EXPORT	boolean read_state_variables(
	const IO_TYPE	*io_type,
	int		n_vars,
	INTERFACE	*restart_intfc,
	INPUT_SOLN	**is,
	int		dim)
{
	int  var;

	debug_print("restart","Entered read_state_variables()\n");

	if ((dim == 1) &&
	    read_state_variables1d(io_type,n_vars,restart_intfc,is))
		return FUNCTION_SUCCEEDED;

	for (var = 0; var < n_vars; ++var)
	{
	    if (!read_states(io_type,restart_intfc,var,is[var],dim)) 
	    {
	    	(void) printf("WARNING in read_state_variables(), "
		              "read_states() failed\n");
	    	return FUNCTION_FAILED;
	    }
	}

	debug_print("restart","Left read_state_variables()\n");
	return FUNCTION_SUCCEEDED;
}		/*end read_state_variables*/

/*ARGSUSED*/
LOCAL	boolean read_state_variables1d(
	const IO_TYPE *io_type,
	int	      n_vars,
	INTERFACE     *restart_intfc,
	INPUT_SOLN    **is)
{
	FILE	          *file = io_type->file;
	INTERFACE         *intfc;
	RECT_GRID         *cgr = computational_grid(restart_intfc);
	POINT	          *p;
	const char        *line;
	char	          Line[2048];
	char              *name;
	double	          x;
	double	          *states1d;
	double	          sl, sr;
	double             dummy;
	int               np_vars;
	int	          num_points;
	int	          i, ix;
	int	          var;
	int	          ch;
	int	          xmax;
	static const char *blanks = " \t\n";

	line = next_output_line_containing_string(file,
				"Printout of one dimensional states");
	if (line == NULL)
	{
	    (void) printf("WARNING in read_state_variables1d(), states "
	                  "printout for Printout of one dimensional states "
			  "not found\n");
	    return FUNCTION_FAILED;
	}
	read_grid_fit_and_smoothness(file,&cgr->Remap,&is[0]->grid,
	                             &is[0]->fit,&is[0]->smoothness,1);
	if ((is[0]->fit | is[0]->smoothness) != (LINEAR_FIT | SINGULAR))
	{
	    (void) printf("WARNING in read_state_variables1d(), "
	                  "only linear-singular printout supported\n");
	    return FUNCTION_FAILED;
	}
	for (var = 1; var < n_vars; ++var)
	{
	    set_rect_grid(is[0]->grid.L,is[0]->grid.U,is[0]->grid.L,
			  is[0]->grid.U,NOBUF,NOBUF,is[0]->grid.gmax,
			  1,&computational_grid(restart_intfc)->Remap,
			  &is[var]->grid);
	    is[var]->fit = is[0]->fit;
	    is[var]->smoothness = is[0]->smoothness;
	}
	for (var = 0; var < n_vars; ++var)
	{
	    if (allocate_input_soln(restart_intfc,is[var],1) == FUNCTION_FAILED)
	    {
	    	(void) printf("WARNING in read_state_variables1d(), "
	    	              "no more storage for input solutions\n");
	    	return FUNCTION_FAILED;
	    }
	}

	line = next_output_line_containing_string(file,"POSITION");
	if (line == NULL)
	{
	    (void) printf("WARNING in read_state_variables1d(), "
	                  "POSITION not found\n");
	    return FUNCTION_FAILED;
	}

	/* Count number of printed variables */
	(void) strcpy(Line,line);
	(void) strtok(Line,blanks);
	for (np_vars = 0, name = strtok(NULL,blanks); name != NULL;
	     name = strtok(NULL,blanks), ++np_vars);

	xmax = is[0]->grid.gmax[0];
	ch = getc(file);
	(void) ungetc(ch,file);
	if (ch != '\f') 	/* NOBINARY */
	{
	    for (ix = 0; ix < xmax; ++ix)
	    {
	        (void) fscan_float(file,&x);
	        for (var = 0; var < n_vars; ++var)
		{	        
	            states1d = is[var]->states1d;
	            (void) fscan_float(file,&states1d[ix]);
	        }
		while ((ch = getc(file)) != '\n');/*Clear line*/
	    }
	}
	else 			/* BINARY */
	{
	    for (ix = 0; ix < xmax; ++ix)
	    {
		(void) getc(file);		/* "\f" */
	        (void) getc(file);		/* "%c" */
		(void) read_binary_real_array(&x,1,io_type);
	        for (var = 0; var < n_vars; ++var)
	        {
	            states1d = is[var]->states1d;
		    (void) read_binary_real_array(&states1d[ix],1,io_type);
	        }
		for (; var < np_vars; ++var)
		    (void) read_binary_real_array(&dummy,1,io_type);
	    }
	}

	if (restart_intfc == NULL)
	    return FUNCTION_SUCCEEDED;

	for (var = 0; var < n_vars; ++var)
	{
	    intfc = is[var]->intfc;
	    if (intfc == NULL)
	    {
	    	(void) printf("WARNING in read_state_variables1d(), "
			      "NULL intfc\n");
	    	return FUNCTION_FAILED;
	    }
	    if (size_of_state(intfc) < 1)
	    {
	    	(void) printf("WARNING in read_state_variables1d(), "
	    	              "size_of_state(intfc) < 1\n");
	    	return FUNCTION_FAILED;
	    }
	    if (is[var]->set_intfc_states == NULL)
	    {
	    	(void) printf("WARNING in read_state_variables1d(), "
	    	              "is[%d]->set_intfc_states == NULL\n",var);
	    	return FUNCTION_FAILED;
	    }
	}
	if (fgetstring(file,"Start of intfc states") == FUNCTION_FAILED)
	{
	    (void) printf("WARNING in read_state_variables1d(), "
	                  "intfc states not supplied\n");
	    return FUNCTION_FAILED;
	}
	(void) fgets(Line,2046,file);/*Clear line*/
	(void) fgets(Line,2046,file);/*Get header line*/
	Line[strlen(Line)-1] = '\0';
	for (np_vars = 0, name = strtok(Line,blanks); name != NULL;
	     name = strtok(NULL,blanks), ++np_vars);
	num_points = restart_intfc->num_points;
	ch = getc(file);
	(void) ungetc(ch,file);
	if (ch != '\f') 	/* NOBINARY */
	{
	    for (i = 0; i < num_points; ++i)
	    {
		(void) fscan_float(file,&x);
		for (var = 0; var < n_vars; ++var)
		{
		    p = is[var]->intfc->points[i];
		    (void) fscan_float(file,&sl);
		    (void) fscan_float(file,&sr);
		    (*is[var]->set_intfc_states)(&sl,&sr,var,p,NULL,
						 Hyper_surf(p));
		}
		(void) fgets(Line,2046,file);/*Clear line*/
	    }
	}
	else			/* BINARY */
	{
	    int nread;
	    for (i = 0; i < num_points; ++i)
	    {
		(void) getc(file);		/* "\f" */
		(void) getc(file);		/* "%c" */
		(void) read_binary_real_array(&x,1,io_type);
		nread = 1;
		for (var = 0; var < n_vars; ++var)
		{
		    (void) read_binary_real_array(&sl,1,io_type);
		    (void) read_binary_real_array(&sr,1,io_type);
		    nread += 2;
		    p = is[var]->intfc->points[i];
		    (*is[var]->set_intfc_states)(&sl,&sr,var,p,NULL,
						 Hyper_surf(p));
		}
		for (; nread < np_vars; ++nread)
		    (void) read_binary_real_array(&dummy,1,io_type);
	    }
	}
	(void) fgetstring(file,"End of intfc states\n");
	return FUNCTION_SUCCEEDED;
}		/*end read_state_variables1d*/


LOCAL	boolean read_states(
	const IO_TYPE	*io_type,
	INTERFACE	*restart_intfc,
	int		var,
	INPUT_SOLN	*input_soln,
	int		dim)
{
	FILE	   *file = io_type->file;
	RECT_GRID  *cgr = computational_grid(restart_intfc);
	const char *line;

	debug_print("restart","Entered read_states()\n");

	line = next_output_line_containing_string(file,input_soln->name);
	if (line == NULL)
	{
	    (void) printf("WARNING in read_states(), states printout for "
	                  "%s not found\n",input_soln->name);
	    return FUNCTION_FAILED;
	}

	switch (read_version)
	{
	default:
	case 0:
	    /* input_soln->fit and input_soln->smoothness supplied */
	    set_grid(line,restart_intfc,input_soln->fit,
		     input_soln->smoothness,&input_soln->grid,dim);
	    break;
	case 1:
	    read_grid_fit_and_smoothness(file,&cgr->Remap,&input_soln->grid,
			                 &input_soln->fit,
					 &input_soln->smoothness,dim);
	    break;
	}

	if (allocate_input_soln(restart_intfc,input_soln,dim) ==
							FUNCTION_FAILED)
	{
	    (void) printf("WARNING in read_states(), "
	                  "no more storage for input solutions\n");
	    return FUNCTION_FAILED;
	}
	if (read_printout_of_states(io_type,var,input_soln) == FUNCTION_FAILED)
	{
	    (void) printf("WARNING in read_states(), "
	                  "read_printout_of_states failed\n");
	    return FUNCTION_FAILED;
	}

	debug_print("restart","Left read_states()\n");
	return FUNCTION_SUCCEEDED;
}		/*end read_states*/



LOCAL	void set_grid(
	const char *line,
	INTERFACE  *intfc,
	int	   fit,
	int	   smoothness,
	RECT_GRID  *gr,
	int	   dim)
{
	int		i;

	if (intfc != NULL)
	    copy_rect_grid(gr,computational_grid(intfc));
	else
	{
	    static	double L[3] = {0.0, 0.0, 0.0}, U[3] = {1.0, 1.0, 1.0};
	    switch (dim)
	    {
	    case 1:
	        (void) sscanf(line,"%*s%d",&gr->gmax[0]);
	        break;
	    case 2:
	        (void) sscanf(line,"%*s%d%d",&gr->gmax[1],&gr->gmax[0]);
	        break;
	    case 3:
	        (void) sscanf(line,"%*s%d%d%d",
			      &gr->gmax[2],&gr->gmax[1],&gr->gmax[0]);
	        break;
	    }
#if defined(STATES_AT_NODES)
	    if (fit == LINEAR_FIT)
	    {
#else /* defined(STATES_AT_NODES) */
	    if (fit == LINEAR_FIT && smoothness == SMOOTH)
	    {
#endif /* defined(STATES_AT_NODES) */
	    	for (i = 0; i < dim; ++i)
		    --gr->gmax[i];
	    }
	    set_rect_grid(L,U,L,U,NOBUF,NOBUF,gr->gmax,dim,remap_info(),gr);
	}
}		/*end set_grid*/

LOCAL const char *READ_STATES_FORMAT1 = "%*s%*s%lf%*s%*s%lf%*s%*s%d";
LOCAL const char *READ_STATES_FORMAT2 =
    "%*s%*s%lf%*s%*s%lf%*s%*s%lf%*s%*s%lf%*s%*s%d%*s%*s%d";
LOCAL const char *READ_STATES_FORMAT3 =
    "%*s%*s%lf%*s%*s%lf%*s%*s%lf%*s%*s%lf"
    "%*s%*s%lf%*s%*s%lf%*s%*s%d%*s%*s%d%*s%*s%d";

LOCAL	void read_grid_fit_and_smoothness(
	FILE	  *file,
	REMAP     *remap,
	RECT_GRID *gr,
	int	  *fit,
	int	  *smoothness,
	int	  dim)
{
	char	Line[2048], s_fit[120], s_smoothness[120];
	int	i;

	gr->dim = dim;
	(void) fgets(Line,2046,file);
	switch (dim)
	{
	case 1:
	    (void) sscanf(Line,READ_STATES_FORMAT1,&gr->L[0],&gr->U[0],
			  &gr->gmax[0]);
		break;
	case 2:
	    (void) sscanf(Line,READ_STATES_FORMAT2,&gr->L[0],&gr->L[1],
			  &gr->U[0],&gr->U[1],&gr->gmax[0],&gr->gmax[1]);
		break;
	case 3:
	    (void) sscanf(Line,READ_STATES_FORMAT3,&gr->L[0],&gr->L[1],
			  &gr->L[2],&gr->U[0],&gr->U[1],&gr->U[2],
			  &gr->gmax[0],&gr->gmax[1],&gr->gmax[2]);
	    break;
	}

	(void) fgets(Line,2046,file);
	(void) sscanf(Line,"%s%s",s_fit,s_smoothness);

	switch (s_fit[0])
	{
	default:
	case 'c':
		*fit = CONSTANT_FIT;
		break;
	case 'l':
		*fit = LINEAR_FIT;
		break;
	}

	switch (s_smoothness[1])
	{
	default:
	case 'm':
		*smoothness = SMOOTH;
		break;
	case 'i':
		*smoothness = SINGULAR;
		break;
	}

	if (*fit == LINEAR_FIT)
	{
#if !defined(STATES_AT_NODES)
	    if (*smoothness == SMOOTH)
#endif /* !defined(STATES_AT_NODES) */
	    	for (i = 0; i < dim; ++i)
		    gr->gmax[i]--;
	}
	set_rect_grid(gr->L,gr->U,gr->L,gr->U,NOBUF,NOBUF,gr->gmax,dim,
		      remap,gr);
}		/*end read_grid_fit_and_smoothness*/



LOCAL boolean allocate_input_soln(
	INTERFACE	*restart_intfc,
	INPUT_SOLN	*is,
	int		dim)
{
	int		i;
	int		gmax[MAXD];

	for (i = 0; i < dim; ++i)
	    gmax[i] = is->grid.gmax[i];

	if (is->fit == LINEAR_FIT) 
	{
	    for (i = 0; i < dim; ++i)
	        ++gmax[i];
	}

	is->states1d = NULL;
	is->states2d = NULL;
	is->states3d = NULL;
	switch(dim)
	{
	case 1:

	    uni_array(&is->states1d,gmax[0],FLOAT);
	    if (is->states1d == NULL)
	        return FUNCTION_FAILED;

	    break;
	case 2:

	    bi_array(&is->states2d,gmax[1],gmax[0],FLOAT);
	    if (is->states2d == NULL)
	        return FUNCTION_FAILED;

	    break;
	case 3:

	    tri_array(&is->states3d,gmax[2],gmax[1],gmax[0],FLOAT);
	    if (is->states3d == NULL)
	        return FUNCTION_FAILED;

	    break;
	}

	if (is->smoothness == SMOOTH)
	    return FUNCTION_SUCCEEDED;

	is->intfc = restart_intfc;

	return FUNCTION_SUCCEEDED;
}		/*end allocate_input_soln*/



EXPORT void free_input_soln(
	INPUT_SOLN	*is)
{
	if (is->states1d != NULL) { free(is->states1d); is->states1d = NULL; }
	if (is->states2d != NULL) { free(is->states2d); is->states2d = NULL; }
	if (is->states3d != NULL) { free(is->states3d); is->states3d = NULL; }
	if (is->intfc  != NULL) is->intfc = NULL;
}		/*end free_input_soln*/


LOCAL boolean read_printout_of_states(
	const IO_TYPE *io_type,
	int	      var,
	INPUT_SOLN    *is)
{
	RECT_GRID	*gr = &is->grid;
	int		dim = gr->dim;
	int		i, gmax[MAXD];

	switch (is->fit | is->smoothness)
	{
	default:
	case CONSTANT_FIT | SMOOTH:
	    (*read_array_of_states[dim])(io_type,gr->gmax,is);
		break;
	case LINEAR_FIT | SMOOTH:
	    for (i = 0; i < dim; ++i)
		gmax[i] = gr->gmax[i] + 1;
	    (*read_array_of_states[dim])(io_type,gmax,is);
		break;
	case LINEAR_FIT | SINGULAR:
#if defined(STATES_AT_NODES)
	    for (i = 0; i < dim; ++i)
		gmax[i] = gr->gmax[i] + 1;
	    (*read_array_of_states[dim])(io_type,gmax[0],is);
#else /* defined(STATES_AT_NODES) */
	    (*read_array_of_states[dim])(io_type,gr->gmax,is);
#endif /* defined(STATES_AT_NODES) */
	    if (read_intfc_states(io_type,var,is) == FUNCTION_FAILED)
	    {
	    	(void) printf("WARNING in read_printout_of_states(), "
	    	              "read_intfc_states() failed\n");
	    	return FUNCTION_FAILED;
	    }
	    break;
	}

	return FUNCTION_SUCCEEDED;
}		/*end read_printout_of_states*/


/*ARGSUSED*/
LOCAL void read_array_of_states1d(
	const IO_TYPE *io_type,
	int	      *gmax,
	INPUT_SOLN    *is)
{
	FILE		*file = io_type->file;
	double		*states1d = is->states1d;
	double		x;
	int		ch;
	int		ix;
	int		xmax;

	if (next_output_line_containing_string(file,is->name) == NULL)
	{
	    screen("ERROR in read_array_of_states1d(), "
	           "state printout for %s not found\n",is->name);
	    clean_up(ERROR);
	}

	xmax = gmax[0];
	while ((ch = getc(file)) == '\n')
		;
	(void) ungetc(ch,file);
	if (ch != '\f') 	/* NOBINARY */
	{
	    for (ix = 0; ix < xmax; ++ix)
	    {
	    	(void) fscan_float(file,&x);
	    	(void) fscan_float(file,&states1d[ix]);
	    }
	}
	else 			/* BINARY */
	{
	    (void) getc(file);
	    (void) getc(file);		/* "\f%c" */
	    for (ix = 0; ix < xmax; ++ix)
	    {
		(void) read_binary_real_array(&x,1,io_type);
		(void) read_binary_real_array(&states1d[ix],1,io_type);
	    }
	}
}		/*end read_array_of_states1d*/


/*ARGSUSED*/
LOCAL	void read_array_of_states2d(
	const IO_TYPE *io_type,
	int	      *gmax,
	INPUT_SOLN    *is)
{
	FILE		*file = io_type->file;
	double		**states2d = is->states2d;
	int		ch;
	int		ix,iy;
	int		xmax, ymax;

	xmax = gmax[0];	ymax = gmax[1];
	while ((ch = getc(file)) == '\n')
		;
	(void) ungetc(ch,file);

	if (ch != '\f') 	/* NOBINARY */
	{
	    for (iy = ymax - 1; iy >= 0; --iy)
	    {
	    	for (ix = 0; ix < xmax; ++ix)
	    	{
	    	    (void) fscan_float(file,&states2d[iy][ix]);
	    	}
	    }
	}
	else 			/* BINARY */
	{
	    for (iy = ymax - 1; iy >= 0; --iy)
	    {
	    	(void) getc(file);
	    	(void) getc(file);	/* "\f%c" */
		(void) read_binary_real_array(states2d[iy],xmax,io_type);
	    }
	}
}		/*end read_array_of_states2d*/

/*ARGSUSED*/
LOCAL	void read_array_of_states3d(
	const IO_TYPE *io_type,
	int	      *gmax,
	INPUT_SOLN    *is)
{
	FILE		*file = io_type->file;
	double		***states3d = is->states3d;
	int		ch;
	int		ix,iy,iz;
	int		xmax, ymax, zmax;

	xmax = gmax[0];	ymax = gmax[1]; zmax = gmax[2];
	while ((ch = getc(file)) == '\n')
		;
	(void) ungetc(ch,file);
	if (ch != '\f') 	/* NOBINARY */
	{
	    for (iz = 0; iz < zmax; ++iz)
	    {
	        for (iy = ymax - 1; iy >= 0; --iy)
	        {
	            for (ix = 0; ix < xmax; ++ix)
	    	    {
			(void) fscan_float(file,&states3d[iz][iy][ix]);
		    }
		}
	    }
	}
	else 			/* BINARY */
	{
	    for (iz = 0; iz < zmax; ++iz)
	    {
	        for (iy = ymax - 1; iy >= 0; --iy)
	        {
	            /* "\f%c" */
	            (void) getc(file);
	            (void) getc(file);
		    (void) read_binary_real_array(states3d[iz][iy],
		                                  xmax,io_type);
	        }
	    }
	}
}		/*end read_array_of_states3d*/


#define fscan_states(file,sl,sr) (void) fscanf(file,"%lf%lf",(sl),(sr))

/*ARGSUSED*/
LOCAL boolean read_intfc_states(
	const IO_TYPE *io_type,
	int	      var,
	INPUT_SOLN    *is)
{
	FILE		   *file = io_type->file;
	INTERFACE	   *intfc = is->intfc;
	char		   Line[2048];
	int		   i, dim = intfc->dim, ch;
	boolean		   binary = NO;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF	   *hs, *hslast;
	POINT		   *p;
	double		   sl, sr;
	char		   hslabel;
	static const char  *dhslabel = "pcs";

	hslabel = dhslabel[intfc->dim-1];
	if (intfc == NULL || size_of_state(intfc) < 1)
	{
	    (void) printf("WARNING in read_intfc_states(), NULL intfc\n");
	    return FUNCTION_FAILED;
	}
	if (is->set_intfc_states == NULL)
	{
	    (void) printf("WARNING in read_intfc_states(), "
	                  "is->set_intfc_states == NULL\n");
	    return FUNCTION_FAILED;
	}
	if (fgetstring(file,"Start of intfc states") == FUNCTION_FAILED)
	{
	    (void) printf("WARNING in read_intfc_states(), "
	                  "intfc states not supplied\n");
	    return FUNCTION_FAILED;
	}

	hslast = NULL;
	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (hs != hslast)
	    {
	    	hslast = hs;
	       	(void) fgets(Line,2046,file);/*Clear line*/
	    	(void) fgets(Line,2046,file);
	    	if (Line[0] != hslabel)
		{
		    (void) printf("WARNING in read_intfc_states(), "
		                  "Line[0] != hslabel\n");
		    (void) printf("hslabel = %c, Line = %s\n",hslabel,Line);
		    return FUNCTION_FAILED;
	    	}
	    	binary = ((ch = getc(file)) == '\f') ? YES : NO;
		(void) ungetc(ch,file);
	    }
	    if (binary == NO)
	    {
	    	for (i = 0; i < dim; ++i)
		    (void) fscanf(file,"%*f");
	    	fscan_states(file,&sl,&sr);
		(*is->set_intfc_states)(&sl,&sr,var,p,hse,hs);
	    }
	    else
	    {
	    	(void) getc(file);
	    	(void) getc(file);
	    	(void) read_binary_real_array(&sl,1,io_type);
	    	(void) read_binary_real_array(&sr,1,io_type);
	    	(*is->set_intfc_states)(&sl,&sr,var,p,hse,hs);
	    }
	}
	(void) fgetstring(file,"End of intfc states\n");

	return FUNCTION_SUCCEEDED;
}		/*end read_intfc_states*/
#undef fscan_states



EXPORT	void print_INPUT_SOLN_structure(
	INPUT_SOLN	*in)
{
	(void) printf("\n\n\n\t\tINPUT_SOLN %p structure\n",(POINTER)in);
	if (in == NULL)
	{
	    (void) printf("\t\tstructure not yet allocated\n");
	    (void) printf("\n\t\tEnd INPUT_SOLN %p structure\n\n",(POINTER)in);
	    return;
	}

	(void) printf("name %s\tfit %d\tsmoothness %d\n",in->name,in->fit,
		      in->smoothness);
	(void) printf("\t[fits: CONSTANT_FIT %d LINEAR_FIT %d]\n",
		      CONSTANT_FIT,LINEAR_FIT);
	(void) printf("\t[smoothness: SMOOTH %d SINGULAR %d]\n",
		      SMOOTH,SINGULAR);
	(void) printf("states1d %p\n",in->states1d);
	(void) printf("states2d %p\n",in->states2d);
	(void) printf("states3d %p\n",in->states3d);
	(void) printf("intfc %llu (SINGULAR smoothness only)\n",
		      interface_number(in->intfc));
	(void) printf("set_intfc_states  %p\n",in->set_intfc_states);
	print_RECT_GRID_structure(&in->grid);

	(void) printf("\n\t\tEnd INPUT_SOLN %p structure\n\n",(POINTER)in);
}		/*end print_INPUT_SOLN_structure*/


EXPORT	void print_OUTPUT_SOLN_structure(
	OUTPUT_SOLN	*out)
{
	(void) printf("\n\n\n\t\tOUTPUT_SOLN %p structure\n",(POINTER)out);
	if (out == NULL)
	{
	   (void) printf("\t\tstructure not yet allocated\n");
	   (void) printf("\n\t\tEnd OUTPUT_SOLN %p structure\n\n",(POINTER)out);
	   return;
	}

	(void) printf("name= %s\tfit= %d\tsmoothness= %d\n",out->name,out->fit,
		      out->smoothness);
	(void) printf("\t[fits: CONSTANT_FIT= %d LINEAR_FIT= %d]\n",
		      CONSTANT_FIT,LINEAR_FIT);
	(void) printf("\t[smoothness: SMOOTH= %d SINGULAR= %d]\n",
		      SMOOTH,SINGULAR);
	(void) printf("solution var= %d\t\t",out->var);
	(void) printf("solution()= %p\n",out->solution);
	(void) printf("intfc_solution()= %p (LINEAR_FIT,SINGULAR only)\n",
		      out->intfc_solution);
	(void) printf("intfc= %llu (SINGULAR smoothness only)\n",
		      interface_number(out->intfc));
	(void) printf("RECT_GRID below needed only if intfc == NULL\n");
	print_RECT_GRID_structure(out->grid);
	(void) printf("extra[]= %p\n",(POINTER)out->extra);

	(void) printf("\n\t\tEnd OUTPUT_SOLN %p structure\n\n",(POINTER)out);
}		/*end print_OUTPUT_SOLN_structure*/
