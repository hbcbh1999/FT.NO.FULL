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
*				gitabregion.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Enables state lookup from tables.  The tables correspond to a
*	rectangular region R with upper and lower corners given by the 
*	input coordinates (UX, UY, UZ) and (LX, LY, LZ).  Data in the
*	table correspond to cell centers of a regular rectangle lattice
*	over R with mesh spacings dx = (UX - LX)/nx, dy = (UY - LY)/ny
*	and dz = (UZ - LZ)/nz where nx, ny, nz are the number of mesh zones
*	in the respective coordinate directions.  Data in the table
*	is stored as a linear array of the form
*
*			val(i,j,k) = entry(i+nx*(j+ny*k))
*
*	where i,j, and k are the the indices in the x, y, and z directions
*	respectively. Geometrically we can interpret the data as a collection
*	of nz nx X ny matrices where the rows correspond to values in the
*	x direction,  columns to values in the y direction,  and successive
*	matrices to data in planes parallel to the x-y plane.  Note the
*	where viewed as a matrix the y coordinates increase with increasing
*	row number so the the data would appear in the reverse of the normal
*	plane view is printed directions,  ie the data layout looks like:
*
*				------> x
*				|
*				|               in plane 0
*				|
*			       \ /
*				y
*
*				------> x
*				|
*				|               in plane 1
*				|
*			       \ /
*				y
*
*				.
*				.
*				.
*	
*	For one and two dimensional runs the layout reduces in the natural
*	way,  for example two dimensional data looks like a single plane
*	and one dimension data a single row.  The tables may be printed in
*	either ASCII characters corresponding to decimal numbers separated
*	by white space (spaces, tabs, or newlines) or in unformated binary.
*	The binary format is equivalent to that produced by the C code
*	fragment:
*
*		REAL    *data;
*		int N = nx*ny*nz;
*		(void) fprintf(file,"\f%c",N);
*		(void) fwrite((const void*)data,sizeof(REAL),N,file);
*
*	where REAL is either double or double.  Note that the code currently
*	only supports binary reads for formats that contain the same precision
*	as the running code,  so if you compiled FronTier with
*	PRECISION=double,  the table file must have been printed using doubles.
*	Similarlly for single precion.  Also note that endian restrictions
*	apply so that binary files are not portable from big to little endian
*	machines or the reverse.  ASCII file formats are portable to all
*	machines but may contain round-off error.
*
*	Once the data has been read it is interpolated bilinearly to yield
*	a continuous flow field on the specificed region.  This in turn
*	is used to initialize the FronTier states overlapping that region.
*	Thus there is no need for the table grid sizes to agree with the
*	computational grid being used in FronTier.
*
*	SUGGESTED UPGRADE:
*		Use HDF for binary input to remove precision and
*		portability restrictions.
*/

#include <ginit/ginit.h>

/* LOCAL Function Prototypes */
LOCAL	void	read_table_data(const char*,double*,int);
LOCAL	void	free_tabulated_region_comp_type(COMP_TYPE*);
LOCAL	void	get_state_tabulated_region(double*,Locstate,COMP_TYPE*,
                                           HYPER_SURF*,INTERFACE*,INIT_DATA*,
					   int);

/*ARGSUSED*/
EXPORT	void	set_tabulated_region_comp_type(
	COMP_TYPE	*ctype,
	Front		*front)
{
	const COMPONENT       *comps;
	TABULATED_REGION_DATA *table = NULL;
	int                   i, N;

	if (ctype->type == TABULATED_REGION) /*ALREADY SET*/
	    return;

	if (ctype->free_comp_type_extra != NULL)
	    (*ctype->free_comp_type_extra)(ctype);

	comps = comps_with_type(TABULATED_REGION,&N);
	if (comps != NULL)
	{
	    COMPONENT comp;
	    char      s[Gets_BUF_SIZE];

	    screen("%d component%s, comp%s",N,(N>1)?"s":"",(N>1)?"s":"");
	    for (i = 0; i < N; ++i)
	        screen(" %d",comps[i]);
	    screen(", %s type TABULATED_REGION\n",(N>1)?"have":"has");
	    screen("Enter a component number to share table data: ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	        if (sscanf(s,"%d",&comp) == 1)
		{
		    for (i = 0; i < N; ++i)
		    {
		        if (comps[i] == comp)
			{
			    COMP_TYPE *ct = comp_type(comp);
	                    table = Tabulated_region_data(ct);
			    screen("Component %d shares component %d's "
			           "TABULATED_REGION_DATA table\n",
			           ctype->comp,comp);
			    break;
			}
		    }
		}
	    }
	}
	
	ctype->type = TABULATED_REGION;
	if (table == NULL)
	{
	    scalar(&table,sizeof(TABULATED_REGION_DATA));
	    table->allocated_from = ctype;
	    table->tab_region_set = NO;
	}
	ctype->extra = (POINTER)table;

	ctype->_get_state = get_state_tabulated_region;
	ctype->free_comp_type_extra = free_tabulated_region_comp_type;

}		/*end set_tabulated_region_comp_type*/

EXPORT void set_up_read_state_from_file_region(
	COMP_TYPE  *ctype,
	const char *message,
	Front	   *front)
{       
	TABULATED_REGION_DATA *table;
	RECT_GRID             *gr = front->rect_grid;
	Gas_param             *params = ctype->params;
	char                  s[Gets_BUF_SIZE];
	int                   i, dim = front->rect_grid->dim;
	int                   len;
	double                 min_pressure = params->min_pressure;
	double                 vacuum_dens = params->vacuum_dens;
	const char            **dnm = gr->Remap.dnm;
	const char            **Dnm = gr->Remap.Dnm;

	table = Tabulated_region_data(ctype);
	if (table->tab_region_set)
	    return;

	for (i = 0; i < dim; ++i)
	{
	    table->GL[i] = gr->GL[i];
	    table->GU[i] = gr->GU[i];
	    table->n[i] = gr->gmax[i];
	    table->h[i] = gr->h[i];
	}
	for (i = dim; i < 3; ++i)
	{
	    table->n[i] = 1;
	    table->GL[i] =  HUGE_VAL;
	    table->GU[i] = -HUGE_VAL;
	}
	screen("Initialize the tabulated region state%s -\n",message);
	screen("\nEnter the table dimensions n%s",dnm[0]);
	if (dim > 1)
	{
	    for (i = 1; i < (dim-1); ++i)
	        screen(", n%s",dnm[i]);
	    screen(", and n%s",dnm[i]);
	}
	screen("\n\t(default =");
	for (i = 0; i < dim; ++i)
	    screen(" %d",table->n[i]);
	screen("): ");
	(void) Gets(s);
	if (s[0] != 0)
	{
	    if (sscanf(s,"%d %d %d\n",table->n,table->n+1,table->n+2) != dim)
	    {
	        screen("ERROR in set_up_read_state_from_file_region(), "
	               "table dimension entered incorrectly\n");
	        clean_up(ERROR);
	    }
	}
	for (len = 1, i = 0; i < dim; ++i)
	    len *= table->n[i];
	table->len = len;
	uni_array(&table->rho,len,FLOAT);
	uni_array(&table->p,len,FLOAT);
	for (i = 0; i < dim; ++i)
	    uni_array(table->v+i,len,FLOAT);
	screen("Enter the tabulated region bounding box");
	for (i = 0; i < dim; ++i)
	    screen(" L%s",Dnm[i]);
	for (i = 0; i < dim; ++i)
	    screen(" U%s",Dnm[i]);
	screen("\n\t(default =");
	for (i = 0; i < dim; ++i)
	    screen(" %g",table->GL[i]);
	for (i = 0; i < dim; ++i)
	    screen(" %g",table->GU[i]);
	screen("): ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    double bb[6];
	    const char *fmt = "%lf %lf %lf %lf %lf %lf";

	    if (sscanf(s,fmt,bb,bb+1,bb+2,bb+3,bb+4,bb+5) != 2*dim)
	    {
	        screen("ERROR in set_up_read_state_from_file_region(), "
	               "bounding box entered incorrectly\n");
	        clean_up(ERROR);
	    }
	    for (i = 0; i < dim; ++i)
	    {
	        table->GL[i] = bb[i];
		table->GU[i] = bb[dim+i];
	    }
	    (void) screen("Bounding box = ");
	    for (i = 0; i < dim; ++i)
	        screen(" %g",table->GL[i]);
	    for (i = 0; i < dim; ++i)
	        screen(" %g",table->GU[i]);
	    screen("\n");
	}
	for (i = 0; i < dim; ++i)
	{
	    table->h[i] = (table->GU[i] - table->GL[i])/table->n[i];
	    table->GL[i] += 0.5*table->h[i];
	    table->GU[i] -= 0.5*table->h[i];
	}

	screen("Enter the file name with the density data: ");
	(void) Gets(s);
	read_table_data(s,table->rho,len);
	screen("Enter the file name with the pressure data: ");
	(void) Gets(s);
	read_table_data(s,table->p,len);
	for (i = 0; i < dim; ++i)
	{
	   screen("Enter the file name with the %s component of "
	          "velocity data: ",dnm[i]);
	   (void) Gets(s);
	   read_table_data(s,table->v[i],len);
	}
	for (i = 0; i < len; ++i)
	{
	    if (table->rho[i] < vacuum_dens)
	        table->rho[i] = vacuum_dens;
	    if (table->p[i] < min_pressure)
	        table->p[i] = min_pressure;
	}
	table->tab_region_set = YES;
}		/*end set_up_read_state_from_file_region*/


/*ARGSUSED*/
LOCAL	void	get_state_tabulated_region(
	double		*coords,
	Locstate	s,
	COMP_TYPE	*ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		stype)
{
	TABULATED_REGION_DATA *table;
	double                 a[3], f[8];
	double                 *GL, *GU, *h; 
	double                 *rho, *p, **v;
	int                   indx[8];
	int                   *n;
	int                   i, j, dim = ct->params->dim;
	int                   li[3], ui[3];

	debug_print("init_states","Entered get_state_tabulated_region()\n");
	table = Tabulated_region_data(ct);
	GL = table->GL;
	GU = table->GU;
	h = table->h;
	n = table->n;
	a[0] = a[1] = a[2] = 0.0;
	li[0] = li[1] = li[2] = 0;
	ui[0] = ui[1] = ui[2] = 0;
	for (i = 0; i < dim; ++i)
	{
	    if (coords[i] < GL[i])
	    {
		li[i] = 0;
	        a[i] = 0.0;
	    }
	    else if (coords[i] > GU[i])
	    {
		li[i] = n[i] - 1;
	        a[i] = 1.0;
	    }
	    else
	    {
	        li[i] = (int) ((coords[i] - GL[i])/h[i]);
	        a[i] = (coords[i] - (GL[i]+li[i]*h[i]))/h[i];
	    }
	    ui[i] = li[i]+1;
	}
	f[0] = (1.0-a[0])*(1.0-a[1])*(1.0-a[2]);
	f[1] =      a[0] *(1.0-a[1])*(1.0-a[2]);
	f[2] = (1.0-a[0])*     a[1] *(1.0-a[2]);
	f[3] =      a[0] *     a[1] *(1.0-a[2]);
	f[4] = (1.0-a[0])*(1.0-a[1])*     a[2] ;
	f[5] =      a[0] *(1.0-a[1])*     a[2] ;
	f[6] = (1.0-a[0])*     a[1] *     a[2] ;
	f[7] =      a[0] *     a[1] *     a[2] ;
	indx[0] = li[0] + n[0]*li[1] + n[0]*n[1]*li[2];
	indx[1] = ui[0] + n[0]*li[1] + n[0]*n[1]*li[2];
	indx[2] = li[0] + n[0]*ui[1] + n[0]*n[1]*li[2];
	indx[3] = ui[0] + n[0]*ui[1] + n[0]*n[1]*li[2];
	indx[4] = li[0] + n[0]*li[1] + n[0]*n[1]*ui[2];
	indx[5] = ui[0] + n[0]*li[1] + n[0]*n[1]*ui[2];
	indx[6] = li[0] + n[0]*ui[1] + n[0]*n[1]*ui[2];
	indx[7] = ui[0] + n[0]*ui[1] + n[0]*n[1]*ui[2];

	set_type_of_state(s,TGAS_STATE);
	rho = table->rho;
	p = table->p;
	for (Dens(s) = 0.0, Press(s) = 0.0, i = 0; i < 8; ++i)
	{
	    Dens(s) += f[i]*rho[indx[i]];
	    Press(s) += f[i]*p[indx[i]];
	}
	v = table->v;
	for (j = 0; j < dim; ++j)
	    for (Vel(s)[j] = 0.0, i = 0; i < 8; ++i)
	        Vel(s)[j] += f[i]*v[j][indx[i]];
        reset_gamma(s);

	Init_params(s,ct->params);
	set_state(s,stype,s);
}		/*end get_state_tabulated_region*/

LOCAL	void	free_tabulated_region_comp_type(
	COMP_TYPE	*comp_type)
{
	COMPONENT             comp;
	TABULATED_REGION_DATA *table;
	int                   i;

	if (comp_type->type != TABULATED_REGION)
	    return;

	table = Tabulated_region_data(comp_type);
	if ((table == NULL) || (table->allocated_from != comp_type))
	    return;

	if (table->rho != NULL)
	    free(table->rho);
	if (table->p != NULL)
	    free(table->p);
	for (i = 0; i < 3; ++i)
	    if (table->v[i] != NULL)
	        free(table->v[i]);

	free(table);
	comp = comp_type->comp;
	zero_scalar(comp_type,sizeof(COMP_TYPE));
	comp_type->type = UNSET_COMP_TYPE;
	comp_type->comp = comp;
}		/*end free_trans_layer_comp_type*/

LOCAL	void	read_table_data(
	const char *tablename,
	double      *data,
	int        len)
{
	IO_TYPE IO_type;
	FILE    *file;
	int     c;
	int     i;

	if ((file = fopen(tablename,"r")) == NULL)
	{
	    screen("ERROR in read_table_data(), can't open %s for reading\n",
	           tablename);
	    clean_up(ERROR);
	}
	determine_io_type(file,&IO_type);
	if ((c = getc(file)) != '\f') /* NOBINARY */
	{
	    (void) ungetc(c,file);
	    for (i = 0; i < len; ++i)
	    {
	        if (fscan_float(file,data+i) != 1)
	        {
		    screen("ERROR in read_table_data(), invalid table %s\n",
		           tablename);
		    clean_up(ERROR);
	        }
	    }
	}
	else
	{
	    if (read_binary_real_array(data,len,&IO_type) != len)
	    {
		screen("ERROR in read_table_data(), invalid table %s\n",
		       tablename);
		clean_up(ERROR);
	    }
	}
	(void) fclose(file);
}		/*end read_table_data*/
