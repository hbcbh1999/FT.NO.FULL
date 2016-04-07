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
*                                gintext.c:
*
*       Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*       Contains diagnostic routines for interface extrema
*
*       TODO: Generalize to more than 2 distinct materials
*/


#include <gprt/glayer.h>

	/* LOCAL Function Prototypes */

LOCAL	boolean	accumulate_fractions_in_layer(Front*,double,
					      double*,double,const double*);
LOCAL	boolean	height_at_fraction(Front*,double*,double,int,
				   double,double*,int,double,const double*);
LOCAL	boolean	is_on_local_grid(const RECT_GRID*,const double*);
LOCAL	void 	add_state_to_totals(const Locstate,Big_State*,
				    Turbulence_moments*);
LOCAL 	void 	convert_to_spherical(double*,const double*,const double*,int);
LOCAL	void 	copy_Big_State(const Big_State*,Big_State*);
LOCAL	void 	copy_state_to_Big_State(const Locstate,Big_State*);
LOCAL   void 	print_data_at_height(const Grid*,const Front*,OUTPUT_DATA*,
				     const Big_State*,const Big_State*,
				     FILE*,FILE*,FILE*,double,double,boolean);
LOCAL   void    print_data_at_height_headers(FILE*,FILE*,FILE*,int);
LOCAL   void    print_intfc_extrema(const Grid*,const Front*,
	                            OUTPUT_DATA*,const Big_State[2],
				    const Big_State[2],boolean,boolean);
LOCAL   void    print_intfc_extrema_headers(const Intfc_extrema*,int);
LOCAL	void 	record_intfc_extrema(Grid*,Wave*,Front*,Printplot*,
				     OUTPUT_DATA*,boolean);
LOCAL	void    spherical_unit_vectors(double*,double*,double*,
				       const double*,const double*,int);
LOCAL	void	zero_state_totals(Big_State*,Turbulence_moments*,int);


LOCAL   void    init_vol_frac (int  rw, int *do_vol_frac, 
			      int *print_freq, char *fn);
LOCAL   void    init_print_allstates(int  rw, int *print_allstates, 
				    int *print_freq, char *fn);
LOCAL   void    init_print_effective_atwood_num(int  rw, int step, 
						Front *front, int *do_spikes, 
						int *print_freq, char *b_fn);
LOCAL   void    print_effective_atwood_num(int *do_spikes, char *b_fn,
					   Grid*, Front*, Wave*, double percent,
					   double bub_height,double spike_height, 
					   double rfactor, const double *origin, 
					   double mix_zone_frac);
LOCAL	void    old_print_vol_frac(Wave *wave, Front *front,
				   double h_min, double h_max,
				   int dir, int index, double rfactor,
				   const double *origin, char *b_fn);


LOCAL   double   find_intersection_by_bisection(Front *front, double *coords_a,
					       double *coords_b, int dir);
LOCAL   double   find_particular_fluid_cell_volume(Front *front, 
						  COMPONENT vol_comp,
						  double *crds, double *dh);
LOCAL	void 	print_vol_frac(Wave *wave, Front *front, 
			       double h_min, double h_max,
			       COMPONENT comp_for_vol_frac, double rfactor,
			       const double *origin, char *b_fn);
LOCAL	void 	new_accumulate_fractions_in_layer(Front *front, double height,
						  double *frac, double rfactor,
						  const double *origin,
						  COMPONENT comp_for_vol_frac);
LOCAL	double 	new_height_at_fraction(Front *front, double h0, double dh,
				       int dir, double fraction, double *frac,
				       COMPONENT vol_comp, double rfactor,
				       const double *origin);

LOCAL  int cmp_for_qsort_dec( const void *vp, const void *vq);
LOCAL  int cmp_for_qsort_inc( const void *vp, const void *vq);
LOCAL  double compute_intfc_area(INTERFACE*);
LOCAL  double compute_tri_area(TRI*);


EXPORT	void init_intfc_extrema(
	INIT_DATA	*init,
	Front		*front,
	Grid		*grid,
	Printplot	*prt)
{
	Gas_param       **prms_list;
	Intfc_extrema	*iext;
	const int	dim = front->rect_grid->dim;
	char 		s[Gets_BUF_SIZE];
	char		*fname, *dname;
	int 		result;
	int             i, nprms;

	if (dim == 1)
	    return;

	screen("Type 'y' to request interface extrema data: ");
	(void) Gets(s);
	if ((s[0] != 'Y') && (s[0] != 'y'))
	    return;

	scalar(&iext,sizeof(Intfc_extrema));

        screen("Type 'y' to request surface area calculation: ");
        (void) Gets(s);
        if ((s[0] != 'Y') && (s[0] != 'y'))
	{
            prt->surface_area = NO;
	}
        else
	{
            prt->surface_area = YES;
        }

	/* defaults */
	Output_mode(&iext->odata) = MESH_TIME;
	Output_step_freq(&iext->odata) = 1;
	Output_start_step(&iext->odata) = 0;
	Output_in_binary(&iext->odata) = NO;

	init_output_data(init,&iext->odata,grid,prt,NULL,NO,NO,YES);
	add_user_output_function(record_intfc_extrema,&iext->odata,prt);

	/* Added 24 Oct. 2003 - egeorge */
	init_print_effective_atwood_num(0, restart_time_step(init), 
					front, NULL, NULL, NULL);
	init_vol_frac(0, NULL, NULL, NULL);
	init_print_allstates(0, NULL, NULL, NULL);

	iext->geom = planar;	/* default */
	screen("Compute interface extrema for planar ('p'%s) or"
	       " radial geometry ('r'%s): ",
	       (iext->geom == planar)    ? ", default" : "",
	       (iext->geom == spherical) ? ", default" : "");
	(void) Gets(s);
	if (s[0] == 'R' || s[0] == 'r')
	    iext->geom = spherical;

	if (iext->geom == spherical)
	{
	    double *origin = iext->origin;
	    origin[0] = origin[1] = origin[0] = 0.0;
	    sprint_general_vector(s,"(default = ",origin,dim,")");
	    screen("Enter the coordinates of the origin %s: ",s);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
		const char *fmt = "%lf %lf %lf";
	        result = sscanf(s,fmt,origin,origin+1,origin+2);
	        if (result < dim)
	        {
		    screen("ERROR in init_intfc_extrema(), "
		           "Insufficient number of coordinates.\n");
		    clean_up(ERROR);
	        }
	    }
	}

	iext->rfactor = 2.0;		/* default */
	screen("Enter a sub-grid refinement factor for the averaging of the\n"
	       " ambient state at the interface extrema"
	       " (default = %g): ",iext->rfactor);
	(void) Gets(s);
	if (s[0] !='\0')
	    (void) sscan_float(s,&iext->rfactor);

	screen("Current gas param list\n");
	nprms = return_params_list(&prms_list);
	screen("Number of params = %d\n\n",nprms);
	for (i = 0; i < nprms; i++)
	{
	    IMPORT boolean suppress_prompts;
	    screen("Param[%d]\n",i);
	    fprint_Gas_param(stdout,prms_list[i]);
	    if (suppress_prompts == NO)
	        fprint_Gas_param(stderr,prms_list[i]);
	    screen("\n");
	}
	screen("Enter the EOS indices of the %s and %s materials, "
	       "respectively: ",
	       (iext->geom == spherical ? "inner" : "lower"),
	       (iext->geom == spherical ? "outer" : "upper"));
	Scanf("%d %d\n",&iext->index_at_max,&iext->index_at_min);

	set_default_data_file_names("intfc_extrema",".min",&dname,&fname,init);
	iext->out_min = open_data_file(front,"interface minimum",YES,YES,dname,
				       &iext->min_dname,fname,&iext->min_fname);

	set_default_data_file_names("intfc_extrema",".max",&dname,&fname,init);
	iext->out_max = open_data_file(front,"interface maximum",YES,YES,dname,
				       &iext->max_dname,fname,&iext->max_fname);

	set_default_data_file_names("intfc_extrema",".amp",&dname,&fname,init);
	iext->out_amp = open_data_file(front,"interface amplitude",YES,YES,
	                               dname,&iext->amp_dname,
	                               fname,&iext->amp_fname);

	screen("Type 'y' to get data for 1%%-99%% levels: ");
	(void) Gets(s);
	if (s[0] == 'Y' || s[0] == 'y')
	{
	    char fname[Gets_BUF_SIZE];
	    char *tmpdname, *tmpfname;

	    iext->do_01 = YES;

	    (void) sprintf(fname,"%s%s",
			   (iext->min_fname != NULL)?iext->min_fname:"","01");
	    iext->out_min1 = open_data_file(front,"interface 1% minimum",
	                                    YES,YES,iext->min_dname,
	                                    &tmpdname,fname,&tmpfname);

	    (void) sprintf(fname,"%s%s",
			   (iext->max_fname != NULL)?iext->max_fname:"","01");
	    iext->out_max1 = open_data_file(front,"interface 1% maximum",
	                                    YES,YES,iext->max_dname,
	                                    &tmpdname,fname,&tmpfname);

	    (void) sprintf(fname,"%s%s",
			   (iext->amp_fname != NULL)?iext->amp_fname:"","01");
	    iext->out_amp1 = open_data_file(front,"interface 1% amplitude",
	                                    YES,YES,iext->amp_dname,
	                                    &tmpdname,fname,&tmpfname);
	}
	screen("Type 'y' to get data for 5%%-95%% levels: ");
	(void) Gets(s);
	if (s[0] == 'Y' || s[0] == 'y')
	{
	    char fname[Gets_BUF_SIZE];
	    char *tmpdname, *tmpfname;

	    iext->do_05 = YES;

	    strcpy(fname,iext->min_fname);
	    strcat(fname,"05");
	    iext->out_min5 = open_data_file(front,"interface 5% minimum",
	                                    YES,YES,iext->min_dname,
	                                    &tmpdname,fname,&tmpfname);
	    strcpy(fname,iext->max_fname);
	    strcat(fname,"05");
	    iext->out_max5 = open_data_file(front,"interface 5% maximum",
	                                    YES,YES,iext->max_dname,
	                                    &tmpdname,fname,&tmpfname);
	    strcpy(fname,iext->amp_fname);
	    strcat(fname,"05");
	    iext->out_amp5 = open_data_file(front,"interface 5% amplitude",
	                                    YES,YES,iext->amp_dname,
	                                    &tmpdname,fname,&tmpfname);
	}
}		/*end init_intfc_extrema*/


LOCAL void convert_to_spherical(
	double        *p_sph,
	const double  *p_rect,
	const double  *origin,
	int          dim)
{
	int   i;
	double pr[3];

	for (i = 0; i < dim; i++)
	    pr[i] = p_rect[i];
	if (origin != NULL)
	    for (i = 0; i < dim; i++)
	        pr[i] -= origin[i];
	for (; i < 3; i++)
	    pr[i] = 0.0;

	p_sph[0] = mag_vector(pr,dim);
	if (dim == 3)
	    p_sph[1] = acos(pr[2]/p_sph[0]);
	p_sph[dim-1] = atan2(pr[1],pr[0]);
}		/*end convert_to_spherical*/

/*
*			spherical_unit_vectors():
*
*	Compute rectangular components of standard orthonormal uni_arrays ("r hat",
*	"theta hat", "phi hat") for polar (dim == 2) or spherical (dim == 3)
*	coordinates given a point p_rect in rectangular coordinates and an
*	origin.
*/

LOCAL	void	spherical_unit_vectors(
	double       *rhat,
	double       *thhat,
	double       *phhat,
	const double *p_rect,
	const double *origin,
	int         dim)
{
	double p_sph[3], snph, csph;

	convert_to_spherical(p_sph,p_rect,origin,dim);
	switch(dim)
	{
	case 1:
	    screen("ERROR in spherical_unit_vectors(), 1D not supported\n");
	    clean_up(ERROR);
	    break;
#if defined(TWOD)
	case 2:  /* polar coordinates */
	    snph = sin(p_sph[1]);
	    csph = cos(p_sph[1]);
	     rhat[0] =  csph;  rhat[1] = snph;
	    phhat[0] = -snph; phhat[1] = csph;
	    break;
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:
	    {
	        double snth, csth;
	        snph = sin(p_sph[2]);
	        csph = cos(p_sph[2]);
	        snth = sin(p_sph[1]);
	        csth = cos(p_sph[1]);
	         rhat[0] = snth*csph;  rhat[1] = snth*snph;  rhat[2] =  csth;
	        thhat[0] = csth*csph; thhat[1] = csth*snph; thhat[2] = -snth;
	        phhat[0] = -snph;     phhat[1] = csph;      phhat[2] =  0.0;
	    }
	    break;
#endif /* defined(THREED) */
	}
}		/*end spherical_unit_vectors*/

LOCAL	boolean	is_on_local_grid(
	const RECT_GRID *rgrid,
	const double     *coords)
{
	int   i;
	const int   dim = rgrid->dim;
	const double *L = rgrid->L, *U = rgrid->U;
	for (i = 0; i < dim; i++)
	    if (coords[i] < L[i] || coords[i] > U[i])
		return NO;
	return YES;
}		/*end is_on_local_grid*/

LOCAL	void 	copy_state_to_Big_State(
	const Locstate	st,
	Big_State	*bst)
{
	int	i;

	bst->d = Dens(st);
	for (i = 0; i < MAXD; i++)
	    bst->v[i] = vel(i,st);
	bst->p = pressure(st);
	bst->k = kinetic_energy(st);
	bst->e = internal_energy(st);
}		/*end copy_state_to_Big_State*/

LOCAL	void 	copy_Big_State(
	const Big_State	*bst0,
	Big_State	*bst)
{
	int	i;

	bst->d = bst0->d;
	for (i = 0; i < MAXD; i++)
	    bst->v[i] = bst0->v[i];
	bst->p = bst0->p;
	bst->k = bst0->k;
	bst->e = bst0->e;
	bst->frac = bst0->frac;
	bst->count = bst0->count;
}		/*end copy_Big_State*/

LOCAL	boolean accumulate_fractions_in_layer(
	Front	    *front,
	double	    height,
	double	    *frac,
	double	    rfactor,
	const double *origin)
{
	INTERFACE		  *intfc = front->interf;
	const RECT_GRID		  *rgrid = front->rect_grid;
	const enum intfc_geometry geom = (origin == NULL ? planar : spherical);
	const int		  dim = rgrid->dim;
	const int		  n_params = num_gas_params(intfc);
	const int		  nn = pp_numnodes();
	const int		  zdir = dim - 1;
	double			  dx, dphi;
	double			  nor_fac, coords[3];
	double			  x_left, phi;
	int			  Nx, Nphi;
	int			  p;
	register int 		  i;
	COMPONENT		  comp;
	const double		  dh = (dim == 2) ?
				      min(rgrid->h[0],rgrid->h[1])/rfactor
				      : min(min(rgrid->h[0],rgrid->h[1]),
					    rgrid->h[2])/rfactor;

	static long int		*count = NULL;

	debug_print("glayer","Entered accumulate_fractions_in_layer(), h = %g\n",
		       height);

	if (count == NULL)
	    uni_array(&count,n_params,sizeof(long int));

	zero_scalar(count,n_params*sizeof(long int));

	for (i = 0; i < 3; i++)
	    coords[i] = 0.0;

	switch(geom)
	{
	case planar:
	    coords[zdir] = height;
	    if (height < rgrid->L[zdir] || height >= rgrid->U[zdir])
	      break;
	    Nx = irint(rgrid->gmax[0]*rfactor);
	    dx = (rgrid->U[0] - rgrid->L[0])/Nx;
	    x_left = rgrid->L[0] + 0.5*dx;
	    switch(dim)
	    {
	    case 1:
	        screen("ERROR in accumulate_fractions_in_layer(), "
		       "1D not supported\n");
	        clean_up(ERROR);
	        break;
#if defined(TWOD)
	    case 2:
	        for (i = 0; i < Nx; i++)
		{
		    coords[0] = x_left + i*dx;
		    comp = component(coords,intfc);
		    count[index_of_Gas_param(
				      gas_params_for_comp(comp,intfc))]++;
		}
		break;
#endif /* defined(TWOD) */
#if defined(THREED)
	    case 3:
		{
		    double dy, y_left;
		    int j, Ny;
		    Ny = irint(rgrid->gmax[1]*rfactor);
		    dy = (rgrid->U[1] - rgrid->L[1])/Ny;
		    y_left = rgrid->L[1] + 0.5*dy;
		    for (j = 0; j < Ny; j++)
		    {
		        coords[1] = y_left + j*dy;
		        for (i = 0; i < Nx; i++)
		        {
			    coords[0] = x_left + i*dx;
			    comp = component(coords,intfc);
			    count[index_of_Gas_param(
					  gas_params_for_comp(comp,intfc))]++;
		        }
		    }
		}
		break;
#endif /* defined(THREED) */
	    }
	    break;

	case spherical:
	    switch(dim)
	    {
	    case 1:
	        screen("ERROR in accumulate_fractions_in_layer(), "
		       "1D not supported\n");
	        clean_up(ERROR);
	        break;
#if defined(TWOD)
	    case 2:
	        Nphi = irint(2.0*PI*height/dh);
		dphi = 2.0*PI/Nphi;
		for (i = 0; i < Nphi; i++)
		{
		    phi = i*dphi;
		    coords[0] = height*cos(phi) + origin[0];
		    coords[1] = height*sin(phi) + origin[1];
		    if (is_on_local_grid(rgrid,coords) != YES)
		      continue;
		    comp = component(coords,intfc);
		    count[index_of_Gas_param(
				      gas_params_for_comp(comp,intfc))]++;
		}
	        break;
#endif /* defined(TWOD) */
#if defined(THREED)
	    case 3:
		{
	            double dth, th, rpolar;
		    int   j, Nth;
		    Nth = irint(PI*height/dh);
		    dth = PI/Nth;
		    for (j = 0; j < Nth; j++)
		    {
		        th = j*dth;
		        rpolar = height*sin(th);
		        Nphi = irint(2.0*PI*rpolar/dh);
		        dphi = 2.0*PI/Nphi;
		        for (i = 0; i < Nphi; i++)
		        {
			    phi = i*dphi;
			    coords[0] = rpolar*cos(phi) + origin[0];
			    coords[1] = rpolar*sin(phi) + origin[1];
			    coords[2] = height*cos(th) + origin[2];
			    if (is_on_local_grid(rgrid,coords) != YES)
			        continue;
			    comp = component(coords,intfc);
			    count[index_of_Gas_param(
				      gas_params_for_comp(comp,intfc))]++;
		        }
		    }
		}
		break;
#endif /* defined(THREED) */
	    }
	    break;
	}

	if (debugging("glayer"))
	{
	    (void) printf("\tLocal accumulation of fractions at "
			  "h = %g completed.\n",height);
	    for (p = 0; p < n_params; p++)
	        (void) printf("\t\tcount[%d] = %ld\t",p,count[p]);
	    (void) printf("\n");
	}

	if (nn > 1)
	    pp_global_lsum(count,(long int)n_params);

	if (debugging("glayer"))
	{
	    (void) printf("\tGlobal accumulation of fractions at "
			  "h = %g completed.\n",height);
	    for (p = 0; p < n_params; p++)
	        (void) printf("\t\tcount[%d] = %ld\t",p,count[p]);
	    (void) printf("\n");
	}

	for (p = 0; p < n_params && count[p] == 0; p++);
	if (p == n_params)
	{
	    (void) printf("WARNING in accumulate_fractions_in_layer(), "
			  "No points found in layer at h = %g\n",height);
	    return FUNCTION_FAILED;
	}

	for (p = 0, nor_fac = 0.0; p < n_params; p++)
	    nor_fac += (double) count[p];
	for (p = 0; p < n_params; p++)
	    frac[p] = ((double) count[p])/nor_fac;

	debug_print("glayer","Left accumulate_fractions_in_layer()\n");
	
	return FUNCTION_SUCCEEDED;
}		/*end accumulate_fractions_in_layer*/

/*
*			height_at_fraction():
*
*	Find the height (z coordinate or radius) closest to initial height
*	h0 in the given direction dir (+/- 1) where the given layer fraction
*	occurs.  This fraction corresponds to the material whose params
*	has the given index.
*/

LOCAL 	boolean	height_at_fraction(
	Front	    *front,
	double	    *h0,
	double	    dh,
	int	    dir,
	double	    fraction,
	double	    *frac,
	int	    index,
	double	    rfactor,
	const double *origin)
{
	boolean             result;
	const int        n_params = num_gas_params(front->interf);
	double            height = *h0, hlo, hhi;
	int              p, iter2;
	static const int max_num_iters = 3;
	static double     *new_frac = NULL, *frac_hi = NULL, *frac_lo = NULL;

	debug_print("gfrac","Entered height_at_fraction()\n");

	if (new_frac == NULL)
	{
	    uni_array(&new_frac,n_params,FLOAT);
	    uni_array(&frac_hi,n_params,FLOAT);
	    uni_array(&frac_lo,n_params,FLOAT);
	}

	(void) accumulate_fractions_in_layer(front,height,frac,rfactor,origin);
	if (frac[index] >= fraction)
	{
	    height += dir*dh*fraction/frac[index];
	    *h0 = height;
	    (void) printf("WARNING in height_at_fraction(), "
			  "frac >= %g occurred at initial height of %g.\n",
			  fraction,height);
	    debug_print("gfrac","Left height_at_fraction()\n");
	    return FUNCTION_SUCCEEDED;
	}

	if (debugging("gfrac"))
	{
	    (void) printf("\theight = %g, fraction = %g, index = %d\n",
			  height,fraction,index);
	    (void) printf("\tAfter initialization, "
			  "frac[%d] = %g, frac[%d] = %g\n",
			  index,frac[index],1-index,frac[1-index]);
	}

	do
	{
	    result = accumulate_fractions_in_layer(front,height+dir*dh,
						   new_frac,rfactor,origin);
	    if (result == FUNCTION_SUCCEEDED)
	    {
	        if ((  (frac[index] <= fraction) &&
		       (fraction <= new_frac[index])) ||
		     ( (new_frac[index] <= fraction) &&
		       (fraction <= frac[index])))
	        {   /* we've bracketed the height, so interpolate */
	            height += dir*dh*(fraction-frac[index]) /
			             (new_frac[index]-frac[index]);
		    break;
	        }
	        else
	        {
	            for (p = 0; p < n_params; p++)
		        frac[p] = new_frac[p];
		    height += dir*dh;
	        }
	    }
	} while (result == FUNCTION_SUCCEEDED);

	if (result == FUNCTION_FAILED)
	{
	    (void) printf("WARNING in height_at_fraction(), "
	                  "unable to locate height at which frac[%d] = %g\n",
			  index,fraction);
	    (void) printf("\toriginal height = %g, direction = %d\n",*h0,dir);
	    (void) printf("\tcurrent height = %g\n",height);
	    debug_print("gfrac","Left height_at_fraction()\n");
	    return FUNCTION_FAILED;
	}

	if ((frac[index] <= fraction) && (fraction <= new_frac[index]))
	{
	    hlo = height - dir*dh;
	    hhi = height;
	    for (p = 0; p < n_params; p++)
	    {
	        frac_lo[p] = frac[p];
		frac_hi[p] = new_frac[p];
	    }
	}
	else
	{
	    hlo = height;
	    hhi = height - dir*dh;
	    for (p = 0; p < n_params; p++)
	    {
	        frac_lo[p] = new_frac[p];
		frac_hi[p] = frac[p];
	    }
	}

	if (debugging("gfrac"))
	{
	    (void) printf("\theight at frac = %g has been framed.\n",fraction);
	    (void) printf("\thlo = %g; frac_lo[%d] = %g\n",
			  hlo,index,frac_lo[index]);
	    (void) printf("\thhi = %g; frac_hi[%d] = %g\n",
			  hhi,index,frac_hi[index]);
	}

	for (iter2 = 0; iter2 < max_num_iters; iter2++)
	{
	    height = 0.5*(hlo + hhi);
	    result = accumulate_fractions_in_layer(front,height,
						   new_frac,rfactor,origin);
	    if (result == FUNCTION_FAILED)
	    {
	        screen("ERROR in height_at_fraction(), "
		       "\taccumulate_fractions_in_layer() failed during"
		       "bisection loop.\n");
		clean_up(ERROR);
	    }
	    if (fraction <= new_frac[index])
	    {
	        hhi = height;
		for (p = 0; p < n_params; p++)
		    frac_hi[p] = new_frac[p];
	    }
	    else
	    {
	        hlo = height;
		for (p = 0; p < n_params; p++)
		    frac_lo[p] = new_frac[p];
	    }
	}

	height = ((fraction-frac_lo[index])*hhi+(frac_hi[index]-fraction)*hlo)/
		 (frac_hi[index]-frac_lo[index]);

	(void) accumulate_fractions_in_layer(front,height,frac,rfactor,origin);

	if (frac[index] == 0.0 || frac[index] == 1.0)
	{
	    (void) printf("WARNING in height_at_fraction(), "
	                  "Interpolated height of %g has fraction = %g\n",
			  height,frac[index]);
	    if (frac[index] == 0.0)
	    {
	        height = hhi;
		for (p = 0; p < n_params; p++)
		    frac[p] = frac_hi[p];
		(void) printf("\tWill use h = %g, frac = %g instead.\n",
			      hhi, frac[index]);
	    }
	    else
	    {
	        height = hlo;
		for (p = 0; p < n_params; p++)
		    frac[p] = frac_lo[p];
		(void) printf("\tWill use h = %g, frac = %g instead.\n",
			      hlo, frac[index]);
	    }
	}

	if (debugging("gfrac"))
	{
	    (void) printf("\tAfter %d bisection iterations,\n",max_num_iters);
	    (void) printf("\theight = %0.8e; frac[%d] = %0.8e\n",
		          height,index,frac[index]);
	    (void) printf("\thlo = %0.8e; frac_lo[%d] = %0.8e\n",
		          hlo,index,frac_lo[index]);
	    (void) printf("\thhi = %0.8e; frac_hi[%d] = %0.8e\n",
		          hhi,index,frac_hi[index]);
	}

	*h0 = height;
	debug_print("gfrac","Left height_at_fraction()\n");
	return FUNCTION_SUCCEEDED;
}		/*end height_at_fraction*/

/*			record_intfc_extrema()
*
*	Find the extremal points of all contact surfaces, the states at these
*	points, and the ambient states (averaged across the layer at the same
*	height or radius), and print the results.
*
*	NOTE:  The performance of this function in situations where one or more
*	extremal points occur on pieces of the contact surface that have broken
*	off of the main branch is untested.
*/

/*ARGSUSED*/
LOCAL	void record_intfc_extrema(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*out,
	boolean		about_to_stop)
{
	Big_State	   bst[2], amb_st[2];
	HYPER_SURF	   *hs, *hs_max, *hs_min;
	HYPER_SURF_ELEMENT *hse, *hse_max, *hse_min;
	INTERFACE	   *intfc = front->interf;
	Intfc_extrema 	   *iext = Intfc_extrema_data(out);
	Locstate	   sminl, sminr, smaxl, smaxr;
	POINT		   *p, *p_max, *p_min;
	const double	   *origin = (iext->geom == spherical) ?
				      iext->origin : NULL;
	const int	   dim = front->rect_grid->dim;
	const int 	   myid = pp_mynode();
	const int	   nn = pp_numnodes();
	const int	   zdir = dim - 1;
	const double        *L = front->rect_grid->L;
	const double        *U = front->rect_grid->U;
	double		   g_h_at_min, g_h_at_max;
	double		   h, h_at_min, h_at_max;
	double		   nor_at_max[MAXD], nor_at_min[MAXD];
	int                index_at_max = iext->index_at_max;
	int                index_at_min = iext->index_at_min;
	int		   i;
	static const size_t size_bst = sizeof(Big_State);
	static double 	    **min_by_node = NULL;
	static double 	    **max_by_node = NULL;

	/* The following type declarations added 24 Oct. 2003 - egeorge */
	static boolean	    first = YES;
	static long         general_counter = 0; 
	static int          eff_at_print_freq;
	static int          do_spikes = NO;
	static char         b_fn[256]; 
	static int          do_vol_frac = NO;
	static int          vol_frac_print_freq;
	static char         vf_fn[256];

	static int          print_allstates = NO;
	static int          allstates_print_freq;
	static char         allstates_fn[256];

	double               percent;
	COMPONENT	    comp_at_max, comp_at_min;
        boolean  surface_area;


	debug_print("gintext","Entered record_intfc_extrema()\n");
	start_clock("record_intfc_extrema");
	
        /*printf("SURFACE AREA ON THE INTERFACE : %f\n", compute_intfc_area(intfc)); */


	if (first == YES)
	{
	    first = NO;

	    init_print_effective_atwood_num(1, 0, NULL, &do_spikes, 
					    &eff_at_print_freq, b_fn);
	    init_vol_frac(1, &do_vol_frac, &vol_frac_print_freq, vf_fn);
	    init_print_allstates(1, &print_allstates, &allstates_print_freq, 
				 allstates_fn);
	}



	if (min_by_node == NULL)
	{
	    bi_array(&min_by_node,nn,MAXD,FLOAT);
	    bi_array(&max_by_node,nn,MAXD,FLOAT);
	    if (is_io_node(myid))
	    {
	        print_intfc_extrema_headers(iext,dim);
		if (iext->do_01 == YES)
		    print_data_at_height_headers(iext->out_min1,iext->out_max1,
						 iext->out_amp1,dim);
		if (iext->do_05 == YES)
		    print_data_at_height_headers(iext->out_min5,iext->out_max5,
						 iext->out_amp5,dim);
	    }
	}

	h_at_min =  HUGE_VAL;
	h_at_max = -HUGE_VAL;
	for (i = 0; i < dim; i++)
	{
	    iext->pos[0][i] = HUGE_VAL;
	    iext->pos[1][i] = -HUGE_VAL;
	}

	/*
	*	Loop over all points in every interface, to find the
	*	extremal points of the contact surfaces.
	*
	* 	The method used here is a little wasteful, as it loops
	*	over all points of irrelevant interfaces.
	*/

	start_clock("find_interface_extrema");

	p_max = p_min = NULL;
	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (!is_scalar_wave(wave_type(hs)))
		continue;

	    /* Make sure p is on the actual grid */
	    for (i = 0; i < dim; i++)
		if (Coords(p)[i] < L[i] || Coords(p)[i] >= U[i])
		    break;
	    if (i < dim)
		continue;

	    h = (iext->geom == planar) ? Coords(p)[zdir] :
	                                 distance_between_positions(Coords(p),
								    origin,dim);

	    if (h < h_at_min)
	    {
	        h_at_min = h;
	        p_min = p;
	        hs_min = hs;
	        hse_min = hse;
	    }
	    if (h > h_at_max)
	    {
	        h_at_max = h;
	        p_max = p;
	        hs_max = hs;
	        hse_max = hse;
	    }
	}

	/*  Now find normal vector and states at min & max */

	if ((p_max != NULL) && (p_min != NULL))
	{
	    normal(p_max,hse_max,hs_max,nor_at_max,front);
	    normal(p_min,hse_min,hs_min,nor_at_min,front);
	    slsr(p_max,hse_max,hs_max,&smaxl,&smaxr);
	    slsr(p_min,hse_min,hs_min,&sminl,&sminr);

	    for (i = 0; i < dim; i++)
	    {
		iext->pos[1][i] = Coords(p_max)[i];
		iext->pos[0][i] = Coords(p_min)[i];
	    }

	    if (index_of_Gas_param(gas_params_for_comp(hs_max->neg_comp,intfc))
			    == index_at_max)
	    {
	        copy_state_to_Big_State(smaxl,&bst[1]);
		comp_at_max = hs_max->neg_comp;
	    }
	    else
	    {
		copy_state_to_Big_State(smaxr,&bst[1]);
		comp_at_max = hs_max->pos_comp;
	    }

	    if (index_of_Gas_param(gas_params_for_comp(hs_min->pos_comp,intfc))
		            == index_at_min)
	    {
		copy_state_to_Big_State(sminr,&bst[0]);
		comp_at_min = hs_max->pos_comp;
	    }
	    else
	    {
		copy_state_to_Big_State(sminl,&bst[0]);
		comp_at_min = hs_max->neg_comp;
	    }
	}  /* if ((p_max != NULL) and ... ) */

	/* Compute global max and min */

	if (nn > 1)
	{
	    int i_min, i_max;

	    for (i = 0; i < dim; i++)
	    {
	        min_by_node[myid][i] = iext->pos[0][i];
		max_by_node[myid][i] = iext->pos[1][i];
	    }

	    /* I/O node collects each of the other nodes'
	       minimum and maximum position. */

	    if (is_io_node(myid))
	    {
	        for (i = 0; i < nn; i++)
		{
		    if (i == myid)
			continue;
		    pp_recv(L_MIN_MAX_ID+i,i,min_by_node[i],MAXD*FLOAT);
		    pp_recv(L_MIN_MAX_ID+nn+i,i,max_by_node[i],MAXD*FLOAT);
		}
	    }
	    else
	    {
	        pp_send(L_MIN_MAX_ID+myid,min_by_node[myid],
			MAXD*FLOAT,IO_NODE_ID);
	        pp_send(L_MIN_MAX_ID+nn+myid,max_by_node[myid],
			MAXD*FLOAT,IO_NODE_ID);
	    }

	    if (debugging("gintext"))
	        (void) printf("Local extrema have been communicated.\n");

	    /* I/O node has a list of all the local extrema, arranged 
	       by node index, so now it can find the global extrema. */

	    if (is_io_node(myid))
	    {
	        switch(iext->geom)
		{
		case planar:
		    g_h_at_min = min_by_node[myid][zdir];
		    g_h_at_max = max_by_node[myid][zdir];
		    break;
		case spherical:
		    g_h_at_min = distance_between_positions(min_by_node[myid],
							    origin,dim);
		    g_h_at_max = distance_between_positions(max_by_node[myid],
							    origin,dim);
		    break;
		}

		i_max = i_min = myid;

		for (i = 0; i < nn; i++)
		{
		    switch(iext->geom)
		    {
		    case planar:
			h_at_min = min_by_node[i][zdir];
			h_at_max = max_by_node[i][zdir];
			break;
		    case spherical:
			h_at_min = distance_between_positions(min_by_node[i],
							      origin,dim);
			h_at_max = distance_between_positions(max_by_node[i],
							      origin,dim);
			break;
		    }
		    if (h_at_min < g_h_at_min)
		    {
			g_h_at_min = h_at_min;
			i_min = i;
		    }
		    if (h_at_max > g_h_at_max)
		    {
			g_h_at_max = h_at_max;
			i_max = i;
		    }
		}
		pp_send_all(G_IMIN_ID,&i_min,INT);
		pp_send_all(G_IMAX_ID,&i_max,INT);
		pp_send_all(G_HMIN_ID,&g_h_at_min,FLOAT);
		pp_send_all(G_HMAX_ID,&g_h_at_max,FLOAT);
	    }
	    else
	    {
	        pp_recv(G_IMIN_ID,IO_NODE_ID,&i_min,INT);
		pp_recv(G_IMAX_ID,IO_NODE_ID,&i_max,INT);
		pp_recv(G_HMIN_ID,IO_NODE_ID,&g_h_at_min,FLOAT);
		pp_recv(G_HMAX_ID,IO_NODE_ID,&g_h_at_max,FLOAT);
	    }

	    /* See if this node has a global extremum, and if so,
	       send the associated information to the I/O node. */

	    if (myid == i_min)
	    {
	        if (debugging("gintext"))
		    (void) printf("I have the global minimum.\n");
		pp_send_all(MIN_POS_ID,iext->pos[0],dim*FLOAT);
		pp_send_all(MIN_STATE_ID,&bst[0],size_bst);
	    }
	    else
	    {
	        if (debugging("gintext"))
		    (void) printf("Node %d has the global minimum.\n",i_min);
		pp_recv(MIN_POS_ID,i_min,iext->pos[0],dim*FLOAT);
		pp_recv(MIN_STATE_ID,i_min,&bst[0],size_bst);
	    }
	    if (myid == i_max)
	    {
	        if (debugging("gintext"))
		    (void) printf("I have the global maximum.\n");
		pp_send_all(MAX_POS_ID,iext->pos[1],dim*FLOAT);
		pp_send_all(MAX_STATE_ID,&bst[1],size_bst);
	    }
	    else
	    {
	        if (debugging("gintext"))
		    (void) printf("Node %d has the global maximum.\n",i_max);
		pp_recv(MAX_POS_ID,i_max,iext->pos[1],dim*FLOAT);
		pp_recv(MAX_STATE_ID,i_max,&bst[1],size_bst);
	    }

	}  /* if (nn > 1) */
	else
	{
	    g_h_at_min = h_at_min;
	    g_h_at_max = h_at_max;
	}

	/* Get the ambient states by a global layer average. */
	
	accumulate_state_in_layer(wave,front,g_h_at_min,&amb_st[0],NULL,
				  iext->rfactor,origin,YES,YES,YES,NO);
	accumulate_state_in_layer(wave,front,g_h_at_max,&amb_st[1],NULL,
				  iext->rfactor,origin,YES,YES,YES,NO);

	stop_clock("find_interface_extrema");
	

	/* egeorge */
	if (do_vol_frac == YES && (general_counter % vol_frac_print_freq == 0)) 
	{ 
	    /* I intend to eventually make volume fraction calculation for 
	      spherical geometry revert back to the old method. Note however 
	      that the spherical case is not yet implemented.
	    */
	    if (iext->geom == planar)
	    {
		old_print_vol_frac(wave, front, g_h_at_min, g_h_at_max, 1, 
		   index_of_Gas_param(gas_params_for_comp(comp_at_min,intfc)),
				   iext->rfactor,origin,vf_fn);
	       
		
		print_vol_frac(wave, front, g_h_at_min, g_h_at_max,
			       comp_at_min, iext->rfactor, origin,vf_fn); 
	    }
	    else if (iext->geom == spherical)
	    {
		old_print_vol_frac(wave, front, g_h_at_min, g_h_at_max, 1,
		   index_of_Gas_param(gas_params_for_comp(comp_at_min,intfc)),
				   iext->rfactor,origin,vf_fn);
	    }
	}


	if (print_allstates == YES && (general_counter % allstates_print_freq == 0)) 
	{ 
	    /* TODO - the spherical case (it is not yet fully implemented
	       as of 16 Apr. 2004).
	    */

	    FILE       *ALL_STATES;
	    double      frac[2]; /* assumes only 2 gas parameters */
	    double      allstates_height;
	    double      allstates_dh = (g_h_at_max - g_h_at_min)/100.0;
	    char       filename[256];
	    int        start_index;
	    Big_State  allstates_bst[2];


	    sprintf(filename,"%s-allstates-ts%d-t%07.2f",
		    allstates_fn,front->step,front->time);

	    if (myid == 0)
	    {
		ALL_STATES = fopen(filename,"w");
		
		fprintf(ALL_STATES,"%12s %12s %12s %12s %38s %38s %12s %12s ", 
		  "P_h","P_l","D_h","D_l","V_spike","V_bub","Beta_l","Height");
		fprintf(ALL_STATES,"%12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
		  "Spk_edge", "Bub_edge","Time","D3_h","D3_l","D5_h","D5_l","D6_h","D6_l");
	    }



	    if (iext->geom == planar)
	    {
	        if (g_h_at_min >= front->rect_grid->GL[zdir] &&
		    g_h_at_max <= front->rect_grid->GU[zdir])
		{
		    start_index = ((g_h_at_min - front->rect_grid->GL[zdir] > 
				    allstates_dh) ? 
				   ((int) ((g_h_at_min - front->rect_grid->GL[zdir])
					   / allstates_dh) - 1)
				   : 0);

		    allstates_height = front->rect_grid->GL[zdir] + 
		      start_index*allstates_dh;

		    while (allstates_height < g_h_at_min)
		        allstates_height +=allstates_dh;

		    while (allstates_height < g_h_at_max)
		    {
		        accumulate_state_in_layer(wave,front,allstates_height,
						  allstates_bst,NULL,
						  iext->rfactor,
						  origin,NO,YES,YES,YES);

			/* Until the way in which accumulate_state_in_layer() 
			   averages states is made more compatible with 
			   the way new_accumulate_fractions_in_layer() works, 
			   we will use the old version of that latter 
			   function. 3 May 2004.

			
			   new_accumulate_fractions_in_layer(front,
			   allstates_height,&frac[1],iext->rfactor,
			   origin, comp_at_min);
			
			   frac[0] = 1.0 - frac[1];  
			*/
			
			accumulate_fractions_in_layer(front,allstates_height,
						      frac,iext->rfactor,
						      origin);
			
			if (myid == 0)
			{	
			    fprintf(ALL_STATES,
			      "%12g %12g %12g %12g %12g %12g %12g %12g %12g ",
			      allstates_bst[0].p, allstates_bst[1].p, 
			      allstates_bst[0].d, allstates_bst[1].d, 
			      allstates_bst[0].v[0], allstates_bst[0].v[1],
			      allstates_bst[0].v[2], allstates_bst[1].v[0],
			      allstates_bst[1].v[1]);
			    fprintf(ALL_STATES,
			      "%12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12\n",	
			      allstates_bst[1].v[2], frac[1], allstates_height,
			      g_h_at_max, g_h_at_min, front->time,
			      allstates_bst[0].global_max_d,
			      allstates_bst[1].global_min_d,
			      allstates_bst[0].partial_d,
			      allstates_bst[1].partial_d,
			      allstates_bst[0].median_d,
			      allstates_bst[1].median_d);
			}

			allstates_height += allstates_dh;
		    }

		}
	    }
	    else if (iext->geom == spherical)
	    {
	        (void) printf("WARNING - spherical geometry printing of\n");
		(void) printf("state data at fixed radius and time not\n");
		(void) printf("yet implemented.\n");
	    }

	    if (myid == 0)
	        fclose(ALL_STATES);
	}



	if (eff_at_print_freq > 0 && (general_counter % eff_at_print_freq == 0))
	{
	    percent = 0.0;
	    print_effective_atwood_num(&do_spikes, b_fn, grid, front, wave, 
				       percent, g_h_at_min, g_h_at_max,
				       iext->rfactor, origin, 1.0/3.0);
	    print_effective_atwood_num(&do_spikes, b_fn, grid, front, wave, 
				       percent, g_h_at_min, g_h_at_max, 
				       iext->rfactor, origin, 0.5);
	    print_effective_atwood_num(&do_spikes, b_fn, grid, front, wave, 
				       percent, g_h_at_min,  g_h_at_max,
				       iext->rfactor, origin, 2.0/3.0);
	    print_effective_atwood_num(&do_spikes, b_fn, grid, front, wave, 
				       percent, g_h_at_min,  g_h_at_max,
				       iext->rfactor, origin, 1.0);
	}

        surface_area = NO;
        if(prt->surface_area == YES)
           surface_area = YES;

	if (is_io_node(myid))
	    print_intfc_extrema(grid,front,out,bst,amb_st,about_to_stop,surface_area);

	/* Now that we're finished with intfc extrema, we work on the
	   volume fraction levels, if any were specified.

	   We need to find the height NEAREST THE EDGE corresponding to a given
	   volume fraction.  Bisection is not robust because the volume
	   fraction may not be monotone; instead we need to start at the edge
	   and then move into the mixing zone. */

	if (iext->do_01 == YES)     /* 1% specified? */
	{
	    Big_State tmp_bst[2];
	    boolean   min_result, max_result;
	    double     dh = (g_h_at_max - g_h_at_min)/100.0;
	    double     frac_at_min[2], frac_at_max[2];

	    start_clock("find_1%%_levels");

	    h_at_min = g_h_at_min - dh;
	    h_at_max = g_h_at_max + dh;

	    if (iext->geom == planar)
	    {	
		h_at_min = new_height_at_fraction(front,g_h_at_min,dh,1,0.01,
						  &frac_at_min[index_at_min],
						  comp_at_min,iext->rfactor,
						  origin);
		h_at_max = new_height_at_fraction(front,g_h_at_max,dh,-1,0.01,
						  &frac_at_max[index_at_max],
						  comp_at_max,iext->rfactor,
						  origin);

		/* The following added 11 Nov. 2003 to go with the if  
		   statement below controlling print_data_at_height(). Ensures 
		   that printing goes ahead even if new_height_at_fraction()
		   is used as opposed to height_at_fraction() - egeorge.
		*/
		min_result = FUNCTION_SUCCEEDED;
		max_result = FUNCTION_SUCCEEDED;
	    }
	    else if (iext->geom == spherical)
	    {
		min_result = height_at_fraction(front,&h_at_min,dh,1,0.01,
						frac_at_min,index_at_min,
						iext->rfactor,origin);
		max_result = height_at_fraction(front,&h_at_max,dh,-1,0.01,
						frac_at_max,index_at_max,
						iext->rfactor,origin);
	    }

	    if (debugging("gfrac"))
	    {
	        (void) printf("After calls to height_at_fraction():\n");
		(void) printf("\th_at_min = %g; frac[%d] = %g\n",
		              h_at_min,index_at_min,frac_at_min[index_at_min]);
		(void) printf("\th_at_max = %g; frac[%d] = %g\n",
		              h_at_max,index_at_max,frac_at_max[index_at_max]);
	    }

	    accumulate_state_in_layer(wave,front,h_at_min,tmp_bst,NULL,
				      iext->rfactor,origin,NO,YES,YES,NO);

	    copy_Big_State(&tmp_bst[index_at_min],&bst[0]);
	    copy_Big_State(&tmp_bst[1-index_at_min],&amb_st[0]);

	    accumulate_state_in_layer(wave,front,h_at_max,tmp_bst,NULL,
				      iext->rfactor,origin,NO,YES,YES,NO);

	    copy_Big_State(&tmp_bst[index_at_max],&bst[1]);
	    copy_Big_State(&tmp_bst[1-index_at_max],&amb_st[1]);

	    if (debugging("gfrac"))
	    {
	        (void) printf("After calls to accumulate_state_in_layer():\n");
		(void) printf("\t\tAt min, frac[%d] = %g\n",
			      index_at_min,bst[0].frac);
		(void) printf("\t\tAt max, frac[%d] = %g\n",
			      index_at_max,bst[1].frac);
	    }

	    /* Following if statement added 24 Oct. 2003 - egeorge */
	    if (eff_at_print_freq > 0 && (general_counter % eff_at_print_freq == 0))
	    {	    
	        percent = 1.0;
		print_effective_atwood_num(&do_spikes, b_fn, grid, front, 
					   wave, percent, h_at_min, h_at_max,
					   iext->rfactor, origin, 1.0/3.0);
		print_effective_atwood_num(&do_spikes, b_fn, grid, front, 
					   wave, percent, h_at_min, h_at_max,
					   iext->rfactor, origin, 0.5);
		print_effective_atwood_num(&do_spikes, b_fn, grid, front, 
					   wave, percent, h_at_min, h_at_max, 
					   iext->rfactor, origin, 2.0/3.0);
		print_effective_atwood_num(&do_spikes, b_fn, grid, front,
					   wave, percent, h_at_min, h_at_max,
					   iext->rfactor, origin, 1.0);
	    }

	    if (is_io_node(myid) &&
		(min_result == FUNCTION_SUCCEEDED) &&
		(max_result == FUNCTION_SUCCEEDED) &&
		(bst[0].count > 0) &&
		(bst[1].count > 0)) 
	        print_data_at_height(grid,front,out,bst,amb_st,
				     iext->out_min1,iext->out_max1,
				     iext->out_amp1,h_at_min,
				     h_at_max,about_to_stop);

	    stop_clock("find_1%%_levels");

	} /* if (iext->do_01 == YES) */

	if (iext->do_05 == YES)     /* 5% specified? */
	{
	    double dh = (g_h_at_max - g_h_at_min)/100.0;
	    double frac_at_min[2], frac_at_max[2];
	    Big_State tmp_bst[2];
	    boolean min_result, max_result;

	    start_clock("find_5%%_levels");

	    h_at_min = g_h_at_min - dh;
	    h_at_max = g_h_at_max + dh;

	    if (iext->geom == planar)
	    {	
		h_at_min = new_height_at_fraction(front,g_h_at_min,dh,1,0.05,
						  &frac_at_min[index_at_min],
						  comp_at_min,iext->rfactor,
						  origin);
		h_at_max = new_height_at_fraction(front,g_h_at_max,dh,-1,0.05,
						  &frac_at_max[index_at_max],
						  comp_at_max,iext->rfactor,
						  origin);

		/* The following added 11 Nov. 2003 to go with the if 
		   statement below controlling print_data_at_height(). Ensures 
		   that printing goes ahead even if new_height_at_fraction()
		   is used as opposed to height_at_fraction() - egeorge.
		*/
		min_result = FUNCTION_SUCCEEDED;
		max_result = FUNCTION_SUCCEEDED;
	    }
	    else if (iext->geom == spherical)
	    {
		min_result = height_at_fraction(front, &h_at_min, dh, 1, 0.05,
						frac_at_min, index_at_min,
						iext->rfactor,origin);
		max_result = height_at_fraction(front, &h_at_max, dh,-1, 0.05,
						frac_at_max, index_at_max,
						iext->rfactor,origin);
	    }

	    if (debugging("gfrac"))
	    {
	        printf("After calls to height_at_fraction():\n");
		printf("\th_at_min = %g; frac[%d] = %g\n",
		       h_at_min,index_at_min,frac_at_min[index_at_min]);
		printf("\th_at_max = %g; frac[%d] = %g\n",
		       h_at_max,index_at_max,frac_at_max[index_at_max]);
	    }

	    accumulate_state_in_layer(wave,front,h_at_min,
				      tmp_bst,NULL,
				      iext->rfactor,origin,
				      NO,YES,YES,NO);

	    copy_Big_State(&tmp_bst[index_at_min],&bst[0]);
	    copy_Big_State(&tmp_bst[1-index_at_min],&amb_st[0]);

	    accumulate_state_in_layer(wave,front,h_at_max,
				      tmp_bst,NULL,
				      iext->rfactor,origin,
				      NO,YES,YES,NO);

	    copy_Big_State(&tmp_bst[index_at_max],&bst[1]);
	    copy_Big_State(&tmp_bst[1-index_at_max],&amb_st[1]);

	    if (debugging("gfrac"))
	    {
	        (void) printf("After calls to accumulate_state_in_layer():\n");
		(void) printf("\t\tAt min, frac[%d] = %g\n",
			      index_at_min,bst[0].frac);
		(void) printf("\t\tAt max, frac[%d] = %g\n",
			      index_at_max,bst[1].frac);
	    }


	    /* Following if statement added 24 Oct. 2003 - egeorge */
	    if (eff_at_print_freq > 0 && (general_counter % eff_at_print_freq == 0))
	    {	    
	        percent = 5.0;
		print_effective_atwood_num(&do_spikes, b_fn, grid, front,
					   wave, percent, h_at_min, h_at_max,
					   iext->rfactor, origin, 1.0/3.0);
		print_effective_atwood_num(&do_spikes, b_fn, grid, front,
					   wave, percent, h_at_min, h_at_max,
					   iext->rfactor, origin, 0.5);
		print_effective_atwood_num(&do_spikes, b_fn, grid, front,
					   wave,percent, h_at_min, h_at_max,
					   iext->rfactor, origin, 2.0/3.0);
		print_effective_atwood_num(&do_spikes, b_fn, grid, front,
					   wave, percent, h_at_min, h_at_max,
					   iext->rfactor, origin, 1.0);
	    }




	    if (is_io_node(myid) &&
		(min_result == FUNCTION_SUCCEEDED) &&
		(max_result == FUNCTION_SUCCEEDED) &&
		(bst[0].count > 0) &&
		(bst[1].count > 0))
	        print_data_at_height(grid,front,out,bst,amb_st,
			       iext->out_min5,iext->out_max5,
			       iext->out_amp5,h_at_min,
			       h_at_max,about_to_stop);


	    stop_clock("find_5%%_levels");

	} /* if (iext->do_05 == YES) */

	/* Incrementing of general_counter added 24 Oct. 2003 - egeorge */
	general_counter++;

	stop_clock("record_intfc_extrema");
	debug_print("gintext","Left record_intfc_extrema()\n");
}		/*end record_intfc_extrema*/

LOCAL   void print_intfc_extrema(
	const Grid      *grid,
	const Front     *front,
	OUTPUT_DATA     *out,
	const           Big_State bst[2],
	const           Big_State amb_st[2],
	boolean            about_to_stop,
        boolean            surface_area)
{
	const Intfc_extrema *iext = (Intfc_extrema*)out;
	const int 	    dim = front->rect_grid->dim;
	double		    vel_at_max[MAXD], vel_at_min[MAXD];
	double		    amb_vel_at_max[3], amb_vel_at_min[3];
	double		    sph_pos_at_min[3], sph_pos_at_max[3];
	double		    sph_vel_at_min[3], sph_vel_at_max[3];
	double		    rhat_min[3], thhat_min[3], phhat_min[3];
	double		    rhat_max[3], thhat_max[3], phhat_max[3];
	register int	    i;
	static int	    n_min, n_max, n_amp;
	static double	    *fmax = NULL, *fmin = NULL, *famp = NULL;
	FILE                *omin = iext->out_min;
	FILE                *omax = iext->out_max;
	FILE                *oamp = iext->out_amp;

	debug_print("gintext","Entered print_intfc_extrema()\n");
	start_clock("print_intfc_extrema");

	if (fmin == NULL)
	{
	    n_min = n_max = 10;
	    n_amp = 3+4*dim;
	    uni_array(&fmin,n_min,FLOAT);
	    uni_array(&fmax,n_max,FLOAT);
	    uni_array(&famp,n_amp,FLOAT);
	}

	for (i = 0; i < dim; i++)
	{
	    vel_at_min[i] = bst[0].v[i];
	    vel_at_max[i] = bst[1].v[i];
	    amb_vel_at_min[i] = amb_st[0].v[i];
	    amb_vel_at_max[i] = amb_st[1].v[i];
	}

	if (iext->geom == spherical)
	{
	    convert_to_spherical(sph_pos_at_min,iext->pos[0],iext->origin,dim);
	    convert_to_spherical(sph_pos_at_max,iext->pos[1],iext->origin,dim);
	}

	fmin[0] = fmax[0] = famp[0] = grid->time;

	switch(iext->geom)
	{
	case planar:
	    fmin[1] = iext->pos[0][dim-1];
	    fmax[1] = iext->pos[1][dim-1];
  	    fmin[2] = vel_at_min[dim-1];
	    fmax[2] = vel_at_max[dim-1];
	    fmin[3] = amb_vel_at_min[dim-1];
	    fmax[3] = amb_vel_at_max[dim-1];
	    famp[1] = 0.5*(iext->pos[1][dim-1] - iext->pos[0][dim-1]);
	    famp[2] = 0.5*(vel_at_max[dim-1] - vel_at_min[dim-1]);

	    for (i = 1; i <= dim; i++)
	    {
	        famp[2+i] = iext->pos[0][i-1];
		famp[2+dim+i] = iext->pos[1][i-1];
		famp[2+2*dim+i] = vel_at_min[i-1];
		famp[2+3*dim+i] = vel_at_max[i-1];
	    }
	    break;

	case spherical:
	    spherical_unit_vectors(rhat_min,thhat_min,phhat_min,
				   iext->pos[0],iext->origin,dim);
	    spherical_unit_vectors(rhat_max,thhat_max,phhat_max,
				   iext->pos[1],iext->origin,dim);
	    sph_vel_at_min[0] = scalar_product(vel_at_min,rhat_min,dim);
	    sph_vel_at_max[0] = scalar_product(vel_at_max,rhat_max,dim);
	    sph_vel_at_min[1] = scalar_product(vel_at_min,phhat_min,dim);
	    sph_vel_at_max[1] = scalar_product(vel_at_max,phhat_max,dim);
	    if (dim == 3)
	    {
	        sph_vel_at_min[2] = scalar_product(vel_at_min,thhat_min,dim);
	        sph_vel_at_max[2] = scalar_product(vel_at_max,thhat_max,dim);
	    }
	    fmin[1] = sph_pos_at_min[0];
	    fmax[1] = sph_pos_at_max[0];
	    fmin[2] = sph_vel_at_min[0];
	    fmax[2] = sph_vel_at_max[0];
	    fmin[3] = amb_vel_at_min[0];
	    fmax[3] = amb_vel_at_max[0];
	    famp[1] = 0.5*(sph_pos_at_max[0] - sph_pos_at_min[0]);
	    famp[2] = 0.5*(fmax[2] - fmin[2]);

	    for (i = 1; i <= dim; i++)
	    {
	        famp[2+i] = sph_pos_at_min[i-1];
		famp[2+dim+i] = sph_pos_at_max[i-1];
		famp[2+2*dim+i] = sph_vel_at_min[i-1];
		famp[2+3*dim+i] = sph_vel_at_max[i-1];
	    }
	    break;
	}

	fmin[4] = bst[0].d;
	fmax[4] = bst[1].d;
	fmin[5] = amb_st[0].d;
	fmax[5] = amb_st[1].d;
	fmin[6] = bst[0].p;
	fmax[6] = bst[1].p;
	fmin[7] = amb_st[0].p;
	fmax[7] = amb_st[1].p;
	fmin[8] = bst[0].e;
	fmax[8] = bst[1].e;
	fmin[9] = amb_st[0].e;
	fmax[9] = amb_st[1].e;

	for (i = 0; i < n_min; i++)
	    (void) fprintf(omin,"%15.5e",fmin[i]);
	if (debugging("gintext"))
	    (void) fprintf(omin,"%10d",amb_st[0].count);
	(void) fprintf(omin,"\n");

	for (i = 0; i < n_max; i++)
	    (void) fprintf(omax,"%15.5e",fmax[i]);
	if (debugging("gintext"))
	    (void) fprintf(omax,"%10d",amb_st[1].count);
	(void) fprintf(omax,"\n");

	for (i = 0; i < n_amp; i++)
	    (void) fprintf(oamp,"%15.5e",famp[i]);

        if(surface_area == YES)
	   (void) fprintf(oamp,"%15.5e",compute_intfc_area(front->interf));

	(void) fprintf(oamp,"\n");

	if (about_to_stop == YES)
	{
	    (void) fprintf(omin,"\n");
	    (void) fprintf(omax,"\n");
	    (void) fprintf(oamp,"\n");
	}

	stop_clock("print_intfc_extrema");
	debug_print("gintext","Left print_intfc_extrema()\n");
}		/*end print_intfc_extrema*/

LOCAL   void    print_intfc_extrema_headers(
	const Intfc_extrema*    iext,
	int                     dim)
{
	FILE *omax = iext->out_max;
	FILE *omin = iext->out_min;
	FILE *oamp = iext->out_amp;
	char xyz[3] = {'x', 'y', 'z'};
	char rtp[3] = {'r', 't', 'p'};
	int  i;

	(void) foutput(omin);
	(void) foutput(omax);
	(void) foutput(oamp);
	(void) fprintf(omin,"\n# %13s","time");
	(void) fprintf(omax,"\n# %13s","time");
	(void) fprintf(oamp,"\n# %13s","time");

	(void) fprintf(oamp," %14s %14s", "a", "adot");

	switch(iext->geom)
	{
	case planar:
	    for (i = 0; i < dim; i++)
	        (void) fprintf(oamp,"          %c_%3s",xyz[i],"min");
	    for (i = 0; i < dim; i++)
	        (void) fprintf(oamp,"          %c_%3s",xyz[i],"max");
	    for (i = 0; i < dim; i++)
	        (void) fprintf(oamp,"         v%c_%3s",xyz[i],"min");
	    for (i = 0; i < dim; i++)
	        (void) fprintf(oamp,"         v%c_%3s",xyz[i],"max");
	    break;
	case spherical:
	    for (i = 0; i < dim; i++)
	        (void) fprintf(oamp,"          %c_%3s",rtp[i],"min");
	    for (i = 0; i < dim; i++)
	        (void) fprintf(oamp,"          %c_%3s",rtp[i],"max");
	    for (i = 0; i < dim; i++)
	        (void) fprintf(oamp,"         v%c_%3s",rtp[i],"min");
	    for (i = 0; i < dim; i++)
	        (void) fprintf(oamp,"         v%c_%3s",rtp[i],"max");
	    break;
	}

	fprintf(oamp,"   intf_surf_area");

	(void) fprintf(omin,"          %c_%3s",'h',"min");
	(void) fprintf(omax,"          %c_%3s",'h',"max");

	(void) fprintf(omin," %14s %14s %14s %14s",
		       "vel", "amb_vel", "den", "amb_den");
	(void) fprintf(omin," %14s %14s", "pre", "amb_pre");
	(void) fprintf(omax," %14s %14s %14s %14s",
		       "vel", "amb_vel", "den", "amb_den");
	(void) fprintf(omax," %14s %14s", "pre", "amb_pre");

	(void) fprintf(omin," %14s %14s", "ien", "amb_ien");
	(void) fprintf(omax," %14s %14s", "ien", "amb_ien");

	(void) fprintf(omin,"\n");
	(void) fprintf(omax,"\n");
	(void) fprintf(oamp,"\n");
}		/*end print_intfc_extrema_headers*/

LOCAL   void print_data_at_height(
	const Grid      *grid,
	const Front     *front,
	OUTPUT_DATA     *out,
	const Big_State bst[2],
	const Big_State amb_st[2],
	FILE		*omin,
	FILE		*omax,
	FILE		*oamp,
	double		h_at_min,
	double		h_at_max,
	boolean         about_to_stop)
{
	const Intfc_extrema *iext = (Intfc_extrema*)out;
	const int 	dim = front->rect_grid->dim;
	double		vel_at_max[3], vel_at_min[3];
	double		amb_vel_at_max[3], amb_vel_at_min[3];
	register int	i;
	static int	n_min, n_max, n_amp;
	static double	*fmax = NULL, *fmin = NULL, *famp = NULL;

	debug_print("gintext","Entered print_data_at_height()\n");
	start_clock("print_data_at_height");

	if (fmin == NULL)
	{
	    n_min = n_max = 10;
	    n_amp = 5 + 2*dim;
	    uni_array(&fmin,n_min,FLOAT);
	    uni_array(&fmax,n_max,FLOAT);
	    uni_array(&famp,n_amp,FLOAT);
	}

	for (i = 0; i < dim; i++)
	{
	    vel_at_min[i] = bst[0].v[i];
	    vel_at_max[i] = bst[1].v[i];
	    amb_vel_at_min[i] = amb_st[0].v[i];
	    amb_vel_at_max[i] = amb_st[1].v[i];
	}

	fmin[0] = fmax[0] = famp[0] = grid->time;
	fmin[1] = h_at_min;
	fmax[1] = h_at_max;

	switch(iext->geom)
	{
	case planar:
  	    fmin[2] = vel_at_min[dim-1];
	    fmax[2] = vel_at_max[dim-1];
	    fmin[3] = amb_vel_at_min[dim-1];
	    fmax[3] = amb_vel_at_max[dim-1];
	    break;

	case spherical:
	    fmin[2] = vel_at_min[0];
	    fmax[2] = vel_at_max[0];
	    fmin[3] = amb_vel_at_min[0];
	    fmax[3] = amb_vel_at_max[0];
	    break;
	}

	famp[1] = 0.5*(h_at_max - h_at_min);
	famp[2] = 0.5*(fmax[2] - fmin[2]);

	fmin[4] = bst[0].d;
	fmax[4] = bst[1].d;
	fmin[5] = amb_st[0].d;
	fmax[5] = amb_st[1].d;
	fmin[6] = bst[0].p;
	fmax[6] = bst[1].p;
	fmin[7] = amb_st[0].p;
	fmax[7] = amb_st[1].p;
	fmin[8] = bst[0].e;
	fmax[8] = bst[1].e;
	fmin[9] = amb_st[0].e;
	fmax[9] = amb_st[1].e;

	for (i = 1; i <= dim; i++)
	{
	    famp[2+i] = vel_at_min[i-1];
	    famp[2+dim+i] = vel_at_max[i-1];
	}
	famp[3+2*dim] = bst[0].frac;
	famp[4+2*dim] = bst[1].frac;

	for (i = 0; i < n_min; i++)
	    (void) fprintf(omin,"%15.5e",fmin[i]);
	if (debugging("gintext"))
	    (void) fprintf(omin,"%15.5e",bst[0].frac);
	(void) fprintf(omin,"\n");

	for (i = 0; i < n_max; i++)
	    (void) fprintf(omax,"%15.5e",fmax[i]);
	if (debugging("gintext"))
	    (void) fprintf(omax,"%15.5e",bst[1].frac);
	(void) fprintf(omax,"\n");

	for (i = 0; i < n_amp; i++)
	    (void) fprintf(oamp,"%15.5e",famp[i]);
	(void) fprintf(oamp,"\n");

	if (about_to_stop == YES)
	{
	    (void) fprintf(omin,"\n");
	    (void) fprintf(omax,"\n");
	    (void) fprintf(oamp,"\n");
	}

	stop_clock("print_data_at_height");
	debug_print("gintext","Left print_data_at_height()\n");
}		/*end print_data_at_height*/


LOCAL   void    print_data_at_height_headers(
	FILE *omin,
	FILE *omax,
	FILE *oamp,
	int  dim)
{
	int i;

	char xyz[3] = {'x', 'y', 'z'};

	(void) foutput(omin);
	(void) foutput(omax);
	(void) foutput(oamp);

	(void) fprintf(omin,"\n# %13s","time");
	(void) fprintf(omax,"\n# %13s","time");
	(void) fprintf(oamp,"\n# %13s","time");

	(void) fprintf(omin,"          %c_%3s",'h',"min");
	(void) fprintf(omax,"          %c_%3s",'h',"max");

	(void) fprintf(omin," %14s %14s %14s %14s",
		       "vel", "amb_vel", "den", "amb_den");
	(void) fprintf(omin," %14s %14s", "pre", "amb_pre");
	(void) fprintf(omax," %14s %14s %14s %14s",
		       "vel", "amb_vel", "den", "amb_den");
	(void) fprintf(omax," %14s %14s", "pre", "amb_pre");

	(void) fprintf(oamp," %14s %14s", "a", "adot");

	for (i = 0; i < dim; i++)
	    (void) fprintf(oamp,"         v%c_%3s",xyz[i],"min");
	for (i = 0; i < dim; i++)
	    (void) fprintf(oamp,"         v%c_%3s",xyz[i],"max");

	(void) fprintf(oamp," %14s %14s", "lower_frac", "upper_frac");

	(void) fprintf(omin," %14s %14s", "ien", "amb_ien");
	(void) fprintf(omax," %14s %14s", "ien", "amb_ien");

	(void) fprintf(omin,"\n");
	(void) fprintf(omax,"\n");
	(void) fprintf(oamp,"\n");
}		/*end print_data_at_height_headers*/

LOCAL	void	zero_state_totals(
	Big_State          *bst,
	Turbulence_moments *turb,
	int                n_params)
{
	const size_t 	size_bst = sizeof(Big_State);
	const size_t	size_turb = sizeof(Turbulence_moments);
	int 		p;

	if (bst != NULL)
	    for (p = 0; p < n_params; p++)
	        zero_scalar(&bst[p],size_bst);
	if (turb != NULL)
	    for (p = 0; p < n_params; p++)
	        zero_scalar(&turb[p],size_turb);
}		/*end zero_state_totals*/

LOCAL	void 	add_state_to_totals(
	const Locstate     st,
	Big_State          *bst,
	Turbulence_moments *turb)
{
	int	i, j;
	double	d, p, k, e;
	double	v[MAXD];

	bst->count++;
	if (is_obstacle_state(st))
	    return;

	d = Dens(st);
	p = pressure(st);
	k = kinetic_energy(st);
	e = internal_energy(st);
	bst->d += d;
	for (i = 0; i < MAXD; i++)
	{
	    v[i] = vel(i,st);
	    bst->v[i] += v[i];
	}
	bst->p += p;
	bst->k += k;
	bst->e += e;
	if (turb != NULL)
	{
	    turb->dd += d*d;
	    turb->de += d*e;
	    for (i = 0; i < MAXD; i++)
	    {
	        turb->dv[i] += d*v[i];
		turb->dev[i] += d*e*v[i];
		turb->dkv[i] += d*k*v[i];
		turb->pv[i] += p*v[i];
		for (j = i; j < MAXD; j++)
		{
		    turb->dvv[sym_index(i,j)] += d*v[i]*v[j];
		    turb->vv[sym_index(i,j)] += v[i]*v[j];
		}
	    }
	}
}		/*end add_state_to_totals*/

EXPORT	void 	accumulate_state_totals(
	const Big_State          *bst1,
	const Turbulence_moments *turb1,
	Big_State                *bst2,
	Turbulence_moments       *turb2)
{
	int	i, j;

	if (bst2 != NULL && bst1 != NULL)
	{
	    bst2->count += bst1->count;
	    bst2->d += bst1->d;
	    for (i = 0; i < MAXD; i++)
	        bst2->v[i] += bst1->v[i];
	    bst2->p += bst1->p;
	    bst2->k += bst1->k;
	    bst2->e += bst1->e;
	}

	if (turb2 != NULL && turb1 != NULL)
	{
	    turb2->dd += turb1->dd;
	    turb2->de += turb1->de;
	    for (i = 0; i < MAXD; i++)
	    {
	        turb2->dv[i] += turb1->dv[i];
		turb2->dev[i] += turb1->dev[i];
		turb2->dkv[i] += turb1->dkv[i];
		turb2->pv[i] += turb1->pv[i];
		for (j = i; j < MAXD; j++)
		{
		    turb2->vv[sym_index(i,j)] += turb1->vv[sym_index(i,j)];
		    turb2->dvv[sym_index(i,j)] += turb1->dvv[sym_index(i,j)];
		}
	    }
	}
}		/*end accumulate_state_totals*/

EXPORT	void	normalize_state_totals(
	Big_State          *bst,
	Turbulence_moments *turb,
	int                n_points)
{
	int i, j;
	double nf;

	if (bst->count == 0)
	    return;

	nf = (double)bst->count;

	bst->d /= nf;
	bst->p /= nf;
	bst->e /= nf;
	bst->k /= nf;
	for (i = 0; i < MAXD; i++)
	    bst->v[i] /= nf;
	bst->frac = nf/n_points;

	if (turb != NULL)
	{
	    turb->dd /= nf;
	    turb->de /= nf;
	    for (i = 0; i < MAXD; i++)
	    {
	        turb->dv[i] /= nf;
		turb->dkv[i] /= nf;
		turb->dev[i] /= nf;
		turb->pv[i] /= nf;
		for (j = i; j < MAXD; j++)
		{
		    turb->vv[sym_index(i,j)] /= nf;
		    turb->dvv[sym_index(i,j)] /= nf;
		}
	    }
	}
}		/*end normalize_state_totals*/

/*
*			accumulate_state_in_layer():
*
*   Sum the gas states in the planar or spherical layer at the given height.
*
*   origin == NULL means that the layer is a plane of given height; otherwise
*   the layer is a spherical shell whose radius (height) is measured with
*   respect to origin[].
*
*   ignore_params == YES means to combine the sums for the different materials
*   into a single global sum (otherwise sums are kept for the individual
*   materials).
*
*   normalize == YES means convert the raw totals to averages prior to exiting.
*
*   communicate == YES means that, in the case of parallel computation, to
*   combine the totals from the different nodes.
*/

EXPORT	void 	accumulate_state_in_layer(
	Wave			*wave,
	Front			*front,
	double			height,
	Big_State		*bst,
	Turbulence_moments	*turb,
	double			rfactor,
	const double		*origin,
	boolean			ignore_params,
	boolean			normalize,
	boolean			communicate,
	boolean			do_atwood)
{
	INTERFACE		*intfc = front->interf;
	const RECT_GRID		*rgrid = front->rect_grid;
	const int		n_params = num_gas_params(intfc);
	const int		nn = pp_numnodes();
	const int		myid = pp_mynode();
	const int		dim = rgrid->dim;
	const int		zdir = dim - 1;
	const size_t		size_bst = sizeof(Big_State);
	const size_t		size_turb = sizeof(Turbulence_moments);
	const enum intfc_geometry  geom = (origin == NULL ? planar : spherical);
	const double		dh = (dim == 2 ?
				      min(rgrid->h[0],rgrid->h[1])/rfactor
				      : min(min(rgrid->h[0],rgrid->h[1]),
					    rgrid->h[2])/rfactor);
	int			Nx, Ny, Nphi;
	double			dx, dy, dphi;
	double			x_left, y_left, phi;
	register int 		i, j;
	int			k, n, p, n_points;
	double			coords[3], vtmp[3];
	double			rhat[MAXD],  phhat[MAXD];
	COMPONENT		comp;

	static Locstate		st = NULL;
	static Big_State	*commbst = NULL, *tmpbst = NULL;
	static Turbulence_moments *commturb = NULL, *tmpturb = NULL;


	/* Following two variables for computation of atwood3. 
	   NB - Assumes only 2 parameters. Used instead of the 
	   corresponding Big_State variables in the MPI_Allreduce() 
	   function calls nested in pp_global_min() and 
	   pp_global_max(). - egeorge
	*/
	double		        *heavy_dens, *light_dens;

	/* Remaining type declarations to help compute atwood5 
	   and atwood6 - egeorge */
	long                    *at5_count;
	long                    **vector_global_at5_count;
	double                   *sort_dens_l = NULL;
	double                   *sort_dens_h = NULL; 
	double                   *big_sort_dens_l = NULL; 
	double                   *big_sort_dens_h = NULL;
	long                    ii, jj;


	debug_print("glayer","Entered accumulate_state_in_layer(), h = %g\n",height);

	if (st == NULL)
	{
	    alloc_state(intfc,&st,front->sizest);
	    uni_array(&commbst,n_params,size_bst);
	    uni_array(&tmpbst,n_params,size_bst);
	    uni_array(&commturb,n_params,size_turb);
	    uni_array(&tmpturb,n_params,size_turb);
	}

	if (do_atwood == YES)
	{
	    uni_array(&heavy_dens,n_params,FLOAT);
	    uni_array(&light_dens,n_params,FLOAT);
	    uni_array(&at5_count,n_params,sizeof(long));
	    bi_array(&vector_global_at5_count,n_params,nn,sizeof(long));
	}


	if ((normalize == YES) && (communicate == NO))
	{
	    (void) printf("WARNING in accumulate_state_in_layer(), "
	                  "Normalization with no communication will give\n"
	                  "\tincorrect results in parallel computation.\n");
	}

	zero_state_totals(tmpbst,tmpturb,n_params);
	n_points = 0;

	/* Following for statement added 21 Aug. 2003 - egeorge */
	for (p = 0; p < n_params; p++)
	{
	    tmpbst[p].min_d = tmpbst[p].global_min_d = 1.e+6;
	    tmpbst[p].max_d = tmpbst[p].global_max_d = -1.e+6;
	    
	    bst[p].partial_d = tmpbst[p].partial_d = 0.0;
	    bst[p].median_d = 0.0;

	    if (do_atwood == YES)
	        at5_count[p] = 0;
	}


	switch(geom)
	{
	case planar:
	    coords[zdir] = height;
	    if (height < rgrid->L[zdir] || height >= rgrid->U[zdir])
	        break;
	    Nx = irint(rgrid->gmax[0]*rfactor);
	    dx = (rgrid->U[0] - rgrid->L[0])/Nx;
	    x_left = rgrid->L[0] + 0.5*dx;
	    switch(dim)
	    {
	    case 1:
	        screen("ERROR in accumulate_state_in_layer(), "
		       "1D not supported\n");
	        clean_up(ERROR);
	        break;
#if defined(TWOD)
	    case 2:
	        /* 
	       TODO test 2d case of atwood5 computation and do spherical case.

	       NOTE - see the corresponding lines of the 3d case below for 
	       comments. 6 Apr. 2004 egeorge.
		*/
		if (do_atwood == YES)
		{
		    for (i = 0; i < Nx; i++)
		    {
		        coords[0] = x_left + i*dx;
			comp = component(coords,intfc);
			hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,
				     front,wave,st,NULL);
			if (index_of_Gas_param(Params(st)) == 0) 
			  (at5_count[0])++;  /* heavy fluid */ 
			else if  (index_of_Gas_param(Params(st)) == 1) 
			  (at5_count[1])++;  /* light fluid */ 
		    }

		    uni_array(&sort_dens_l,at5_count[1],FLOAT);
		    uni_array(&sort_dens_h,at5_count[0],FLOAT);

		    for (i = 0; i < at5_count[1]; i++)
		        sort_dens_l[i] = HUGE;
		    for (i = 0; i < at5_count[0]; i++)
		        sort_dens_h[i] = -HUGE;
		    ii = jj = 0;
		}


		for (i = 0; i < Nx; i++)
		{
		    n_points++;
		    coords[0] = x_left + i*dx;
		    comp = component(coords,intfc);
		    hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,
				 front,wave,st,NULL);


		    if (do_atwood == YES)
		    {
		        tmpbst[index_of_Gas_param(Params(st))].min_d = 
			  min(Dens(st),
			      tmpbst[index_of_Gas_param(Params(st))].min_d);
			tmpbst[index_of_Gas_param(Params(st))].max_d = 
			  max(Dens(st),
			      tmpbst[index_of_Gas_param(Params(st))].max_d);

			/* heavy fluid */
			if (index_of_Gas_param(Params(st)) == 0) 
			{
			    sort_dens_h[ii] = Dens(st);
			    ii++;
			}
			/* light fluid */
			else if  (index_of_Gas_param(Params(st)) == 1) 
			{ 
			    sort_dens_l[jj] = Dens(st);
			    jj++;
			}
		    }

		    add_state_to_totals(st,
				   &tmpbst[index_of_Gas_param(Params(st))],
				   &tmpturb[index_of_Gas_param(Params(st))]);
		}


		if (do_atwood == YES)
		{
		    if ((communicate == YES) && (nn > 1))
		    {
		        for (p = 0; p < n_params; p++)
			{
			    pp_l_allgather(&at5_count[p], 1, 
					   vector_global_at5_count[p], 1);

			    pp_global_lsum(&at5_count[p], 1);
			}
		    }

		    uni_array(&big_sort_dens_l,at5_count[1],FLOAT);
		    uni_array(&big_sort_dens_h,at5_count[0],FLOAT);
	
		    if ((communicate == YES) && (nn > 1))
		    {
		        pp_f_allgatherv(sort_dens_l, jj, big_sort_dens_l, 
					vector_global_at5_count[1]);
			pp_f_allgatherv(sort_dens_h, ii, big_sort_dens_h, 
					vector_global_at5_count[0]);
		    }
		    else if ((communicate == NO) && (nn == 1)) /*serial case*/
		    {
		        for (ii = 0; ii < at5_count[0]; ii++)
			    big_sort_dens_h[ii] = sort_dens_h[ii];
			for (jj = 0; jj < at5_count[1]; jj++)
			    big_sort_dens_l[jj] = sort_dens_l[jj];
		    }
		}
		break;
#endif /* defined(TWOD) */
#if defined(THREED)
	    case 3:
		Ny = irint(rgrid->gmax[1]*rfactor);
		dy = (rgrid->U[1] - rgrid->L[1])/Ny;
		y_left = rgrid->L[1] + 0.5*dy;


		if (do_atwood == YES)
		{
		    for (j = 0; j < Ny; j++)
		    {
		        coords[1] = y_left + j*dy;
			for (i = 0; i < Nx; i++)
			{
			    coords[0] = x_left + i*dx;
			    comp = component(coords,intfc);
			    hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,
					 front,wave,st,NULL);
			    if (index_of_Gas_param(Params(st)) == 0) 
			        (at5_count[0])++;	/* heavy fluid */  
			    else if (index_of_Gas_param(Params(st)) == 1) 
			        (at5_count[1])++;	/* light fluid */ 
			}
		    }  
	   
		    uni_array(&sort_dens_l,at5_count[1],FLOAT);
		    uni_array(&sort_dens_h,at5_count[0],FLOAT);

		    for (i = 0; i < at5_count[1]; i++)
		        sort_dens_l[i] = HUGE;
		    for (i = 0; i < at5_count[0]; i++)
		        sort_dens_h[i] = -HUGE;
		    ii = jj = 0;
		}


		for (j = 0; j < Ny; j++)
		{
		    coords[1] = y_left + j*dy;
		    for (i = 0; i < Nx; i++)
		    {
		        n_points++;
			coords[0] = x_left + i*dx;
			comp = component(coords,intfc);
			hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,
				     front,wave,st,NULL);
		    

			/* The assumption below with the min./max. layer 
			   density computation is that the heavy fluid index 
			   is always 0 and the light fluid index is always 1. 
			   So this assumes only two layers in the R-T problem.
			   21 Aug. 2003 - egeorge
			*/
			if (do_atwood == YES)
			{
			    tmpbst[index_of_Gas_param(Params(st))].min_d = 
			      min(Dens(st),
			      tmpbst[index_of_Gas_param(Params(st))].min_d);
			    tmpbst[index_of_Gas_param(Params(st))].max_d = 
			      max(Dens(st),
			      tmpbst[index_of_Gas_param(Params(st))].max_d);
			
			    /* heavy fluid */
			    if (index_of_Gas_param(Params(st)) == 0) 
			    {
			        sort_dens_h[ii] = Dens(st);
				ii++;			  
			    }
			    /* light fluid */ 
			    else if (index_of_Gas_param(Params(st)) == 1) 
			    {
			        sort_dens_l[jj] = Dens(st);
				jj++;	
			    }
			}


			add_state_to_totals(st,
				    &tmpbst[index_of_Gas_param(Params(st))],
				    &tmpturb[index_of_Gas_param(Params(st))]);
		    }
		}

		/* Following for atwood5 and atwood6 - egeorge */
		if (do_atwood == YES)
		{
		    if ((communicate == YES) && (nn > 1))
		    {
		        for (p = 0; p < n_params; p++)
			{
			    pp_l_allgather(&at5_count[p], 1, 
					   vector_global_at5_count[p], 1);
      
			    pp_global_lsum(&at5_count[p], 1);
			}
		    }

		    uni_array(&big_sort_dens_l,at5_count[1],FLOAT);
		    uni_array(&big_sort_dens_h,at5_count[0],FLOAT);


		    if ((communicate == YES) && (nn > 1))
		    {	
		        pp_f_allgatherv(sort_dens_l, jj, big_sort_dens_l, 
					vector_global_at5_count[1]);
			pp_f_allgatherv(sort_dens_h, ii, big_sort_dens_h, 
					vector_global_at5_count[0]);
		    }
		    else if ((communicate == NO) && (nn == 1)) /*serial case*/
		    {
		        for (ii = 0; ii < at5_count[0]; ii++)
			    big_sort_dens_h[ii] = sort_dens_h[ii];
			for (jj = 0; jj < at5_count[1]; jj++)
			    big_sort_dens_l[jj] = sort_dens_l[jj];
		    }
		}
		break;
#endif /* defined(THREED) */
	    }
	    break;

	case spherical:
	    switch(dim)
	    {
	    case 1:
	        screen("ERROR in accumulate_state_in_layer(), "
		       "1D not supported\n");
	        clean_up(ERROR);
	        break;
#if defined(TWOD)
	    case 2:
	        Nphi = irint(2.0*PI*height/dh);
		dphi = 2.0*PI/Nphi;
		for (i = 0; i < Nphi; i++)
		{
		    phi = i*dphi;
		    coords[0] = height*cos(phi) + origin[0];
		    coords[1] = height*sin(phi) + origin[1];
		    if (is_on_local_grid(rgrid,coords) != YES)
		      continue;
		    n_points++;
		    comp = component(coords,intfc);
		    hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,
				 front,wave,st,NULL);

		    /* Get components of velocity vector along
		       polar unit directions */

		    set_state(st,TGAS_STATE,st);
		    spherical_unit_vectors(rhat,NULL,phhat,coords,origin,dim);
		    vtmp[0] = scalar_product(Vel(st),rhat,dim);
		    vtmp[1] = scalar_product(Vel(st),phhat,dim);
		    for (k = 0; k < dim; k++)
		        Vel(st)[k] = vtmp[k];

		    add_state_to_totals(st,
				    &tmpbst[index_of_Gas_param(Params(st))],
				    &tmpturb[index_of_Gas_param(Params(st))]);
		}
	    break;
#endif /* defined(TWOD) */
#if defined(THREED)
	    case 3:
		{
		    double th, dth, rpolar, thhat[3];
	            int j, Nth;
		    Nth = irint(PI*height/dh);
		    dth = PI/Nth;
		    for (j = 0; j < Nth; j++)
		    {
		        th = j*dth;
		        rpolar = height*sin(th);
		        Nphi = irint(2.0*PI*rpolar/dh);
		        dphi = 2.0*PI/Nphi;
		        for (i = 0; i < Nphi; i++)
		        {
			    phi = i*dphi;
			    coords[0] = rpolar*cos(phi) + origin[0];
			    coords[1] = rpolar*sin(phi) + origin[1];
			    coords[2] = height*cos(th) + origin[2];
			    if (is_on_local_grid(rgrid,coords) != YES)
			        continue;
			    n_points++;
			    comp = component(coords,intfc);
			    hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,
				         front,wave,st,NULL);

			    /* Get components of velocity vector along
			       spherical unit directions */

			    set_state(st,TGAS_STATE,st);
			    spherical_unit_vectors(rhat,thhat,phhat,
					           coords,origin,dim);
			    vtmp[0] = scalar_product(Vel(st),rhat,dim);
			    vtmp[1] = scalar_product(Vel(st),thhat,dim);
			    vtmp[2] = scalar_product(Vel(st),phhat,dim);
			    for (k = 0; k < dim; k++)
			        Vel(st)[k] = vtmp[k];

			    add_state_to_totals(st,
				    &tmpbst[index_of_Gas_param(Params(st))],
				    &tmpturb[index_of_Gas_param(Params(st))]);
		        }
		    }
		}
		break;
#endif /* defined(THREED) */
	    }
	    break;
	}

	if (debugging("glayer"))
	{
	    (void) printf("\tI found %d points in this layer.\n"
	                  "\th = %g; geom = %s\n",n_points,height,
		          (geom == spherical ? "spherical" : "planar"));
	}


	if (do_atwood == YES)
	{
	    if ((communicate == YES) && (nn > 1)) 
	    {
	        for (p = 0; p < n_params; p++)
		{	
		    /* Following added 21 Aug. 2003 - egeorge */
		    light_dens[p] = tmpbst[p].min_d;
		    heavy_dens[p] = tmpbst[p].max_d;

		    pp_global_min(&light_dens[p], 1); 		
		    pp_global_max(&heavy_dens[p], 1); 
		}
	    }
	}



	if ((communicate == YES) && (nn > 1))
	{
	    if (myid == IO_NODE_ID)
	    {
	        for (n = 0; n < nn; n++)
		{
		    if (n != myid)
		    {
		        pp_recv(LAYER_SUM_ID,n,(POINTER)commbst,
				n_params*size_bst);
			pp_recv(LAYER_SUM_ID+1,n,(POINTER)commturb,
				n_params*size_turb);
			for (p = 0; p < n_params; p++)
			    accumulate_state_totals(&commbst[p],&commturb[p],
						    &tmpbst[p],&tmpturb[p]);
		    }
		}
		pp_send_all(LAYER_SUM_ID+2,(POINTER)tmpbst,n_params*size_bst);
		pp_send_all(LAYER_SUM_ID+3,(POINTER)tmpturb,n_params*size_turb);
	    }
	    else
	    {
	        pp_send(LAYER_SUM_ID,(POINTER)tmpbst,
			n_params*size_bst,IO_NODE_ID);
		pp_send(LAYER_SUM_ID+1,(POINTER)tmpturb,
			n_params*size_turb,IO_NODE_ID);
		pp_recv(LAYER_SUM_ID+2,IO_NODE_ID,(POINTER)tmpbst,
			n_params*size_bst);
		pp_recv(LAYER_SUM_ID+3,IO_NODE_ID,(POINTER)tmpturb,
			n_params*size_turb);
	    }

	    pp_global_isum(&n_points,1L);
	}

	if (ignore_params == YES)
	{
	    zero_state_totals(bst,turb,1);
	    for (p = 0; p < n_params; p++)
	        accumulate_state_totals(&tmpbst[p],&tmpturb[p],bst,turb);
	    if (normalize == YES)
	        normalize_state_totals(bst,turb,n_points);
	}
	else
	{
	    zero_state_totals(bst,turb,n_params);
	    if (turb == NULL)
	    {
		for (p = 0; p < n_params; p++)
		{
		    accumulate_state_totals(&tmpbst[p],NULL,&bst[p],NULL);
		    if (normalize == YES)
		        normalize_state_totals(&bst[p],NULL,n_points);
		}
	    }
	    else
	    {
		for (p = 0; p < n_params; p++)
		{
		    accumulate_state_totals(&tmpbst[p],&tmpturb[p],
					    &bst[p],&turb[p]);
		    if (normalize == YES)
		        normalize_state_totals(&bst[p],&turb[p],n_points);
		}
	    }
	}

	/* The following ft_assignments to bst[p] members must be made here or 
	   near the end of this function due to calls to zero_state_totals()
	   in the preceding few lines. 19 Nov. 2003 - egeorge 
	*/
	if (do_atwood == YES)
	{
	    if (ignore_params == NO)
	    {
	        for (p = 0; p < n_params; p++)
		{
		    bst[p].global_min_d = light_dens[p];
		    bst[p].global_max_d = heavy_dens[p];
		    /* DEBUGGING 23/8/2003 - egeorge   
		       (void)printf("GLOBAL MAX/MIN 1 den %g %g height %g\n",
		       bst[p].global_max_d,bst[p].global_min_d,height);	
		    */
		}
	    }
	
	    /* Following 17 Nov. 2003 for atwood5 computation - egeorge.
	       Sorting the heavy and light density arrays merged from the 
	       different processors then determining the average of the 
	       upper and lower (respectively) 50% of values.
	    */
	    qsort(big_sort_dens_l, at5_count[1], sizeof(double), 
		  cmp_for_qsort_inc);
	    qsort(big_sort_dens_h, at5_count[0], sizeof(double), 
		  cmp_for_qsort_dec);

	    for (ii = 0; ii < ((at5_count[0] + 1)/2); ii++)
	        bst[0].partial_d += big_sort_dens_h[ii];

	    for (ii = 0; ii < ((at5_count[1] + 1)/2); ii++)
	        bst[1].partial_d += big_sort_dens_l[ii];


	    for (p = 0; p < n_params; p++)
	        bst[p].partial_d /= ( (double)( (at5_count[p] + 1)/2) );

	    /* Now find the median density values in each phase - needed 
	       to compute atwood6.  Again, the assumption is that there are 
	       only two fluid phases. 
	       TODO - rewrite this function to remove that assumption. 
	       6 Apr. 2004 - egeorge.
	    */
	    if ((at5_count[0] % 2 == 0) && (at5_count[0] != 0))
	    {
	        bst[0].median_d = 0.5 * (big_sort_dens_h[(at5_count[0]/2)] 
				   + big_sort_dens_h[((at5_count[0]/2) - 1)]);
	    }
	    else
	    {
	        bst[0].median_d = big_sort_dens_h[((at5_count[0] - 1) / 2)];
	    }

	    if ((at5_count[1] % 2 == 0) && (at5_count[1] != 0))
	    {
	        bst[1].median_d = 0.5 * (big_sort_dens_l[(at5_count[1]/2)] 
				   + big_sort_dens_l[((at5_count[1]/2) - 1)]);
	    }
	    else
	    {
	        bst[1].median_d = big_sort_dens_l[((at5_count[1] - 1) / 2)];
	    }


	    free(heavy_dens);
	    free(light_dens);
	    free(at5_count);
	    free(sort_dens_l); 
	    free(sort_dens_h);
	    free(big_sort_dens_l); 
	    free(big_sort_dens_h);
	    free(vector_global_at5_count);
	} /* end of 'if (do_atwood == YES)' */

	debug_print("glayer","Left accumulate_state_in_layer()\n");
}			/*end accumulate_state_in_layer*/




LOCAL   void   init_vol_frac (
	int	              rw,
	int	              *do_vol_frac,
	int                   *print_freq,
	char                  *fn)
{
        static int    s_do_vol_frac = NO;
	static int    s_print_freq = 100;
	static char   s_fn[256];
	char          s[Gets_BUF_SIZE];

	if (rw == 0)
	{
	    sprintf(s_fn, "FronTier");
      

	    screen("Type 'y' to print volume fractions: ");
	    (void) Gets(s);
	    if (s[0] == 'Y' || s[0] == 'y')
	    {	    
	        s_do_vol_frac = YES;
	        screen("The volume fractions will be printed at a\n");
		screen("\tfixed integer multiple of  the interface extrema data\n");
		screen("\tprinting frequency.\n");
		screen("\tEnter a POSITIVE INTEGER choice for this multiple\n"
		       "\t\t(default = %d): ", s_print_freq);
		(void) Gets(s);
		if (s[0] != '\0')
		    (void) sscanf(s,"%d", &s_print_freq);

		screen("All volume fraction output files will be stored\n");
		screen("\tin the same directory under the names, \n");
		screen("\t'[basename]-vol_frac-tXXXX.YY',\n");
		screen("\twhere XXXX.YY is the simulation time.\n");
		screen("\tEnter a choice for the basename of the output files\n");
		screen("\t(give the full directory path AND create any directories\n");
		screen("\tin this path which do not already exist)\n"
		       "\t\t(default = %s): ", s_fn);
		(void) Gets(s);
		if (s[0] != '\0')
		    (void) strcpy(s_fn, s);
		screen("\n");
		{
		    FILE *ifp;
		    char tmp_filename[28];
		    
		    sprintf(tmp_filename,"%s-dtp-tMP_file-xyz.%d",
			    s_fn,pp_mynode());
		    if ((ifp = fopen(tmp_filename,"w")) == NULL)
		    {
		        fclose(ifp);
			remove(tmp_filename);
		        screen("ERROR in init_vol_frac() "
			       "Cannot open file with base name %s.\n",s_fn);
			clean_up(ERROR);
		    }
		    else
		    {
		        fclose(ifp);
			remove(tmp_filename);
		    }
		}
	    }
	}
	else if (rw == 1)
	{
	    *do_vol_frac = s_do_vol_frac;
	    *print_freq = s_print_freq;
	    sprintf(fn, s_fn);
	}

}  /*end init_vol_frac */


LOCAL   void   init_print_allstates (
	int	              rw,
	int	              *print_allstates,
	int                   *print_freq,
	char                  *fn)
{
        static int    s_print_allstates = NO;
	static int    s_print_freq = 100;
	static char   s_fn[256];
	char          s[Gets_BUF_SIZE];

	if (rw == 0)
	{
	    sprintf(s_fn, "FronTier");
      

	    screen("Type 'y' to print\n");
	    screen("averaged mixing zone state data as a function of Z and t: ");
	    (void) Gets(s);
	    if (s[0] == 'Y' || s[0] == 'y')
	    {
	        s_print_allstates = YES;
	        screen("The state data will be printed at a\n");
		screen("\tfixed integer multiple of  the interface extrema data\n");
		screen("\tprinting frequency.\n");
		screen("\tEnter a POSITIVE INTEGER choice for this multiple\n"
		       "\t\t(default = %d): ", s_print_freq);
		(void) Gets(s);
		if (s[0] != '\0')
		    (void) sscanf(s,"%d", &s_print_freq);

		screen("All state data output files will be stored\n");
		screen("\tin the same directory and will be named, \n");
		screen("\t'[basename]-allstates-tsX-tY',\n");
		screen("\twhere X is the printing timestep and\n");
		screen("\tY is the printing time.\n");
		screen("\tEnter a choice for the basename of the output files\n");
		screen("\t(give the full directory path AND create any directories\n");
		screen("\tin this path which do not already exist)\n"
		       "\t\t(default = %s): ", s_fn);
		(void) Gets(s);
		if (s[0] != '\0')
		    (void) strcpy(s_fn, s);
		screen("\n\n");
		{
		    FILE *ifp;
		    char tmp_filename[28];
		    
		    sprintf(tmp_filename,"%s-h57-TmP_fIlE-xyz.%d",s_fn,pp_mynode());
		    if ((ifp = fopen(tmp_filename,"w")) == NULL)
		    {
		        fclose(ifp);
		        remove(tmp_filename);
		        screen("ERROR in init_print_allstates() "
			       "Cannot open file with base name %s.\n",s_fn);
			clean_up(ERROR);
		    }
		    else
		    {
		        fclose(ifp);
		        remove(tmp_filename);
		    }
		}
	    }
	}
	else if (rw == 1)
	{
	    *print_allstates = s_print_allstates;
	    *print_freq = s_print_freq;
	    sprintf(fn, s_fn);
	}

}  /*end init_print_allstates */


/* 
*			old_print_vol_frac():
*
*       Computes and prints volume fraction data as a function of domain 
*       height at pre-determined time step intervals. Much of the details 
*       of the algorithms used are contained in Chapter 4 of Erwin George's 
*       doctoral dissertation. The key difference between this function and 
*       print_vol_frac() is that this function uses a more discrete (and 
*       less accurate) algorithm for determining the volume fractions.
*
*       This function is presently only defined for 2d and 3d planar
*       geometry.
*
*       TODO: Modify following function and/or print_vol_frac() to deal 
*       with spherical geometry - 19 Apr. 2004. 
*/

LOCAL	void 	old_print_vol_frac(
	Wave			*wave,
	Front			*front,	
	double			h_min,
	double			h_max,
	int			dir,
	int			index,
	double			rfactor,
	const double		*origin,
	char                    *b_fn)
{
        FILE                    *OLD_VOLFRAC;
	double                   height;
	double                   frac[2];
	double                   dh = (h_max - h_min)/100.0;
	char                    filename[256];
	int                     vert_dir = front->rect_grid->dim - 1;
	int                     my_id = pp_mynode();
	const enum intfc_geometry  geom = (origin == NULL ? planar : spherical);



	
	sprintf(filename,"%s-OLDvol_frac-t%07.2f",b_fn,front->time);
	if (my_id == 0)
	    OLD_VOLFRAC = fopen(filename, "w");


	if (geom == planar)
	{
	    height = front->rect_grid->GL[vert_dir] + dh;
	    
	    sprintf(filename,"%s-old_allstates-ts%d-t%07.2f",
		    b_fn,front->step,front->time);

	    while (height < front->rect_grid->GU[vert_dir])
	    {
		accumulate_fractions_in_layer(front,height,frac,
					      rfactor,origin);

		if (my_id == 0)
		    fprintf(OLD_VOLFRAC,"%lg  %lg\n",height,frac[1]);
		
		height += dir*dh;
	    }
	}
	else /* if geom == spherical */
	{
	    return;
	}

	
	if (my_id == 0)
	{
	    fclose(OLD_VOLFRAC);
	}
	return;
}  /*end old_print_vol_frac*/



LOCAL   void    init_print_effective_atwood_num(int   rw, 
						int   step,
						Front *front,
						int   *do_spikes,
						int   *print_freq,
						char  *b_fn)
{
	char           s[Gets_BUF_SIZE];
        static char    s_fn[256];
	static int     s_print_freq = 100;
	static int     s_do_spikes = NO;

	if (rw == 0)
	{
	    sprintf(s_fn, "FronTier");
      
	    screen("Type 'y' to also print Effective Atwood Number data: ");
	    (void) Gets(s);
	    if (s[0] == 'Y' || s[0] == 'y')
	    {
	        screen("The Effective Atwood Number data will be printed at a\n");
		screen("\tfixed integer multiple of  the interface extrema data\n");
		screen("\tprinting frequency.\n");
		screen("\tEnter a POSITIVE INTEGER choice for this multiple\n"
		       "\t\t(default = %d): ", s_print_freq);
		(void) Gets(s);
		if (s[0] != '\0')
		    (void) sscanf(s,"%d", &s_print_freq);

		screen("All Effective Atwood Number output files will be stored\n");
		screen("\tin the same directory and will be named\n"); 
		screen("\t'[basename]-EFF-ATWOOD-0[015].00_topX.YY',\n");
		screen("\twhere X.YY = fraction of outer portion of bubble mixing\n"); 
		screen("\tzone used to compute the Atwood number.\n");
		screen("\tEnter a choice for the basename of the output files\n");
		screen("\t(give the full directory path AND create any directories\n");
		screen("\tin this path which do not already exist)\n"
		       "\t\t(default = %s): ", s_fn);
		(void) Gets(s);
		if (s[0] != '\0')
		    (void) strcpy(s_fn, s);
		{
		    FILE *ifp;
		    char tmp_filename[28];
		    
		    sprintf(tmp_filename,"%s-A2c-tMp_file-xyz.%d",s_fn,pp_mynode());
		    if ((ifp = fopen(tmp_filename,"w")) == NULL)
		    {
		        fclose(ifp);
		        remove(tmp_filename);
		        screen("ERROR in print_effective_atwood_num() "
			       "Cannot open file with base name %s.\n",s_fn);
			clean_up(ERROR);
		    }
		    else
		    {
		        fclose(ifp);
		        remove(tmp_filename);
		    }
		}
		screen("Effective Atwood number data for the bubbles will");
		screen(" be printed.\n");
		screen("\tType 'y' to also print data for the spikes\n"
		       "\t(these files will be similarly named as their "
		       "bubble counterparts,\n"
		       "\twith 'spike' preceding EFF-ATWOOD).\n"
		       "\t\t(default = no): ");
		(void) Gets(s);
		if (s[0] == 'Y' || s[0] == 'y')
		    s_do_spikes = YES;

		/* If this is a restart, the mean initial height of the 
		   interface must be obtained here, not in 
		   init_random_surface() which is not called on restart. */
		if (step > 0)
		{
		    double mean_height;

		    mean_height = 0.5* (
			   front->rect_grid->GU[front->interf->dim - 1] - 
			   front->rect_grid->GL[front->interf->dim - 1]);
		    screen("Since this is a restart, enter here the mean\n");
		    screen("\tposition of the front above L[%d].\n",
			   front->interf->dim - 1);
		    screen("\t\t(default = %lg): ", mean_height);
		    (void) Gets(s);
		    if (s[0] != '\0')
		        (void) sscanf(s,"%lg", &mean_height);
		    mean_height += front->rect_grid->GL[front->interf->dim-1];
		    (void) get_mean_initial_interface_height(0, &mean_height);
		}

		screen("\n");
	    }
	    /* 
	       if Effective Atwood Number data is NOT to be printed, the 
	       following else statement provides a convenient way of 
	       passing this information back to the calling function without 
	       having to introduce another static variable and/or another 
	       argument to print_effective_atwood_number().
	    */
	    else 
	        s_print_freq = -1;
	} /* end of if (rw == 0) */
	else if (rw != 0) /* pass the default or user-input data 
			     to the calling function */
	{
	    *print_freq = s_print_freq;
	    *do_spikes = s_do_spikes;
	    sprintf(b_fn, s_fn);
	}

}  /*end init_print_effective_atwood_num */



/*
                    print_effective_atwood_num():
* 
*   Computes and prints time dependent effective Atwood number data for 
*   the bubble portion, and, as an option, the spike portion of the mixing
*   zone.  NOTE that to avoid the use of an extra static variable in 
*   record_intfc_extrema(), the print interval is set to -1 during the 
*   first call to print_effective_atwood_num() IF printing of effective 
*   Atwood number data is turned OFF. This provides an easy check on 
*   whether to print Effective Atwood Number data in record_intfc_extrema(), 
*   using the single static variable which stores print frequency.
*   
*   There are currently (23 Apr. 2004) 5 versions of the effective Atwood
*   number computed and printed. An integer index (or lack thereof)   
*   immediately after the 'EFF-ATWOOD-' in the file name indicates the 
*   version present in that file. See also Erwin George's doctoral 
*   dissertation and other references cited below for more on this topic. 
*   The five versions of the effective Atwood number are:
*
*   1) atwood - no integer index in the file name. This was the original and 
*      simply uses the max. and min. planar density values to compute 
*      a time and Z dependent Atwood number, which is then averaged over 
*      the outer mix_zone_frac fraction of the mixing zone to produce an
*      Atwood number dependent on time only. This version was used in 
*      "A comparison of experimental, theoretical, and numerical simulation 
*      Rayleigh-Taylor mixing rates", E. George, J. Glimm, X. L. Li, and
*      Z. L. Xu, Proc. National Acad. of Sci. USA, 99:2587-2592, 2002.
*
*   In computing all subsequent effective Atwood numbers, one stays strictly
*   within the relevant fluid component when figuring out heavy fluid and 
*   light fluid densities to use in computing the local time and Z dependent
*   Atwood numbers (the averaging over height to produce a time dependent 
*   Atwood number remains the same as in 1) above).  Staying within the 
*   relevant components is important for capturing the attenuating effect of 
*   the increased density stratification which accompanies more compressible
*   simulations under isothermal initialization.
*
*   2) atwood2 - integer index 2 in the file name. Average of all densities 
*      in the heavy (light) fluid component is used as the heavy (light) fluid 
*      density in computing the local time and space dependent Atwood number.
*
*   3) atwood3 - integer index 3 in the file name. The maximum (minimum)  
*      density in the heavy (light) fluid component is used as the heavy 
*      (light) fluid density in computing the local time and space dependent 
*      Atwood number. (Essentially equivalent to 1) at low compressibility).
*
*   4) atwood5 - integer index 5 in the file name. The average of the largest 
*      (smallest) 50 percent of densities in the heavy (light) fluid
*      component is used as the heavy (light) fluid density in computing 
*      the local time and space dependent Atwood number. This version was 
*      used in "Self similarity of Rayleigh-Taylor mixing rates", E. 
*      George and J. Glimm, submitted to Phys. of Fluids, 2004.
*
*   5) atwood6 - integer index 6 in the file name. The median density value 
*      in the heavy (light) fluid component is used as the heavy (light) fluid 
*      density in computing the local time and space dependent Atwood number.
*
*   You may wonder why no atwood4. This was implemented in the TVD code using 
*   the respective averages of only those heavy and light fluid densities near 
*   to the interface, but it was not (as of 30 Apr. 2004) implemented in 
*   FronTier.
*
*   This function uses the initial print time and a multiple of the print 
*   interval entered in the input file for "interface extrema data". 
*   1 Feb. 2002. 
*   NOTE this function assumes no cutting in the vertical (Z for 3d, Y for 
*   2d) direction when decomposing the domain. - 11 Sep. 2001.
*
*   TODO - as of 23 Apr. 2004, the spherical geometry case. 
*/

LOCAL   void    print_effective_atwood_num(
					   int         *do_spikes,
					   char        *b_fn,
					   Grid*       grid,
					   Front*      front,
					   Wave*       wave,
					   double       percent,
					   double       bub_height,
					   double       spike_height,
					   double       rfactor,
					   const double *origin,
					   double       mix_zone_frac)
{
	double	         *L = front->rect_grid->GL, *U = front->rect_grid->GU;
	double            *coords;
	
	double            dens_minmax[2], atwood_sum;
	double            atwood;
	double            mix_zone_frac_outerbound;
	double            adjusted_bub_or_spike_height;
	double            mean_initial_interface_height;
	
	int              first; /* Tells whether the Effective Atwood Number 
				   output files have already been opened once, 
				   and hence should only be opened for 
				   "appending". This also has the effect of 
				   ensuring that the old Effective Atwood 
				   Number output files are not overwritten 
				   if the run is a restart.
				*/
	int              counter;
	int              atwood_counteradjuster; /* So that mixing zone planes 
						    in which there is no light 
						    fluid are not counted. */
	int              icoords[MAXD];
	int              h_extremum_index;
	int              i;
	int              ix, iy, iz;
	char             filename[256];
	FILE             *ATWOOD;	
	COMPONENT        comp;
	Locstate         state;
	
	double            atwood_sum2, atwood2;
	FILE             *ATWOOD2;
	char             filename2[256];
	Big_State        bst[2];

	double            atwood_sum3, atwood3;
	FILE             *ATWOOD3;
	char             filename3[256];

	double            atwood_sum5, atwood5;
	FILE             *ATWOOD5;
	char             filename5[256];

	double            atwood_sum6, atwood6;
	FILE             *ATWOOD6;
	char             filename6[256];

	const RECT_GRID *rgrid = wave->rect_grid;
	const double     *dh = rgrid->h;
	const int       dim = front->rect_grid->dim;
	int             vert_dir = dim -1;
	int             xmax = rgrid->gmax[0];
	int             ymax = rgrid->gmax[1];
	int             zmax = rgrid->gmax[2];



	/* _________________________________________________________ */
  

	if ((dim != 3 && dim != 2) || (origin != NULL))
	{
	    (void) printf("\nWARNING - print_effective_atwood_num() is ");
	    (void) printf("only implemented for planar 2d and 3d runs.\n"); 
	    return;
	}

	/* initializations */
	first = YES;
	atwood_counteradjuster = 0;
	counter = 0;
	dens_minmax[0] = HUGE;      
	dens_minmax[1] = -HUGE;
	atwood_sum = 0.0;
	atwood_sum2 = 0.0;
	atwood_sum3 = 0.0;
	atwood_sum5 = 0.0; 
	atwood_sum6 = 0.0;
	    
	mean_initial_interface_height = get_mean_initial_interface_height(
								    1, NULL);
	(void) printf("In print_effective_atwood_num() mean initial ");
	(void) printf("interface height = %lg\n", 
		      mean_initial_interface_height);

	/* The fraction of the outer portion of the (bubble) mixing 
	   zone used for averaging to determine the time dependent 
	   effective Atwood number is determined by the fifth argument 
	   to this function, mix_zone_frac. This mix_zone_frac is also  
	   interjected into the name of the specific Effective Atwood 
	   Number output file.
	   
	   NOTE that the light fluid is on top of the heavy fluid in 
	   FronTier Rayleigh-Taylor simulations, and this must be 
	   accounted for in computing the "top" mix_zone_frac of the 
	   mixing zone. 
	*/
	
	adjusted_bub_or_spike_height = 2.0*mean_initial_interface_height - bub_height;
	mix_zone_frac_outerbound = mean_initial_interface_height - 
	             (1.0 - mix_zone_frac) * 
	             (mean_initial_interface_height -  bub_height);
	    
	    
	/* Note that it is the centres of the cell blocks and not the 
	   edges which are being used in obtaining densities - 
	   egeorge 10/9/2001.
	*/
	if (dim == 3)
	    h_extremum_index = (int) ((zmax * (bub_height-(L[vert_dir] + 
						      0.5*dh[vert_dir]))) /
				 (U[vert_dir]-L[vert_dir]));
	else
	    h_extremum_index = (int) ((ymax * (bub_height-(L[vert_dir] + 
						      0.5*dh[vert_dir]))) /
				 (U[vert_dir]-L[vert_dir]));
	    

	    
	icoords[0] = icoords[1] = 0;
	icoords[vert_dir] = h_extremum_index;
	coords = Rect_coords(icoords,wave);
	while (coords[vert_dir] < bub_height)
	{
	    h_extremum_index++;
	    icoords[vert_dir] = h_extremum_index;
	    coords = Rect_coords(icoords,wave);
	}
	    
	    

	if (pp_mynode() == 0)
	{
	    if (front->time > 0 || front->step > 0)
	        first = NO;
      
	    sprintf(filename,"%s-EFF-ATWOOD-%05.2f_top%3.2f",
		    b_fn,percent,mix_zone_frac);
	    sprintf(filename2,"%s-EFF-ATWOOD-2-%05.2f_top%3.2f",
		    b_fn,percent,mix_zone_frac);
	    sprintf(filename3,"%s-EFF-ATWOOD-3-%05.2f_top%3.2f",
		    b_fn,percent,mix_zone_frac);
	    sprintf(filename5,"%s-EFF-ATWOOD-5-%05.2f_top%3.2f",
		    b_fn,percent,mix_zone_frac);
	    sprintf(filename6,"%s-EFF-ATWOOD-6-%05.2f_top%3.2f",
		    b_fn,percent,mix_zone_frac);     
	    if(first)
	    {
	        ATWOOD = fopen(filename, "w");
		ATWOOD2 = fopen(filename2, "w");
		ATWOOD3 = fopen(filename3, "w");
		ATWOOD5 = fopen(filename5, "w");
		ATWOOD6 = fopen(filename6, "w");
	    }
	    else
	    {
	        ATWOOD = fopen(filename, "a"); 
		ATWOOD2 = fopen(filename2, "a");
		ATWOOD3 = fopen(filename3, "a");
		ATWOOD5 = fopen(filename5, "a");
		ATWOOD6 = fopen(filename6, "a");
	    }
	    
	    
	    fprintf(ATWOOD,"  TIME    ");
	    fprintf(ATWOOD,"    Z        (adj. Z)   ");
	    fprintf(ATWOOD,"  EFF_AT_NO \n"); 
	    
	    fprintf(ATWOOD2,"  TIME    ");
	    fprintf(ATWOOD2,"    Z        (adj. Z)   ");
	    fprintf(ATWOOD2,"  EFF_AT_NO \n"); 

	    fprintf(ATWOOD3,"  TIME    ");
	    fprintf(ATWOOD3,"    Z        (adj. Z)   ");
	    fprintf(ATWOOD3,"  EFF_AT_NO \n"); 
	    
	    fprintf(ATWOOD5,"  TIME    ");
	    fprintf(ATWOOD5,"    Z        (adj. Z)   ");
	    fprintf(ATWOOD5,"  EFF_AT_NO \n"); 
	    
	    fprintf(ATWOOD6,"  TIME    ");
	    fprintf(ATWOOD6,"    Z        (adj. Z)   ");
	    fprintf(ATWOOD6,"  EFF_AT_NO \n");
	}
	
	    
	if (dim == 3)
	    iz = h_extremum_index - 1;
	else 
	    iy = h_extremum_index - 1;
	while (coords[vert_dir] <= mix_zone_frac_outerbound)
	{
	    dens_minmax[0] = HUGE;
	    dens_minmax[1] = -HUGE;
	    
	    counter++;
		
	    if (dim == 3)
	    {
	        iz++;  
		icoords[vert_dir] = iz;
		for (iy = 0; iy < ymax; iy++)
		{
		    icoords[1] = iy;
		    for (ix = 0; ix < xmax; ix++)
		    {
		        icoords[0] = ix;
			coords = Rect_coords(icoords,wave);
			comp = Rect_comp(icoords,wave);
			state = Rect_state(icoords,wave);
			dens_minmax[0] = min(Dens(state), dens_minmax[0]);
			dens_minmax[1] = max(Dens(state), dens_minmax[1]);
		    }
		}
	    }
	    else /* if (dim == 2) */
	    {
	        iy++;  
		icoords[vert_dir] = iy;
		
		for (ix = 0; ix < xmax; ix++)
		{
		    icoords[0] = ix;
		    coords = Rect_coords(icoords,wave);
		    comp = Rect_comp(icoords,wave);
		    state = Rect_state(icoords,wave);
		    dens_minmax[0] = min(Dens(state), dens_minmax[0]);
		    dens_minmax[1] = max(Dens(state), dens_minmax[1]);
		}
	    }
	    
  
	    accumulate_state_in_layer(wave,front,coords[vert_dir],bst,
				      NULL,rfactor,origin,NO,YES,YES,YES);

	    pp_gsync();
      
	    pp_global_min(&dens_minmax[0], 1);
	    pp_global_max(&dens_minmax[1], 1);




	    atwood = (dens_minmax[1] - dens_minmax[0]) /
	             (dens_minmax[1] + dens_minmax[0]);

	    /* The assumption below is that the heavy fluid index is 
	       always 0 and the light fluid index is always 1. So this 
	       assumes only two layers in the R-T problem. */
	    atwood2 = (bst[0].d - bst[1].d) / (bst[0].d + bst[1].d) ;
	    atwood3 = (bst[0].global_max_d - bst[1].global_min_d) / 
	              (bst[0].global_max_d + bst[1].global_min_d) ;      
	    atwood5 = (bst[0].partial_d - bst[1].partial_d) / 
	              (bst[0].partial_d + bst[1].partial_d);
	    atwood6 = (bst[0].median_d - bst[1].median_d) / 
	              (bst[0].median_d + bst[1].median_d);

	    if (coords[vert_dir] <= mix_zone_frac_outerbound)
	    {
	        atwood_sum += atwood;
		if (bst[1].d == 0)
		    atwood_counteradjuster++;
		else
		{
		    atwood_sum2 += atwood2;
		    atwood_sum3 += atwood3;
		    atwood_sum5 += atwood5;
		    atwood_sum6 += atwood6;
		}
	    }
  
  
	    if (pp_mynode() == 0)
	    {
	        fprintf(ATWOOD,"%-10g  ", front->time);
		fprintf(ATWOOD,"%-10g %-10g   ", coords[vert_dir], 
			2.0*mean_initial_interface_height - 
			coords[vert_dir]);
		fprintf(ATWOOD,"%-10g  \n", atwood);
		
		fprintf(ATWOOD,"\tdens min=%g, dens max=%g,",
			dens_minmax[0],dens_minmax[1]);
		fprintf(ATWOOD," atwood_sum=%g, count=%d\n\n", 
			atwood_sum,counter);    
		
		

		fprintf(ATWOOD2,"%-10g  ", front->time);
		fprintf(ATWOOD2,"%-10g %-10g   ", coords[vert_dir], 
			2.0*mean_initial_interface_height - 
			coords[vert_dir]);
		fprintf(ATWOOD2,"%-10g  \n", atwood2);
		fprintf(ATWOOD2,"\tdens l=%g, dens h=%g,",
			bst[1].d,bst[0].d);
		fprintf(ATWOOD2," atwood_sum2=%g, count=%d\n\n", 
			atwood_sum2,counter);
		


		fprintf(ATWOOD3,"%-10g  ", front->time);
		fprintf(ATWOOD3,"%-10g %-10g   ", coords[vert_dir], 
			2.0*mean_initial_interface_height - 
			coords[vert_dir]);
		fprintf(ATWOOD3,"%-10g  \n", atwood3);
		fprintf(ATWOOD3,"\tdens l=%g, dens h=%g,", 
			bst[1].global_min_d,bst[0].global_max_d);
		fprintf(ATWOOD3," atwood_sum3=%g, count=%d\n\n", 
			atwood_sum3,counter);
		
		
		
		fprintf(ATWOOD5,"%-10g  ", front->time);
		fprintf(ATWOOD5,"%-10g %-10g   ", coords[vert_dir], 
			2.0*mean_initial_interface_height - 
			coords[vert_dir]);
		fprintf(ATWOOD5,"%-10g  \n", atwood5);
		fprintf(ATWOOD5,"\tdens l=%g, dens h=%g,",
			bst[1].partial_d,bst[0].partial_d);
		fprintf(ATWOOD5," atwood_sum5=%g, count=%d\n\n", 
			atwood_sum5,counter);
		
		
		
		fprintf(ATWOOD6,"%-10g  ", front->time);
		fprintf(ATWOOD6,"%-10g %-10g   ", coords[vert_dir], 
			2.0*mean_initial_interface_height - 
			coords[vert_dir]);
		fprintf(ATWOOD6,"%-10g  \n", atwood6);
		fprintf(ATWOOD6,"\tdens l=%g, dens h=%g,",
			bst[1].median_d,bst[0].median_d);
		fprintf(ATWOOD6," atwood_sum6=%g, count=%d\n\n", 
			atwood_sum6,counter);
	    }
	    
	} /* end of   'while (coords[vert_dir] <= mix_zone_frac_outerbound)' loop */
	


	if (pp_mynode() == 0)
	{  
	    counter--; /* Done since counter is incremented once 
			  more than needed in the above while loop.
		       */
	    atwood = atwood_sum/counter;
	    atwood2 = atwood_sum2 / (counter - atwood_counteradjuster);
	    atwood3 = atwood_sum3 / (counter - atwood_counteradjuster);
	    atwood5 = atwood_sum5 / (counter - atwood_counteradjuster);
	    atwood6 = atwood_sum6 / (counter - atwood_counteradjuster);
	    
	    fprintf(ATWOOD,"\nAVERAGE EFFECTIVE ATWOOD NUMBER = %g,", atwood);
	    fprintf(ATWOOD," time = %g, bub. height = %g\n", 
		    front->time, adjusted_bub_or_spike_height);
	    fprintf(ATWOOD,"\nMix zone lower(actual) bd used  %lg ",
		    (2.0*mean_initial_interface_height-mix_zone_frac_outerbound));
	    fprintf(ATWOOD,"(%lg), real bub h %lg, counter  %d\n\n", 
		    mix_zone_frac_outerbound,bub_height,counter);
		
		
	    fprintf(ATWOOD2,"\nAVERAGE EFFECTIVE ATWOOD NUMBER = %g,", 
		    atwood2);
	    fprintf(ATWOOD2," time = %g, bub. height = %g\n", 
		    front->time, adjusted_bub_or_spike_height);
	    fprintf(ATWOOD2,"\nMix zone lower(actual) bd used  %lg ",
		    (2.0*mean_initial_interface_height-mix_zone_frac_outerbound));
	    fprintf(ATWOOD2,"(%lg), real bub h %lg, counter  %d ",
		    mix_zone_frac_outerbound,bub_height,counter);
	    fprintf(ATWOOD2,"counter_adj %d\n\n", atwood_counteradjuster); 
	    
	    
	    fprintf(ATWOOD3,"\nAVERAGE EFFECTIVE ATWOOD NUMBER = %g,",
		    atwood3);
	    fprintf(ATWOOD3," time = %g, bub. height = %g\n", 
		    front->time, adjusted_bub_or_spike_height);
	    fprintf(ATWOOD3,"\nMix zone lower(actual) bd used  %lg",
		    (2.0*mean_initial_interface_height-mix_zone_frac_outerbound));
	    fprintf(ATWOOD3," (%lg), real bub h %lg, counter  %d ",
		    mix_zone_frac_outerbound,bub_height,counter);
	    fprintf(ATWOOD3,"counter_adj %d\n\n", atwood_counteradjuster); 
		
		
	    fprintf(ATWOOD5,"\nAVERAGE EFFECTIVE ATWOOD NUMBER = %g,",
		    atwood5);
	    fprintf(ATWOOD5," time = %g, bub. height = %g\n", 
		    front->time, adjusted_bub_or_spike_height);
	    fprintf(ATWOOD5,"\nMix zone lower(actual) bd used  %lg ",
		    (2.0*mean_initial_interface_height-mix_zone_frac_outerbound));
	    fprintf(ATWOOD5,"(%lg), real bub h %lg, counter  %d ",
		    mix_zone_frac_outerbound,bub_height,counter);
	    fprintf(ATWOOD5,"counter_adj %d\n\n", atwood_counteradjuster);
	    

	    fprintf(ATWOOD6,"\nAVERAGE EFFECTIVE ATWOOD NUMBER = %g,",
		    atwood6);
	    fprintf(ATWOOD6," time = %g, bub. height = %g\n", 
		    front->time, adjusted_bub_or_spike_height);
	    fprintf(ATWOOD6,"\nMix zone lower(actual) bd used  %lg ",
		    (2.0*mean_initial_interface_height-mix_zone_frac_outerbound));
	    fprintf(ATWOOD6,"(%lg), real bub h %lg, counter  %d ",
		    mix_zone_frac_outerbound,bub_height,counter);
	    fprintf(ATWOOD6,"counter_adj %d\n\n", atwood_counteradjuster);
	    
		
	    fclose(ATWOOD);
	    fclose(ATWOOD2);
	    fclose(ATWOOD3);
	    fclose(ATWOOD5);
	    fclose(ATWOOD6);
	} /* end of  'if (pp_mynode() == 0)' */
	    






	/* _________________________________________________________ */
  

	/* If requested, collect and print all of the effective Atwood 
	   number data for the spike region*/
	
	if (*do_spikes == YES)
	{
	    /* initializations */
	    if (front->time == 0 || front->step == 0)
	        first = YES;

	    atwood_counteradjuster = 0;
	    counter = 0;
	    dens_minmax[0] = HUGE;      
	    dens_minmax[1] = -HUGE;
	    atwood_sum = 0.0;
	    atwood_sum2 = 0.0;
	    atwood_sum3 = 0.0;
	    atwood_sum5 = 0.0; 
	    atwood_sum6 = 0.0;
	    

	    /* The fraction of the outer portion of the (spike) mixing 
	       zone used for averaging to determine the time dependent 
	       effective Atwood number is determined by the fifth argument 
	       to this function, mix_zone_frac. This mix_zone_frac is also  
	       interjected into the name of the specific Effective Atwood 
	       Number output file.
	       
	       NOTE that the light fluid is on top of the heavy fluid in 
	       FronTier Rayleigh-Taylor simulations, and this must be 
	       accounted for in computing the "bottom" mix_zone_frac of the 
	       mixing zone. 
	    */
	    
	    adjusted_bub_or_spike_height = 2.0*mean_initial_interface_height - 
	                                   spike_height;
	    mix_zone_frac_outerbound = mean_initial_interface_height - 
	                     (1.0 - mix_zone_frac) * 
	                     (mean_initial_interface_height -  spike_height);
	    
	    
	    /* Note that it is the centres of the cell blocks and not the 
	       edges which are being used in obtaining densities - 
	       egeorge 10/9/2001.
	    */
	    if (dim == 3)
	        h_extremum_index = (int) ((zmax * (spike_height-(L[vert_dir] + 
						    0.5*dh[vert_dir]))) /
				     (U[vert_dir]-L[vert_dir]));
	    else
	        h_extremum_index = (int) ((ymax * (spike_height-(L[vert_dir] + 
						      0.5*dh[vert_dir]))) /
				 (U[vert_dir]-L[vert_dir]));
	    

	    
	    icoords[0] = icoords[1] = 0;
	    icoords[vert_dir] = h_extremum_index;
	    coords = Rect_coords(icoords,wave);
	    while (coords[vert_dir] < spike_height)
	    {
	        h_extremum_index++;
		icoords[vert_dir] = h_extremum_index;
		coords = Rect_coords(icoords,wave);
	    }
	    while (coords[vert_dir] > spike_height)
	    {
	        h_extremum_index--;
		icoords[vert_dir] = h_extremum_index;
		coords = Rect_coords(icoords,wave);
	    }
	    
	    

	    if (pp_mynode() == 0)
	    {      
	        sprintf(filename,"%s-spikeEFF-ATWOOD-%05.2f_top%3.2f",
			b_fn,percent,mix_zone_frac);
		sprintf(filename2,"%s-spikeEFF-ATWOOD-2-%05.2f_top%3.2f",
			b_fn,percent,mix_zone_frac);
		sprintf(filename3,"%s-spikeEFF-ATWOOD-3-%05.2f_top%3.2f",
			b_fn,percent,mix_zone_frac);
		sprintf(filename5,"%s-spikeEFF-ATWOOD-5-%05.2f_top%3.2f",
			b_fn,percent,mix_zone_frac);
		sprintf(filename6,"%s-spikeEFF-ATWOOD-6-%05.2f_top%3.2f",
			b_fn,percent,mix_zone_frac);     
		if(first)
		{
		    ATWOOD = fopen(filename, "w");
		    ATWOOD2 = fopen(filename2, "w");
		    ATWOOD3 = fopen(filename3, "w");
		    ATWOOD5 = fopen(filename5, "w");
		    ATWOOD6 = fopen(filename6, "w");
		}
		else
		{
		    ATWOOD = fopen(filename, "a"); 
		    ATWOOD2 = fopen(filename2, "a");
		    ATWOOD3 = fopen(filename3, "a");
		    ATWOOD5 = fopen(filename5, "a");
		    ATWOOD6 = fopen(filename6, "a");
		}
	    
	    
		fprintf(ATWOOD,"  TIME    ");
		fprintf(ATWOOD,"    Z        (adj. Z)   ");
		fprintf(ATWOOD,"  EFF_AT_NO \n"); 
	    
		fprintf(ATWOOD2,"  TIME    ");
		fprintf(ATWOOD2,"    Z        (adj. Z)   ");
		fprintf(ATWOOD2,"  EFF_AT_NO \n"); 

		fprintf(ATWOOD3,"  TIME    ");
		fprintf(ATWOOD3,"    Z        (adj. Z)   ");
		fprintf(ATWOOD3,"  EFF_AT_NO \n"); 
	    
		fprintf(ATWOOD5,"  TIME    ");
		fprintf(ATWOOD5,"    Z        (adj. Z)   ");
		fprintf(ATWOOD5,"  EFF_AT_NO \n"); 
	    
		fprintf(ATWOOD6,"  TIME    ");
		fprintf(ATWOOD6,"    Z        (adj. Z)   ");
		fprintf(ATWOOD6,"  EFF_AT_NO \n");
	    }
	
	    
	    if (dim == 3)
	        iz = h_extremum_index + 1;
	    else 
	        iy = h_extremum_index + 1;
	    while (coords[vert_dir] >= mix_zone_frac_outerbound)
	    {
	        dens_minmax[0] = HUGE;
		dens_minmax[1] = -HUGE;
	    
		counter++;
	    
		if (dim == 3)
		{
		    iz--;  
		    icoords[vert_dir] = iz;
		    for (iy = 0; iy < ymax; iy++)
		    {
		        icoords[1] = iy;
			for (ix = 0; ix < xmax; ix++)
			{
			    icoords[0] = ix;
			    coords = Rect_coords(icoords,wave);
			    comp = Rect_comp(icoords,wave);
			    state = Rect_state(icoords,wave);
			    dens_minmax[0] = min(Dens(state), dens_minmax[0]);
			    dens_minmax[1] = max(Dens(state), dens_minmax[1]);
			}
		    }
		}
		else /* if (dim == 2) */
		{
		    iy--;  
		    icoords[vert_dir] = iy;
		
		    for (ix = 0; ix < xmax; ix++)
		    {
		        icoords[0] = ix;
			coords = Rect_coords(icoords,wave);
			comp = Rect_comp(icoords,wave);
			state = Rect_state(icoords,wave);
			dens_minmax[0] = min(Dens(state), dens_minmax[0]);
			dens_minmax[1] = max(Dens(state), dens_minmax[1]);
		    }
		}
		
  
		accumulate_state_in_layer(wave,front,coords[vert_dir],bst,
					  NULL,rfactor,origin,NO,YES,YES,YES);

		pp_gsync();
      
		pp_global_min(&dens_minmax[0], 1);
		pp_global_max(&dens_minmax[1], 1);




		atwood = (dens_minmax[1] - dens_minmax[0]) /
		         (dens_minmax[1] + dens_minmax[0]);

		/* The assumption below is that the heavy fluid index is 
		   always 0 and the light fluid index is always 1. So this 
		   assumes only two layers in the R-T problem. */
		atwood2 = (bst[0].d - bst[1].d) / (bst[0].d + bst[1].d) ;
		atwood3 = (bst[0].global_max_d - bst[1].global_min_d) / 
		          (bst[0].global_max_d + bst[1].global_min_d) ;      
		atwood5 = (bst[0].partial_d - bst[1].partial_d) / 
	                  (bst[0].partial_d + bst[1].partial_d);
		atwood6 = (bst[0].median_d - bst[1].median_d) / 
	                  (bst[0].median_d + bst[1].median_d);

		if (coords[vert_dir] >= mix_zone_frac_outerbound)
		{
		    atwood_sum += atwood;
		    if (bst[1].d == 0)
		        atwood_counteradjuster++;
		    else
		    {
		        atwood_sum2 += atwood2;
			atwood_sum3 += atwood3;
			atwood_sum5 += atwood5;
			atwood_sum6 += atwood6;
		    }
		}
  
  
		if (pp_mynode() == 0)
		{
		    fprintf(ATWOOD,"%-10g  ", front->time);
		    fprintf(ATWOOD,"%-10g %-10g   ", coords[vert_dir], 
			    2.0*mean_initial_interface_height - 
			    coords[vert_dir]);
		    fprintf(ATWOOD,"%-10g  \n", atwood);
		
		    fprintf(ATWOOD,"\tdens min=%g, dens max=%g,",
			    dens_minmax[0],dens_minmax[1]);
		    fprintf(ATWOOD," atwood_sum=%g, count=%d\n\n", 
			    atwood_sum,counter);    
		
		

		    fprintf(ATWOOD2,"%-10g  ", front->time);
		    fprintf(ATWOOD2,"%-10g %-10g   ", coords[vert_dir], 
			    2.0*mean_initial_interface_height - 
			    coords[vert_dir]);
		    fprintf(ATWOOD2,"%-10g  \n", atwood2);
		    fprintf(ATWOOD2,"\tdens l=%g, dens h=%g,",
			    bst[1].d,bst[0].d);
		    fprintf(ATWOOD2," atwood_sum2=%g, count=%d\n\n", 
			    atwood_sum2,counter);
		    


		    fprintf(ATWOOD3,"%-10g  ", front->time);
		    fprintf(ATWOOD3,"%-10g %-10g   ", coords[vert_dir], 
			    2.0*mean_initial_interface_height - 
			    coords[vert_dir]);
		    fprintf(ATWOOD3,"%-10g  \n", atwood3);
		    fprintf(ATWOOD3,"\tdens l=%g, dens h=%g,", 
			    bst[1].global_min_d,bst[0].global_max_d);
		    fprintf(ATWOOD3," atwood_sum3=%g, count=%d\n\n", 
			    atwood_sum3,counter);
		    
		
		
		    fprintf(ATWOOD5,"%-10g  ", front->time);
		    fprintf(ATWOOD5,"%-10g %-10g   ", coords[vert_dir], 
			    2.0*mean_initial_interface_height - 
			    coords[vert_dir]);
		    fprintf(ATWOOD5,"%-10g  \n", atwood5);
		    fprintf(ATWOOD5,"\tdens l=%g, dens h=%g,",
			    bst[1].partial_d,bst[0].partial_d);
		    fprintf(ATWOOD5," atwood_sum5=%g, count=%d\n\n", 
			    atwood_sum5,counter);
		
		
		
		    fprintf(ATWOOD6,"%-10g  ", front->time);
		    fprintf(ATWOOD6,"%-10g %-10g   ", coords[vert_dir], 
			    2.0*mean_initial_interface_height - 
			    coords[vert_dir]);
		    fprintf(ATWOOD6,"%-10g  \n", atwood6);
		    fprintf(ATWOOD6,"\tdens l=%g, dens h=%g,",
			    bst[1].median_d,bst[0].median_d);
		    fprintf(ATWOOD6," atwood_sum6=%g, count=%d\n\n", 
			    atwood_sum6,counter);
		}
	    
	    } /*end of 'while (coords[vert_dir] <= mix_zone_frac_outerbound)'*/
	


	    if (pp_mynode() == 0)
	    {  
	        counter--; /* Done since counter is incremented once 
			      more than needed in the above while loop.
			   */
		atwood = atwood_sum/counter;
		atwood2 = atwood_sum2 / (counter - atwood_counteradjuster);
		atwood3 = atwood_sum3 / (counter - atwood_counteradjuster);
		atwood5 = atwood_sum5 / (counter - atwood_counteradjuster);
		atwood6 = atwood_sum6 / (counter - atwood_counteradjuster);
		
		fprintf(ATWOOD,"\nAVERAGE EFFECTIVE ATWOOD NUMBER = %g,", 
			atwood);
		fprintf(ATWOOD," time = %g, spike height = %g\n", 
			front->time, adjusted_bub_or_spike_height);
		fprintf(ATWOOD,"\nMix zone upper(actual) bd used  %lg ",
		 (2.0*mean_initial_interface_height-mix_zone_frac_outerbound));
		fprintf(ATWOOD,"(%lg), real spike h %lg, counter  %d\n\n", 
			mix_zone_frac_outerbound,spike_height,counter);
		
		
		fprintf(ATWOOD2,"\nAVERAGE EFFECTIVE ATWOOD NUMBER = %g,", 
			atwood2);
		fprintf(ATWOOD2," time = %g, spike height = %g\n", 
			front->time, adjusted_bub_or_spike_height);
		fprintf(ATWOOD2,"\nMix zone upper(actual) bd used  %lg ",
		 (2.0*mean_initial_interface_height-mix_zone_frac_outerbound));
		fprintf(ATWOOD2,"(%lg), real spike h %lg, counter  %d ",
			mix_zone_frac_outerbound,spike_height,counter);
		fprintf(ATWOOD2,"counter_adj %d\n\n", atwood_counteradjuster); 
	    
	    
		fprintf(ATWOOD3,"\nAVERAGE EFFECTIVE ATWOOD NUMBER = %g,",
			atwood3);
		fprintf(ATWOOD3," time = %g, spike height = %g\n", 
			front->time, adjusted_bub_or_spike_height);
		fprintf(ATWOOD3,"\nMix zone upper(actual) bd used  %lg",
		 (2.0*mean_initial_interface_height-mix_zone_frac_outerbound));
		fprintf(ATWOOD3," (%lg), real spike h %lg, counter  %d ",
			mix_zone_frac_outerbound,spike_height,counter);
		fprintf(ATWOOD3,"counter_adj %d\n\n", atwood_counteradjuster); 
		
		
		fprintf(ATWOOD5,"\nAVERAGE EFFECTIVE ATWOOD NUMBER = %g,",
			atwood5);
		fprintf(ATWOOD5," time = %g, spike height = %g\n", 
			front->time, adjusted_bub_or_spike_height);
		fprintf(ATWOOD5,"\nMix zone upper(actual) bd used  %lg ",
		 (2.0*mean_initial_interface_height-mix_zone_frac_outerbound));
		fprintf(ATWOOD5,"(%lg), real spike h %lg, counter  %d ",
			mix_zone_frac_outerbound,spike_height,counter);
		fprintf(ATWOOD5,"counter_adj %d\n\n", atwood_counteradjuster);
	    

		fprintf(ATWOOD6,"\nAVERAGE EFFECTIVE ATWOOD NUMBER = %g,",
			atwood6);
		fprintf(ATWOOD6," time = %g, spike height = %g\n", 
			front->time, adjusted_bub_or_spike_height);
		fprintf(ATWOOD6,"\nMix zone upper(actual) bd used  %lg ",
		 (2.0*mean_initial_interface_height-mix_zone_frac_outerbound));
		fprintf(ATWOOD6,"(%lg), real spike h %lg, counter  %d ",
			mix_zone_frac_outerbound,spike_height,counter);
		fprintf(ATWOOD6,"counter_adj %d\n\n", atwood_counteradjuster);
	    
		
		fclose(ATWOOD);
		fclose(ATWOOD2);
		fclose(ATWOOD3);
		fclose(ATWOOD5);
		fclose(ATWOOD6);
	    } /* end of  'if (pp_mynode() == 0)' */
	} /* end of 'if (do_spikes == YES)' */    


	return;	  
}		/*end print_effective_atwood_num*/



/*                       find_intersection_by_bisection():
*
*   The following attempts to approximate the intersection of the 
*   interface and the line joining coords_a and coords_b by bisection 
*   in direction dir, performing a maximum of 10 bisections. This function 
*   is only to be called if the component of coords_a and coords_b differ.
*/
LOCAL   double    find_intersection_by_bisection(
	Front			*front,
	double                   *coords_a,
	double                   *coords_b,
	int                     dir)
{
        double                 ans;
	double                 crds_a[MAXD], crds_b[MAXD], cross_crds[MAXD];
	INTERFACE	      *intfc = front->interf;
	COMPONENT	      comp, comp_a, comp_b;
	int                   i;
	int                   max_num_bisections = 10;



	for (i = 0; i < MAXD; i++)
	{
	    crds_a[i] = cross_crds[i] = coords_a[i];
	    crds_b[i] = coords_b[i];
	}

	/* Now arrange it so that the crds_a[dir] < crds_b[dir] */

	if (coords_a[dir] > coords_b[dir])
	{
	    for (i = 0; i < MAXD; i++)
	    {
		crds_a[i] = crds_b[i];
		crds_b[i] = cross_crds[i];
	    }
	}
	  

	/* First check if coords_a or coords_b happen to be ONFRONT */
	if (component(crds_a, intfc) == ONFRONT)
	{
	    ans = crds_a[dir];
	}
	else if (component(crds_b, intfc) == ONFRONT)
	{
	    ans = crds_b[dir];
	}
	else
	{
	    for (i = 0; i < max_num_bisections; i++)
	    {      	
		cross_crds[dir] = 0.5*(crds_a[dir] + crds_b[dir]);
		comp = component(cross_crds,intfc);

		if (comp == ONFRONT)
		    break;		
		comp_a = component(crds_a,intfc);
		comp_b = component(crds_b,intfc);
	
		if (comp == comp_a)
		    crds_a[dir] = cross_crds[dir];
		else
		    crds_b[dir] = cross_crds[dir];
	    }
	    ans = cross_crds[dir];
	}
   
	return ans;
}		/*end find_intersection_by_bisection*/ 
        



/*                find_particular_fluid_cell_volume():
*
*   The following takes the lower left hand corner of a 1-d or 2-d (for now)
*   rectangular grid cell and computes the volume (i.e., length for 1d cell 
*   and area for 2d cell) of fluid with component label vol_comp in that cell.
*/

LOCAL   double    find_particular_fluid_cell_volume(Front *front, 
						   COMPONENT vol_comp,
						   double *crds,
						   double *dh)
{
        INTERFACE		*intfc = front->interf;
        double                   ans = 0.0;
	double                   previous_coords[MAXD], coords[MAXD];
	double                   *intersection_coord;
	/* The following, intersection_length, is by default computed
	   as the length of the line segment from the LOWER x or y bound 
	   (respectively) of a grid cell line to where the interface cuts 
	   that x or y (resp.) grid cell line. The values are 
	   adjusted appropriately in case it is the length on 
	   the other side of the intersection which is needed.
	*/
	double                   *intersection_length;
	COMPONENT		*comp;
	int                     vol_comp_count, num_points;
	int                     i;


	uni_array(&intersection_coord, (int)pow(2,MAXD-1), sizeof(double));
	uni_array(&intersection_length, (int)pow(2,MAXD-1), sizeof(double));
	uni_array(&comp, (int)pow(2,MAXD-1), sizeof(int));

	for (i = 0; i < MAXD; i++)
	    previous_coords[i] = coords[i] = crds[i];
	/* Note in 2d, we really only need intersection_coord and 
	   intersection_length to be scalars.
	*/
	for (i = 0; i < (int)pow(2,MAXD-1); i++)
	    intersection_coord[i] = intersection_length[i] = -HUGE;


	/* If the dimension of the interface is not 2 or 3, return 
	   the initial value of ans  (0.0).
	*/
	if (front->rect_grid->dim == 2)
	{
	    coords[0] += dh[0]; 
	    comp[0] = component(previous_coords, intfc);
	    comp[1] = component(coords, intfc);
		  
	    if (comp[0] != comp[1])
	    {
		intersection_coord[0] = find_intersection_by_bisection(
					     front,previous_coords,coords,0);
		intersection_length[0] = intersection_coord[0] - 
		                         previous_coords[0];
	    }

	    if (comp[0] == vol_comp && comp[1] == vol_comp)
	        ans = dh[0];
	    else if (comp[0] == vol_comp && comp[1] != vol_comp)
	        ans = intersection_length[0];
	    else if (comp[0] != vol_comp && comp[1] == vol_comp)
	        ans = dh[0] - intersection_length[0];
	    /* else, return 0.0 */
	}



	else if (front->rect_grid->dim == 3)
	{
	    comp[0] = component(coords,intfc);


	    coords[0] += dh[0];
	    comp[1] = component(coords,intfc);
	    if (comp[1] != comp[0])
	    {
		intersection_coord[0] = find_intersection_by_bisection(
					 front,previous_coords, coords, 0);
		intersection_length[0] = intersection_coord[0] - 
		                         previous_coords[0];
	    }
	    previous_coords[0] = coords[0];


	    coords[1] += dh[1]; 
	    comp[2] = component(coords,intfc);
	    if (comp[2] != comp[1])
	    {
		intersection_coord[1] = find_intersection_by_bisection(
					  front,previous_coords, coords, 1);
		intersection_length[1] = intersection_coord[1] - 
		                         previous_coords[1];
	    }
	    previous_coords[1] = coords[1];
		    

	    coords[0] -= dh[0];
	    comp[3] = component(coords,intfc);
	    if (comp[3] != comp[2])
	    {
		intersection_coord[2] = find_intersection_by_bisection(
					  front,previous_coords, coords, 0);
		intersection_length[2] = intersection_coord[2] - coords[0];
	    }
	    previous_coords[0] = coords[0];


	    coords[1] -= dh[1];
	    if (comp[0] != comp[3])
	    {
		intersection_coord[3] = find_intersection_by_bisection(
					  front,previous_coords, coords, 1);
		intersection_length[3] = intersection_coord[3] - coords[1];
	    }




	    vol_comp_count = 0;
	    num_points = 4;
	    for (i = 0; i < num_points; i++)
	    {
		if (comp[i] == vol_comp)
		    vol_comp_count++;
	    }
	    
	    /* There are sixteen basic cases to consider, which 
	       I have amalgamated into 3 main cases.  
	    */
	    
	    /* CASE 1: No cell edges are crossed by the interface. 
	       This accounts for 2 of the basic cases.
	    */
	    if (vol_comp_count == 0 || vol_comp_count == 4)
	    {
		if (comp[0] == vol_comp)
		    ans = dh[0]*dh[1];
		else
		    ans = 0.0;
	    }
	    
	    
	    /* CASE 2: Two adjacent edges are crossed by the interface
	       and no other edges are crossed. This contains four 
	       sub-cases. In total, this main case takes care of 
	       all eight of the basic cases in which only one of the 
	       cell nodes differs in component value from the other three.
	    */
	    else if (vol_comp_count == 1 || vol_comp_count == 3)
	    {
		if ((comp[0] != comp[1]) && (comp[1] != comp[2]))
		{
		    ans = 0.5 * (dh[0] - intersection_length[0]) * 
		      (intersection_length[1]);
		    if (comp[1] != vol_comp)
		    {
			ans *= -1.0;
			ans += dh[0]*dh[1];
		    }
		}
		
		
		else if ((comp[1] != comp[2]) && (comp[2] != comp[3]))
		{
		    ans = 0.5 * (dh[1]-intersection_length[1]) * 
		      (dh[0]-intersection_length[2]);
		    if (comp[2] != vol_comp)
		    {
			ans *= -1.0;
			ans += dh[0]*dh[1];
		    }
		}
		
		
		else if ((comp[2] != comp[3]) && (comp[3] != comp[0]))
		{
		    ans = 0.5 * (intersection_length[2]) * 
		      (dh[1]-intersection_length[3]);
		    if (comp[3] != vol_comp)
		    {
			ans *= -1.0;
			ans += dh[0]*dh[1];
		    }
		}
		
		else
		{
		    ans = 0.5 * (intersection_length[3]) * 
		      (intersection_length[0]);
		    if (comp[0] != vol_comp)
		    {
			ans *= -1.0;
			ans += dh[0]*dh[1];
		    }
		}
	    }
	    


	    /* CASE 3: Two non-adjacent edges are crossed by the interface
	       and no other edges are crossed or all four edges are crossed 
	       by the interface. This contains the remaining six basic cases.
	    */
	    else if (vol_comp_count == 2)
	    {
		/* First the four cases in which two adjacent cell nodes 
		   have the same component value
		*/
		if (comp[0] == comp[1])
		{
		    ans = dh[0] * 0.5 * (intersection_length[1] + 
					 intersection_length[3]); 
		    if (comp[0] != vol_comp)
		    {
			ans *= -1.0;
			ans += dh[0]*dh[1];
		    }
		}
		
		else if (comp[1] == comp[2])
		{
		    ans = dh[1] * 0.5* (intersection_length[0] + 
					intersection_length[2]); 
		    if (comp[1] == vol_comp)
		    {
			ans *= -1.0;
			ans += dh[0]*dh[1];
		    }  
		}
		
		/* Finally, there are the two tricky cases in which 
		   the interface cuts each of the four cell edges. 
		   There are two possible solutions to each case.
		   I break ties here by determining the component 
		   in the cell centre.
		*/
		else 
		{
		    coords[0] = crds[0] + 0.5*dh[0];
		    coords[1] = crds[1] + 0.5*dh[1];
		    
		    if (comp[0] == vol_comp)
		    {
			if (component(coords,intfc) != vol_comp)
			{
			    ans = 0.5 * (intersection_length[0]) * 
			          intersection_length[3];
			    ans += 0.5 * (dh[0] - intersection_length[2]) *
			          (dh[1]-intersection_length[1]);
			}
			else
			{
			    ans = 0.5 * (dh[0] - intersection_length[0]) *
			          intersection_length[1];
			    ans += 0.5 * intersection_length[2] * 
			           (dh[1] - intersection_length[3]);
			    ans *= -1.0;
			    ans += dh[0]*dh[1];
			}
		    }
		    else /* if comp[1] == vol_comp */
		    {
			if (component(coords,intfc) != vol_comp)
			{
			    ans = 0.5 * (dh[0]-intersection_length[0]) * 
			          intersection_length[1];
			    ans += 0.5 * (intersection_length[2]) * 
			           (dh[1]-intersection_length[3]);
			}
			else
			{
			    ans = 0.5 * (dh[1]-intersection_length[1]) * 
			          (dh[0]-intersection_length[2]);
			    ans += 0.5 * (intersection_length[0] * 
					  intersection_length[3]);
			    ans *= -1.0;
			    ans += dh[0]*dh[1];
			}
		    }
		}
	    }
	  
	    
	    /* Finally, if none of the four cell nodes has the same 
	       component value as vol_comp. This may seem to be a sub-case 
	       of CASE 1 above, but it in fact also takes care of the 
	       situation in which at least one of the cell nodes has 
	       NULL component value - which seems to happen within about 
	       3 cell blocks of the top and bottom of the global domain 
	       at time = 0. 4 Apr. 2002 - egeorge.
	    */
	    if ((comp[0] != vol_comp) && (comp[1] != vol_comp))
	    {
		if ((comp[2] != vol_comp) && (comp[3] != vol_comp))
		{
		    ans = 0.0;
		}
	    }
	} /* end of else if (front->rect_grid->dim == 3) */


	free(intersection_coord);
	free(intersection_length);
	free(comp);
  
	return ans;
}		/*end find_particular_fluid_cell_volume*/




LOCAL	void 	new_accumulate_fractions_in_layer(
	Front			*front,
	double			height,
	double			*frac,
	double			rfactor,
	const double		*origin,
	COMPONENT               comp_for_vol_frac)
{
	const RECT_GRID		*rgrid = front->rect_grid;
	const int		dim = rgrid->dim;
	const int		zdir = dim - 1;
	const enum intfc_geometry  geom = (origin == NULL ? planar : spherical);
	int			Nx;
	int 		        i, j;

	double                   dh[MAXD - 1];
	double                   coords[MAXD];
	const double             *L = rgrid->L;
	const double             *U = rgrid->U;
	const double             *GL = rgrid->GL;
	const double             *GU = rgrid->GU;
	double                   light_fluid_volume = 0.0;
	double                   ans;

	debug_print("glayer","Entered new_accumulate_fractions_in_layer(), h = %g\n",
		       height);

	for (i = 0; i < MAXD; i++)
	    coords[i] = 0.0;


	switch(geom)
	{
	case planar:
	    coords[zdir] = height;
	    if (height < L[zdir] || height >= U[zdir])
	        break;
	    Nx = irint(rgrid->gmax[0]*rfactor);
	    dh[0] = (U[0] - L[0])/Nx;

	    switch(dim)
	    {
#if defined(ONED)
	    case 1:
	        screen("ERROR in new_accumulate_fractions_in_layer(), "
		       "1D not supported\n");
	        clean_up(ERROR);
	        break;
#endif /* defined(ONED) */
#if defined(TWOD)
	    case 2:
	        for (i = 0; i < Nx; i++)
		{
		    coords[0] = L[0] + i*dh[0];
		  
		    light_fluid_volume += find_particular_fluid_cell_volume(
					front, comp_for_vol_frac, coords, dh);
		}

		if (debugging("glayer"))
		{
		    (void) printf("Accumulation of layer fractions "
				  "at h = %g completed.\n",height);
		    (void) printf("Local volume for component %d fluid = %lg\n",
				  comp_for_vol_frac, light_fluid_volume);
		}


		if (pp_numnodes() > 1)
		{
		    pp_global_sum(&light_fluid_volume, 1);
		}
		    /* The following assumes that horizontal slices take up 
		       the entire horizontal dimensions of the global 
		       domain. This could, I guess, not be the case, 
		       though I cannot think of a reason why such a volume 
		       fraction would be desired - 5 Apr. 2002 egeorge.
		    */
		ans = light_fluid_volume/(GU[0] - GL[0]);
		*frac = ans;





		break;
#endif /* defined(TWOD) */
#if defined(THREED)
	    case 3:
		{
		    double Ny = rgrid->gmax[1]*rfactor;    
		    dh[1] = (U[1] - L[1])/Ny;
		    for (j = 0; j < Ny; j++)
		    {
		        coords[1] = L[1] + j*dh[1];
			for (i = 0; i < Nx; i++)
		        {			    
			    coords[0] = L[0] + i*dh[0];
			    
			    light_fluid_volume += 
			      find_particular_fluid_cell_volume(front, 
					     comp_for_vol_frac, coords, dh);
			}
		    }
		}

		if (debugging("glayer"))
		{
		    (void) printf("Accumulation of layer fractions "
				  "at h = %g completed.\n",height);
		    (void) printf("Local volume for component %d fluid = %lg\n",
				  comp_for_vol_frac, light_fluid_volume);
		}

		if (pp_numnodes() > 1)
		{
		    pp_global_sum(&light_fluid_volume, 1);
		    pp_gsync();
		}
		    /* The following assumes that horizontal slices take up 
		       the entire horizontal dimensions of the global 
		       domain. This could, I guess, not be the case, 
		       though I cannot think of a reason why such a volume 
		       fraction would be desired - 5 Apr. 2002 egeorge.
		    */

		ans = light_fluid_volume/
		  ((GU[0] - GL[0]) * (GU[1] - GL[1]));
		*frac = ans;
		
		break;
#endif /* defined(THREED) */
	    }
	    




	    break;


	case spherical:
	    switch(dim)
	    {
#if defined(ONED)
	    case 1:
	      screen("ERROR in new_accumulate_fractions_in_layer(),\n"); 
	      screen("spherical 1D not supported. In gas/gintext.c, change\n");
	      screen("new_accumulate_fraction_in_layer() function calls\n");
	      screen("to accumulate_fraction_in_layer() function calls\n");
	        clean_up(ERROR);
	        break;
#endif /* defined(ONED) */
#if defined(TWOD)
	    case 2:
	      screen("ERROR in new_accumulate_fractions_in_layer(), \n");
	      screen("spherical 2D not supported. In gas/gintext.c, change\n"); 
	      screen("new_accumulate_fraction_in_layer() function calls\n"); 
              screen("to accumulate_fraction_in_layer() function calls\n");
	        break;
#endif /* defined(TWOD) */
#if defined(THREED)
	    case 3:
	      screen("ERROR in new_accumulate_fractions_in_layer(), \n");
	      screen ("spherical 3D not supported. In gas/gintext.c, change\n");
	      screen("new_accumulate_fraction_in_layer() function calls\n");
              screen("to accumulate_fraction_in_layer() function calls\n");
		break;
#endif /* defined(THREED) */
	    }
	    break;
	}
	

	if (debugging("glayer"))
	  {
	    (void) printf("Global volume for component %d fluid = %lg\n",
			  comp_for_vol_frac, light_fluid_volume);
	    (void) printf("\n");
	  }



	debug_print("glayer","Leaving new_accumulate_fractions_in_layer()\n");
	return;
}		/*end new_accumulate_fractions_in_layer*/



/*
*			new_height_at_fraction():
*
*	Find the height (z coordinate or radius) closest to initial height
*	h0 in the given direction dir (+/- 1) where the given layer fraction
*	occurs.  This fraction corresponds to the material whose component 
*	value is vol_comp. The main difference between this function and 
*       height_at_fraction() is the use here of 
*       new_accumulate_fractions_in_layer() instead of 
*       accumulate_fractions_in_layer(), in order to get a more accurate 
*       volume fraction. 
*
*       As of 19 Apr. 2004, this function is implemented and used only for 
*       planar geometry and not spherical geometry. TODO - fix this.
*/

LOCAL double 	new_height_at_fraction(
        Front			*front,
        double			h0,
        double			dh,
        int			dir,
        double			fraction,
        double			*frac,
        COMPONENT		vol_comp,
        double			rfactor,
        const double		*origin)
{	
        boolean                  stop;
	double                   height = h0;
	static double            old_frac, new_frac;


	old_frac = *frac;
	new_accumulate_fractions_in_layer(front,height,&old_frac,
					  rfactor,origin,vol_comp);
	if (old_frac >= fraction)
	{
	    height += dir*dh*fraction/old_frac;
	    if (debugging("warning"))
	    	(void) printf("WARNING in new_height_at_fraction(), "
			  "frac >= %g occurred at initial height of %g.\n",
			  fraction,height);
	    return height;
	}


	stop = NO;
	do 
	{
	    new_accumulate_fractions_in_layer(front,height+dir*dh,&new_frac,
					      rfactor,origin,vol_comp);
	    if ( ((old_frac <= fraction) && (fraction <= new_frac)) || 
		 ((new_frac <= fraction) && (fraction <= old_frac)) )
	    {   /* we've bracketed the height, so interpolate */
	        height += dir*dh*(fraction-old_frac) /
			  (new_frac-old_frac);
		stop = YES;
	    }
	    else
	    {
		old_frac = new_frac;
		height += dir*dh;
	    }
	}
	while (stop == NO && 
	       height >= front->rect_grid->GL[front->rect_grid->dim - 1] && 
	       height <= front->rect_grid->GU[front->rect_grid->dim - 1] );

	new_accumulate_fractions_in_layer(front,height,&old_frac,
					  rfactor,origin,vol_comp);
	*frac = old_frac;
	    
	return height;
}		/*end new_height_at_fraction*/

/* 
*			print_vol_frac():
*
*       Computes and prints volume fraction data as a function of domain 
*       height at pre-determined time step intervals. Much of the details 
*       of the algorithms used are contained in Chapter 4 and Appendix A 
*       of Erwin George's doctoral dissertation. The key difference 
*       between this function and old_print_vol_frac() is that the latter 
*       uses a more discrete (and less accurate) algorithm for determining
*       the volume fractions.
*
*       This function is presently only defined for 2d and 3d planar
*       geometry.
*
*       TODO: Modify following function and/or print_vol_frac() to deal 
*       with spherical geometry - 19 Apr. 2004. 
*/

LOCAL	void 	print_vol_frac(
	Wave			*wave,
	Front			*front,
	double			h_min,
	double			h_max,
	COMPONENT		comp_for_vol_frac,
	double			rfactor,
	const double		*origin,
	char                    *b_fn)
{
        const RECT_GRID		*rgrid = front->rect_grid;
        FILE                    *VOLFRAC;
	const double             *GL = rgrid->GL;
	const double             *GU = rgrid->GU;
	double                   dh = (h_max - h_min)/100.0; 
     
	double                   height, vol_frac;
	int                     k;
	int                     my_id = pp_mynode();
	const int               vert_dir = rgrid->dim - 1;
	char                    filename[256];


    


	if (rgrid->dim != 2 && rgrid->dim != 3)
	{
	    (void) printf("\nWARNING - entire domain volume fractions ");
	    (void) printf("currently only printed for 2d and 3d runs\n");
	    return;
	}


	sprintf(filename,"%s-vol_frac-t%07.2f",b_fn,front->time);

	if (my_id == 0)
	    VOLFRAC = fopen(filename, "w");


	height = GL[vert_dir] + dh;
	while (height < GU[vert_dir])
	{	  
	    pp_gsync(); 
	    new_accumulate_fractions_in_layer(front,height,&vol_frac,rfactor,
					      origin,comp_for_vol_frac);
	    if (my_id == 0)
	        fprintf(VOLFRAC,"%lg  %lg\n", height, vol_frac);	
	    
	    height += dh;
	}

	if (my_id == 0)
	    fclose(VOLFRAC);

	return;
}		/*end print_vol_frac*/


/*
* 			cmp_for_qsort_dec():
*
*      Used with the C standard library function qsort() to sort the entries  
*      of an array into decreasing order. Originally (17 Nov. 2003) implemented 
*      to help compute atwood5 and atwood6 (see also 
*      print_effective_atwood_num() in this file, gintext.c).
*/
LOCAL  int cmp_for_qsort_dec( 
	const void *vp, 
	const void *vq)
{
        const double *p = (const double *)vp;
        const double *q = (const double *)vq;
	double diff = *p - *q;

	return ((diff >= 0.0) ? ((diff > 0.0) ? -1 : 0) : +1);
}		/* end cmp_for_qsort_dec */

/*
* 			cmp_for_qsort_inc():
*
*      Used with the C standard library function qsort() to sort the entries  
*      of an array into increasing order. Originally (17 Nov. 2003) implemented 
*      to help compute atwood5 and atwood6 (see also 
*      print_effective_atwood_num() in this file, gintext.c).
*/
LOCAL int cmp_for_qsort_inc( 
        const void *vp, 
	const void *vq)
{
        const double *p = (const double*)vp;
	const double *q = (const double*)vq;
	double diff = *p - *q;

	return ((diff >= 0.0) ? ((diff > 0.0) ? +1 : 0) : -1);
}		/* end cmp_for_qsort_inc */


/*
* 			get_mean_initial_interface_height():
*
*      If version == 0, this function obtains the initial mean height 
*      of the interface for the random surface problem from 
*      gas/ginit/gipert.c's init_random_surface() or, if the run is 
*      a restart, from an added prompt in gas/gprt/gintext.c's 
*      init_print_effective_atwood_num().
*
*      If version != 0, this function passes that mean height to  
*      gas/gprt/gintext.c's print_effective_atwood_num(), where is it 
*      used in computing time dependent Atwood numbers.
*/
EXPORT double get_mean_initial_interface_height( 
					       int version,
					       double *mean_height)
{
        static double ans;

        if (version == 0)
	    ans = *mean_height;

	return ans;
}		/* end get_mean_initial_interface_height */

/* 
 * Following function will compute the interface area 
 **/
LOCAL double compute_intfc_area(INTERFACE *intfc)
{
        SURFACE **surfs;
        RECT_GRID *gr=computational_grid(intfc);
       
        TRI             *t;
    
        double *L = gr->L;      /* Lower corner of global grid */
        double *U = gr->U;      /* Upper corner of global grid */
        int   i,j,count;
        double p_coord,surf_area = 0.0 ;
        static double centroid[3];

        for (surfs = intfc->surfaces; surfs && *surfs; ++surfs)
        {
            if (wave_type(*surfs) < FIRST_PHYSICS_WAVE_TYPE)
	        continue;
	    for (t = first_tri(*surfs); !at_end_of_tri_list(t,*surfs); t = t->next)
	    {
	        centroid[0] = centroid[1] = centroid[2] = 0.0;
	        for (j = 0; j < 3; ++j)
	        {   
	            for (i = 0; i < 3; ++i)
		    {   
		        p_coord = Coords(Point_of_tri(t)[j])[i];
                        centroid[i] += p_coord;
                    }
	       }
	       count = 0;
	       for (i = 0; i < 3; ++i)
	       {
	           centroid[i] /= 3.0;
		   if (centroid[i] < L[i] || centroid[i] > U[i])
		       count++;
	       }
	       if (count == 0)
	       {
	           /*The centroid of the tri is inside the processor domain */
		   surf_area += compute_tri_area(t);
	       }
	    } 
       }
		
       if (pp_numnodes() > 1)
       {
           pp_global_sum(&surf_area, 1);
       }
       return surf_area;

} /*end compute_intfc_area*/


LOCAL  double compute_tri_area(
	TRI *tri)
{
        const double  *p[3];
        double        tri_area, s[3][3];
        static double n[3];
        int   i;

        p[0] = Coords(Point_of_tri(tri)[0]);
        p[1] = Coords(Point_of_tri(tri)[1]);
        p[2] = Coords(Point_of_tri(tri)[2]);

        for (i = 0; i < 3; ++i)
        {
            s[i][0] = p[Next_m3(i)][0] - p[i][0];
            s[i][1] = p[Next_m3(i)][1] - p[i][1];
            s[i][2] = p[Next_m3(i)][2] - p[i][2];
        }

        n[0] = s[2][1]*s[0][2] - s[2][2]*s[0][1];
        n[1] = s[2][2]*s[0][0] - s[2][0]*s[0][2];
        n[2] = s[2][0]*s[0][1] - s[2][1]*s[0][0];
	
        tri_area = 0.5*Mag3d(n);
	return tri_area;

}

