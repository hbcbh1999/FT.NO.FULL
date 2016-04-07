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
*				giglobs.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains global variables for the use in the gas dynamics code.
*	All globals in this file are local to the file, but can be
*	evaluated or set by calls to the functions in this file.
*
*
*	TODO:  Introduce a physics structure to pass these parameters
*
*	Contains
*
*		g_init_composition_type()
*		init_gravity()
*		max_num_comps()
*		g_prompt_for_maximum_number_of_components()
*		comp_type()
*		new_comp_type()
*		check_int_input()
*		check_float_input()
*/

#include <ginit/ginit.h>

LOCAL	int max_n_comps = 100;
LOCAL	COMP_TYPE	**ct = NULL;


	/* LOCAL Function Prototypes*/
LOCAL	COMP_TYPE	*new_comp_type(COMPONENT);
LOCAL	void	read_gravity_data(FILE*,int,double**,int);
#if defined(ONED)
LOCAL	void	set_oned_state_from_interface(double*,Locstate,COMP_TYPE*,
					      ONED_OVERLAY*);
#endif /* defined(ONED) */

EXPORT	void	g_prompt_for_composition_type(
	INIT_DATA	*init)
{
	int	i;
	int	dim = i_intfc(init)->dim;
	char	s[Gets_BUF_SIZE];
	static	Prompt_type ctypes[] = {
	    {"PURE_NON_REACTIVE","PNR",3,{PURE_NON_REACTIVE} },
#if defined(MULTI_COMPONENT)
	    {"MULTI_COMP_NON_REACTIVE","MCNR",4,{MULTI_COMP_NON_REACTIVE}},
#endif /* defined(MULTI_COMPONENT) */
#if defined(COMBUSTION_CODE)
	    {"PTFLAME","PTF",3,{PTFLAME}},
	    {"ZND","ZND",3,ZND},
	    {"TWO_CONSTITUENT_REACTIVE","TCR",3,{TWO_CONSTITUENT_REACTIVE}},
	    {"THINFLAME","THF",3,{THINFLAME}},
#endif /* defined(COMBUSTION_CODE) */
	    {NULL, NULL, 0, {ERROR} }
	};

	material_composition_type(init) = PURE_NON_REACTIVE;
	screen("\nRequest composition type of materials. "
	       "Available types are\n");
	for (i = 0; ctypes[i].prompt != NULL; ++i)
	{
	    screen("\t\t%s (%s",ctypes[i].prompt,ctypes[i].select);
	    if (ctypes[i].type.itype == material_composition_type(init))
	    	screen(", default");
	    screen(")\n");
	}
	screen("\tEnter choice here: ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    for (i = 0; ctypes[i].prompt != NULL; ++i)
	    {
	        if (strncasecmp(s,ctypes[i].select,ctypes[i].ncmp) == 0)
	        {
	            material_composition_type(init) = ctypes[i].type.itype;
	            break;
	        }
	    }
	}

        /* NOTE: NumberFloats(init) will be modified in g_compute_sizest() **/
        if(MULTI_COMP_NON_REACTIVE == material_composition_type(init))
        {
            (void) printf("Request for maximum number of composed materials: ");
            (void) Scanf("%d\n",&NumberFloats(init));
	    (void) printf("The maximum number of composed materials is %d\n",
                       NumberFloats(init));
        }
	g_compute_sizest(material_composition_type(init),
			 &StateSize(init),&NumberFloats(init),dim);
        /* This change should be universal */
        if(MULTI_COMP_NON_REACTIVE == material_composition_type(init))
        {
            g_set_sizeof_state(NULL,StateSize(init),NumberFloats(init));
            set_composition_type(material_composition_type(init));
        }
}		/*end g_prompt_for_composition_type*/

/*
*			g_init_composition_type():
*/


/*ARGSUSED*/
EXPORT int g_init_composition_type(
	INIT_PHYSICS	*ip,
	INIT_DATA	*init,
	size_t		*sizest,
	int		*nfloats)
{
	*sizest = StateSize(init);
	*nfloats = NumberFloats(init);
	return material_composition_type(init);
}		/*end g_init_composition_type*/


EXPORT	void	g_compute_sizest(
	int	composition_type,
	size_t	*sizest,
	int	*nfloats,
	int	dim)
{
	size_t	sizest_not_returned;
	int	nfloats_not_returned;

	if (sizest == NULL)
	    sizest = &sizest_not_returned;
	if (nfloats == NULL)
	    nfloats = &nfloats_not_returned;
	switch (composition_type)
	{
	    case PURE_NON_REACTIVE:
	    	*sizest = sizeof(Gas);
	    	*nfloats = dim + 2;
	    	break;
            case MULTI_COMP_NON_REACTIVE:
                *sizest = sizeof(MGas) + ((*nfloats)-1)*sizeof(double);
                *nfloats = dim + 2 + (*nfloats);
                break;
#if defined(COMBUSTION_CODE)
	    case PTFLAME:
	    	*sizest = sizeof(MGas);
	    	/* the 5th variable = 1. if BURNED, 0. if UNBURNED */
	    	*nfloats = dim + 3;
	    	break;
	    case ZND:
	    	*sizest = sizeof(MGas);
	    	*nfloats = dim + 3;
	    	break;
	    case TWO_CONSTITUENT_REACTIVE:
	    	*sizest = sizeof(MGas) + sizeof(double);
	    	*nfloats = dim + 4;
	    	break;
	    case THINFLAME:
                *sizest = sizeof(Gas);
                *nfloats = dim + 2;
                break;
#endif /* defined(COMBUSTION_CODE) */
	    default:
	    	screen("ERROR in g_compute_sizest(), "
		       "unknown composition type\n");
	    	clean_up(ERROR);
	}
}		/*end g_compute_sizest*/

/*
*			prompt_for_gravity():
*/

EXPORT void prompt_for_gravity(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	RECT_GRID	*gr = &Comp_grid(init);
	double		*L, *U;
	char		s[Gets_BUF_SIZE];
	double		**m;
	int		c, len;
	static char	dname[3][2] = {"x", "y", "z"};
	double		*g;
	GRAVITY		*grav_data;
	int		i, dim;

	gravity_data(init) = NULL;
	screen("The following choices are available for a gravitational "
	       "acceleration\n"
	       "\tNo gravity (N or default)\n"
	       "\tConstant gravity (C or Y)\n"
	       "\tTime dependent gravity (T)\n"
	       "\tAstrophysical (central force) gravity (A)\n"
	       "\tGeneralized Astrophysical gravity (G)\n"
	       "\tRadial gravity with constant magnitude (R)\n");
	screen("Enter choice: ");
	(void) Gets(s);

	scalar(&grav_data,sizeof(GRAVITY));
	gravity_data(init) = grav_data;

	dim = gr->dim;
	grav_data->dim = dim;
	switch (s[0])
	{
	case 'N':
	case 'n':
	case '\0':
	    grav_data->type = NO_GRAVITY;
	    return;
	case 'C':
	case 'c':
	case 'Y': /*For historical compatibility of input files*/
	case 'y': /*For historical compatibility of input files*/
	    grav_data->type = CONSTANT_GRAVITY;
	    g = grav_data->g;
	    g[0] = g[1] = g[2] = 0.0;
	    for (i = 0; i < dim; ++i)
	    {
	        screen("\tEnter %s component of gravity (dflt = 0): ",dname[i]);
	        (void) Gets(s);
	        if (s[0] != '\0')
	    	    (void) sscan_float(s,g+i);
	    }
	    return;
	case 'T':
	case 't':
	    screen("Gravity data consists of an array of uni_arrays (time, g)\n");
	    screen("This data can be entered directly or input from a file\n");
	    screen("Enter either a filename or the number of data points\n");
	    screen("Enter filename or integer: ");
	    (void) Gets(s);
	    if (s[0] == '\0')
	    {
	        grav_data->type = NO_GRAVITY;
	        return;
	    }
	    grav_data->type = TIME_DEPENDENT_GRAVITY;
	    /*Was an integer input?*/
	    len = (int) strlen(s);
	    for (i = 0; i < len; ++i)
	    {
		c = s[i];
		if (!isdigit(c))
		    break;
	    }
	    if (i == len)
	    {
		(void) sscanf(s,"%d",&grav_data->num_time_points);
		bi_array(&m,grav_data->num_time_points,dim+1,FLOAT);
		grav_data->g_of_t = m;
		screen("Enter %d data points consisting of a monotonically "
		       "increasing time value\n",grav_data->num_time_points);
		screen("\tfollowed by a %d vector of gravity "
		       "values for that time\n",dim);
		read_gravity_data(stdin,dim,m,grav_data->num_time_points);
	    }
	    else
	    {
	        FILE	*file = fopen(s,"r");
		int	nf;
		double	x;
		if (file == NULL)
		{
		    screen("ERROR in prompt_for_gravity(), can't open %s",s);
		    clean_up(ERROR);
		}
		/*count entries in file*/
		for (nf = 0; ((c = getc(file)) != EOF); ++nf)
		{
		    (void)ungetc(c,file);
		    if (fscan_float(file,&x) != 1)
			break;
		}
		if ((nf%(dim+1)) != 0)
		{
		    screen("ERROR in prompt_for_gravity(), invalid file\n");
		    clean_up(ERROR);
		}
		grav_data->num_time_points = nf/(dim+1);
		rewind(file);
		bi_array(&m,grav_data->num_time_points,dim+1,FLOAT);
		grav_data->g_of_t = m;
		read_gravity_data(file,dim,m,grav_data->num_time_points);
	    }
	    /*Check for monotonicity of time values*/
	    for (i = 1; i < grav_data->num_time_points; ++i)
	    {
		if (m[i][0] <= m[i-1][0])
		{
		    screen("ERROR in prompt_for_gravity(), "
			   "invalid time values\n");
		    clean_up(ERROR);
		}
	    }
	    return;
	case 'A':
	case 'a':
	    grav_data->type = ASTROPHYSICAL_GRAVITY;
	    grav_data->G = 0.0;
	    grav_data->M = 0.0;
	    L = gr->GL;
	    U = gr->GU;
	    for (i = 0; i < dim; ++i)
		grav_data->center[i] = 0.5*(L[i]+U[i]);
	    screen("Enter the coordinates of the gravity center (dflt =");
	    for (i = 0; i < dim; ++i)
		screen(" %g",grav_data->center[i]);
	    screen("): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
		double	*cen = grav_data->center;
		static const char *fmt = "%lf %lf %lf";
		if (sscanf(s,fmt,cen,cen+1,cen+2) != dim)
		{
		    screen("ERROR in prompt_for_gravity(), "
			   "invalid coordinate for gravity center\n");
		    clean_up(ERROR);
		}
	    }
	    screen("Enter the gravitational constant (dflt = %g): ",
		   grav_data->G);
	    (void) Gets(s);
	    if (s[0] != '\0')
		sscan_float(s,&grav_data->G);
	    screen("Enter the mass of the gravity center (dflt = %g): ",
		   grav_data->M);
	    (void) Gets(s);
	    if (s[0] != '\0')
		sscan_float(s,&grav_data->M);
	    return;
	case 'R':
	case 'r':
	    grav_data->type = RADIAL_GRAVITY;
	    grav_data->G = 0.0;
	    grav_data->M = 0.0;
	    L = gr->GL;
	    U = gr->GU;
	    for (i = 0; i < dim; ++i)
		grav_data->center[i] = 0.5*(L[i]+U[i]);
	    screen("Enter the coordinates of the gravity center (dflt =");
	    for (i = 0; i < dim; ++i)
		screen(" %g",grav_data->center[i]);
	    screen("): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
		double	*cen = grav_data->center;
		static const char *fmt = "%lf %lf %lf";
		if (sscanf(s,fmt,cen,cen+1,cen+2) != dim)
		{
		    screen("ERROR in prompt_for_gravity(), "
			   "invalid coordinate for gravity center\n");
		    clean_up(ERROR);
		}
	    }
	    screen("Enter the magnitude of the gravitational acceleration "
		   "(dflt = %g): ",grav_data->G);
	    (void) Gets(s);
	    if (s[0] != '\0')
		sscan_float(s,&grav_data->G);
	    return;
	case 'G':
	case 'g':
	    grav_data->type = GENERALIZED_ASTROPHYSICAL_GRAVITY;
            grav_data->roots = &ip->root;
            grav_data->nroots = 1;
	    grav_data->G = 0.0;
	    grav_data->M = 0.0;
	    grav_data->N = 100;
	    for (i = 0; i < dim; ++i)
		grav_data->center[i] = 0.0;
	    screen("Enter the coordinates of the gravity center (dflt =");
	    for (i = 0; i < dim; ++i)
		screen(" %g",grav_data->center[i]);
	    screen("): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
		double	*cen = grav_data->center;
		static const char *fmt = "%lf %lf %lf";

		if (sscanf(s,fmt,cen,cen+1,cen+2) != dim)
		{
		    screen("ERROR in prompt_for_gravity(), "
			   "invalid coordinate for gravity center\n");
		    clean_up(ERROR);
		}
	    }
	    screen("Enter the gravitational constant (dflt = %g): ",
		   grav_data->G);
	    (void) Gets(s);
	    if (s[0] != '\0')
		sscan_float(s,&grav_data->G);

            for(i = 0; i < dim; ++i)
            {
                grav_data->GL[i] = gr->GL[i];
                grav_data->GU[i] = gr->GU[i];
                grav_data->gmax[i] = gr->gmax[i];
            }

	    screen("Enter the number of rings the whole domain"
		   "should be divided into (dflt = %d): ", grav_data->N);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscanf(s,"%d",&grav_data->N);
	    /*grav_data->Mp = (double*)malloc(grav_data->N*sizeof(double)); */
	    return;
	default:
	    screen("ERROR in prompt_for_gravity(), "
		   "invalid choice %s\n",s);
	    clean_up(ERROR);
	}
}		/*end prompt_for_gravity*/


LOCAL	void	read_gravity_data(
	FILE	*file,
	int	dim,
	double	**m,
	int	npts)
{
	int i;
	static const char *fmts[] = {"","%lf %lf","%lf %lf %lf","%lf %lf %lf"};
	const char *fmt = fmts[dim];

	switch (dim)
	{
	case 1:
	    for (i = 0; i < npts; ++i)
	    {
		if (fscanf(file,fmt,&m[i][0],&m[i][1]) != 2)
		{
		    screen("ERROR in read_gravity_data(), "
			   "invalid input of gravity data\n");
		    clean_up(ERROR);
		}
	    }
	    return;
	case 2:
	    for (i = 0; i < npts; ++i)
	    {
	        if (fscanf(file,fmt,&m[i][0],&m[i][1],&m[i][2]) != 3)
	        {
	            screen("ERROR in prompt_for_gravity(), "
	                   "invalid input of gravity data\n");
	            clean_up(ERROR);
	        }
	    }
	    return;
	case 3:
	    for (i = 0; i < npts; ++i)
	    {
	        if (fscanf(file,fmt,&m[i][0],&m[i][1],&m[i][2],&m[i][3]) != 4)
	        {
	            screen("ERROR in prompt_for_gravity(), "
	                   "invalid input of gravity data\n");
	            clean_up(ERROR);
	        }
	    }
	    return;
	}
}		/*end read_gravity_data*/

 
EXPORT int max_num_comps(void)
{
	return max_n_comps;
}		/*end max_num_comps*/
 

EXPORT void g_prompt_for_maximum_number_of_components(void)
{
	char		s[Gets_BUF_SIZE];
 
	screen("Enter an upper bound for the "
	       "number of components (default = %d): ",max_n_comps);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%d",&max_n_comps);
	max_n_comps += FIRST_DYNAMIC_COMPONENT;
	check_int_input("for the number of components",max_n_comps,1,1,GE_);
}		/*end g_prompt_for_maximum_number_of_components*/
 
EXPORT	const COMPONENT *comps_with_type(
	COMP_TYPE_TYPE type,
	int            *num)
{
	static COMPONENT *compstore = NULL;
	static int       len_compstore = 0;
	int              i, N;

	if (len_compstore <= max_n_comps)
	{
	    if (compstore != NULL)
	        free(compstore);
	    len_compstore = max_n_comps+1;
	    uni_array(&compstore,len_compstore,sizeof(COMPONENT));
	}
	for (i = 0; i <= max_n_comps; ++i)
	    compstore[i] = NO_COMP;
	for (N = 0, i = 0; i <= max_n_comps; ++i)
	{
	    if ((ct[i] != NULL) && (ct[i]->type == type))
	        compstore[N++] = ct[i]->comp;
	}
	if (num != NULL)
	    *num = N;
	return (N > 0) ? compstore : NULL;
}		/*end comps_with_type*/

#if defined(ONED)
EXPORT	boolean	overlays_exist(void)
{
	int i;
	for (i = 0; i <= max_n_comps; ++i)
	{
	    if ((ct[i] != NULL) && (ct[i]->type == ONE_DIMENSIONAL_OVERLAY))
	    {
		ONED_OVERLAY *olay = One_d_overlay(ct[i]);
		if (olay->overlay_type != OVERLAY_TYPE_UNSET)
		    return YES;
	    }
	}
	return NO;
}		/*end overlays_exist*/
#endif /* defined(ONED) */


EXPORT	COMP_TYPE *comp_type(
	COMPONENT	comp)
{
	static	int		maxcomp = -1;

	if (ct == NULL)
	{
	    int i;
	    uni_array(&ct,max_n_comps+1,sizeof(COMP_TYPE*));
	    for (i = 0; i <= max_n_comps; ++i)
	    	ct[i] = NULL;
	}

	if (comp < 0)
	    comp = maxcomp + 1;

	if (comp > max_n_comps)
	{
	    COMP_TYPE	**new_ct;
	    int		i, new_max_n_comps = max(comp,2*max_n_comps);

	    uni_array(&new_ct,new_max_n_comps+1,sizeof(COMP_TYPE*));
	    for (i = 0; i <= max_n_comps; ++i)
	    	new_ct[i] = ct[i];
	    for (;i <= new_max_n_comps; ++i)
	    	new_ct[i] = NULL;
	    free(ct);
	    ct = new_ct;
	    max_n_comps = new_max_n_comps;
	}

	if (maxcomp < comp)
	    maxcomp = comp;

	if (ct[comp] == NULL)
	    ct[comp] = new_comp_type(comp);
	else if (ct[comp]->comp != comp)
	{
	    screen("\nERROR in comp_type(), Inconsistent component number!\n");
	    clean_up(ERROR);
	}

	return	ct[comp];
}		/*end comp_type*/

EXPORT	void	free_comp_types(void)
{
	int	i;

	for (i = 0; i <= max_n_comps; ++i)
	{
	    if (ct[i] != NULL)
	    {
	    	if (ct[i]->free_comp_type_extra != NULL)
	    	    (*ct[i]->free_comp_type_extra)(ct[i]);
	    	free(ct[i]);
	    }
	}
	free(ct);
	ct = NULL;
	max_n_comps = -1;
}		/*end free_comp_types*/


LOCAL	COMP_TYPE	*new_comp_type(
	COMPONENT	comp)
{
	COMP_TYPE	*newct;

	scalar(&newct, sizeof(COMP_TYPE));
	newct->type = UNSET_COMP_TYPE;

	newct->comp = comp;
	newct->params = NULL;
	return	newct;
}		/*end new_comp_type*/


/*ARGSUSED*/
EXPORT	void	get_state_ambient(
	double		*coords,
	Locstate	s,
	COMP_TYPE	*ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		stype)
{
	Locstate	ct_state = Ambient(ct);
        int             save_stype = stype;

	debug_print("init_states","Entered get_state_ambient()\n");
	if ( ct->type != AMBIENT )
	{
	    screen("\nERROR in get_state_ambient(), "
		   "inconsistent comp_type->type\n");
	    clean_up(ERROR);
	}

	if ( ct_state != NULL )
        {
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                stype = TGAS_STATE;
            set_state(s,stype,ct_state);
        }
	else
	{
	    screen("ERROR in get_state_ambient(), NULL state encountered\n");
	    clean_up(ERROR);
	}

        /* This is for the oned sod-tube test problem */
        /* init partial density */
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(debugging("sod_tube"))
            {
                Gas_param *params = Params(s);
                int    num_comps, i;
                double   *rho;
                if((num_comps = params->n_comps) != 1)
                {
                    rho = pdens(s);
                    if(coords[0] >= 0.0 && coords[0] <= 0.5)
                        rho[0] = 0.8*Dens(s);
                    else if(coords[0] > 0.5 && coords[0] <= 0.7)
                        rho[0] = 0.3*Dens(s);
                    else
                        rho[0] = 0.1*Dens(s);
                    rho[1] = 0.19*sqr(sin(20*3.1415927*coords[0]))*Dens(s);
                    rho[2] = Dens(s) - rho[1] - rho[0];
                    if(rho[0] < 0.0 || rho[1] < 0.0 || rho[2] < 0.0)
                    {
                        printf("ERROR: get_state_ambient, partial den < 0.0\n");
                        clean_up(ERROR);
                    }
                }
            }
            else if(debugging("blast_wave"))
            {
                Gas_param *params = Params(s);
                int    num_comps, i;
                double   *rho;
                if((num_comps = params->n_comps) != 1)
                {
                    rho = pdens(s);
                    rho[0] = 0.5*sqr(coords[0])*Dens(s);
                    rho[1] = 0.5*sqr(sin(20*coords[0]))*Dens(s);
                    rho[2] = Dens(s) - rho[1] - rho[0];
                    if(rho[0] < 0.0 || rho[1] < 0.0 || rho[2] < 0.0)
                    {
                        printf("ERROR: get_state_ambient, partial den < 0.0\n");
                        clean_up(ERROR);
                    }
                }
            }
            set_state(s,save_stype,s);
            /* verbose_print_state("get_state_ambient", s); */
        }
}		/*end get_state_ambient*/

/*ARGSUSED*/
EXPORT	void	get_state_2d_riemann(
	double		*coords,
	Locstate	s,
	COMP_TYPE	*ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		stype)
{
	int		i,num_sections;
	double		angle,x,y,xm,ym;
	double		*start_angle;
	Locstate	*states;
	RECT_GRID	*gr = computational_grid(intfc);

	debug_print("init_states","Entered get_state_2d_riemann()\n");

	num_sections = Riemann_2d(ct)->num_sections;
	start_angle = Riemann_2d(ct)->start_angle;
	states = Riemann_2d(ct)->states;

	xm = 0.5*(gr->GL[0] + gr->GU[0]);
	ym = 0.5*(gr->GL[1] + gr->GU[1]);
	x = coords[0] - xm;
	y = coords[1] - ym;

	if (x == 0.0)
	    angle = (x > 0.0) ? 0.0 : PI;
	else
	{
	    angle = atan(fabs(y/x));
	    if (x < 0.0 && y > 0.0) 
		angle = PI - angle;
	    else if (x < 0.0 && y < 0.0) 
		angle = PI + angle;
	    else if (x > 0.0 && y < 0.0) 
		angle = 2.0*PI - angle;
	}

	if (angle < start_angle[0] || angle >= start_angle[num_sections-1])
	    set_state(s,stype,states[num_sections-1]);
	else
	{
	    for (i = 0; i < num_sections-1; ++i)
	    {
		if (angle >= start_angle[i] && angle < start_angle[i+1])
	    	    set_state(s,stype,states[i]);
	    }
	}
}	/*end get_state_2d_riemann*/

#if defined(ONED)
/*ARGSUSED*/
EXPORT	void	get_state_1d_overlay(
	double		*coords,
	Locstate	s,
	COMP_TYPE	*ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		stype)
{
	ONED_OVERLAY	*olay = (ONED_OVERLAY*)ct->extra;
	INPUT_SOLN	**is = olay->is;
	RECT_GRID	*gr = &is[0]->grid;
	INTERFACE	*intfc1d = olay->intfc1d;
	double		x, coords1d[3], m, sp;
	double		dcoords[3];
	double		c0[3], c1[3];
	double		alpha;
	int		i0[3], i1[3], i, dim = ct->params->dim;
	int		nvars = olay->prt->n_restart_vars;
	static Locstate	s0 = 0, s1 = 0;

	debug_print("init_states","Entered get_state_1d_overlay()\n");
	if (s0 == NULL)
	{
	    (*ct->params->_alloc_state)(&s0,ct->params->sizest);
	    (*ct->params->_alloc_state)(&s1,ct->params->sizest);
	}

	for (i = 0; i < dim; ++i)
	    dcoords[i] = coords[i] - olay->origin[i];
	switch (olay->overlay_type)
	{
	case RADIAL_OVERLAY:
	    x = mag_vector(dcoords,dim);
	    break;
	case CYLINDRICAL_OVERLAY:
	    sp = scalar_product(dcoords,olay->direction,dim);
	    for (i = 0; i < dim; ++i)
	    	dcoords[i] -= sp*olay->direction[i];
	    x = mag_vector(dcoords,dim);
	    break;
	case RECTANGULAR_OVERLAY:
	    x = scalar_product(dcoords,olay->direction,dim);
	    break;
	case OVERLAY_TYPE_UNSET:
	default:
	    x = -HUGE_VAL;
	    screen("ERROR in get_state_1d_overlay(), unset overlay type\n");
	    clean_up(ERROR);
	    break;
	}
	if (x > gr->U[0])
	    x = gr->U[0];
	if (x < gr->L[0])
	    x = gr->L[0];
	coords1d[0] = x;
	if (rect_in_which(coords1d,i0,gr) == FUNCTION_FAILED)
	{
	    screen("ERROR in get_state_1d_overlay(), rect_in_which() failed\n");
	    print_general_vector("dcoords = ",dcoords,dim,"\n");
	    (void) printf("overlay type = %d, x = %g\n",olay->overlay_type,x);
	    (void) printf("rectangular grid gr\n");
	    print_rectangular_grid(gr);
	    (void) printf("One dimensional interface\n");
	    print_interface(intfc1d);
	    clean_up(ERROR);
	}
	c0[0] = cell_center(i0[0],0,gr);
	if (x < c0[0])
	{
	    i1[0] = i0[0];
	    i0[0]--;
	    if (i0[0] < 0)
		i0[0] = 0;
	    c1[0] = c0[0];
	    c0[0] = cell_center(i0[0],0,gr);
	}
	else
	{
	    i1[0] = i0[0] + 1;
	    if (i1[0] >= gr->gmax[0])
	    	i1[0] = i0[0];
	    c1[0] = cell_center(i1[0],0,gr);
	}
	if (component(c0,intfc1d) != ct->comp)
	{
	    set_oned_state_from_interface(c0,s0,ct,olay);
	}
	else
	{
	    g_restart_initializer(i0,nvars,ct->comp,is,s0,init);
	    Init_params(s0,ct->params);
	}
	if (component(c1,intfc1d) != ct->comp)
	{
	    set_oned_state_from_interface(c1,s1,ct,olay);
	}
	else
	{
	    g_restart_initializer(i1,nvars,ct->comp,is,s1,init);
	    Init_params(s1,ct->params);
	}
	alpha = (c1[0] > c0[0]) ? (x - c0[0])/(c1[0] - c0[0]) : 0.5;
	ct->params->dim = 1;
	interpolate_states(olay->front,1.0-alpha,alpha,c0,s0,c1,s1,s);
	ct->params->dim = dim;

	m = Mom(s)[0];
	switch (olay->overlay_type)
	{
	case RADIAL_OVERLAY:
	case CYLINDRICAL_OVERLAY:
	    if (x > 0.0)
	    {
	        for (i = 0; i < dim; ++i)
	    	    Mom(s)[i] = m*dcoords[i]/x;
	    }
	    else
	    {
	        for (i = 0; i < dim; ++i)
	    	    Mom(s)[i] = 0.0;
	    }
	    break;
	case RECTANGULAR_OVERLAY:
	    for (i = 0; i < dim; ++i)
	    	Mom(s)[i] = m*olay->direction[i];
	    break;
	case OVERLAY_TYPE_UNSET:
	default:
	    screen("ERROR in get_state_1d_overlay(), unset overlay type\n");
	    clean_up(ERROR);
	    break;
	}
	set_state(s,stype,s);
}		/*end get_state_1d_overlay*/

LOCAL	void	set_oned_state_from_interface(
	double	     *coords,
	Locstate     s,
	COMP_TYPE    *ct,
	ONED_OVERLAY *olay)
{
	HYPER_SURF		*hs;
	HYPER_SURF_ELEMENT	*hse;
	POINT			*p;
	double			coords_on[MAXD];

	if (nearest_interface_point(coords,ct->comp,olay->intfc1d,
				    INCLUDE_BOUNDARIES,NULL,
				    coords_on,NULL,&hse,&hs) != YES)
	{
	    screen("ERROR in set_oned_state_from_interface(), "
	           "nearest_interface_point() failed\n");
	    clean_up(ERROR);
	}
	coords[0] = coords_on[0];
	p = Point_of_hs(hs);
	if (ct->comp == positive_component(hs))
	    ft_assign(s,right_state(p),ct->params->sizest);
	else if (ct->comp == negative_component(hs))
	    ft_assign(s,left_state(p),ct->params->sizest);
	else
	{
	    screen("ERROR in set_oned_state_from_interface(), "
	           "ct->comp not on interface\n");
	    clean_up(ERROR);
	}
	Init_params(s,ct->params);
}		/*end set_oned_state_from_interface*/
#endif /* defined(ONED) */


/*
*		check_int_input():
*		check_float_input();
*
*	Produces a fatal error if i is not in the desired range and defined
*	by il, iu,  and opt.  This fatal error will be produced if any of
*	the following conditions are TRUE.
*
*	opt = GE_AND_LE and i < il or i > iu.
*	opt = LE_OR_GE and il < i < iu.
*	opt = GE_ and i < il.
*	opt = _LE and i > iu.
*/

EXPORT	void	check_int_input(
	const char *mesg,
	int	   i,
	int	   il,
	int	   iu,
	int	   opt)
{
	int		wrong = 0;
	char		s[100];
 
	switch ( opt )
	{
	case GE_AND_LE:
	    if (i < il || i > iu)
		wrong = 1;
	    (void) sprintf(s,">= %d and <= %d",il,iu);
	    break;
	case LE_OR_GE:
	    if (i > il && i < iu)
		wrong = 1;
	    (void) sprintf(s,"<= %d, or >= %d",il,iu);
	    break;
	case GE_:
	    if (i < il)
		wrong = 1;
	    (void) sprintf(s,">= %d",il);
	    break;
	case _LE:
	    if (i > iu)
		wrong = 1;
	    (void) sprintf(s,"<= %d",iu);
	    break;
	default:
	    break;
	}
 
	if ( wrong )
	{
	    screen("\nERROR in check_int_input(), "
	           "The input value %d %s is out of range!\n"
	           "It should be %s.\n",i,mesg,s);
	    clean_up(ERROR);
	}
}		/*end check_int_input*/
 

EXPORT	void	check_float_input(
	const char *mesg,
	double	   i,
	double	   il,
	double	   iu,
	int	   opt)
{
	boolean		wrong = NO;
	char		s[100];
 
	switch (opt)
	{
	case GE_AND_LE:
	    if (i < il || i > iu)
		wrong = YES;
	    (void) sprintf(s,">= %g and <= %g",il,iu);
	    break;
	case LE_OR_GE:
	    if (i > il && i < iu)
		wrong = YES;
	    (void) sprintf(s,"<= %g, or >= %g",il,iu);
	    break;
	case GE_:
	    if (i < il)
		wrong = YES;
	    (void) sprintf(s,">= %g",il);
	    break;
	case _LE:
	    if (i > iu)
		wrong = YES;
	    (void) sprintf(s,"<= %g",iu);
	    break;
	default:
	    break;
	}
 
	if (wrong == YES)
	{
	    screen("\nERROR in check_float_input(), "
	           "The input value %g %s is out of range!\n"
	           "It should be %s.\n",i,mesg,s);
	    clean_up(ERROR);
	}
}		/*end check_float_input*/
