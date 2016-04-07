/*
*			gcauchy.c:
*
*	Contains initialization routines for cauchy depositions
*/

#include <ghyp/ghyp.h>
#define MAX_DEPOSITION_NUMBER 10

LOCAL	void	add_ext_deposition(Locstate,COMPONENT,double*,int);

int	dep_number;

enum DEPOSITION_REGION {
	DEP_RECT,
	DEP_ELLIP
};

enum DEPOSITION_PROFILE {
	GAUSSIAN,
	POWER,
	SINUSOIDAL,
	LINEAR,
	CONSTANT_PROFILE
};

enum DEPOSITION_TYPE {
	EXT_ENERGY,
	EXT_PRESSURE,
	EXT_TEMPER,
	EXT_VX,
	EXT_VY,
	EXT_VZ,
	EXT_DENSITY
};
const char dep_name[] = "EPTXYZD";

struct	DEPOSITION_PARAMS { /* dim=0: time, dim>0: space */
	int			comp;
	enum DEPOSITION_TYPE	type;
	double			base; /* base deposition */
	double			amp; /* amplitude */
	enum DEPOSITION_REGION	region;
	double			bd[MAXD+1][2]; /* RECT: 0=lower, 1=upper; ELLIP: 0=center, 1=radii */
	enum DEPOSITION_PROFILE	prof[MAXD+1]; /* profile */
	double			orig[MAXD+1]; /* origin of profile */
	double			param[MAXD+1][3];
	boolean			flag[MAXD+1][2];
} deparam[MAX_DEPOSITION_NUMBER];


EXPORT	void	g_cauchy_deposition(
	Wave	*wave,
	Front	*front)
			
{				
	char    axis[]={'t','x','y','z'}, s[Gets_BUF_SIZE];
	int     i, d, dim = wave->rect_grid->dim;
	int	gridmax, vgmax[MAXD], tmp, ic[MAXD];

	debug_print("deposition","Entered g_cauchy_deposition()\n");

	/* prompt for Cauchy deposition */
	screen("Type \'y\' to set cauchy deposition: ");
	Gets(s);
	if (s[0] != 'Y' && s[0] != 'y') return;

	screen("Enter the number of depositions (1 to %d): ",MAX_DEPOSITION_NUMBER);
	Scanf("%d\n", &dep_number);
	if (dep_number < 1 || dep_number > MAX_DEPOSITION_NUMBER)
	{
	    screen("Illegal number of depositions.\n");
	    clean_up(ERROR);
	}
	
	for (i = 0; i < dep_number; i++)
	{
	    screen("For deposition %d\n", i);
    	    
	    screen("Enter the component number for the deposition: ");
	    Scanf("%d\n",&deparam[i].comp);

	    screen("Choose the deposition type from\n");
	    screen("\tEnergy(dflt), Pressure, Temperature, X-Velocity,\n");
	    screen("\tY-Velocity, Z-Velocity, Density\n");
	    screen("Enter the deposition type: ");
	    Gets(s);
	    if (s[0] == 'P' || s[0] == 'p')
		deparam[i].type = EXT_PRESSURE;
	    else if (s[0] == 'T' || s[0] == 't')
		deparam[i].type = EXT_TEMPER;
	    else if (s[0] == 'X' || s[0] == 'x')
		deparam[i].type = EXT_VX;
	    else if (s[0] == 'Y' || s[0] == 'y')
		deparam[i].type = EXT_VY;
	    else if (s[0] == 'Z' || s[0] == 'z')
		deparam[i].type = EXT_VZ;
	    else if (s[0] == 'D' || s[0] == 'd')
		deparam[i].type = EXT_DENSITY;
	    else
		deparam[i].type = EXT_ENERGY;
	    printf("deparam[%d].type = %d\n",i,deparam[i].type);
	    
	    screen("Enter the base value of deposition: ");
	    Scanf("%f\n",&deparam[i].base);
	    screen("Enter the amplitude of deposition: ");
	    Scanf("%f\n",&deparam[i].amp);

	    screen("Enter the shape of deposition region (Rect or Ellip): ");
	    Gets(s);
	    if (s[0] == 'E' || s[0] == 'e')
		deparam[i].region = DEP_ELLIP;
	    else
		deparam[i].region = DEP_RECT;
	    
	    for (d = 0; d < dim+1; d++)
	    {
		screen("Along %c direction\n",axis[d]);
		if (d == 0 || deparam[i].region == DEP_RECT)
		{
		    deparam[i].bd[d][0] = -HUGE_VAL;
		    deparam[i].bd[d][1] = HUGE_VAL;
		    screen("Enter the lower and upper boundary: ");
		}
		else
		{
		    deparam[i].bd[d][0] = 0;
		    deparam[i].bd[d][1] = HUGE_VAL;
		    screen("Enter the center and radius: ");
		}

		Scanf("%f %f\n",deparam[i].bd[d],deparam[i].bd[d]+1);
		if (deparam[i].bd[d][0] == deparam[i].bd[d][1])
		{
		    deparam[i].prof[d] = CONSTANT_PROFILE;
		    continue; /* depo layer length[d] = 0: CONSTANT */
		}

		screen("Choose deposition profile from\n");
		screen("\tGaussian, Power, Sinusoidal, Constant(dflt)\n"); 
		screen("Enter the deposition profile: ");
		Gets(s);
		if (s[0] == 'g' || s[0] == 'G')
		{
		    deparam[i].prof[d] = GAUSSIAN;

		    screen("Enter the origin of the profile: ");
		    Scanf("%f\n",&deparam[i].orig[d]);
		    screen("Enter the decay length: ");
		    Scanf("%f\n",deparam[i].param[d]);
		    screen("Enter if it is exponential decay: ");
		    Gets(s);
		    if (s[0] == 'y' || s[0] == 'Y')
			deparam[i].flag[d][0] = TRUE;
		    else
			deparam[i].flag[d][0] = FALSE;
		}
		else if (s[0] == 'p' || s[0] == 'P')
		{
		    deparam[i].prof[d] = POWER;

		    screen("Enter the origin of the profile: ");
		    Scanf("%f\n",&deparam[i].orig[d]);
		    screen("Enter the decay exponent: ");
		    Scanf("%f\n",deparam[i].param[d]);
		}
		else if (s[0] == 's' || s[0] == 'S')
		{
		    deparam[i].prof[d] = SINUSOIDAL;

		    screen("Enter the origin of the profile: ");
		    Scanf("%f\n",&deparam[i].orig[d]);
		    screen("Enter the period: ");
		    Scanf("%f\n",deparam[i].param[d]);
		    screen("Enter the phase at origin in degrees: ");
		    Scanf("%f\n",deparam[i].param[d]+1);
		} 
		else
		{
		    deparam[i].prof[d] = CONSTANT_PROFILE;
		}
	    }
	}

	return;

	
	/* Modify front states  */
	/*
	HYPER_SURF		*hs;
	HYPER_SURF_ELEMENT 	*hse;
	POINT			*pt;
	Locstate	 	stl, str;

	debug_print("deposition","Front deposition\n");
	next_point(front->interf,NULL,NULL,NULL);
	while (next_point(front->interf,&pt,&hse,&hs))
	{	
	    slsr(pt,hse,hs,&stl,&str);
	    add_ext_deposition(str,positive_component(hs),Coords(pt),dim);
	    add_ext_deposition(stl,negative_component(hs),Coords(pt),dim);
	}
	*/

	printf("#modify interior states.\n");
	
	/* Modify interior states  */

	debug_print("deposition","\nInterior deposition\n");
	
	gridmax = 1;
	for (d = 0; d < dim; d++)
	{
	    vgmax[d] = wave->rect_grid->gmax[d] + wave->rect_grid->lbuf[d]
			+  wave->rect_grid->ubuf[d];
	    gridmax *= vgmax[d];
	}
	for (i = 0; i < gridmax; i++)
	{
	    tmp = i;
	    for (d = 0; d < dim; d++)
	    {
		ic[d] = tmp % vgmax[d] - wave->rect_grid->lbuf[d];
		tmp /= vgmax[d];	
	    }

	    add_ext_deposition(Rect_state(ic,wave),Rect_comp(ic,wave),Rect_coords(ic,wave),dim);
	}
	
	debug_print("deposition","\nLeft g_cauchy_deposition()\n");
}	/* end of set_ext_deposition */

LOCAL	void		add_ext_deposition(
	Locstate	state,
	COMPONENT	comp,
	double		*coords,
	int		dim)
{
	double	dep, tcrds[MAXD+1];
	int	i, d;

	if (is_obstacle_state(state)) return;

	/* set tcrds */
	tcrds[0] = 0;
	for (d = 0; d < dim; d++)  tcrds[d+1] = coords[d];

	for (i = 0; i < dep_number; i++)
	{
	    if (comp != -1 && comp != deparam[i].comp) continue; /* comp = -1 for all components */
	    
	    /* check if within deposition boundary */
	    if ( tcrds[0] < deparam[i].bd[0][0] || tcrds[0] > deparam[i].bd[0][1] )
		continue; /* check time */
	    if (deparam[i].region == DEP_RECT)
	    {
		for (d = 1; d < dim+1; d++)
		    if ( tcrds[d] < deparam[i].bd[d][0] || tcrds[d] > deparam[i].bd[d][1] )
			    break;
		if (d < dim+1)  continue;
	    }
	    else /* DEP_ELLIP */
	    {
		double	sum = 0;
		for (d = 1; d < dim+1; d++)
		    sum += sqr((tcrds[d]-deparam[i].bd[d][0])/deparam[i].bd[d][1]);
		if (sum > 1) continue;
	    }

	    dep = deparam[i].amp;
	    
	    for (d = 0; d < dim+1; d++)
	    {
		if (deparam[i].prof[d] == GAUSSIAN)
		{   /* param[d][0] = decay length */
		    if (deparam[i].flag[d][0] == FALSE) /* Gaussian */
			dep /= exp(sqr((tcrds[d]-deparam[i].orig[d])/deparam[i].param[d][0]));
		    else /* exponential */
			dep /= exp((tcrds[d]-deparam[i].orig[d])/deparam[i].param[d][0]);
		}
		else if (deparam[i].prof[d] == POWER)
		{   /* param[d][0] = decay exponent */
		    dep /= pow(fabs(tcrds[d]-deparam[i].orig[d]),deparam[i].param[d][0]);
		}
		else if (deparam[i].prof[d] == SINUSOIDAL)
		{   /* param[d][0] = wavelength, params[d][1] = phase */
		    dep *= sin(2*PI*( (tcrds[d] - deparam[i].orig[d])
			/ deparam[i].param[d][0] + deparam[i].param[d][1]/360 ));
		}
	    }

	    dep += deparam[i].base;
	    switch(deparam[i].type)
	    {
		case EXT_ENERGY: /* total energy density */
		    Energy(state) = dep;
		    break;
		case EXT_PRESSURE:
		    set_state(state,TGAS_STATE,state);
		    Press(state) = dep;
		    set_state(state,GAS_STATE,state);
		    break;
		case EXT_TEMPER: /* keep density fixed */
		    set_state(state,FGAS_STATE,state);
		    Temperature(state) = dep;
		    set_state(state,GAS_STATE,state);
		    break;
		case EXT_VX:
		    set_state(state,TGAS_STATE,state);
		    Vel(state)[0] = dep;
		    set_state(state,GAS_STATE,state);
		    break;
		case EXT_VY:
		    set_state(state,TGAS_STATE,state);
		    Vel(state)[1] = dep;
		    set_state(state,GAS_STATE,state);
		    break;
		case EXT_VZ:
		    set_state(state,TGAS_STATE,state);
		    Vel(state)[2] = dep;
		    set_state(state,GAS_STATE,state);
		    break;
		case EXT_DENSITY:
		    Dens(state) = dep;
		    break;
	    }
	    if (debugging("deposition"))
	    {
		printf("{%c}dep[",dep_name[deparam[i].type]);
		for (d = 0;d < dim-1;d++)  printf("%g,",coords[d]);
		printf("%g]=%g, ",coords[dim-1],dep);
	    }
	}
}	/* end add_ext_deposition() */

/* temporary for sinusoidal incident wave */
/* Returns YES if the state is determined here, otherwise
 * return NO for default reflection boundary state. */
/*
EXPORT	int	boundary_deposition(
	Locstate	state,
	INTERFACE	*intfc,
	double		*pt,
	double		*p,
	double		*nor)
{
	double	dep, crds[MAXD], c;
	int	i, d;
	int	dim = intfc->dim;
	boolean	closed_end = TRUE;
*/
	/* use Flow_Through for outer boundaries in pellet ablation */
/*
	if (debugging("dep2ft"))
	{
	    if(nor[0] == 1 || nor[1] == 1)
		return NO;
	    else
		return YES;
	}
*/
	/* deposition only at two ends along y-axis */
/*
	if (nor[1] != -1 && nor[1] != 1) return NO;
	if (nor[1] == 1 && closed_end) return NO;
*/
	/* flow_through at ymin */
/*
	if (nor[1] == 1)  return YES;
*/
	/* set coords of deposition */
/*
	for(d=0;d<dim;d++) crds[d] = pt[d];
	crds[1] += 2*(p[1] - pt[1]);
	c = sound_speed(state);

	for (i = 0; i < dep_number; i++)
	{
	    dep = deparam[i].amp;
	    
	    for (d = 1; d <= dim; d++)
	    {
		if (deparam[i].prof[d] == SINUSOIDAL) 
		{   *//* param[d][0] = wavelength, params[d][1] = phase */
/*
		    dep *= sin(2*PI*( (crds[d-1] - deparam[i].orig[d] + c*cur_time)
			/ deparam[i].param[d][0] + deparam[i].param[d][1]/360 ));
		}
	    }

	    dep += deparam[i].base;
	    switch(deparam[i].type)
	    {
		case EXT_PRESSURE:
		    set_state(state,TGAS_STATE,state);
		    Press(state) = dep;
*/
		    /*Press(state) = 2*dep - Press(state);  forcing */
/*
		    set_state(state,GAS_STATE,state);
		    break;
		case EXT_VY:
		    set_state(state,TGAS_STATE,state);
		    Vel(state)[1] = dep;
		    set_state(state,GAS_STATE,state);
		    break;
	    }
	    if (debugging("bdry_dep") && crds[0] < 0 && crds[2] < 0)
	    {
		printf("t=%g, c=%g, bdry_dep[",cur_time,c);
		for (d = 0;d < dim-1;d++)  printf("%g,",crds[d]);
		printf("%g] = %g(type=%d)\n",crds[dim-1],dep,deparam[i].type);
	    }
	}
	return YES;
}
*/
