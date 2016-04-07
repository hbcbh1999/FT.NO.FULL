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
*				gprint.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains printing routines for gas dynamics that are specific
*	to dynamic simulations.
*
*	g_show_front_states() is set to a function pointer in the
*	Front structure.
*
*/



#include <gdecs/gdecs.h>

        /* LOCAL Function Declarations */
LOCAL	void	g_fprint_combustion_params(FILE*,Gas_param*);
LOCAL	void	fprint_viscosity_params(FILE*,AVISC*,boolean);
LOCAL	void	g_FPrintUSWSSten2d(FILE*,USWSSten2d*);
LOCAL	void	print_max_min_val(FILE*,const char*,double,double*,double,double*,
                                  double,double*,double,double*,int);
LOCAL	void	print_probe_field(double,Locstate,Locstate,Locstate,double,
                                  FILE*,boolean,int);
LOCAL	void	print_state_variable(int,int,Locstate);
LOCAL  	void 	print_comp_params(const char*,int);

#if defined(ONED)
LOCAL  void	g1d_accumulate_conserved_vars(FILE*,Front*,Wave*);
#endif /*defined(ONED)*/

 
/*test functions */
#include "gprtst.c"

EXPORT	void test_print_Tan_stencil(
	Tan_stencil	*sten)
{
	int	   i;
	int	   nsts = sten->npts/2+1;
	Locstate   st;

	printf("#npts = %d\n", sten->npts);
	
	for(i = -nsts+1; i < nsts; i++)
	{
	    if(i == 0)
	        continue;
	    printf("%d  ", sten->hs[i]);
	}
	printf("\n");

	for(i=-nsts+1; i<nsts; i++)
	{
	    if(i == 0)
	        continue;
	    printf("%d  ", Params(sten->leftst[i]));
	}
	printf("\n");
	for(i=-nsts+1; i<nsts; i++)
	{
	    if(i == 0)
	        continue;
	    printf("%d  ", Params(sten->rightst[i]));
	}
	printf("\n");
}

EXPORT	void	g_FPrintWSSten(
	FILE	*file,
	WSSten	*sten)
{
	INTERFACE	*intfc = current_interface();
	double   v[3];
	char	s[80];
	int	i, dim, nsts;

	(void) fprintf(file,"WSSten 0x%p\n",(POINTER)sten);
	if (sten == NULL)
	    return;

	dim = intfc->dim;

	nsts = sten->nsts;
	(void) fprintf(file,"sten->nsts = %d\n",nsts);
	(void) fprintf(file,"sten->ncomp = %d, sten->pcomp = %d\n",
	                    sten->ncomp,sten->pcomp);
	(void) printf("w_type = %d, ",sten->w_type);
	fprint_wave_type(file,"",sten->w_type,"\n",intfc);
	(void) fprintf(file,"sten->hs = %llu\n",hypersurface_number(sten->hs));
	fprint_general_vector(file,"sten->nor = ",sten->nor,dim,"\n");
	(void) fprintf(file,"sten->dn = %g\n",sten->dn);
	(void) fprintf(file,"sten->dt = %g\n",sten->dt);
	(void) fprintf(file,"sten->pjump = %g\n",sten->pjump);
	(void) fprintf(file,"sten->front = 0x%p\n",(POINTER)sten->front);
	(void) fprintf(file,"sten->wave = 0x%p\n",(POINTER)sten->wave);
	(void) fprintf(file,"_ws_interpolate = 0x%p\n",sten->_ws_interpolate);
	(void) fprintf(file,"_set_ws_slopes = 0x%p\n",sten->_set_ws_slopes);
	(void) fprintf(file,"_ClearWSStenData = 0x%p\n",sten->_ClearWSStenData);
	(void) fprintf(file,"_FPrintWSSten = 0x%p\n",sten->_FPrintWSSten);
	fprint_general_vector(file,"sten->coords = ",sten->coords,dim,"\n");
	for (i = nsts-1; i >= 0; --i)
	{
	    (void) sprintf(s,"sten->lcrds[%d] = ",i);
	    fprint_general_vector(file,s,sten->lcrds[i],dim," ");
	    (void) VelocityVector(sten->sl[i],v);
	    (void) sprintf(s,"sten->sl[%d], normal velocity = %g",i,
			   scalar_product(v,sten->nor,dim));
	    verbose_fprint_state(file,s,sten->sl[i]);
	    (void) VelocityVector(sten->tsl[i],v);
	    (void) sprintf(s,"sten->tsl[%d], normal velocity = %g",i,
			   scalar_product(v,sten->nor,dim));
	    verbose_fprint_state(file,s,sten->tsl[i]);
	    (void) fprintf(file,"sten->dsl[%d]\n",i);
	    fprint_raw_gas_data(file,sten->dsl[i],dim);
	}
	for (i = 0; i < nsts; ++i)
	{
	    (void) sprintf(s,"sten->rcrds[%d] = ",i);
	    fprint_general_vector(file,s,sten->rcrds[i],dim," ");
	    (void) VelocityVector(sten->sr[i],v);
	    (void) sprintf(s,"sten->sr[%d], normal velocity = %g",i,
			   scalar_product(v,sten->nor,dim));
	    verbose_fprint_state(file,s,sten->sr[i]);
	    (void) VelocityVector(sten->tsr[i],v);
	    (void) sprintf(s,"sten->tsr[%d], normal velocity = %g",i,
			   scalar_product(v,sten->nor,dim));
	    verbose_fprint_state(file,s,sten->tsr[i]);
	    (void) fprintf(file,"sten->dsr[%d]\n",i);
	    fprint_raw_gas_data(file,sten->dsr[i],dim);
	}
}		/*end g_FPrintWSSten*/

EXPORT	void print_WSStenData(
	WSSten			*sten)
{
	Front	*fr = sten->front;
	char s[120], s1[120], s2[120];
	Locstate sl, sr;
	int	i, dim = fr->rect_grid->dim;
	double   v[3];

	sl = sten->sl[0];
	sr = sten->sr[0];
	print_point_propagate_data(sten->p,sten->hse,sten->hs,dim);
	print_general_vector("sten->p = ",Coords(sten->p),dim,"\n");
	print_general_vector("sten->nor = ",sten->nor,dim,"\n");
	printf("sten->pjump = %24.16e\n", sten->pjump);
	if(sten->V != NULL)
	    print_general_vector("sten->V = ",sten->V, dim,"\n");

	for (i = sten->nsts-1; i >= 0; i--)
	{
	    (void) sprintf(s1,"sten->sl[%d]",i);
	    sprint_general_vector(s2,"coords = ",sten->lcrds[i],dim,"");
	    (void) sprintf(s,"%s, %s",s1,s2);
	    if (debugging("bad_state") &&
		    is_bad_state(sten->sl[i],YES,"print_WSStenData"))
	    {
	        screen("ERROR in print_WSStenData(), %s is bad\n",s);
	        fprint_raw_gas_data(stdout,sten->sl[i],current_interface()->dim);
	        clean_up(ERROR);
	    }
	    else
	    {
	        (void) VelocityVector(sten->sl[i],v);
		(void) printf("Normal velocity = %g\n",
			      scalar_product(v,sten->nor,dim));
	        verbose_print_state(s,sten->sl[i]);
	    }
	}
	for (i = 0; i < sten->nsts; ++i)
	{
	    (void) sprintf(s1,"sten->sr[%d]",i);
	    sprint_general_vector(s2,"coords = ",sten->rcrds[i],dim,"");
	    (void) sprintf(s,"%s, %s",s1,s2);
	    if (debugging("bad_state") &&
		    is_bad_state(sten->sr[i],YES,"print_WSStenData"))
	    {
	        screen("ERROR in print_WSStenData(), %s is bad\n",s);
	        fprint_raw_gas_data(stdout,sten->sr[i],current_interface()->dim);
	        clean_up(ERROR);
	    }
	    else
	    {
	        (void) VelocityVector(sten->sr[i],v);
		(void) printf("Normal velocity = %g\n",
			      scalar_product(v,sten->nor,dim));
 	        verbose_print_state(s,sten->sr[i]);
	    }
	}
	for (i = sten->nsts-1; i >= 0; i--)
	{
	    if (Different_params(sl,sten->sl[i])) 
	    {
	    	screen("ERROR in print_WSStenData(), "
		       "params sten->sl[%d] != params sl",i);
	    	clean_up(ERROR);
	    }
	}
	for (i = 0; i < sten->nsts; ++i)
	{
	    if (Different_params(sr,sten->sr[i])) 
	    {
	        screen("ERROR in print_WSStenData(), "
		       "params sten->sr[%d] != params sr",i);
	    	clean_up(ERROR);
	    }
	}
}		/*end print_WSStenData*/


EXPORT	void	g_PrintUSWSSten2d(
	USWSSten2d	*sten)
{
	g_FPrintUSWSSten2d(stdout,sten);
}		/*end g_PrintUSWSSten2d*/

LOCAL	void	g_FPrintUSWSSten2d(
	FILE		*file,
	USWSSten2d	*sten)
{
	int	i, dim;
	int	nor_rad, tan_rad;

	(void) fprintf(file,"USWSSten2d structure %p\n",(POINTER)sten);

	if (sten == NULL)
	{
	    (void) fprintf(file,"\tNULL USWSSten2d structure\n");
	    (void) fprintf(file,"End USWSSten2d structure %p\n",(POINTER)sten);
	    return;
	}

	if (sten->fr == NULL)
	{
	    (void) fprintf(file,"\tsten->fr == NULL, can't print \n");
	    (void) fprintf(file,"End USWSSten2d structure %p\n",(POINTER)sten);
	    return;
	}

	dim = sten->fr->rect_grid->dim;

	(void) fprintf(file,"sten->fr = %p, sten->wave = %p\n",
	               (POINTER)sten->fr,(POINTER)sten->wave);
	tan_rad = sten->tan_rad;
	nor_rad = sten->nor_rad;
	(void) fprintf(file,"sten->nor_rad = %d, sten->tan_rad = %d\n",
	               sten->nor_rad,sten->tan_rad);
	for (i = -tan_rad; i <= tan_rad; ++i)
	{
	    char label[21];
	    (void) sprintf(label,"nor[%d] = ",i);
	    fprint_general_vector(file,label,nor_vec(sten)[i],dim,"\n");
	    (void) fprintf(file,"\n");
	}
	fprint_general_vector(file,"sten->tngt = ",sten->tngt,dim,"\n");
	(void) fprintf(file,"sten->ds = %g, sten->dn = %g, sten->dt = %g\n",
	               sten->ds,sten->dn,sten->dt);
	(void) fprintf(file,"sten->hs\n");
	fprint_hypersurface(file,sten->hs);

	for (i = -nor_rad+1; i < nor_rad; ++i)
	{
	    (void) fprintf(file,"tan_fr(sten)[%d]\n",i);
	    fprint_Tan_stencil(file,sten->fr,tan_fr(sten)[i]);
	}
	for (i = -tan_rad; i <= tan_rad; ++i)
	{
	    (void) fprintf(file,"nor_fr(sten)[%d]\n",i);
	    FPrintWSSten(file,nor_fr(sten)[i]);
	}

	(void) fprintf(file,"sten->nor_ans\n");
	fprint_Tan_stencil(file,sten->fr,sten->nor_ans);

	(void) fprintf(file,"sten->tan_ans\n");
	FPrintWSSten(file,sten->tan_ans);

	(void) fprintf(file,"End USWSSten2d structure %p\n",(POINTER)sten);
}		/*end g_FPrintUSWSSten2d*/

/*
*			g_show_wave_states():
*
*	Prints the interior states in the form of arrays suitable for
*	use with "remaplot".  May be used for wave->show_wave_states.
*	      (see init_physics).
*/

EXPORT void g_show_wave_states(
	Wave		*wave)
{
	HYPER_SURF	   *hs, *hslast;
	HYPER_SURF_ELEMENT *hse;
	INTERFACE	   *intfc = wave_tri_soln(wave)->intfc;
	Locstate	   sl, sr, state;
	POINT		   *p;
	RECT_GRID	   *gr = wave->rect_grid;
	const char	   *hsname;
	const char	   *varname[20];
	double		   *VL = gr->VL, *VU = gr->VU;
	int		   i, j, ix, iy, iz, icoords[3];
	int		   xmax, ymax, zmax;
	int		   xmin, ymin, zmin;
	int		   dim = gr->dim;
	int		   gmax[MAXD], *lbuf = gr->lbuf, *ubuf = gr->ubuf;
	static const char         *Vname[] = {"X","Y","Z"};
	static const char         *vname[3] = {"x","y","z"};
#if defined(COMBUSTION_CODE)
	int		   composition_type = g_composition_type();
#endif /* defined(COMBUSTION_CODE) */

	(void) output();
	(void) printf("WAVE RECTANGULAR STATES\n");

	i = 0;
	varname[i++] = "DENSITY";
	varname[i++] = "ENERGY";
	varname[i++] = "X-MOMENTUM";
	switch (dim)
	{
#if defined(ONED)
 	case 1:
	    hsname = "point";
	    gmax[0] = gr->gmax[0] + lbuf[0] + ubuf[0];
	    xmax = gr->gmax[0] + ubuf[0];	xmin = -lbuf[0];
	    ymax = 1;    ymin = 0;
	    zmax = 1;    zmin = 0;
	    break;
#endif /* defined(ONED) */
#if defined(TWOD)
	case 2:
	    hsname = "curve";
	    varname[i++] = "Y-MOMENTUM";
	    gmax[0] = gr->gmax[0] + lbuf[0] + ubuf[0];
	    gmax[1] = gr->gmax[1] + lbuf[1] + ubuf[1];
	    xmax = gr->gmax[0] + ubuf[0];    xmin = -lbuf[0];
	    ymax = gr->gmax[1] + ubuf[1];    ymin = -lbuf[1];
	    zmax = 1;    zmin = 0;
	    break;
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:
	    hsname = "surface";
	    varname[i++] = "Y-MOMENTUM";
	    varname[i++] = "Z-MOMENTUM";
	    gmax[0] = gr->gmax[0] + lbuf[0] + ubuf[0];
	    gmax[1] = gr->gmax[1] + lbuf[1] + ubuf[1];
	    gmax[2] = gr->gmax[2] + lbuf[2] + ubuf[2];
	    xmax = gr->gmax[0] + ubuf[0];    xmin = -lbuf[0];
	    ymax = gr->gmax[1] + ubuf[1];    ymin = -lbuf[1];
	    zmax = gr->gmax[2] + ubuf[2];    zmin = -lbuf[2];
	    break;
#endif /* defined(THREED) */
	}
#if defined(COMBUSTION_CODE)
	if (composition_type == PTFLAME)
	    varname[i++] = "REACTION_PROGRESS";
	else if (composition_type == ZND)
	    varname[i++] = "PRODUCT_DENSITY";
	varname[i++] = "DENSITY_1";
#endif /* defined(COMBUSTION_CODE) */

	for (i = 0; i < wave->nfloats; ++i) 
	{
	    (void) printf("\n\n");
	    (void) output();
	    (void) printf("\t\t\t%s\n",varname[i]);
	    for (j = 0; j < dim; ++j)
	        (void) printf("%sL = %-6g ",Vname[j],VL[j]);
	    for (j = 0; j < dim; ++j)
	        (void) printf("%sU = %-6g ",Vname[j],VU[j]);
	    for (j = 0; j < dim; ++j)
	        (void) printf("%smax = %-6d%s",vname[j],gmax[j],
	        	      (j == dim-1) ? "\n" : " ");
	    (void) printf("linear singular\n\n");
	    for (iz = zmin; iz < zmax; ++iz)
	    {
	        icoords[2] = iz;
	        for (iy = ymax - 1; iy >= ymin; --iy) 
	        {
	            icoords[1] = iy;
	            for (ix = xmin; ix < xmax; ++ix) 
	            {
	                icoords[0] = ix;
	                print_state_variable(i,dim,Rect_state(icoords,wave));
	            }
	            (void) printf("\n");
	        }
	        (void) printf("\n\n");
	    }
	    (void) printf("\n\n");
	    (void) printf("\nStart of intfc states\n");
/* TEMP add */
            if(dim == 2)
            {
                if(intfc == NULL)
                {
                   (void) printf("\nEnd of intfc states, intfc is NULL\n");
                   continue;
                }
                if(intfc->curves == NULL)
                {
                   (void) printf("\nEnd of intfc states, intfc is NULL\n");
                   continue;
                }
            }
/* End of TEMP add */

	    hslast = NULL;
	    (void) next_point(intfc,NULL,NULL,NULL);
	    while (next_point(intfc,&p,&hse,&hs))
	    {
	        if (hs != hslast)
	        {
	            hslast = hs;
	            (void) printf("%s %llu\n",hsname,hypersurface_number(hs));
	        }
	        slsr(p,hse,hs,&sl,&sr);
	        for (j = 0; j < dim; ++j)
	            (void) printf("%-10g ",Coords(p)[j]);
	        print_state_variable(i,dim,sl);
	        print_state_variable(i,dim,sr);
	        (void) printf("\n");
	    }
	    (void) printf("End of intfc states\n\n");
	}
	(void) printf("\n\n");
	(void) output();
	(void) printf("\t\t\tGAS_PARAMS\n");
	for (j = 0; j < dim; ++j)
	    (void) printf("%sL = %-6g ",Vname[j],VL[j]);
	for (j = 0; j < dim; ++j)
	    (void) printf("%sU = %-6g ",Vname[j],VU[j]);
	for (j = 0; j < dim; ++j)
	    (void) printf("%smax = %-6d%s",vname[j],gmax[j],
	            (j == dim-1) ? "\n" : " ");
	(void) printf("linear singular\n\n");
	for (iz = zmin; iz < zmax; ++iz)
	{
	    icoords[2] = iz;
	    for (iy = ymax - 1; iy >= ymin; --iy) 
	    {
	        icoords[1] = iy;
	        for (ix = xmin; ix < xmax; ++ix) 
	        {
	            icoords[0] = ix;
	            state = Rect_state(icoords,wave);
	            (void) printf("%-10llu ",gas_param_number(Params(state)));
	        }
	        (void) printf("\n");
	    }
	    (void) printf("\n");
	}
	(void) printf("\n\n");
	(void) printf("\nStart of intfc states\n");
/* TEMP add */
            if(dim == 2)
            {
                if(intfc == NULL)
                {
                   (void) printf("\nEnd of intfc states, intfc is NULL\n");
                   goto end_this_function;
                }
                if(intfc->curves == NULL)
                {
                   (void) printf("\nEnd of intfc states, intfc is NULL\n");
                   goto end_this_function;
                }
            }
/* End of TEMP add */

	hslast = NULL;
	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (hs != hslast)
	    {
	        hslast = hs;
	        (void) printf("%s %llu\n",hsname,hypersurface_number(hs));
	    }
	    slsr(p,hse,hs,&sl,&sr);
	    for (j = 0; j < dim; ++j)
	        (void) printf("%-10g ",Coords(p)[j]);
	    (void) printf("%-10llu %-10llu\n",
	           gas_param_number(Params(sl)),
	           gas_param_number(Params(sr)));
	}
	(void) printf("End of intfc states\n\n");
/* TEMP add */
end_this_function:
/* End of TEMP add */
	(void) printf("\n\n");
	(void) output();
	(void) printf("END WAVE RECTANGULAR STATES\n");
}		/*end g_show_wave_states*/


LOCAL	void	print_state_variable(
	int		i,
	int		dim,
	Locstate	state)
{
	double		variable;

	if (i == 0)
	    variable = Dens(state);
	else if (i == 1)
	    variable = Energy(state);
	else if (i-2 < dim && i > 1)
	    variable = Mom(state)[i-2];
	else
	{
	    switch(i-2-dim)
	    {
#if defined(COMBUSTION_CODE)
	    case 0:
	    	if ((Composition_type(state) == PTFLAME) ||
		    (Composition_type(state) == THINFLAME))
	    	{
	    		variable = (Burned(state)) ? 1.0 : 0.0;
	    	}
	    	else
	    		variable = Prod(state);
	    	break;
	     case 1:
	    	variable = Dens1(state);
	    	break;
#endif /* defined(COMBUSTION_CODE) */
	     default:
	    	screen("ERROR in g_print_state_variable(), unknown variable\n");
	    	clean_up(ERROR);
	     }
	}
	(void) printf("%-10g ",variable);
}		/*end print_state_variable*/

/*
*		g_print_Gas_param_list():
*
*/

EXPORT void g_fprint_Gas_param_list(
	FILE		*file,
	INTERFACE	*intfc)
{
	COMPONENT	comp;
	Gas_param	*params;
	Gas_param	**prms_list;
	int		i, nprms;

	nprms = return_params_list(&prms_list);
	(void) foutput(file);
	(void) fprintf(file,"\tEquation of State Params List\n");
	(void) fprintf(file,"Number of params = %d\n\n",nprms);
	for (i = 0; i < nprms; ++i)
	{
	    fprint_Gas_param(file,prms_list[i]);
	    (void) fprintf(file,"\n");
	}
	(void) foutput(file);
	(void) fprintf(file,"\tEnd Equation of State Params List\n\n");

	(void) foutput(file);
	(void) fprintf(file,"\tComponents and Gas_param of Interface %llu\n",
	               interface_number(intfc));
	(void) fprintf(file,"\n");
       
        if(intfc->dim == 3)
	{
	    /*#bjet2 */
            (void) fprintf(file,"min_component = %d   max_component = %d\n",
                min_component(intfc), max_component(intfc));
            (void) fprintf(file,"\n");
	}

	for (comp = min_component(intfc); comp <= max_component(intfc); ++comp)
	{
	    params = gas_params_for_comp(comp,intfc);
	    (void) fprintf(file,"Component = %d Gas_param = %llu",comp,
	                   gas_param_number(params));
#if defined(COMBUSTION_CODE)
	    if (params && params->composition_type != PURE_NON_REACTIVE)
	    {
	    	(void) fprintf(file," other_params = %llu",
	    		       gas_param_number(params->other_params));
	    }
#endif /* defined(COMBUSTION_CODE) */
	    (void) fprintf(file,"\n\n");
	}
	(void) foutput(file);
	(void) fprintf(file,
	               "\tEnd Components and Gas_param of Interface %llu\n",
	               interface_number(intfc));
	(void) fprintf(file,"\n\n");
}		/*end g_print_Gas_param_list*/


/*
*			g_fprint_Gas_param():
*/

EXPORT void g_fprint_Gas_param(
	FILE		*file,
	Gas_param	*params)
{
	if (params == NULL)
	{
	    (void) fprintf(file,"Gas_param = 0%s\n",avisc_print_style());
	    (void) fprintf(file,"\tEquation of state = %d",OBSTACLE_EOS);
	    (void) fprintf(file," OBSTACLE\n\n");
	    return;
	}
	(void) fprintf(file,"%s = %llu%s\n",params->prt_label,
	               gas_param_number(params),
		       avisc_print_style());
	fprint_EOS_params(file,params);
	fprint_avisc_structure(file,&params->avisc,NO);
	g_fprint_combustion_params(file,params);
	fprint_thermodynamic_restrictions(file,params);
	(void) fprintf(file,"\n");
}		/*end fprint_Gas_param*/

EXPORT	void	fprint_thermodynamic_restrictions(
	FILE		*file,
	Gas_param	*params)
{
	(void) fprintf(file,"\tmin_energy = ");
	fprint_vector_of_floats(file,1,&params->min_energy);
	(void) fprintf(file,"\tmin_pressure = ");
	fprint_vector_of_floats(file,1,&params->min_pressure);
	(void) fprintf(file,"\tvacuum_dens = ");
	fprint_vector_of_floats(file,1,&params->vacuum_dens);
	(void) fprintf(file,"\traref_press = ");
	fprint_vector_of_floats(file,1,&params->raref_press);
#if defined(COMBUSTION_CODE)
	if (params->composition_type != PURE_NON_REACTIVE)
	{
	    (void) fprintf(file,"\ttol_alpha = ");
	    fprint_vector_of_floats(file,1,&params->tol_alpha);
	    (void) fprintf(file,"\ttol_press = ");
	    fprint_vector_of_floats(file,1,&params->tol_press);
	}
#endif /* defined(COMBUSTION_CODE) */
}		/*end fprint_thermodynamic_restrictions*/

/*ARGSUSED*/
LOCAL	void g_fprint_combustion_params(
	FILE		*file,
	Gas_param	*params)
{
#if !defined(COMBUSTION_CODE)
	(void) fprintf(file,"\tcomposition_type = %d %s\n",
	               PURE_NON_REACTIVE,"PURE_NON_REACTIVE");
#else /* !defined(COMBUSTION_CODE) */
	(void) fprintf(file,"\tcomposition_type = %d",params->composition_type);
	switch(params->composition_type)
	{
	case PURE_NON_REACTIVE:
	    (void) fprintf(file," PURE_NON_REACTIVE\n");
	    return;
	case PTFLAME:
	    (void) fprintf(file," PTFLAME\n");
	    break;
	case ZND:
	    (void) fprintf(file," ZND\n");
	    break;
	case TWO_CONSTITUENT_REACTIVE:
	    (void) fprintf(file," TWO_CONSTITUENT_REACTIVE\n");
	    break;
	case THINFLAME:
            (void) fprintf(file," THINFLAME\n");
            break;
	}
	(void) fprintf(file,"\tcritical_temperature = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void*) &params->critical_temperature,
	    	          FLOAT,1,file);
	    (void) fprintf(file,"\n");
	}
	else
	    (void) fprintf(file,"%"FFMT"\n",params->critical_temperature);
	(void) fprintf(file,"\tburned = %d q = ",params->burned);
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void*) &params->q,FLOAT,1,file);
	    (void) fprintf(file,"\n");
	}
	else
	    (void) fprintf(file,"%"FFMT"\n",params->q);
	(void) fprintf(file,"\trate_mult = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void*) &params->rate_mult,FLOAT,1,file);
	}
	else
	    (void) fprintf(file,"%"FFMT" ",params->rate_mult);
	(void) fprintf(file,"activ_en = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void*) &params->activ_en,FLOAT,1,file);
	    (void) fprintf(file,"\n");
	}
	else
	    (void) fprintf(file,"%"FFMT"\n",params->activ_en);

	if(params->composition_type == THINFLAME)
	{
	    int i;

	    (void) fprintf(file,"\tflame_velocity = ");
	    (void) fprintf(file,"%c %d ", params->code, params->ncoeffs);
	    if(is_binary_output() == YES)
	    {
		/*(void) fprintf(file,"\f%c",1); */
		/*(void) fwrite((const void*) &params->code,CHAR,1,file); */
		/*(void) fprintf(file,"\f%c",1); */
		/*(void) fwrite((const void*) &params->ncoeffs,INT,1,file); */
		for(i=0; i<params->ncoeffs; i++)
		{
		    (void) fprintf(file,"\f%c",1);
		    (void) fwrite((const void*) &params->flame_coeff[i],
				  FLOAT,1,file);
		}
	    }
	    else
	    {
		for(i=0; i<params->ncoeffs; i++)
		    (void) fprintf(file,"%"FFMT" ",params->flame_coeff[i]);
	    }
	    (void) fprintf(file,"\n");
	}

	    /* What about other params? */
#endif /* !defined(COMBUSTION_CODE) */
}		/*end g_fprint_combustion_params*/

EXPORT	Gas_param *gas_params_for_comp(
	COMPONENT	comp,
	INTERFACE	*intfc)
{
	Gas_param	*params = NULL;
	HYPER_SURF	**hslist;
	HYPER_SURF	*hs, **hss;
	int		dim = intfc->dim;
	int		comp_on_intfc = NO;
	SIDE		side;

	if (is_exterior_comp(comp,intfc))
	    return params;
	hslist = hyper_surf_list(intfc);
	hs = NULL;
	for (hss = hslist; hss && *hss; ++hss)
	{
	    /*printf("%d  compin=%d comp = %d  %d\n", *hss, comp, negative_component(*hss), positive_component(*hss)); */

	    if (negative_component(*hss) == comp) 
	    {
	    	comp_on_intfc = YES;
	    	if (wave_type(*hss) != SUBDOMAIN_BOUNDARY)
	    	{
	    	    side = NEGATIVE_SIDE;
	    	    hs = *hss;
	    	    break;
	    	}
	    }
	    if (positive_component(*hss) == comp) 
	    {
	    	comp_on_intfc = YES;
	    	if (wave_type(*hss) != SUBDOMAIN_BOUNDARY)
	    	{
	    	    side = POSITIVE_SIDE;
	    	    hs = *hss;
	    	    break;
	    	}
	    }
	}
	if (hs == NULL)
	{
	    if (comp_on_intfc == YES)
	    {
	    	/* No tracked waves except subdomain boundaries */
	    	CHART		*chart;
	    	int		i;
	    	static boolean	first = YES;
	    	static		int icoords[MAXD];

	    	chart = current_chart();
	    	if (	(chart == NULL) ||
	    		(chart->front == NULL) ||
	    		(chart->wave == NULL) ||
	    		(wave_tri_soln(chart->wave) == NULL) ||
	    		(wave_tri_soln(chart->wave)->tri_grid == NULL)
	    	)
	    	    return params;

	        if (first == YES)
	        {
	            first = NO;
	            for (i = 0; i < dim; ++i)
	                icoords[i] = 0;
	        }
	        params = Params(Rect_state(icoords,chart->wave));
	    }
	    else
	    {/*#bjet2 */
	        if(dim == 3)
		    return param_for_comp(comp);
	    }
	    return params;
	}

	switch (dim)
	{
#if defined(ONED)
	case 1:
	    if (side == NEGATIVE_SIDE)
	        params = Params(left_state(Point_of_hs(hs)));
	    else if (side == POSITIVE_SIDE)
	    	params = Params(right_state(Point_of_hs(hs)));
	    else
	    {
		params = NULL;
		screen("ERROR in gas_params_for_comp(), invalid side\n");
		clean_up(ERROR);
	    }
	    break;
#endif /* defined(ONED) */
#if defined(TWOD)
	case 2:
	    if (side == NEGATIVE_SIDE)
		params = Params(left_start_state(Curve_of_hs(hs)));
	    else if (side == POSITIVE_SIDE)
	        params = Params(right_start_state(Curve_of_hs(hs)));
	    else
	    {
		params = NULL;
		screen("ERROR in gas_params_for_comp(), invalid side\n");
		clean_up(ERROR);
	    }
	    break;
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:
	    {
	        SURFACE		*s = Surface_of_hs(hs);
	        TRI 		*tri = first_tri(s);
	        POINT		*p = Point_of_tri(tri)[0];
	        Locstate	sl, sr;

	        slsr(p,Hyper_surf_element(tri),hs,&sl,&sr);
	        
		if(debugging("tst_params"))
		{
		    printf("#paramcomp\n");
		    print_tri(tri, intfc);
		    printf("#surf %d  %d  %d\n", s, negative_component(s), positive_component(s));
		    printf("#params  %p  %p\n", Params(sl), Params(sr));
		}
		
		if (side == NEGATIVE_SIDE)
		    params = Params(sl);
	        else if (side == POSITIVE_SIDE)
		    params = Params(sr);
	        else
	        {
		    params = NULL;
		    screen("ERROR in gas_params_for_comp(), invalid side\n");
		    clean_up(ERROR);
	        }
	    }
	    break;
#endif /* defined(THREED) */
	}
	return params;
}		/*end gas_params_for_comp*/


/*
*		g_fprint_Dirichlet_bdry_states():
*/

EXPORT	void g_fprint_Dirichlet_bdry_states(
	FILE		*file,
	INTERFACE	*intfc)
{
	HYPER_SURF **hs;
	int	   i, dim = intfc->dim;
	static const char *hsname[] = { "point", "curve", "surface" };

	(void) fprintf(file,"\n\n");
	(void) foutput(file);
	(void) fprintf(file,
	               "Hypersurface Dirichlet boundary state information ");
	(void) fprintf(file,"for interface %llu\n",interface_number(intfc));
	for (i = 0, hs = intfc->hss; hs && *hs; ++i, ++hs)
	{
	    if (wave_type(*hs) != DIRICHLET_BOUNDARY)
	        continue;
	    (void) fprintf(file,"Boundary state index for %s %d = %d\n",
	                   hsname[dim-1],i,bstate_index(*hs));
	}
	(void) foutput(file);
	(void) fprintf(file,
	             "End hypersurface Dirichlet boundary state information ");
	(void) fprintf(file,"for interface %llu\n",interface_number(intfc));
	(void) fprintf(file,"\n\n");
}		/*end g_fprint_Dirichlet_bdry_states*/


EXPORT	void	g_fprint_random_velocity_inlet_data(
	FILE		*file,
	INTERFACE	*intfc,
	BOUNDARY_STATE	*bstate)
{
	RANDOM_STATE	*rstate = (RANDOM_STATE*) bstate->_boundary_state_data;
	RECT_GRID	*grid = &rstate->grid;
	int		gmax[3];
	boolean		bin = is_binary_output();
	int		i, j, k;

	f_fprint_boundary_state_data(file,intfc,bstate);

	(void) fprintf(file,"\nPrintout of RANDOM_STATE structure\n\n");

	fprint_state_type(file,"State type = ",state_type(Mean(rstate)));
	(void) fprintf(file,"Mean state = ");
	fprint_state_data(file,Mean(rstate),intfc);
	(void) fprintf(file,"Sigma state = ");
	fprint_state_data(file,Sigma(rstate),intfc);

	(void) fprintf(file,"Q = ");
	if (bin == YES)
	{
	    (void) fprintf(file,"\f%c",9);
	    for (i = 0; i < 3; ++i)
	    	(void) fwrite((const void*)rstate->Q[i],FLOAT,3,file);
	    (void) fprintf(file,"\n");
	}
	else
	{
	    for (i = 0; i < 3; ++i)
	    {
	    	for (j = 0; j < 3; ++j)
	    	    (void) fprintf(file,"\t%"FFMT" ",rstate->Q[i][j]);
	    	(void) fprintf(file,"\n");
	    }
	}

	(void) fprintf(file,"lambda = ");
	if (bin == YES)
	{
	    (void) fprintf(file,"\f%c",3);
	    (void) fwrite((const void*) rstate->lambda,FLOAT,3,file);
	    (void) fprintf(file,"\n");
	}
	else
	{
	    for (i = 0; i < 3; ++i)
	    	(void) fprintf(file,"\t%"FFMT" ",rstate->lambda[i]);
	    (void) fprintf(file,"\n");
	}

	(void) fprintf(file,"A = ");
	if (bin == YES)
	{
	    (void) fprintf(file,"\f%c",9);
	    for (i = 0; i < 3; ++i)
	    	(void) fwrite((const void*) rstate->A[i],FLOAT,3,file);
	    (void) fprintf(file,"\n");
	}
	else
	{
	    for (i = 0; i < 3; ++i)
	    {
	    	for (j = 0; j < 3; ++j)
	    	    (void) fprintf(file,"\t%"FFMT" ",rstate->A[i][j]);
	        (void) fprintf(file,"\n");
	    }
	}

	(void) fprintf(file,"Number of evaluations per tau or lambda = %d\n",
	               rstate->N);

	(void) fprintf(file,"correlation time = ");
	if (bin == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void*) &rstate->tau,FLOAT,1,file);
	    (void) fprintf(file,"\n");
	}
	else
	    (void) fprintf(file,"%"FFMT"\n",rstate->tau);

	(void) fprintf(file,"Evaluation grid =\n");
	fprint_rectangular_grid(file,&rstate->grid);

	(void) fprintf(file,"time_of_save = ");
	if (bin == YES)
	{
	   (void) fprintf(file,"\f%c",1);
	   (void) fwrite((const void*)&rstate->time_of_save,FLOAT,1,file);
	   (void) fprintf(file,"\n");
	}
	else
	   (void) fprintf(file,"%"FFMT"\n",rstate->time_of_save);

	(void) fprintf(file,"Random number seeds after save = %hu %hu %hu\n",
	               rstate->seed_after_save[0],
	               rstate->seed_after_save[1],
	               rstate->seed_after_save[2]);

	for (i = 0; i < 3; ++i)
	    gmax[i] = grid->gmax[i]+1;

	(void) fprintf(file,
	              "Printout of %d X %d X %d grid of saved random states\n",
	               gmax[0],gmax[1],gmax[2]);

	for (i = 0; i < gmax[0]; ++i)
	for (j = 0; j < gmax[1]; ++j)
	for (k = 0; k < gmax[2]; ++k)
	    fprint_state_data(file,rstate->save_st[i][j][k],intfc);

	(void) fprintf(file,
	        "End Printout of %d X %d X %d grid of saved random states\n",
	        gmax[0],gmax[1],gmax[2]);

	if (axisymmetric_random_region_about_origin(rstate,intfc) == YES)
	{
	    (void) fprintf(file,"Radial velocity decay exponent = ");
	    if (bin == YES)
	    {
	        (void) fprintf(file,"\f%c",1);
	    	(void) fwrite((const void*)&RadialVelocityDecayExponent(rstate),
	                      FLOAT,1,file);
	        (void) fprintf(file,"\n");
	    }
	    else
	    {
	    	(void) fprintf(file,"%"FFMT"\n",
			       RadialVelocityDecayExponent(rstate));
	    }
	    (void) fprintf(file,"Radial velocity decay length scale = ");
	    if (bin == YES)
	    {
	        (void) fprintf(file,"\f%c",1);
	    	(void) fwrite((const void*)&RadialVelocityDecayScale(rstate),
	                      FLOAT,1,file);
	        (void) fprintf(file,"\n");
	    }
	    else
	    {
	    	(void) fprintf(file,"%"FFMT"\n",
			       RadialVelocityDecayScale(rstate));
	    }
	}

	(void) fprintf(file,"\nEnd Printout of RANDOM_STATE structure\n");

}		/*end g_fprint_random_velocity_inlet_data*/

EXPORT	const char *avisc_print_style(void)
{
	return "    ";/*Number of spaces = print style*/
}		/*end avisc_print_style*/

EXPORT	void	fprint_avisc_structure(
	FILE		*file,
	AVISC		*avisc,
	boolean		interactive)
{
	if (avisc == NULL)
	    return;

	(void) fprintf(file,"\tArtificial Viscosities and Heat Conductions");
	fprint_viscosity_params(file,avisc,is_binary_output());
	if (interactive == YES)
	    fprint_viscosity_params(stderr,avisc,NO);
}		/*end fprint_avisc_structure*/

LOCAL	void fprint_viscosity_params(
	FILE  *file,
	AVISC *avisc,
	boolean  bio)
{
	char mesg[1024];

	(void) fprintf(file,"\n");

	(void) sprintf(mesg,"\t%s%s\n\t%s",Use_nlav,
	               y_or_n(use_lapidus_artificial_viscosity(*avisc)),
		       Coef_nlav);
	fwrite_float(file,mesg,avisc->lapidus_visc_coef,bio,"%"FFMT,"\n");

	(void) sprintf(mesg,"\t%s%s\n\t%s",Use_lav,
	               y_or_n(use_linear_artificial_viscosity(*avisc)),
		       Coef_lav);
	fwrite_float(file,mesg,avisc->linear_visc_coef,bio,"%"FFMT,"\n");

	(void) sprintf(mesg,"\t%s%s\n\t%s",Use_uwav,
	               y_or_n(use_upwind_artificial_viscosity(*avisc)),
		       Coef_uwav);
	fwrite_float(file,mesg,avisc->upwind_visc_coef,bio,"%"FFMT,"\n");

	(void) sprintf(mesg,"\t%s%s\n""\t%s",Use_msf,
	               y_or_n(use_muscl_slope_flattening(*avisc)),
		       Coef_msf_ieta);
	fwrite_float(file,mesg,avisc->msf_ieta,bio,"%"FFMT,"\n");

	(void) sprintf(mesg,"\t%s",Coef_msf_ms);
	fwrite_float(file,mesg,avisc->min_shock_jump,bio,"%"FFMT,"\n");

	(void) sprintf(mesg,"\t%s",Coef_msf_msvj);
	fwrite_float(file,mesg,avisc->min_sp_vol_jump,bio,"%"FFMT,"\n");

	(void) sprintf(mesg,"\t%s",Coef_hc);
	fwrite_float(file,mesg,avisc->heat_cond,bio,"%"FFMT,"\n");

	(void) sprintf(mesg,"\t%s",Coef_char_speed_cutoff);
	fwrite_float(file,mesg,avisc->char_speed_cutoff,bio,"%"FFMT,"\n");

	(void) sprintf(mesg,"\t%s",Coef_dst);
	fwrite_float(file,mesg,avisc->dynamic_st,bio,"%"FFMT,"\n");

	(void) sprintf(mesg,"\t%s",Coef_sp);
	fwrite_float(file,mesg,avisc->sp_coef,bio,"%"FFMT,"\n");

        (void) sprintf(mesg,"\t%s",Coef_contact_detector);
        fwrite_float(file,mesg,avisc->contact_detector,bio,"%"FFMT,"\n");
}		/*end fprint_viscosity_params*/


EXPORT	void	set_default_data_file_names(
	const char *subdir,
	const char *name,
	char	   **dname,
	char	   **fname,
	INIT_DATA  *init)
{
	char	*ofname = output_filename(init);
	char	buf[1024];
	static char directory[1024], filename[1024];
	boolean dir_in_fname;

	*dname = *fname = NULL;

	if ((ofname == NULL) || (ofname[0] == '\0'))
	    return;
	if ((name == NULL) || (name[0] == '\0'))
	    return;

	dir_in_fname = (strstr(ofname,"/") != NULL) ? YES : NO;
	if (dir_in_fname == YES)
	{
	    strcpy(directory,ofname);
	    (void)get_dirname(directory);
	    strcpy(buf,basename(ofname));

	    if (subdir != NULL)
	    {
	        (void) strcat(directory,"/");
	        (void) strcat(directory,subdir);
	    }
	    (void) sprintf(filename,"%s/%s%s",directory,buf,name);
	    *dname = directory;
	    *fname = filename;
	}
	else
	{
	    (void) sprintf(filename,"%s%s",ofname,name);
	    *dname = NULL;
	    *fname = filename;
	}
}		/*end set_default_data_file_names*/


/*
*			print_graph_header():
*
*	Prints out the header for a statistical data file.  The column header
*	has already been allocated and assembled elsewhere.  The option to
*	print a time header is used for statistics that print an entire table
*	of data at each output.  The time stamp is printed between tables.
*	This is not necessary if the output at each printing is a single line.
*	The time stamp should be included in the line itself.
*/

EXPORT	void print_graph_header(
	FILE		*file,
	const char	*var_name,
	const char	*column_header,
	boolean		print_time_stamp,
	const Front	*front)
{
	(void) foutput(file);
	(void) fprintf(file,"\t\tBEGIN %s DEPENDENT STATISTICAL DATA\n\n",
	               var_name);

	if (print_time_stamp == YES)
	    (void) fprintf(file,"\t\t\t\tt = %"FFMT"\tj = %d\n\n",
	    	           front->time,front->step);

	(void) foutput(file);
	(void) fprintf(file,"%s",column_header);

}		/*end print_graph_header*/


/*
*			print_graph_footer():
*
*	Prints a footer indicating the end of a table of statistical data.  See
*	notes regarding the time stamp in print_graph_header().
*/

EXPORT void print_graph_footer(
	FILE		*file,
	const char	*var_name,
	boolean		trace)
{
	(void) fprintf(file,"\n\n\n");
	(void) foutput(file);
	(void) fprintf(file,"\t\t\tEND %s DEPENDENT STATISTICAL DATA\n\n",
	               var_name);
	if (trace == YES) trace_foutput(file);
}		/*end print_graph_footer*/

EXPORT	const char *print_wave_family(
	FILE	    *file,
	const char  *mesg,
	WAVE_FAMILY l_or_r,
	const char  *end)
{
	const char *family_name;

	family_name = (l_or_r == LEFT_FAMILY) ? "LEFT_FAMILY" : "RIGHT_FAMILY";
	if (file != NULL)
	{
	    if (mesg != NULL)
	        (void) fprintf(file,"%s",mesg);
	    (void) fprintf(file,"%d %s",l_or_r,family_name);
	    if (end != NULL)
	        (void) fprintf(file,end);
	}
	return family_name;
}		/*end print_wave_family*/

/*ARGSUSED*/
EXPORT	void	g_print_extreme_values(
	FILE      *file,
	CHART     *chart,
	Printplot *prt)
{
	EXTREME_VALUES *front_vals, *wave_vals;
	Front          *front = chart->front;
	Wave           *wave = chart->wave;
	int            i, dim = front->rect_grid->dim;
	char           s[50];

	front_vals = &extreme_front_vals(front);
	wave_vals = &extreme_wave_vals(wave);

	(void) fprintf(file,"\n%-24s %-14s %-14s ","state_variable","maximum",
	               "coords_maximum");
	for (i = 1; i < dim; ++i)
	    (void) fprintf(file," %-16s "," ");
	(void) fprintf(file,"on_front %-14s %-14s ",
	              "minimum","coords_minimum");
	for (i = 1; i < dim; ++i)
	    (void) fprintf(file,"%-16s"," ");
	(void) fprintf(file,"  on_front\n");

	print_max_min_val(file,"density",
	                  front_vals->rho_max,front_vals->coords_rho_max,
	                  wave_vals->rho_max,wave_vals->coords_rho_max,
			  front_vals->rho_min,front_vals->coords_rho_min,
			  wave_vals->rho_min,wave_vals->coords_rho_min,dim);
	print_max_min_val(file,"pressure",
	                  front_vals->p_max,front_vals->coords_p_max,
	                  wave_vals->p_max,wave_vals->coords_p_max,
			  front_vals->p_min,front_vals->coords_p_min,
			  wave_vals->p_min,wave_vals->coords_p_min,dim);
	print_max_min_val(file,"specific_internal_energy",
	                  front_vals->e_max,front_vals->coords_e_max,
	                  wave_vals->e_max,wave_vals->coords_e_max,
			  front_vals->e_min,front_vals->coords_e_min,
			  wave_vals->e_min,wave_vals->coords_e_min,dim);
	print_max_min_val(file,"kinetic_energy",
	                  front_vals->ke_max,front_vals->coords_ke_max,
	                  wave_vals->ke_max,wave_vals->coords_ke_max,
			  front_vals->ke_min,front_vals->coords_ke_min,
			  wave_vals->ke_min,wave_vals->coords_ke_min,dim);
	print_max_min_val(file,"total_energy",
	                  front_vals->E_max,front_vals->coords_E_max,
	                  wave_vals->E_max,wave_vals->coords_E_max,
			  front_vals->E_min,front_vals->coords_E_min,
			  wave_vals->E_min,wave_vals->coords_E_min,dim);
	for (i = 0; i < dim; ++i)
	{
	    (void) sprintf(s,"momentum[%d]",i);
	    print_max_min_val(file,s,
	                      front_vals->m_max[i],front_vals->coords_m_max[i],
	                      wave_vals->m_max[i],wave_vals->coords_m_max[i],
			      front_vals->m_min[i],front_vals->coords_m_min[i],
			      wave_vals->m_min[i],wave_vals->coords_m_min[i],
			      dim);
	}
	for (i = 0; i < dim; ++i)
	{
	    (void) sprintf(s,"velocity[%d]",i);
	    print_max_min_val(file,s,
	                      front_vals->v_max[i],front_vals->coords_v_max[i],
	                      wave_vals->v_max[i],wave_vals->coords_v_max[i],
			      front_vals->v_min[i],front_vals->coords_v_min[i],
			      wave_vals->v_min[i],wave_vals->coords_v_min[i],
			      dim);
	}
	(void) fprintf(file,"\n");
#if defined(ONED)
        if (dim == 1)
            g1d_accumulate_conserved_vars(file,front,wave);
#endif /* defined(ONED) */
}		/*end g_print_extreme_values*/

LOCAL	void	print_max_min_val(
	FILE       *file,
	const char *vname,
	double      front_val_max,
	double      *front_val_max_coords,
	double      wave_val_max,
	double      *wave_val_max_coords,
	double      front_val_min,
	double      *front_val_min_coords,
	double      wave_val_min,
	double      *wave_val_min_coords,
	int        dim)
{
	int i;
	(void) fprintf(file,"%-24s ",vname);
	if (front_val_max > wave_val_max)
	{
	    (void) fprintf(file,"%-14g (%-14g",front_val_max,
	                                       front_val_max_coords[0]);
	    for (i = 1; i < dim; ++i)
	        (void) fprintf(file,", %-14g",front_val_max_coords[i]);
	    (void) fprintf(file,") %-8s ","yes");
	}
	else
	{
	    (void) fprintf(file,"%-14g (%-14g",wave_val_max,
	                                       wave_val_max_coords[0]);
	    for (i = 1; i < dim; ++i)
	        (void) fprintf(file,", %-14g",wave_val_max_coords[i]);
	    (void) fprintf(file,") %-8s ","no");
	}
	if (front_val_min < wave_val_min)
	{
	    (void) fprintf(file,"%-14g (%-14g",front_val_min,
	                                       front_val_min_coords[0]);
	    for (i = 1; i < dim; ++i)
	        (void) fprintf(file,", %-14g",front_val_min_coords[i]);
	    (void) fprintf(file,") %-8s ","yes");
	}
	else
	{
	    (void) fprintf(file,"%-14g (%-14g",wave_val_min,
	                                       wave_val_min_coords[0]);
	    for (i = 1; i < dim; ++i)
	        (void) fprintf(file,", %-14g",wave_val_min_coords[i]);
	    (void) fprintf(file,") %-8s ","no");
	}
	(void) fprintf(file,"\n");
}		/*end print_max_min_val*/

/*ARGSUSED*/
EXPORT	void	evaluate_probe(
	Grid        *grid,
	Wave        *wave,
	Front       *front,
	Printplot   *prt,
	OUTPUT_DATA *data,
	boolean        about_to_stop)
{
	COMPONENT comp;
	FILE      *file;
	INTERFACE *intfc = front->interf;
	Locstate  *st, dst, Ndst;
	PROBE     *probe = (PROBE*)data;
	RECT_GRID *gr = front->rect_grid;
	double     *coords;
	int       *n = probe->n;
	int       dim = front->rect_grid->dim;
	int       i, j, k, l;

	st = probe->st;
	dst = probe->dst;
	Ndst = probe->Ndst;
	zero_scalar(st[2],front->sizest);
	set_type_of_state(st[2],TGAS_STATE);
	for (i = -n[0]; i <= n[0]; ++i)
	for (j = -n[1]; j <= n[1]; ++j)
	for (k = -n[2]; k <= n[2]; ++k)
	{
	    Locstate  s;
	    coords = probe_coords(i,j,k,probe);
	    probe_comp(i,j,k,probe) = comp = component(coords,intfc);
	    s = probe_state(i,j,k,probe);
	    hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,front,wave,
			    		s,NULL);
	    Dens(st[2]) += Dens(s);
	    Press(st[2]) += pressure(s);
	    for (l = 0; l < dim; ++l)
	        Vel(st[2])[l] += vel(l,s);
	}
	Dens(st[2]) /= probe->N;
	Press(st[2]) /= probe->N;
	for (l = 0; l < dim; ++l)
	    Vel(st[2])[l] /= probe->N;
	Set_params(st[2],probe_state(0,0,0,probe));
	set_state(st[2],VGAS_STATE,st[2]);

	file = Output_file(data);
	if (file == NULL)
	{
	    file = Output_file(data) = fopen(Output_filename(data),"w");
	    print_machine_parameters(file);
	    (void) foutput(file);
	    (void) fprintf(file,"%-18s %-18s %-18s %-18s %-18s %-18s",
	                   "Time","Wave_tag",
			   "Density","Pressure","Entropy","Sound_speed");
	    for (l = 0; l < dim; ++l)
	    {
	        char vname[256];
	        (void) sprintf(vname,"%s-velocity",gr->Remap.dnm[l]);
	        (void) fprintf(file," %-18s",vname);
	    }
	    (void) fprintf(file," %-18s %-18s %-18s %-18s",
	                   "dDensity","dPressure","dEntropy","dSound_speed");
	    for (l = 0; l < dim; ++l)
	    {
	        char vname[256];
	        (void) sprintf(vname,"d%s-velocity",gr->Remap.dnm[l]);
	        (void) fprintf(file," %-18s",vname);
	    }
	    (void) fprintf(file," %-18s %-18s %-18s %-18s",
	                   "NdDensity","NdPressure","NdEntropy",
			   "NdSound_speed");
	    for (l = 0; l < dim; ++l)
	    {
	        char vname[256];
	        (void) sprintf(vname,"Nd%s-velocity",gr->Remap.dnm[l]);
	        (void) fprintf(file," %-18s",vname);
	    }
	    fprint_general_vector(file," probe at ",probe->c,dim,"\n");
	}
	probe->time[2] = grid->time;
	if (probe->line > 1)
	{
	    double Dt, dt[2];
	    double dvar[2];
	    dt[0] = probe->time[1] - probe->time[0];
	    dt[1] = probe->time[2] - probe->time[1];
	    Dt = probe->L/sound_speed(st[1]);
	    dvar[0] = (Dens(st[1])-Dens(st[0]))/dt[0];
	    dvar[1] = (Dens(st[2])-Dens(st[1]))/dt[1];
	    Dens(dst) = (dt[0]*dvar[1]+dt[1]*dvar[0])/(dt[0]+dt[1]);
	    Dens(Ndst) = Dt*Dens(dst)/Dens(st[1]);
	    dvar[0] = (Press(st[1])-Press(st[0]))/dt[0];
	    dvar[1] = (Press(st[2])-Press(st[1]))/dt[1];
	    Press(dst) = (dt[0]*dvar[1]+dt[1]*dvar[0])/(dt[0]+dt[1]);
	    Press(Ndst) = Dt*Press(dst)/Press(st[1]);
	    dvar[0] = (Entropy(st[1])-Entropy(st[0]))/dt[0];
	    dvar[1] = (Entropy(st[2])-Entropy(st[1]))/dt[1];
	    Entropy(dst) = (dt[0]*dvar[1]+dt[1]*dvar[0])/(dt[0]+dt[1]);
	    Entropy(Ndst) = Dt*Entropy(dst)*temperature(st[1])/
	                    sqr(Sound_speed(st[1]));
	    dvar[0] = (Sound_speed(st[1])-Sound_speed(st[0]))/dt[0];
	    dvar[1] = (Sound_speed(st[2])-Sound_speed(st[1]))/dt[1];
	    Sound_speed(dst) = (dt[0]*dvar[1]+dt[1]*dvar[0])/(dt[0]+dt[1]);
	    Sound_speed(Ndst) = Dt*Sound_speed(dst)/Sound_speed(st[1]);
            reset_gamma(dst);
            reset_gamma(Ndst);
	    if (probe->wave_tag == 0.0)
	    {
	        if (Press(Ndst) > probe->ptol)
		    probe->wave_tag = 1.0;
	        if (Press(Ndst) < -probe->ptol)
		    probe->wave_tag = -1.0;
	    }
	    else
	    {
	        if (fabs(Press(Ndst)) < probe->ptol)
		    probe->wave_tag = 0.0;
	    }

	    for (l = 0; l < dim; ++l)
	    {
	        dvar[0] = (Vel(st[1])[l]-Vel(st[0])[l])/dt[0];
	        dvar[1] = (Vel(st[2])[l]-Vel(st[1])[l])/dt[1];
		Vel(dst)[l] = (dt[0]*dvar[1]+dt[1]*dvar[0])/(dt[0]+dt[1]);
	        Vel(Ndst)[l] = Dt*Vel(dst)[l]/sound_speed(st[1]);
	    }
	    if (probe->line == 2)
	        print_probe_field(probe->time[0],st[0],dst,Ndst,
		                  probe->wave_tag,file,
				  Output_in_binary(data),dim);
	    print_probe_field(probe->time[1],st[1],dst,Ndst,
	                      probe->wave_tag,file,Output_in_binary(data),dim);
	}
	if (probe->line > 0)
	{
	    probe->time[0] = probe->time[1];
	    copy_state(st[0],st[1]);
	}
	probe->time[1] = probe->time[2];
	copy_state(st[1],st[2]);
	++probe->line;
	if (about_to_stop)
	{
	    print_probe_field(probe->time[2],st[2],dst,Ndst,
	                      probe->wave_tag,file,Output_in_binary(data),dim);
	}
}		/*end evaluate_probe*/

LOCAL	void	print_probe_field(
	double    time,
	Locstate st,
	Locstate dst,
	Locstate Ndst,
	double    wave_tag,
	FILE     *file,
	boolean     binary,
	int      dim)
{
	int l;
	if (binary)
	{
	    (void) fprintf(file,"\f%c",2+3*(4+dim));
	    (void) fwrite((const void *)&time,FLOAT,1,file);
	    (void) fwrite((const void *)&wave_tag,FLOAT,1,file);

	    (void) fwrite((const void *)&Dens(st),FLOAT,1,file);
	    (void) fwrite((const void *)&Press(st),FLOAT,1,file);
	    (void) fwrite((const void *)&Entropy(st),FLOAT,1,file);
	    (void) fwrite((const void *)&Sound_speed(st),FLOAT,1,file);
	    (void) fwrite((const void *)Vel(st),FLOAT,dim,file);

	    (void) fwrite((const void *)&Dens(dst),FLOAT,1,file);
	    (void) fwrite((const void *)&Press(dst),FLOAT,1,file);
	    (void) fwrite((const void *)&Entropy(dst),FLOAT,1,file);
	    (void) fwrite((const void *)&Sound_speed(dst),FLOAT,1,file);
	    (void) fwrite((const void *)Vel(dst),FLOAT,dim,file);

	    (void) fwrite((const void *)&Dens(Ndst),FLOAT,1,file);
	    (void) fwrite((const void *)&Press(Ndst),FLOAT,1,file);
	    (void) fwrite((const void *)&Entropy(Ndst),FLOAT,1,file);
	    (void) fwrite((const void *)&Sound_speed(Ndst),FLOAT,1,file);
	    (void) fwrite((const void *)Vel(Ndst),FLOAT,dim,file);

	}
	else
	{
	    (void) fprintf(file,"%-"FFMT" %-"FFMT,time,wave_tag);

	    (void) fprintf(file," %-"FFMT" %-"FFMT" %-"FFMT" %-"FFMT,
		           Dens(st),Press(st),Entropy(st),Sound_speed(st));
	    for (l = 0; l < dim; ++l)
	        (void) fprintf(file," %-"FFMT,Vel(st)[l]);

	    (void) fprintf(file," %-"FFMT" %-"FFMT" %-"FFMT" %-"FFMT,
		           Dens(dst),Press(dst),Entropy(dst),Sound_speed(dst));
	    for (l = 0; l < dim; ++l)
	        (void) fprintf(file," %-"FFMT,Vel(dst)[l]);

	    (void) fprintf(file," %-"FFMT" %-"FFMT" %-"FFMT" %-"FFMT,
		           Dens(Ndst),Press(Ndst),
			   Entropy(Ndst),Sound_speed(Ndst));
	    for (l = 0; l < dim; ++l)
	        (void) fprintf(file," %-"FFMT,Vel(Ndst)[l]);

	    (void) fprintf(file,"\n");
	}
}		/*end print_probe_field*/

EXPORT  void    g_fprint_front(
        Front   *front,
	FILE    *file)
{
	f_fprint_front(front,file);
	
	(void) foutput(file);
	(void) fprintf(file,"Gas Dynamics Front Parameters\n");

	(void) fprintf(file,"include_wall_normal_velocity = %s\n",
	               y_or_n(include_wall_normal_velocity(front)));
	fprint_float(file,"curvature_factor = ",curvature_factor(front),"\n");

	(void) fprintf(file,"\nEnd Gas Dynamics Front Parameters\n");
}/*end g_fprint_front*/

/*
 *      g_print_time_dependent_boundary_state()
 *      Print time dependent boundary state data info
 */

EXPORT  void    g_print_time_dependent_boundary_state(
        FILE            *file,
        INTERFACE       *intfc,
        BOUNDARY_STATE  *bstate)
{
        TD_BSTATE *tds = (TD_BSTATE*)bstate->_boundary_state_data;
        f_fprint_boundary_state_data(file,intfc,bstate);
        (void) fprintf(file,"Time Dependent Boundary State Data:\n");
        (void) fprintf(file,"  Data file = %s\n",tds->fn);
        (void) fprintf(file,"  Params = %llu\n",gas_param_number(tds->params));
        (void) fprintf(file,"End Time Dependent Boundary State Data\n");
}       /* end g_print_time_dependent_boundary_state */

#if defined(ONED)
LOCAL  void	g1d_accumulate_conserved_vars(
	FILE		*file,
	Front		*fr,
	Wave		*wv)
{
	Gas_param	*params;
	Gas_param	**prms_list;
	int		i, j;
	int		npts, nprms;
	INTERFACE	*intfc = fr->interf;
	RECT_GRID	*gr = computational_grid(intfc);
	POINT		**p;
	double		x,xl,xu;
	int		ix,ixs,ixp;
	int		icoords[3];
	double		dl, dr;
	double		dx;
	double		dxn;
	double	        voll, volr;
	double		vol;
	double		voln;
        int             xmax;
	Locstate	sl,sr;
	Locstate	s;
	static Locstate	*accum=NULL;

	if (intfc->dim != 1)
	    return;
		
	nprms = return_params_list(&prms_list);
	if(accum==NULL)
	{
	    uni_array(&accum,nprms,sizeof(Locstate));
	    for(i = 0; i < nprms; ++i)
	    {
                g_alloc_state(accum+i,fr->sizest);
	    }
	}
	for(i = 0; i < nprms; ++i)
	{
	    g_clear_state(accum[i],fr->sizest);
	}

	npts = intfc->num_points;
	ixs = 0;
        xmax = gr->gmax[0];
	dx = gr->h[0];

	for(p = intfc->points, i = 0; p && *p; ++p, ++i)
	{
	    if(Boundary(*p)) continue;
	    x = Coords(*p)[0];
	    ixp = cell_index(x,0,gr);
	    xl = cell_edge(ixp,0,gr); 
	    xu = xl + dx; 

	    if(i>0)
            {
                dl = x - max(xl,Coords(*(p-1))[0]);
            }
	    else
	        dl = x - xl;
	    if(gr->Remap.remap == CYLINDRICAL_REMAP)
	        voll = 2.0*PI*(x-0.5*dl)*dl;
	    else if(gr->Remap.remap == SPHERICAL_REMAP)
	        voll = (4.0*PI/3.0)*(x*x*x-(x-dl)*(x-dl)*(x-dl));
	    else /* IDENTITY */
	        voll = dl;
	    sr = left_state(*p);/*Right state on the interval being accumulated */
	    if ((sr!=NULL) && (params = Params(sr)) != NULL)
	    {
                if ((ixp > 0))
                {
                    if (i > 0)
                    {
                        double xx = Coords(*(p-1))[0];
	                int ixx = cell_index(xx,0,gr);
                        if (ixx == (ixp-1))
                        {
                            if ((right_state(*(p-1)) != NULL) && (Params(right_state(*(p-1))) != NULL))
                                sl = right_state(*(p-1));
                            else
                            {
	                        icoords[0] = ixp-1;
      	                        sl = Rect_state(icoords,wv);
                            }
                        }
                        else/*no intfc point in previous cell */
                        {
	                    icoords[0] = ixp-1;
      	                    sl = Rect_state(icoords,wv);
                        }
                    }
                }
                else
                {
                    if ((i > 0) && (right_state(*(p-1)) != NULL) && (Params(right_state(*(p-1))) != NULL))
                        sl = right_state(*(p-1)); 
                    else
                        sl = sr;
                }
	        for(j = 0; j < nprms; ++j)
		{
	            if( params == prms_list[j])
		    {
		        Dens(accum[j]) += 0.5*(Dens(sl)+Dens(sr))*voll;
		        Mom(accum[j])[0] += 0.5*(Mom(sl)[0]+Mom(sr)[0])*voll;
			Energy(accum[j]) += 0.5*(Energy(sl)+Energy(sr))*voll;
		        break;
		    }
		}
	    }
            while ((i+1) < npts) /*Multiple points in cell */
            {
                double xx = Coords(*(p+1))[0];
	        int ixx = cell_index(xx,0,gr);
                if (ixx != ixp)
                   break;
		dxn = xx - x;
	        if(gr->Remap.remap == CYLINDRICAL_REMAP)
	            voln = 2.0*PI*(x+0.5*dxn)*dxn;
	        else if(gr->Remap.remap == SPHERICAL_REMAP)
	            voln = (4.0*PI/3.0)*((x+dxn)*(x+dxn)*(x+dxn)-x*x*x);
	        else /* IDENTITY */
	            voln = dxn;
	        sl = right_state(*p);
	        sr = left_state(*(p+1));
	        if ((sl!=NULL) && ((params =Params(sl)) != NULL)  &&
			(sr!=NULL) && ((params =Params(sr)) != NULL))
		{
	            for(j = 0; j < nprms; ++j)
		    {
	                if( params == prms_list[j])
		        {
		            Dens(accum[j]) += 0.5*(Dens(sr)+Dens(sl))*voln;
		            Mom(accum[j])[0] += 0.5*(Mom(sr)[0]+Mom(sl)[0])*voln;
			    Energy(accum[j]) += 0.5*(Energy(sr)+Energy(sl))*voln;
		            break;
		        }
		    }
		}
                ++i;
                ++p;
                x = xx;
            }		

	    if (i+1 < npts)
	        dr = min(xu,Coords(*(p+1))[0])-x;
	    else
	        dr = xu - x;
	
	    if(gr->Remap.remap == CYLINDRICAL_REMAP)
	        volr = 2.0*PI*(x+0.5*dr)*dr;
	    else if(gr->Remap.remap == SPHERICAL_REMAP)
	        volr = (4.0*PI/3.0)*((x+dr)*(x+dr)*(x+dr)-x*x*x);
	    else /* IDENTITY */
	        volr = dr;
	    sl = right_state(*p);/*Left state of the interval being accumulated */
	    if ((sl!=NULL) && ((params =Params(sl)) != NULL))
	    {
                if ((ixp < (xmax-1)))
                {
                    if ((i+1) < npts)/*Interface point follows  */
                    {
                        double xx = Coords(*(p+1))[0];
	                int ixx = cell_index(xx,0,gr);
                        if (ixx == (ixp+1))
                        {
                            if ((left_state(*(p+1)) != NULL) && (Params(left_state(*(p+1))) != NULL))
                                sr = left_state(*(p+1));
                            else
                            {
	                        icoords[0] = ixp+1;
      	                        sr = Rect_state(icoords,wv);
                            }
                        }
                        else
                        {
	                    icoords[0] = ixp+1;
      	                    sr = Rect_state(icoords,wv);
                        }
                    }
                }
                else
                {
                    if ((i > 0) && (left_state(*(p+1)) != NULL) && (Params(left_state(*(p+1))) != NULL))
                        sr = left_state(*(p+1));
                    else
                        sr = sl;
                }
	        for(j = 0; j < nprms; ++j)
		{
	            if( params == prms_list[j])
		    {
		        Dens(accum[j]) += 0.5*(Dens(sr)+Dens(sl))*volr;
		        Mom(accum[j])[0] += 0.5*(Mom(sr)[0]+Mom(sl)[0])*volr;
			Energy(accum[j]) += 0.5*(Energy(sr)+Energy(sl))*volr;
		        break;
		    }
		}
	    }
	    for (ix = ixs; ix < ixp; ++ix)
	    {
	        icoords[0] = ix;
      	        s = Rect_state(icoords,wv);
	        if(s!=NULL)
	        {
	            x = cell_center(ix,0,gr); 
	            if(gr->Remap.remap == CYLINDRICAL_REMAP)
	                vol = 2.0*PI*x*dx;
    	            else if(gr->Remap.remap == SPHERICAL_REMAP)
	                vol = (4.0*PI/3.0)*(3.0*x*x+0.25*dx*dx)*dx;
	            else /* IDENTITY */
	                vol = dx;
	            params = Params(s);
	            for(j = 0; j < nprms; ++j)
		    {
	                if( params == prms_list[j])
		        {
		            Dens(accum[j]) += Dens(s)*vol;
		            Mom(accum[j])[0] += Mom(s)[0]*vol;
			    Energy(accum[j]) += Energy(s)*vol;
		            break;
		        }
	    	    }
	        }
	    }
	    ixs = ixp + 1;
	}

        for(j = 0; j < nprms; ++j)
	{
	    double tmp[3];
	    tmp[0] = Dens(accum[j]);
	    tmp[1] = Mom(accum[j])[0];
	    tmp[2] = Energy(accum[j]);
	    pp_global_sum(tmp,3);
	    Dens(accum[j]) = tmp[0];
	    Mom(accum[j])[0] = tmp[1];
	    Energy(accum[j]) = tmp[2];
	}

	(void) fprintf(file,"                     %"SFMT" %"SFMT" %"SFMT" %"SFMT"\n",
                     				  "Time","Mass","Momentum","Energy");	
        for(j = 0; j < nprms; ++j)
	{
	    (void) fprintf(file,"material totals %2d = ", j);	
	    (void) fprintf(file,"%"FFMT" %"FFMT" %"FFMT" %"FFMT"\n",fr->time,
                           Dens(accum[j]),Mom(accum[j])[0],Energy(accum[j]));
	}
	(void) fprintf(file,"\n");	
} /* end g1d_accumulate_conserved_vars */

/*#bjet2 */
typedef  struct {
	int	num_c_p;
	int  	*params;
	boolean 	is_ready;
}  COMP_PARAMS;

LOCAL   COMP_PARAMS comp_params = {
	0,
	NULL,
	NO
};

LOCAL  void print_comp_params(
	const char *msg, 
	int max_n_comps)
{
	int	comp;

	printf("#%s  is_ready=%d  num_c_p=%d  params=%d  max_n_comp=%d\n", msg, 
	    comp_params.is_ready, comp_params.num_c_p, 
	    comp_params.params, max_n_comps);
	for(comp = 0; comp < max_n_comps; comp++)
	    printf("%d\n", comp_params.params[comp]);
}

/*WARNNING make sure, before we call the function, max_num_comps must contain all the comp */
/*in the interface */
EXPORT  void communicate_comp_params(Front *fr, int max_n_comps)
{
INTERFACE	*intfc = fr->interf;
int		comp;
Gas_param  	*param;

	if(intfc->dim != 3)
	    return;

	if(comp_params.params == NULL)
	    uni_array(&comp_params.params, max_n_comps, INT);
	for(comp=0; comp < max_n_comps; comp++)
	    comp_params.params[comp] = -1;
	
	for(comp = min_component(intfc); comp <= max_component(intfc); comp++)
	{
	    param = gas_params_for_comp(comp, intfc);
	    /*printf("#comp %d  %d\n", comp, param); */
	    comp_params.params[comp] = index_of_Gas_param(param);
	}
	pp_global_imax(comp_params.params, (long)max_n_comps);
	comp_params.is_ready = YES;
	comp_params.num_c_p = max_n_comps;
	
	/*update min max comp */
	comp = min_component(intfc);
	pp_global_imin(&comp,1L);
	min_component(intfc) = comp;

	comp = max_component(intfc);
	pp_global_imax(&comp,1L);
	max_component(intfc) = comp;

	/*DEBUG_TMP printf("communicate_comp_params printout.\n"); */
	/*DEBUG_TMP printf("#min_max comp  %d   %d\n", min_component(intfc),  */
			/*DEBUG_TMP max_component(intfc)); */
	/*DEBUG_TMP print_comp_params("comm comp_params af", max_n_comps); */
}

EXPORT Gas_param*  param_for_comp(
	COMPONENT    comp)
{
	if(!comp_params.is_ready)
	    return NULL;

	if(comp<0  || comp >= comp_params.num_c_p)
	{
	    printf("WARNING  param_for_comp, component=%d is wrong, num_c_p=%d\n", 
	        comp, comp_params.num_c_p);
	    /*clean_up(ERROR); */
	    return NULL;
	}
	return  Params_of_index(comp_params.params[comp]);
}


#endif /*defined(ONED)*/
