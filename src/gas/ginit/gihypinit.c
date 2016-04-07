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
*				gihypinit.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Initialization interface to the gas dynamics hyperbolic solver.
*	Functions in this file set function pointers in the wave
*	and front structure to access the finite difference schemes
*	in this library.
*/


#include <ginit/ginit.h>
#include <ghyp/ghyp.h>
#include <gdecs/vecdecs.h>

	/*LOCAL Function Prototypes */
LOCAL	void	prompt_for_muscl_options(INIT_DATA*);
LOCAL	void	prompt_for_pseudo_unsplit_muscl_options(INIT_DATA*);

struct _IrregularSten_PromptType {
	const char *prompt, *select;
	void (*_npt_solver)(double,double,Locstate,const double*,
	                    int,int*,int*,Stencil*);
};
typedef struct _IrregularSten_PromptType IrregularSten_PromptType;

struct _TanSten_PromptType {
	const char        *prompt, *select;
	int	          npts_sten;
	void (*_one_side_npt_tang_solver)(double,double,Tan_stencil*,Locstate,
					  Front*);
};
typedef struct _TanSten_PromptType TanSten_PromptType;

	/*LOCAL Function Prototypes*/
LOCAL	void	g_set_muscl_opts(INIT_DATA*,Front*,Wave*);
LOCAL	void	print_MUSCL_options(Muscl_Opts*,IrregularSten_PromptType*,
				    TanSten_PromptType*);
LOCAL	void	print_cg_params(const char*,CG_PARAMS*);
LOCAL	void	prompt_for_cg_params(INIT_DATA*);

/*
*			g_set_hyp_solvers():
*
*	Sets the "three_pt" pointers in the Wave and Front structures to
*	the appropriate versions of the finite difference routines.
*/

EXPORT void g_set_hyp_solvers(
	INIT_DATA	*init,
	Wave		*wave,
	Front		*fr)
{
	int		dim = fr->rect_grid->dim;

	switch (dim)
	{
	case 1:
	    break;
	case 2:
	    wave->unsplit.flux = NULL;
	    wave->unsplit.flux_obl = NULL;
	    wave->unsplit.sources = NULL;
	    break;
	case 3:	/* TODO */
	    wave->unsplit.flux = NULL;
	    wave->unsplit.flux_obl = NULL;
	    wave->unsplit.sources = NULL;
		break;
	}

	if (strstr(wave->method,"VECTOR_LAX_WENDROFF") != NULL)
	{
	    wave->_vec_solver = LW_vec;
	    wave->_alloc_phys_vecs = LW_alloc_phys_vecs;
	    wave->_free_phys_vecs = LW_free_phys_vecs;
	    g_wave_load_state_vectors(wave) = g_load_state_vectors;
            g_wave_assign_wave_state_vectors(wave) =
                    g_assign_wave_state_vectors;
	    wave->_npt_solver = LW;
	    wave->npts_sten = 3;
	    fr->_npt_tang_solver = g_two_side_npt_tang_solver;
            fr->_npt_parab_tan_solver2d = parab_tan_solver2d;
            fr->_npt_parab_tan_solver3d = parab_tan_solver3d;
	    fr->_one_side_npt_tang_solver = LWoblique;
	    fr->npts_tan_sten = 3;
	}
	else if (strstr(wave->method,"VECTOR_MUSCL") != NULL)
	{
	    wave->_vec_solver = vector_FD;
	    g_wave_oned_interior_scheme(wave) = oned_MUSCL;
	    g_front_oned_tangential_scheme(fr) = oned_MUSCL;
	    g_set_muscl_opts(init,fr,wave);
	}
	else if (strstr(wave->method,"PLM") != NULL)
	{
	    wave->_vec_solver = vector_FD;
	    g_wave_oned_interior_scheme(wave) = oned_PLM;
	    g_front_oned_tangential_scheme(fr) = oned_PLM;
	    g_set_muscl_opts(init,fr,wave);
	}
	else if (strstr(wave->method,"TVD") != NULL)
	{
	    wave->_vec_solver = vector_FD;
	    g_wave_oned_interior_scheme(wave) = oned_TVD;
	    g_front_oned_tangential_scheme(fr) = oned_TVD;
	    wave->_npt_solver = point_FD;
	    fr->_one_side_npt_tang_solver = oblique_FD;
	    fr->_npt_tang_solver = g_two_side_npt_tang_solver;
            fr->_npt_parab_tan_solver2d = parab_tan_solver2d;
            fr->_npt_parab_tan_solver3d = parab_tan_solver3d;
	    g_wave_load_state_vectors(wave) = g_load_state_vectors;
            g_wave_assign_wave_state_vectors(wave) =
                    g_assign_wave_state_vectors;
	    wave->npts_sten = 7;
	    fr->npts_tan_sten = 7;
	    wave->_alloc_phys_vecs = LW_alloc_phys_vecs;
	    wave->_free_phys_vecs = LW_free_phys_vecs;
	}
	else if (strstr(wave->method,"LAX_WENDROFF") != NULL)
	{
	    wave->_npt_solver = LW;
	    fr->_npt_tang_solver = g_two_side_npt_tang_solver;
            fr->_npt_parab_tan_solver2d = parab_tan_solver2d;
            fr->_npt_parab_tan_solver3d = parab_tan_solver3d;
	    fr->_one_side_npt_tang_solver = LWoblique;
	    fr->npts_tan_sten = 3;
	}
	else if (strstr(wave->method,"LAX_FRIEDRICHS") != NULL)
	{
	    wave->_npt_solver = LF;
	    fr->_npt_tang_solver = g_two_side_npt_tang_solver;
            fr->_npt_parab_tan_solver2d = parab_tan_solver2d;
            fr->_npt_parab_tan_solver3d = parab_tan_solver3d;
	    fr->_one_side_npt_tang_solver = LFoblique;
	    fr->npts_tan_sten = 3;
	}
	else if (strstr(wave->method,"GODUNOV") != NULL)
	{
	    wave->_npt_solver = godunov;
	    fr->_npt_tang_solver = g_two_side_npt_tang_solver;
            fr->_npt_parab_tan_solver2d = parab_tan_solver2d;
            fr->_npt_parab_tan_solver3d = parab_tan_solver3d;
	    fr->_one_side_npt_tang_solver = godunovobl;
	    fr->npts_tan_sten = 3;
	}
	else if (strstr(wave->method,"ADVANCE_FRONTS_ONLY") == NULL)
	{
	    screen("ERROR in g_set_hyp_solvers(), Unknown hyp method\n");
	    clean_up(ERROR);
	}
	if (strstr(wave->method,"UNSPLIT") != NULL)
	    h_set_unsplit_options(&USopts(init));
}		/*end g_set_hyp_solvers*/

LOCAL	void	prompt_for_pseudo_unsplit_muscl_options(
	INIT_DATA	*init)
{
	prompt_for_unsplit_options(init);
	prompt_for_muscl_options(init);
}		/*end prompt_for_pseudo_unsplit_muscl_options*/

LOCAL	void	prompt_for_muscl_options(
	INIT_DATA	*init)
{
	Hyp_method                            *method = hyperbolic_method(init);
	const char                            *mname;
	Muscl_Opts	                      *mopts;
	char		                      s[Gets_BUF_SIZE];
	int		                      i, dim = i_intfc(init)->dim;
	MUSCL_PromptType_Reconstructor        *Sintrp;
	MUSCL_PromptType_Rsolver              *Rsolver;
	MUSCL_PromptType_characteristic_solve *Moc;
	IrregularSten_PromptType              Irrsten[3];
	TanSten_PromptType                    Tansten[5];


	zero_scalar(Irrsten,3*sizeof(IrregularSten_PromptType));
	i = 0;
	Irrsten[i].prompt = "Lax_Wendroff";
	Irrsten[i].select = "l";
	Irrsten[i]._npt_solver = LW;
	Irrsten[++i].prompt = "MUSCL";
	Irrsten[i].select = "m";
	Irrsten[i]._npt_solver = point_FD;

	zero_scalar(Tansten,5*sizeof(TanSten_PromptType));
	i = 0;
	Tansten[i].prompt = "Lax_Wendroff";
	Tansten[i].select = "l";
	Tansten[i].npts_sten = 3;
	Tansten[i]._one_side_npt_tang_solver = LWoblique;

	Tansten[++i].prompt = "Lax-Friedrichs";
	Tansten[i].select = "lf";
	Tansten[i].npts_sten = 3;
	Tansten[i]._one_side_npt_tang_solver = LFoblique;

	Tansten[++i].prompt = "First order godunov";
	Tansten[i].select = "g";
	Tansten[i].npts_sten = 3;
	Tansten[i]._one_side_npt_tang_solver = godunovobl;

	Tansten[++i].prompt = "MUSCL";
	Tansten[i].select = "m";
	Tansten[i].npts_sten = method->npts_sten;
	Tansten[i]._one_side_npt_tang_solver = oblique_FD;

	/* Set default options */

	mopts = &MusclOptions(init);
	mname = hyperbolic_method_name(init);
	if (strstr(mname,"PLM") != NULL)
	    set_plm_default_opts(mopts,method,dim);
	else
	    set_muscl_default_opts(mopts,method,dim); 
        if(material_composition_type(init) == MULTI_COMP_NON_REACTIVE)
        {
            mopts->kmax = NumberFloats(init);
            mopts->nfloats = NumberFloats(init);
        }

	print_MUSCL_options(mopts,Irrsten,Tansten);
	screen("Use all defaults for MUSCL code (dflt = y): ");
	(void) Gets(s);

	if (s[0] != 'n' && s[0] != 'N')
	    return;

	Sintrp = mopts->Sintrp;
	Rsolver = mopts->Rsolver;
	Moc = mopts->Moc;
	if ((Sintrp != NULL) && (Sintrp[0].prompt != NULL) &&
	    (Sintrp[1].prompt != NULL)) /* At least two choices available */
	{
	    screen("Choose the desired type of linear reconstruction, "
	           "Choices are\n");
	    for (i = 0; Sintrp[i].prompt != NULL; ++i)
	    {
	        screen("\t%s (%s",Sintrp[i].prompt,Sintrp[i].select);
	        if (mopts->_reconstructor == Sintrp[i].reconstructor)
	    	    screen(", default)\n");
	        else
	    	    screen(")\n");
	    }
	    screen("Enter choice: ");
	    (void) Gets(s);
	    for (i = 0; Sintrp[i].prompt != NULL; ++i)
	    {
	        if (strcasecmp(s,Sintrp[i].select) == 0)
	        {
	    	    mopts->_reconstructor = Sintrp[i].reconstructor;
		    if (Sintrp[i].half_step != NULL)
		    {
	                mopts->_half_step = Sintrp[i].half_step;
	                mopts->_strong_wave_half_step =
			    Sintrp[i].strong_wave_half_step;
		    }
	    	    break;
	        }
	    }
	}

	if ((mopts->_strong_wave_half_step != NULL) &&
	    (mopts->_half_step != mopts->_strong_wave_half_step))
	{
	    screen("Test for negative density and energies "
	           "at half step (dflt=no): ");
	    (void) Gets(s);
	    if (s[0] == 'Y' || s[0] == 'y')
	        mopts->_half_step = mopts->_strong_wave_half_step;
	}

	if ((Rsolver != NULL) && (Rsolver[0].prompt != NULL) &&
	    (Rsolver[1].prompt != NULL))
	{
	    screen("\nChoose the desired Riemann solver, Choices are\n");
	    for (i = 0; Rsolver[i].prompt != NULL; ++i)
	    {
	        screen("\t%s (%s",Rsolver[i].prompt,Rsolver[i].select);
	        if (mopts->_rsolver == Rsolver[i].rsolver)
	    	    screen(", default)\n");
	        else
	    	    screen(")\n");
	    }
	    screen("Enter choice: ");
	    (void) Gets(s);
	    for (i = 0; Rsolver[i].prompt != NULL; ++i)
	    {
	        if (strcasecmp(s,Rsolver[i].select) == 0)
	        {
	    	    mopts->_rsolver = Rsolver[i].rsolver;
	    	    mopts->_rmidstate = Rsolver[i].rmidstate;
	    	    break;
	        }
	    }
	}

	if (mopts->_rsolver == cg_rsolve)
	    prompt_for_cg_params(init);

	if ((Moc != NULL) && (Moc[0].prompt != NULL) &&
	    (Moc[1].prompt != NULL))
	{
	    screen("\nChoose the desired method of characteristic solver, "
	           "Choices are\n");
	    for (i = 0; Moc[i].prompt != NULL; ++i)
	    {
	        screen("\t%s (%s",Moc[i].prompt,Moc[i].select);
	        if (mopts->_characteristic_solve==Moc[i].characteristic_solve)
	    	    screen(", default)\n");
	        else
	    	    screen(")\n");
	    }
	    screen("Enter choice: ");
	    (void) Gets(s);
	    for (i = 0; Moc[i].prompt != NULL; ++i)
	    {
	        if (strcasecmp(s,Moc[i].select) == 0)
	        {
	    	    mopts->_characteristic_solve = Moc[i].characteristic_solve;
	    	    break;
	        }
	    }
	}

	if (Irrsten[1].prompt != NULL)
	{
	    screen("\nChoose the irregular stencil method, Choices are\n");
	    for (i = 0; Irrsten[i].prompt != NULL; ++i)
	    {
	        screen("\t%s (%s",Irrsten[i].prompt,Irrsten[i].select);
	        if (mopts->_npt_solver == Irrsten[i]._npt_solver)
	    	    screen(", default)\n");
	        else
	    	    screen(")\n");
	    }
	    screen("Enter choice: ");
	    (void) Gets(s);
	    for (i = 0; Irrsten[i].prompt != NULL; ++i)
	    {
	        if (strcasecmp(s,Irrsten[i].select) == 0)
	        {
	    	    mopts->_npt_solver = Irrsten[i]._npt_solver;
	    	    break;
	        }
	    }
	}
	if (mopts->_npt_solver == LW)
	{
	    mopts->_one_side_npt_tang_solver = LWoblique;
	    mopts->_npts_tan_sten = 3;
	}
	else if (mopts->_npt_solver == point_FD)
	{
	    mopts->_one_side_npt_tang_solver = oblique_FD;
	    mopts->_npts_tan_sten = method->npts_sten;
	}

	if ((Tansten[1].prompt != NULL) && (dim > 1))
	{
	    screen("\nChoose the tangential sweep method, Choices are\n");
	    for (i = 0; Tansten[i].prompt != NULL; ++i)
	    {
	        screen("\t%s (%s",Tansten[i].prompt,Tansten[i].select);
	        if (mopts->_one_side_npt_tang_solver ==
	            Tansten[i]._one_side_npt_tang_solver)
		    screen(", default)\n");
	        else
		    screen(")\n");
	    }
	    screen("Enter choice: ");
	    (void) Gets(s);
	    for (i = 0; Tansten[i].prompt != NULL; ++i)
	    {
	        if (strcasecmp(s,Tansten[i].select) == 0)
	        {
		    mopts->_one_side_npt_tang_solver =
		        Tansten[i]._one_side_npt_tang_solver;
		    mopts->_npts_tan_sten = Tansten[i].npts_sten;
		    break;
	        }
	    }
	}
	if (strstr(mname,"PLM") != NULL)
	{
	    screen("Enforce monotone state reconstruction at cell edges "
	           "(dflt = %s): ",y_or_n(mopts->monotone_reconstruction));
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	       if ((s[0] == 'n') || (s[0] == 'N'))
	           mopts->monotone_reconstruction = NO;
	       if ((s[0] == 'y') || (s[0] == 'Y'))
	           mopts->monotone_reconstruction = YES;
	    }
	    screen("Link reconstructions,  zero slope in one field implies "
	           "zero slope in all fields (dflt = %s): ",
		   y_or_n(mopts->link_reconstructions));
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	       if ((s[0] == 'n') || (s[0] == 'N'))
	           mopts->link_reconstructions = NO;
	       if ((s[0] == 'y') || (s[0] == 'Y'))
	           mopts->link_reconstructions = YES;
	    }
	}
}		/*end prompt_for_muscl_options*/

LOCAL	void	prompt_for_cg_params(
	INIT_DATA	*init)
{
	CG_PARAMS   *cg_params = &Cg_params(&MusclOptions(init));
	char        s[Gets_BUF_SIZE];

	print_cg_params("",cg_params);
	screen("Use all defaults for Colella-Glaz Riemann solver (dflt = y): ");
	(void) Gets(s);

	if (s[0] != 'n' && s[0] != 'N')
	    return;

	screen("Enter the number of iterations "
	       "for the Riemann Solver (dflt = %d): ",cg_params->pv_iterations);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%d",&cg_params->pv_iterations);
	if (cg_params->pv_iterations <= 0)
	{
	    screen("ERROR in prompt_for_cg_params(), "
	           "number of iterations must be positive\n");
	    clean_up(ERROR);
	}

	screen("Enter the minimum pressure jump below which the\n\t"
	       "mass flux is replaced by the acoustic impedance (dflt = %g): ",
	       cg_params->min_p_jump);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&cg_params->min_p_jump);
	if (cg_params->min_p_jump < 0.0)
	{
	    screen("ERROR in prompt_for_cg_params(), "
	           "min pressure jump must be nonnegative\n");
	    clean_up(ERROR);
	}
	screen("Enter the velocity convergence factor (dflt = %g): ",
	       cg_params->min_v_jump);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&cg_params->min_v_jump);
	if (cg_params->min_v_jump <= 0.0)
	{
	    screen("ERROR in prompt_for_cg_params(), "
	           "velocity convergence factor must be positive\n");
	    clean_up(ERROR);
	}

	screen("Enter the strong wave tolerance to turn "
	       "on exact solver (dflt = %g): ",cg_params->sw_tol);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&cg_params->sw_tol);
	if (cg_params->sw_tol < 0.0)
	{
	    screen("ERROR in prompt_for_cg_params(), "
	           "strong wave tolerance must be nonnegative\n");
	    clean_up(ERROR);
	}

	screen("Enter the minimum allowed mass flux (dflt = %g): ",
	       cg_params->min_mass_flux);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&cg_params->min_mass_flux);
	if (cg_params->min_mass_flux <= 0.0)
	{
	    screen("ERROR in prompt_for_cg_params(), "
	           "minimum mass flux must be positive\n");
	    clean_up(ERROR);
	}
	cg_params->sqr_min_mass_flux = sqr(cg_params->min_mass_flux);
}		/*end prompt_for_cg_params*/

LOCAL	void g_set_muscl_opts(
	INIT_DATA	*init,
	Front		*fr,
	Wave		*wave)
{
	Muscl_Opts *mopts = muscl_options();

	*mopts = MusclOptions(init);
	wave->_npt_solver = mopts->_npt_solver;
	fr->_one_side_npt_tang_solver = mopts->_one_side_npt_tang_solver;
	fr->npts_tan_sten = mopts->_npts_tan_sten;
	fr->_npt_tang_solver = mopts->_npt_tang_solver;
        fr->_npt_parab_tan_solver2d = mopts->_npt_parab_tan_solver2d;
        fr->_npt_parab_tan_solver3d = mopts->_npt_parab_tan_solver3d;
	wave->_alloc_phys_vecs = mopts->_alloc_phys_vecs;
	wave->_free_phys_vecs = mopts->_free_phys_vecs;
	g_wave_load_state_vectors(wave) = mopts->_load_state_vectors;
	g_wave_assign_wave_state_vectors(wave) =
	    mopts->_assign_wave_state_vectors;

	if (mopts->_rsolver == cg_rsolve)
	    set_cg_params(&Cg_params(mopts));
}		/*end g_set_muscl_opts*/

LOCAL	void	print_MUSCL_options(
	Muscl_Opts		       *mopts,
	IrregularSten_PromptType       *Irrsten,
	TanSten_PromptType	       *Tansten)
{
	MUSCL_PromptType_Reconstructor        *Sintrp = mopts->Sintrp;
	MUSCL_PromptType_Rsolver              *Rsolver = mopts->Rsolver;
	MUSCL_PromptType_characteristic_solve *Moc = mopts->Moc;
	int i;

	if (Sintrp != NULL)
	{
	    screen("\nCurrent values for MUSCL parameters\n"
	           "\tState reconstruction = ");
	    for (i = 0; Sintrp[i].prompt != NULL; ++i)
	    {
	        if (mopts->_reconstructor == Sintrp[i].reconstructor)
		{
	    	    screen("%s\n",Sintrp[i].prompt);
		    break;
		}
	    }
	    if (Sintrp[i].prompt == NULL)
	        screen("\tUnknown state reconstructor\n");
	}
	else if (mopts->_reconstructor != NULL)
	    screen("\tUnnamed state reconstructor\n");
	else
	    screen("\tNo state reconstructor used\n");

	if (Rsolver != NULL)
	{
	    screen("\tRiemann flux solver = ");
	    for (i = 0; Rsolver[i].prompt != NULL; ++i)
	    {
	        if (mopts->_rsolver == Rsolver[i].rsolver)
		{
	    	    screen("%s\n",Rsolver[i].prompt);
		    break;
		}
	    }
	    if (Rsolver[i].prompt == NULL)
	        screen("\tUnknown Riemann flux solver\n");
	}
	else if (mopts->_rsolver != NULL)
	    screen("\tUnnamed Riemann flux solver\n");
	else
	    screen("\tNo Riemann flux solver used\n");

	if (mopts->_rsolver == cg_rsolve)
	    print_cg_params("\t",&Cg_params(mopts));

	if (Moc != NULL)
	{
	    screen("\tMethod of characteristic solver = ");
	    for (i = 0; Moc[i].prompt != NULL; ++i)
	    {
	        if (mopts->_characteristic_solve == Moc[i].characteristic_solve)
		{
	    	    screen("%s\n",Moc[i].prompt);
		    break;
		}
	    }
	    if (Moc[i].prompt == NULL)
	        screen("\tUnknown method of characteristic solver\n");
	}
	else if (mopts->_characteristic_solve != NULL)
	    screen("\tUnnamed Method of characteristic solver\n");
	else
	    screen("\tNo method of characteristic solver used\n");

	screen("\tIrregular hyp stencil method = ");
	for (i = 0; Irrsten[i].prompt != NULL; ++i)
	{
	    if (mopts->_npt_solver == Irrsten[i]._npt_solver)
	    	screen("%s\n",Irrsten[i].prompt);
	}
	screen("\tTangential sweep method = ");
	for (i = 0; Tansten[i].prompt != NULL; ++i)
	{
	    if (mopts->_one_side_npt_tang_solver ==
	        Tansten[i]._one_side_npt_tang_solver)
	    	screen("%s\n",Tansten[i].prompt);
	}
	screen("\t%s for negative density and energies at half step\n",
	       (mopts->_half_step == mopts->_strong_wave_half_step) ?
	       "Test" : "Don't test");
	screen("\tEnforce monotone reconstructions at cell edges = %s\n",
	       y_or_n(mopts->monotone_reconstruction));
	screen("\tLink reconstructions (zero slope in one field implies\n"
	       "\tzero slope in all fields) = %s\n",
	       y_or_n(mopts->link_reconstructions));
	screen("End List of current values for MUSCL parameters\n\n");
}		/*end print_MUSCL_options*/

LOCAL  void print_cg_params(
	const char *indent,
	CG_PARAMS  *cg_params)
{
	if (indent == NULL)
	    indent = "";
	screen("\n%sCurrent values of the Colella-Glaz Riemann solver "
	       "parameters\n",indent);
	screen("%s\t%-55s%d\n",indent,
	       "Maximum number of iterations for Riemann solver",
	       cg_params->pv_iterations);
	screen("%s\t%-55s%g\n%s\t\t%s\n",indent,
	       "Minimum pressure jump below which mass flux",
	       cg_params->min_p_jump,
	       indent,"is replaced by acoustic impedance");
	screen("%s\t%-55s%g\n",indent,"Velocity convergence tolerance",
	       cg_params->min_v_jump);
	screen("%s\t%-55s%g\n%s\t\t%s\n",indent,
	       "The strong wave tolerance above which approximate",
	       cg_params->sw_tol,
	       indent,"Riemann solver is replaced by an exact solver");
	screen("%s\t%-55s%g\n",indent,
	       "Minimum allowed mass flux",
	       cg_params->min_mass_flux);
	screen("%sEnd List of Colella-Glaz Riemann solver parameters\n\n",
	       indent);
}        /*end print_default_cg_params*/

EXPORT	void	g_setup_available_hyperbolic_methods_list(
	INIT_DATA	*init)
{
	int               i, dim = i_intfc(init)->dim;
	static Hyp_method g_Methods[20];
	Hyp_method        *dflt;

	i = 0;
	g_Methods[i].ptype.prompt = (dim > 1) ? "Split Lax-Wendroff" :
	                                        "Lax-Wendroff";
	g_Methods[i].ptype.select = "LWS";
	g_Methods[i].ptype.ncmp = 3;
	g_Methods[i].ptype.type.ctype = "LAX_WENDROFF";
	g_Methods[i].npts_sten = 3;
	g_Methods[i].npts_vsten = 0;
	g_Methods[i].sten_rad = 1;
	g_Methods[i].hyp_driver = hyp_scalar_driver;
	g_Methods[i]._prompt_for_hyp_method_options = NULL;

	++i;
	g_Methods[i].ptype.prompt = (dim > 1) ? "Split Lax-Friedrichs" :
	                                        "Lax-Friedrichs";
	g_Methods[i].ptype.select = "LFS";
	g_Methods[i].ptype.ncmp = 3;
	g_Methods[i].ptype.type.ctype = "LAX_FRIEDRICHS";
	g_Methods[i].npts_sten = 3;
	g_Methods[i].npts_vsten = 0;
	g_Methods[i].sten_rad = 1;
	g_Methods[i].hyp_driver = hyp_scalar_driver;
	g_Methods[i]._prompt_for_hyp_method_options = NULL;

	++i;
	g_Methods[i].ptype.prompt = (dim > 1) ? "Split first order Godunov" :
	                                        "First order Godunov";
	g_Methods[i].ptype.select = "G";
	g_Methods[i].ptype.ncmp = 1;
	g_Methods[i].ptype.type.ctype = "GODUNOV";
	g_Methods[i].npts_sten = 3;
	g_Methods[i].npts_vsten = 0;
	g_Methods[i].sten_rad = 1;
	g_Methods[i].hyp_driver = hyp_scalar_driver;
	g_Methods[i]._prompt_for_hyp_method_options = NULL;

	++i;
	g_Methods[i].ptype.prompt = (dim > 1) ? "Vectorized split Lax-Wendroff":
	                                        "Vectorized Lax-Wendroff";
	g_Methods[i].ptype.select = "VLS";
	g_Methods[i].ptype.ncmp = 3;
	g_Methods[i].ptype.type.ctype = "VECTOR_LAX_WENDROFF";
	g_Methods[i].npts_sten = 3;
	g_Methods[i].npts_vsten = 3;
	g_Methods[i].sten_rad = 1;
	g_Methods[i].hyp_driver = hyp_vector_driver;
	g_Methods[i]._prompt_for_hyp_method_options = NULL;

	++i;
	g_Methods[i].ptype.prompt = (dim > 1) ? 
	    "Five point Vectorized split MUSCL" : "Five point Vectorized MUSCL";
	g_Methods[i].ptype.select = "VM";
	g_Methods[i].ptype.ncmp = 2;
	g_Methods[i].ptype.type.ctype = "VECTOR_MUSCL";
	g_Methods[i].npts_sten = 5;
	g_Methods[i].npts_vsten = 5;
	g_Methods[i].sten_rad = 2;
	g_Methods[i].hyp_driver = hyp_vector_driver;
	g_Methods[i]._prompt_for_hyp_method_options = prompt_for_muscl_options;
	dflt = &g_Methods[i];

	++i;
	g_Methods[i].ptype.prompt = "Colella Piecewise Linear Method";
	g_Methods[i].ptype.select = "PLM";
	g_Methods[i].ptype.ncmp = 3;
	g_Methods[i].ptype.type.ctype = "PLM";
	g_Methods[i].npts_sten = 7;
	g_Methods[i].npts_vsten = 7;
	g_Methods[i].sten_rad = 3;
	g_Methods[i].hyp_driver = hyp_vector_driver;
	g_Methods[i]._prompt_for_hyp_method_options = prompt_for_muscl_options;
	dflt = &g_Methods[i];

	++i;
	g_Methods[i].ptype.prompt = (dim > 1) ? 
	    "Vectorized split TVD" : "Vectorized TVD";
	g_Methods[i].ptype.select = "TVD";
	g_Methods[i].ptype.ncmp = 2;
	g_Methods[i].ptype.type.ctype = "VECTOR_TVD";
	g_Methods[i].npts_sten = 7;
	g_Methods[i].npts_vsten = 7;
	g_Methods[i].sten_rad = 3;
	g_Methods[i].hyp_driver = hyp_vector_driver;
	/*g_Methods[i]._prompt_for_hyp_method_options = prompt_for_muscl_options; */
	dflt = &g_Methods[i];

	if (dim > 1)
	{
	    ++i;
	    g_Methods[i].ptype.prompt =
	        "Vectorized pseudo unsplit Lax-Wendroff";
	    g_Methods[i].ptype.select = "PUSLW";
	    g_Methods[i].ptype.ncmp = 5;
	    g_Methods[i].ptype.type.ctype =
	        "PSEUDO_UNSPLIT_VECTOR_LAX_WENDROFF";
	    g_Methods[i].npts_sten = 3;
	    g_Methods[i].npts_vsten = 3;
	    g_Methods[i].sten_rad = 1;
	    g_Methods[i].hyp_driver = pseudo_unsplit_driver;
	    g_Methods[i]._prompt_for_hyp_method_options =
	        prompt_for_unsplit_options;

	    ++i;
	    g_Methods[i].ptype.prompt = "Vectorized pseudo unsplit MUSCL";
	    g_Methods[i].ptype.select = "PUSM";
	    g_Methods[i].ptype.ncmp = 4;
	    g_Methods[i].ptype.type.ctype = "PSEUDO_UNSPLIT_VECTOR_MUSCL";
	    g_Methods[i].npts_sten = 5;
	    g_Methods[i].npts_vsten = 5;
	    g_Methods[i].sten_rad = 2;
	    g_Methods[i].hyp_driver = pseudo_unsplit_driver;
	    g_Methods[i]._prompt_for_hyp_method_options =
	        prompt_for_pseudo_unsplit_muscl_options;

	    ++i;
	    g_Methods[i].ptype.prompt = "Colella pseudo unsplit "
	                                "Piecewise Linear Method";
	    g_Methods[i].ptype.select = "PUSPLM";
	    g_Methods[i].ptype.ncmp = 6;
	    g_Methods[i].ptype.type.ctype = "PSEUDO_UNSPLIT_PLM";
	    g_Methods[i].npts_sten = 7;
	    g_Methods[i].npts_vsten = 7;
	    g_Methods[i].sten_rad = 3;
	    g_Methods[i].hyp_driver = pseudo_unsplit_driver;
	    g_Methods[i]._prompt_for_hyp_method_options =
	        prompt_for_pseudo_unsplit_muscl_options;
	    dflt = &g_Methods[i];
	}

	++i;
	g_Methods[i].ptype.prompt = NULL;
	g_Methods[i].ptype.select = NULL;
	g_Methods[i].ptype.ncmp = 0;
	g_Methods[i].ptype.type.ctype = NULL;
	g_Methods[i].npts_sten = 0;
	g_Methods[i].npts_vsten = 0;
	g_Methods[i].sten_rad = 0;
	g_Methods[i].hyp_driver = NULL;
	g_Methods[i]._prompt_for_hyp_method_options = NULL;

	available_hyperbolic_methods(init) = g_Methods;
	if (default_hyperbolic_method(init) == NULL)
	    default_hyperbolic_method(init) = dflt;
}		/*end g_setup_available_hyperbolic_methods_list*/
