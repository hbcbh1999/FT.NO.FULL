#include "solver.h"

void poisson_solver3d_P0_vd(
        Front *front,       /* front structure containing interface geometry */
        int ilower,
        int iupper,
        int ***ijk_to_I,
        double *source,     /* source term in poisson's equation */
        double *kk,         /* nonhomogeneous Neumann BC in poisson's equation */
        double *soln,       /* solution of the poisson's equation */
        double *max_soln,   /* maximum value of the solution */
        double *min_soln)   /* minimum value of the solution */
{
        int index,index_nb[6],size,sign;
        double rhs,coeff[6];
        int I,I_nb[6];
        int i,j,k,l,icoords[MAXD];
        int imin,imax,jmin,jmax,kmin,kmax;
        COMPONENT comp;
        double aII;
        int num_nb;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
        boolean use_neumann_solver = YES;
        PetscInt num_iter = 0;
        double rel_residual = 0;
        RECT_GRID *rgr = &topological_grid(front->grid_intfc);
        struct Table *T = table_of_interface(front->grid_intfc);
        COMPONENT *top_comp = T->components;
        int *lbuf = front->rect_grid->lbuf;
        int *ubuf = front->rect_grid->ubuf;
        int *top_gmax = rgr->gmax;
        double *top_h = rgr->h;
        HYPER_SURF *hs;
        int reflect[MAXD];

        if (front->grid_intfc == NULL) clean_up(ERROR);
        PETSc solver;
        solver.Create(ilower, iupper-1, 7, 7);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();
        size = iupper - ilower;
        *max_soln = -HUGE;
        *min_soln = HUGE;

        imin = (lbuf[0] == 0) ? 1 : lbuf[0];
        jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
        kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
        imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
        jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
        kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];

        printf("imin = %d	imax = %d\n", imin, imax);
        printf("jmin = %d       jmax = %d\n", jmin, jmax);
        printf("kmin = %d       kmax = %d\n", kmin, kmax);
        fflush(stdout);

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            comp = top_comp[index];
            I = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;

            num_nb = 0;
            for (l = 0; l < 6; ++l)
            {
                if (I_nb[l] == -1);
                else num_nb++;
                coeff[l] = 1.0/top_h[l/2]/top_h[l/2];
            }

            rhs = source[index];

            aII = 0;
            for (l = 0; l < 6; ++l)
            {
                if (num_nb < 2) break;
                if (I_nb[l] != -1)
                {
                    solver.Set_A(I,I_nb[l],coeff[l]);
                    aII += -coeff[l];
                }
                else
                {
                    hs = FT_HyperSurfAtGridCrossing(front,icoords,dir[l]);
                    if (hs != NULL && wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        if (!boundary_state(hs))
                        {
                            /* The boundary condition of phi at present
                               Dirichlet boundary is set to zero, it then
                               uses non-Neumann poisson solver
                            */
                            use_neumann_solver = NO;
                            aII += -coeff[l];
                        }
                    }
                    else //nonhomogeneous Neumann B.C. for p0 on top/bottom: dp/dn = rho*g
                    {
                        sign = (l%2==0)? 1 : -1;
                        rhs += sign*kk[index]/top_h[l/2];
                    }
                }
            }
            /*
             * This change reflects the need to treat point with only one
             * interior neighbor (a convex point). Not sure why PETSc cannot
             * handle such case. If we have better understanding, this should
             * be changed back.
             */
            if(num_nb >= 2)
            {
                solver.Set_A(I,I,aII);
            }
            else
            {
                solver.Set_A(I,I,1.0);
            }
            solver.Set_b(I,rhs);
        }
        use_neumann_solver = pp_min_status(use_neumann_solver);

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

        start_clock("Before Petsc Solver");
        if (use_neumann_solver)
        {
            printf("\nUsing Neumann Solver for p0!\n");
            solver.Solve_withPureNeumann();

/*
            if (debugging("PETSc"))
            {
                double max, min;
                solver.GetExtremeSingularValues(&max, &min); //max=min=-1, only can be used for GMRES
                (void) printf("The max singular value of A = %lf in poisson_solver3d_P0_vd\n", max);
                (void) printf("The min singular value of A = %lf in poisson_solver3d_P0_vd\n", min);
                if ( min != 0.0)
                    (void) printf("The Cond Num of A = %lf in poisson_solver3d_P0_vd\n", max/min);
            }
*/

            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);
            if(rel_residual > 1e-10)
            {
                printf("\n The solution diverges for p0! The residual "
                       "is %le. Solve again using GMRES!\n",rel_residual);
                solver.Reset_x();
                solver.Solve_withPureNeumann_GMRES();

                if (debugging("PETSc"))
                {
                    double max, min;
                    solver.GetExtremeSingularValues(&max, &min);
                    (void) printf("The max singular value of A = %lf in poisson_solver3d_P0_vd\n", max);
                    (void) printf("The min singular value of A = %lf in poisson_solver3d_P0_vd\n", min);
                    if ( min != 0.0)
                        (void) printf("The Cond Num of A = %lf in poisson_solver3d_P0_vd\n", max/min);
                }

                solver.GetNumIterations(&num_iter);
                solver.GetFinalRelativeResidualNorm(&rel_residual);
            }

        }
        else
        {
            printf("\nUsing non-Neumann Solver for p0!\n");
            solver.Solve();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

            if(rel_residual > 1e-10)
            {
                printf("\n The solution diverges for p0! The residual "
                       "is %le. Solve again using GMRES!\n",rel_residual);
                solver.Reset_x();
                solver.Solve_GMRES();
                solver.GetNumIterations(&num_iter);
                solver.GetFinalRelativeResidualNorm(&rel_residual);
            }

        }
        stop_clock("After Petsc Solver");

        double *x;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("In poisson_solver3d_P0_vd(): "
                        "num_iter = %d, rel_residual = %le \n",
                        num_iter, rel_residual);

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            I = ijk_to_I[i][j][k];
            soln[index] = x[I-ilower];
            if (*max_soln < soln[index]) *max_soln = soln[index];
            if (*min_soln > soln[index]) *min_soln = soln[index];
        }
        reflect[0] = reflect[1] =reflect[2] = YES;
        FT_ParallelExchGridArrayBuffer(soln,front,reflect);
        pp_global_max(max_soln,1);
        pp_global_min(min_soln,1);

        if(debugging("step_size"))
        {
            printf("\nThe max value of p0 is %.16g\n",*max_soln);
            printf("\nThe min value of p0 is %.16g\n",*min_soln);
        }

        FT_FreeThese(1,x);
} /* end poisson_solver3d_P0_vd */


void poisson_solver2d(
	Front *front,	    /* front structure containing interface geometry */
	int ilower,
	int iupper,
	int **ij_to_I,
	double *source,	    /* source term in poisson's equation */
	double *k,	    /* coefficient function in poisson's equation */
	double *soln,	    /* solution of the poisson's equation */
	double *max_soln,   /* maximum value of the solution */
	double *min_soln)   /* minimum value of the solution */
{
	int index,index_nb[4],size;
	double k0,k_nb[4];
	double rhs,coeff[4];
	int I,I_nb[4];
	int i,j,l,icoords[MAXD];
	int imin,imax,jmin,jmax;
	COMPONENT comp;
	double aII;
	int num_nb;
	GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
	boolean use_neumann_solver = YES;
	PetscInt num_iter = 0;
	double rel_residual = 0.0;
	RECT_GRID *rgr = &topological_grid(front->grid_intfc);
	struct Table *T = table_of_interface(front->grid_intfc);
	COMPONENT *top_comp = T->components;
	int *lbuf = front->rect_grid->lbuf;
	int *ubuf = front->rect_grid->ubuf;
	int *top_gmax = rgr->gmax;
	double *top_h = rgr->h;
	HYPER_SURF *hs;
    int reflect[MAXD];

	PETSc solver;
	solver.Create(ilower, iupper-1, 5, 0);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();
	size = iupper - ilower;
	*max_soln = -HUGE;
	*min_soln = HUGE;

	imin = (lbuf[0] == 0) ? 1 : lbuf[0];
        jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
        jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    comp = top_comp[index];
	    I = ij_to_I[i][j];
	    if (I == -1) continue;

	    index_nb[0] = d_index2d(i-1,j,top_gmax);
	    index_nb[1] = d_index2d(i+1,j,top_gmax);
	    index_nb[2] = d_index2d(i,j-1,top_gmax);
	    index_nb[3] = d_index2d(i,j+1,top_gmax);
	    I_nb[0] = ij_to_I[i-1][j];
	    I_nb[1] = ij_to_I[i+1][j];
	    I_nb[2] = ij_to_I[i][j-1];
	    I_nb[3] = ij_to_I[i][j+1];
	    icoords[0] = i;
	    icoords[1] = j;

	    k0 = k[index];
	    num_nb = 0;
	    for (l = 0; l < 4; ++l)
	    {
		if (I_nb[l] == -1)
		    index_nb[l] = index;
		else num_nb++;
		k_nb[l] = 0.5*(k0 + k[index_nb[l]]);
	    	coeff[l] = 1.0/k_nb[l]/(top_h[l/2]*top_h[l/2]);
	    }

	    rhs = source[index];

	    aII = 0.0;
	    for (l = 0; l < 4; ++l)
	    {
		if (num_nb < 2) break;
	    	if (I_nb[l] != -1)
		{
		    solver.Set_A(I,I_nb[l],coeff[l]);
                    aII += -coeff[l];
		}
		else
		{
		    hs = FT_HyperSurfAtGridCrossing(front,icoords,dir[l]);
		    if (hs != NULL && wave_type(hs) == DIRICHLET_BOUNDARY)
		    {
			if (!boundary_state(hs))
			{
			    /* The boundary condition of phi at preset
			       Dirichlet boundary is set to zero, it then
			       uses non-Neumann poisson solver
			    */
			    use_neumann_solver = NO;
                    	    aII += -coeff[l];
			}
		    }
		}
	    }
	    /*
	     * This change reflects the need to treat point with only one
	     * interior neighbor (a convex point). Not sure why PETSc cannot
	     * handle such case. If we have better understanding, this should
	     * be changed back.
	     */
	    if(num_nb >= 2)
	    {
                solver.Set_A(I,I,aII);
	    }
            else
            {

                solver.Set_A(I,I,1.0);
            }
            solver.Set_b(I,rhs);
	}
	use_neumann_solver = pp_min_status(use_neumann_solver);

	solver.SetMaxIter(40000);
	solver.SetTol(1e-14);

	start_clock("Before Petsc Solver");
	if (use_neumann_solver)
	{
	    printf("\nUsing Neumann Solver for phi!\n");
	    solver.Solve_withPureNeumann();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	    if(rel_residual > 1)
	    {
		printf("\n The solution diverges for phi! The residual "
		       "is %le. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_withPureNeumann_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	else
	{
	    printf("\nUsing non-Neumann Solver for phi!\n");
	    solver.Solve();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);

	    if(rel_residual > 1)
	    {
		printf("\n The solution diverges for phi! The residual "
		       "is %le. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	stop_clock("After Petsc Solver");

	double *x;
	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

	if (debugging("PETSc"))
	    (void) printf("In poisson_solver(): "
	       		"num_iter = %d, rel_residual = %le \n",
			num_iter, rel_residual);

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    I = ij_to_I[i][j];
	    soln[index] = x[I-ilower];
	    if (*max_soln < soln[index]) *max_soln = soln[index];
	    if (*min_soln > soln[index]) *min_soln = soln[index];
	}
    reflect[0] = reflect[1] = YES;
	FT_ParallelExchGridArrayBuffer(soln,front, reflect);
	pp_global_max(max_soln,1);
	pp_global_min(min_soln,1);

	if(debugging("step_size"))
	{
	    printf("\nThe max solution value is %.16g\n",*max_soln);
	    printf("\nThe min solution value is %.16g\n",*min_soln);
	}

	FT_FreeThese(1,x);
}	/* end computeProjection2d */


void poisson_solver2d_vd(
        Front *front,       /* front structure containing interface geometry */
        int ilower,
        int iupper,
        int **ij_to_I,
        double *source,     /* source term in poisson's equation */
        double *k,          /* coefficient function in poisson's equation */
        double *k_old,      /* coefficient function in poisson's equation */
        double *soln,       /* solution of the poisson's equation */
        double *max_soln,   /* maximum value of the solution */
        double *min_soln)   /* minimum value of the solution */
{
        int index,index_nb[4],size;
        double k0,k_nb[4],k0_old,k_nb_old[4];
        double rhs,coeff[4];
        int I,I_nb[4];
        int i,j,l,icoords[MAXD];
        int imin,imax,jmin,jmax;
        COMPONENT comp;
        double aII;
        int num_nb;
        GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
        boolean use_neumann_solver = YES;
        PetscInt num_iter = 0;
        double rel_residual = 0.0;
        RECT_GRID *rgr = &topological_grid(front->grid_intfc);
        struct Table *T = table_of_interface(front->grid_intfc);
        COMPONENT *top_comp = T->components;
        int *lbuf = front->rect_grid->lbuf;
        int *ubuf = front->rect_grid->ubuf;
        int *top_gmax = rgr->gmax;
        double *top_h = rgr->h;
        HYPER_SURF *hs;
        int reflect[MAXD];

        PETSc solver;
        solver.Create(ilower, iupper-1, 5, 0);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();
        size = iupper - ilower;
        *max_soln = -HUGE;
        *min_soln = HUGE;

        imin = (lbuf[0] == 0) ? 1 : lbuf[0];
        jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
        imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
        jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            comp = top_comp[index];
            I = ij_to_I[i][j];
            if (I == -1) continue;

            index_nb[0] = d_index2d(i-1,j,top_gmax);
            index_nb[1] = d_index2d(i+1,j,top_gmax);
            index_nb[2] = d_index2d(i,j-1,top_gmax);
            index_nb[3] = d_index2d(i,j+1,top_gmax);
            I_nb[0] = ij_to_I[i-1][j];
            I_nb[1] = ij_to_I[i+1][j];
            I_nb[2] = ij_to_I[i][j-1];
            I_nb[3] = ij_to_I[i][j+1];
            icoords[0] = i;
            icoords[1] = j;

            k0 = k[index];
            k0_old = k_old[index];
            num_nb = 0;
            for (l = 0; l < 4; ++l)
            {
                if (I_nb[l] == -1)
                    index_nb[l] = index;  // homogeneous Neumann B.C.
                else num_nb++;
                k_nb[l] = 0.5*(k0 + k[index_nb[l]]);
                k_nb_old[l] = 0.5*(k0_old + k_old[index_nb[l]]);
                coeff[l] = 2.0/(k_nb[l]+k_nb_old[l])/(top_h[l/2]*top_h[l/2]);
            }

            rhs = source[index];

            aII = 0.0;
            for (l = 0; l < 4; ++l)
            {
                if (num_nb < 2) break;
                if (I_nb[l] != -1)
                {
                    solver.Set_A(I,I_nb[l],coeff[l]);
                    aII += -coeff[l];
                }
                else
                {
                    hs = FT_HyperSurfAtGridCrossing(front,icoords,dir[l]);
                    if (hs != NULL && wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        if (!boundary_state(hs))
                        {
                            /* The boundary condition of phi at preset
                               Dirichlet boundary is set to zero, it then
                               uses non-Neumann poisson solver
                            */
                            use_neumann_solver = NO;
                            aII += -coeff[l];
                        }
                    }
                }
            }
            /*
             * This change reflects the need to treat point with only one
             * interior neighbor (a convex point). Not sure why PETSc cannot
             * handle such case. If we have better understanding, this should
             * be changed back.
             */
            if(num_nb >= 2)
            {
                solver.Set_A(I,I,aII);
            }
            else
            {
                solver.Set_A(I,I,1.0);
            }
            solver.Set_b(I,rhs);
        }
        use_neumann_solver = pp_min_status(use_neumann_solver);

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

        start_clock("Before Petsc Solver");
        if (use_neumann_solver)
        {
            printf("\nUsing Neumann Solver for phi!\n");
            solver.Solve_withPureNeumann();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);
            if(rel_residual > 1e-10)
            {
                printf("\n The solution diverges for phi! The residual "
                       "is %le. Solve again using GMRES!\n",rel_residual);
                solver.Reset_x();
                solver.Solve_withPureNeumann_GMRES();
                solver.GetNumIterations(&num_iter);
                solver.GetFinalRelativeResidualNorm(&rel_residual);
            }

        }
        else
        {
            printf("\nUsing non-Neumann Solver for phi!\n");
            solver.Solve();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

            if(rel_residual > 1e-10)
            {
                printf("\n The solution diverges for phi! The residual "
                       "is %le. Solve again using GMRES!\n",rel_residual);
                solver.Reset_x();
                solver.Solve_GMRES();
                solver.GetNumIterations(&num_iter);
                solver.GetFinalRelativeResidualNorm(&rel_residual);
            }

        }
        stop_clock("After Petsc Solver");

        double *x;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("In poisson_solver_vd(): "
                        "num_iter = %d, rel_residual = %le \n",
                        num_iter, rel_residual);

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index2d(i,j,top_gmax);
            I = ij_to_I[i][j];
            soln[index] = x[I-ilower];
            if (*max_soln < soln[index]) *max_soln = soln[index];
            if (*min_soln > soln[index]) *min_soln = soln[index];
        }
        reflect[0] = reflect[1] = YES;
        FT_ParallelExchGridArrayBuffer(soln,front,reflect);
        pp_global_max(max_soln,1);
        pp_global_min(min_soln,1);

        if(debugging("step_size"))
        {
            printf("\nThe max value of phi is %.16g\n",*max_soln);
            printf("\nThe min value of phi is %.16g\n",*min_soln);
        }

        FT_FreeThese(1,x);
}       /* end poisson_solver2d_vd */


void poisson_solver2d_MacPhi_vd(
        Front *front,       /* front structure containing interface geometry */
        int ilower,
        int iupper,
        int **ij_to_I,
        double *source,     /* source term in poisson's equation */
        double *k,          /* coefficient function in poisson's equation */
        double *soln,       /* solution of the poisson's equation */
        double *max_soln,   /* maximum value of the solution */
        double *min_soln)   /* minimum value of the solution */
{
        int index,index_nb[4],size;
        double k0,k_nb[4];
        double rhs,coeff[4];
        int I,I_nb[4];
        int i,j,l,icoords[MAXD];
        int imin,imax,jmin,jmax;
        COMPONENT comp;
        double aII;
        int num_nb;
        GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
        boolean use_neumann_solver = YES;
        PetscInt num_iter = 0;
        double rel_residual = 0.0;
        RECT_GRID *rgr = &topological_grid(front->grid_intfc);
        struct Table *T = table_of_interface(front->grid_intfc);
        COMPONENT *top_comp = T->components;
        int *lbuf = front->rect_grid->lbuf;
        int *ubuf = front->rect_grid->ubuf;
        int *top_gmax = rgr->gmax;
        double *top_h = rgr->h;
        HYPER_SURF *hs;
        int reflect[MAXD];

        PETSc solver;
        solver.Create(ilower, iupper-1, 5, 0);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();
        size = iupper - ilower;
        *max_soln = -HUGE;
        *min_soln = HUGE;

        imin = (lbuf[0] == 0) ? 1 : lbuf[0];
        jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
        imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
        jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            comp = top_comp[index];
            I = ij_to_I[i][j];
            if (I == -1) continue;

            index_nb[0] = d_index2d(i-1,j,top_gmax);
            index_nb[1] = d_index2d(i+1,j,top_gmax);
            index_nb[2] = d_index2d(i,j-1,top_gmax);
            index_nb[3] = d_index2d(i,j+1,top_gmax);
            I_nb[0] = ij_to_I[i-1][j];
            I_nb[1] = ij_to_I[i+1][j];
            I_nb[2] = ij_to_I[i][j-1];
            I_nb[3] = ij_to_I[i][j+1];
            icoords[0] = i;
            icoords[1] = j;

            k0 = k[index];
            num_nb = 0;
            for (l = 0; l < 4; ++l)
            {
                if (I_nb[l] == -1)
                    index_nb[l] = index;
                else num_nb++;
                k_nb[l] = 0.5*(k0 + k[index_nb[l]]);
                coeff[l] = 1.0/k_nb[l]/(top_h[l/2]*top_h[l/2]);
            }

            rhs = source[index];

            aII = 0.0;
            for (l = 0; l < 4; ++l)
            {
                if (num_nb < 2) break;
                if (I_nb[l] != -1)
                {
                    solver.Set_A(I,I_nb[l],coeff[l]);
                    aII += -coeff[l];
                }
                else
                {
                    hs = FT_HyperSurfAtGridCrossing(front,icoords,dir[l]);
                    if (hs != NULL && wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        if (!boundary_state(hs))
                        {
                            /* The boundary condition of phi at preset
                               Dirichlet boundary is set to zero, it then
                               uses non-Neumann poisson solver
                            */
                            use_neumann_solver = NO;
                            aII += -coeff[l];
                        }
                    }
                }
            }
            /*
             * This change reflects the need to treat point with only one
             * interior neighbor (a convex point). Not sure why PETSc cannot
             * handle such case. If we have better understanding, this should
             * be changed back.
             */
            if(num_nb >= 2)
            {
                solver.Set_A(I,I,aII);
            }
            else
            {
                solver.Set_A(I,I,1.0);
            }
            solver.Set_b(I,rhs);
        }
        use_neumann_solver = pp_min_status(use_neumann_solver);

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

        start_clock("Before Petsc Solver");
        if (use_neumann_solver)
        {
            printf("\nUsing Neumann Solver for phi_mac!\n");
            solver.Solve_withPureNeumann();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);
            if(rel_residual > 1e-10)
            {
                printf("\n The solution diverges for phi_mac! The residual "
                       "is %le. Solve again using GMRES!\n",rel_residual);
                solver.Reset_x();
                solver.Solve_withPureNeumann_GMRES();
                solver.GetNumIterations(&num_iter);
                solver.GetFinalRelativeResidualNorm(&rel_residual);
            }

        }
        else
        {
            printf("\nUsing non-Neumann Solver for phi_mac!\n");
            solver.Solve();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

            if(rel_residual > 1e-10)
            {
                printf("\n The solution diverges for phi_mac! The residual "
                       "is %le. Solve again using GMRES!\n",rel_residual);
                solver.Reset_x();
                solver.Solve_GMRES();
                solver.GetNumIterations(&num_iter);
                solver.GetFinalRelativeResidualNorm(&rel_residual);
            }

        }
        stop_clock("After Petsc Solver");

        double *x;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("In poisson_solver2d_MacPhi_vd(): "
                        "num_iter = %d, rel_residual = %le \n",
                        num_iter, rel_residual);

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index2d(i,j,top_gmax);
            I = ij_to_I[i][j];
            soln[index] = x[I-ilower];
            if (*max_soln < soln[index]) *max_soln = soln[index];
            if (*min_soln > soln[index]) *min_soln = soln[index];
        }
        reflect[0] = YES; reflect[1] = YES;
        FT_ParallelExchGridArrayBuffer(soln,front,reflect);
        pp_global_max(max_soln,1);
        pp_global_min(min_soln,1);

        if(debugging("step_size"))
        {
            printf("\nThe max value of phi_MAC is %.16g\n",*max_soln);
            printf("\nThe min value of phi_MAC is %.16g\n",*min_soln);
        }

        FT_FreeThese(1,x);
}       /* end poisson_solver2d_MacPhi_vd */


void poisson_solver3d_vd(
        Front *front,       /* front structure containing interface geometry */
        int ilower,
        int iupper,
        int ***ijk_to_I,
        double *source,     /* source term in poisson's equation */
        double *kk,          /* coefficient function in poisson's equation */
        double *kk_old,      /* coefficient function in poisson's equation */
        double *soln,       /* solution of the poisson's equation */
        double *max_soln,   /* maximum value of the solution */
        double *min_soln)   /* minimum value of the solution */
{
        int index,index_nb[6],size;
        double k0,k_nb[6],k0_old,k_nb_old[6];
        double rhs,coeff[6];
        int I,I_nb[6];
        int i,j,k,l,icoords[MAXD];
        int imin,imax,jmin,jmax,kmin,kmax;
        COMPONENT comp;
        double aII;
        int num_nb;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
        boolean use_neumann_solver = YES;
        PetscInt num_iter = 0;
        double rel_residual = 0.0;
        RECT_GRID *rgr = &topological_grid(front->grid_intfc);
        struct Table *T = table_of_interface(front->grid_intfc);
        COMPONENT *top_comp = T->components;
        int *lbuf = front->rect_grid->lbuf;
        int *ubuf = front->rect_grid->ubuf;
        int *top_gmax = rgr->gmax;
        double *top_h = rgr->h;
        HYPER_SURF *hs;
        int reflect[MAXD];

        PETSc solver;
        solver.Create(ilower, iupper-1, 7, 7);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();
        size = iupper - ilower;
        *max_soln = -HUGE;
        *min_soln = HUGE;

        imin = (lbuf[0] == 0) ? 1 : lbuf[0];
        jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
        kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
        imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
        jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
        kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            comp = top_comp[index];
            I = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;

            k0 = kk[index];
            k0_old = kk_old[index];
            num_nb = 0;
            for (l = 0; l < 6; ++l)
            {
                if (I_nb[l] == -1)
                    index_nb[l] = index;  //homogeneous Neumann B.C.
                else num_nb++;
                coeff[l] = (1.0/k0+1.0/kk[index_nb[l]]+1.0/k0_old+1.0/kk_old[index_nb[l]])/4.0/(top_h[l/2]*top_h[l/2]);
            }

            rhs = source[index];

            aII = 0.0;
            for (l = 0; l < 6; ++l)
            {
                if (num_nb < 2) break;
                if (I_nb[l] != -1)
                {
                    solver.Set_A(I,I_nb[l],coeff[l]);
                    aII += -coeff[l];
                }
                else
                {
                    hs = FT_HyperSurfAtGridCrossing(front,icoords,dir[l]);
                    if (hs != NULL && wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        if (!boundary_state(hs))
                        {
                            /* The boundary condition of phi at preset
                               Dirichlet boundary is set to zero, it then
                               uses non-Neumann poisson solver
                            */
                            use_neumann_solver = NO;
                            aII += -coeff[l];
                        }
                    }
                    else //homogeneous Neumann B.C., do nothing!
                    {
                    }
                }
            }
            //printf("aII = %e i = %d j = %d k = %d\n", aII, i, j, k);
            /*
             * This change reflects the need to treat point with only one
             * interior neighbor (a convex point). Not sure why PETSc cannot
             * handle such case. If we have better understanding, this should
             * be changed back.
             */
            if(num_nb >= 2)
            {
                solver.Set_A(I,I,aII);
            }
            else
            {
                solver.Set_A(I,I,1.0);
            }
            solver.Set_b(I,rhs);
        }
        use_neumann_solver = pp_min_status(use_neumann_solver);

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

        start_clock("Before Petsc Solver");
        if (use_neumann_solver)
        {
            printf("\nUsing Neumann Solver for phi!\n");
            solver.Solve_withPureNeumann();

/*
            if (debugging("PETSc"))
            {
                double max, min;
                solver.GetExtremeSingularValues(&max, &min); //max=min=-1, only can be used for GMRES
                (void) printf("The max singular value of A = %lf in poisson_solver3d_vd\n", max);
                (void) printf("The min singular value of A = %lf in poisson_solver3d_vd\n", min);
                if ( min != 0.0)
                    (void) printf("The Cond Num of A = %lf in poisson_solver3d_vd\n", max/min);
            }
*/

            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);
            if(rel_residual > 1e-10)
            {
                printf("\n The solution diverges for phi! The residual "
                       "is %le. Solve again using GMRES!\n",rel_residual);
                solver.Reset_x();
                solver.Solve_withPureNeumann_GMRES();

                if (debugging("PETSc"))
                {
                    double max, min;
                    solver.GetExtremeSingularValues(&max, &min);
                    (void) printf("The max singular value of A = %lf in poisson_solver3d_vd\n", max);
                    (void) printf("The min singular value of A = %lf in poisson_solver3d_vd\n", min);
                    if ( min != 0.0)
                        (void) printf("The Cond Num of A = %lf in poisson_solver3d_vd\n", max/min);
                }

                solver.GetNumIterations(&num_iter);
                solver.GetFinalRelativeResidualNorm(&rel_residual);
            }

        }
        else
        {
            printf("\nUsing non-Neumann Solver for phi!\n");
            solver.Solve();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

            if(rel_residual > 1e-10)
            {
                printf("\n The solution diverges for phi! The residual "
                       "is %le. Solve again using GMRES!\n",rel_residual);
                solver.Reset_x();
                solver.Solve_GMRES();
                solver.GetNumIterations(&num_iter);
                solver.GetFinalRelativeResidualNorm(&rel_residual);
            }

        }
        stop_clock("After Petsc Solver");

        double *x;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("In poisson_solver_vd(): "
                        "num_iter = %d, rel_residual = %le \n",
                        num_iter, rel_residual);

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            I = ijk_to_I[i][j][k];
            soln[index] = x[I-ilower];
            if (*max_soln < soln[index]) *max_soln = soln[index];
            if (*min_soln > soln[index]) *min_soln = soln[index];
        }
        reflect[0] = reflect[1] = reflect[2] = YES;
        FT_ParallelExchGridArrayBuffer(soln,front,reflect);
        pp_global_max(max_soln,1);
        pp_global_min(min_soln,1);

        if(debugging("step_size"))
        {
            printf("\nThe max value of phi is %.16g\n",*max_soln);
            printf("\nThe min value of phi is %.16g\n",*min_soln);
        }

        FT_FreeThese(1,x);
} /* end poisson_solver3d_vd */


void poisson_solver3d_MacPhi_vd(
        Front *front,       /* front structure containing interface geometry */
        int ilower,
        int iupper,
        int ***ijk_to_I,
        double *source,     /* source term in poisson's equation */
        double *kk,         /* coefficient function in poisson's equation */
        double *soln,       /* solution of the poisson's equation */
        double *max_soln,   /* maximum value of the solution */
        double *min_soln)   /* minimum value of the solution */
{
        int index,index_nb[6],size;
        double k0,k_nb[6];
        double rhs,coeff[6];
        int I,I_nb[6];
        int i,j,k,l,icoords[MAXD];
        int imin,imax,jmin,jmax,kmin,kmax;
        COMPONENT comp;
        double aII;
        int num_nb;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
        boolean use_neumann_solver = YES;
        PetscInt num_iter = 0;
        double rel_residual = 0;
        RECT_GRID *rgr = &topological_grid(front->grid_intfc);
        struct Table *T = table_of_interface(front->grid_intfc);
        COMPONENT *top_comp = T->components;
        int *lbuf = front->rect_grid->lbuf;
        int *ubuf = front->rect_grid->ubuf;
        int *top_gmax = rgr->gmax;
        double *top_h = rgr->h;
        HYPER_SURF *hs;
        int reflect[MAXD];
        reflect[0] = reflect[1] = reflect[2] = YES;

        if (debugging("storage"))
            printf("enter poisson_solver3d_MacPhi_vd()\n");
        if (debugging("storage"))
        {
            char s[100];
            sprintf(s,"Storage before create PETSc solver");
            print_storage(s,"storage");
        }

        PETSc solver;

        if (debugging("storage"))
        {
            char s[100];
            sprintf(s,"Storage after create PETSc solver");
            print_storage(s,"storage");
        }

        if (debugging("storage"))
            printf("create PETSc solver in poisson_solver3d_MacPhi_vd()\n");

        solver.Create(ilower, iupper-1, 7, 7);

        if (debugging("storage"))
        {
            char s[100];
            sprintf(s,"Storage after create sth");
            print_storage(s,"storage");
        }
        if (debugging("storage"))
            printf("create sth in poisson_solver3d_MacPhi_vd()\n");

        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();
        size = iupper - ilower;
        *max_soln = -HUGE;
        *min_soln = HUGE;

        imin = (lbuf[0] == 0) ? 1 : lbuf[0];
        jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
        kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
        imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
        jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
        kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];

        if (debugging("storage"))
            printf("enter for-loop in poisson_solver3d_MacPhi_vd()\n");

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            comp = top_comp[index];
            I = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;

            k0 = kk[index];
            num_nb = 0;
            for (l = 0; l < 6; ++l)
            {
                if (I_nb[l] == -1)
                    index_nb[l] = index;
                else num_nb++;
                if (k0==0 || kk[ index_nb[l] ]==0)	clean_up(ERROR);
                k_nb[l] = 2.0/(1.0/k0 + 1.0/kk[ index_nb[l] ]);
                coeff[l] = 1.0/k_nb[l]/(top_h[l/2]*top_h[l/2]);
            }
            rhs = source[index];

            aII = 0;
            for (l = 0; l < 6; ++l)
            {
                if (num_nb < 2) break;
                if (I_nb[l] != -1)
                {
                    solver.Set_A(I,I_nb[l],coeff[l]);
                    aII += -coeff[l];
                }
                else
                {
                    hs = FT_HyperSurfAtGridCrossing(front,icoords,dir[l]);
                    if (hs != NULL && wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        if (!boundary_state(hs))
                        {
                            /* The boundary condition of phi at present
                               Dirichlet boundary is set to zero, it then
                               uses non-Neumann poisson solver
                            */
                            use_neumann_solver = NO;
                            aII += -coeff[l];
                            printf(" HZ Oops! This is not supposed to be printed out for NEUMANN BOUNDARY CONDITION.\n");
                        }
                    }
                    else // homogeneous Neumann B.C. for phi_mac, do nothing!
                    {
                        if (wave_type(hs) == NEUMANN_BOUNDARY)
                        {
                            printf("BOUNDARY CONDITION is %d\n", wave_type(hs));
                            printf("NEUMANN_BOUNDARY in %s.\n", __func__);
                        }
                        if (wave_type(hs) == REFLECTION_BOUNDARY)
                        {
                            printf("BOUNDARY CONDITION is %d\n", wave_type(hs));
                            printf("REFLECTION_BOUNDARY in %s.\n", __func__);
                        }
                    }
                }
            }

            /*
             * This change reflects the need to treat point with only one
             * interior neighbor (a convex point). Not sure why PETSc cannot
             * handle such case. If we have better understanding, this should
             * be changed back.
             */
            if(num_nb >= 2)
            {
                solver.Set_A(I,I,aII);
            }
            else
            {
                solver.Set_A(I,I,1.0);
            }
            solver.Set_b(I,rhs);
        }
        use_neumann_solver = pp_min_status(use_neumann_solver);

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

        start_clock("Before Petsc Solver");
        if (use_neumann_solver)
        {
            printf("\nUsing Neumann Solver for phi_mac!\n");
            solver.Solve_withPureNeumann();

/*
            if (debugging("PETSc"))
            {
                double max, min;
                solver.GetExtremeSingularValues(&max, &min); //max=min=-1, only can be used for GMRES
                (void) printf("The max singular value of A = %lf in poisson_solver3d_MacPhi_vd\n", max);
                (void) printf("The min singular value of A = %lf in poisson_solver3d_MacPhi_vd\n", min);
                if ( min != 0.0)
                    (void) printf("The Cond Num of A = %lf in poisson_solver3d_MacPhi_vd\n", max/min);
            }
*/

            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);
            if(rel_residual > 1e-10)
            {
                printf("\n The solution diverges for phi_mac! The residual "
                       "is %le. Solve again using GMRES!\n",rel_residual);
                solver.Reset_x();
                solver.Solve_withPureNeumann_GMRES();

                if (debugging("PETSc"))
                {
                    double max, min;
                    solver.GetExtremeSingularValues(&max, &min);
                    (void) printf("The max singular value of A = %lf in poisson_solver3d_MacPhi_vd\n", max);
                    (void) printf("The min singular value of A = %lf in poisson_solver3d_MacPhi_vd\n", min);
                    if ( min != 0.0)
                        (void) printf("The Cond Num of A = %lf in poisson_solver3d_MacPhi_vd\n", max/min);
                }

                solver.GetNumIterations(&num_iter);
                solver.GetFinalRelativeResidualNorm(&rel_residual);
            }

        }
        else
        {
            printf("\nUsing non-Neumann Solver for phi_mac!\n");
            solver.Solve();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

            if(rel_residual > 1e-10)
            {
                printf("\n The solution diverges for phi_mac! The residual "
                       "is %le. Solve again using GMRES!\n",rel_residual);
                solver.Reset_x();
                solver.Solve_GMRES();
                solver.GetNumIterations(&num_iter);
                solver.GetFinalRelativeResidualNorm(&rel_residual);
            }

        }
        stop_clock("After Petsc Solver");

        double *x;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
        solver.Get_x(x);

        if (debugging("PETSc")) {
            (void) printf("In poisson_solver3d_MacPhi_vd(): "
                          "num_iter = %d, rel_residual = %le \n",
                          num_iter, rel_residual);
	}

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            I = ijk_to_I[i][j][k];
            soln[index] = x[I-ilower];
            if (*max_soln < soln[index]) *max_soln = soln[index];
            if (*min_soln > soln[index]) *min_soln = soln[index];
        }
        FT_ParallelExchGridArrayBuffer(soln,front,reflect);
        pp_global_max(max_soln,1);
        pp_global_min(min_soln,1);

        if (debugging("step_size"))
        {
            printf("\nThe max value of phi_MAC is %.16g\n",*max_soln);
            printf("\nThe min value of phi_MAC is %.16g\n",*min_soln);
        }

        FT_FreeThese(1,x);

        if (debugging("storage"))
        {
            char s[100];
            sprintf(s,"Storage after poisson_solver3d_MacPhi_vd");
            print_storage(s,"storage");
        }
} /* end poisson_solver3d_MacPhi_vd */


void poisson_solver3d_Expand_vd(
        Front *front,       /* front structure containing interface geometry */
        int ilower,
        int iupper,
        int ***ijk_to_I,
        double *source,     /* source term in poisson's equation */
        double *kk,         /* coefficient function in poisson's equation */
        double *kk_old,     /* coefficient function in poisson's equation */
        double *soln,       /* solution of the poisson's equation */
        double *max_soln,   /* maximum value of the solution */
        double *min_soln)   /* minimum value of the solution */
{
        int index,index_nb[12],size;
        double k0,k_nb[6],k0_old,k_nb_old[6];
        double rhs,coeff[6];
        int I,I_nb[12];
        int i,j,k,l,icoords[MAXD];
        int imin,imax,jmin,jmax,kmin,kmax;
        COMPONENT comp;
        double aII;
        int num_nb;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
        boolean use_neumann_solver = YES;
        PetscInt num_iter = 0;
        double rel_residual = 0.0;
        RECT_GRID *rgr = &topological_grid(front->grid_intfc);
        struct Table *T = table_of_interface(front->grid_intfc);
        COMPONENT *top_comp = T->components;
        int *lbuf = front->rect_grid->lbuf;
        int *ubuf = front->rect_grid->ubuf;
        int *top_gmax = rgr->gmax;
        double *top_h = rgr->h;
        HYPER_SURF *hs;
        int reflect[MAXD];

        PETSc solver;
        solver.Create(ilower, iupper-1, 7, 7);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();
        size = iupper - ilower;
        *max_soln = -HUGE;
        *min_soln = HUGE;

        imin = (lbuf[0] == 0) ? 1 : lbuf[0];
        jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
        kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
        imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
        jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
        kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            comp = top_comp[index];
            I = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);
            index_nb[6] = d_index3d(i-2,j,k,top_gmax);
            index_nb[7] = d_index3d(i+2,j,k,top_gmax);
            index_nb[8] = d_index3d(i,j-2,k,top_gmax);
            index_nb[9] = d_index3d(i,j+2,k,top_gmax);
            index_nb[10] = d_index3d(i,j,k-2,top_gmax);
            index_nb[11] = d_index3d(i,j,k+2,top_gmax);

            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];
            I_nb[6] = ijk_to_I[i-2][j][k];
            I_nb[7] = ijk_to_I[i+2][j][k];
            I_nb[8] = ijk_to_I[i][j-2][k];
            I_nb[9] = ijk_to_I[i][j+2][k];
            I_nb[10] = ijk_to_I[i][j][k-2];
            I_nb[11] = ijk_to_I[i][j][k+2];

            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;

            k0 = kk[index];
            k0_old = kk_old[index];

            rhs = source[index];
            aII = 0.0;

            // WEST, EAST, SOUTH, and NORTH
            for (l = 6; l < 10; ++l)
            {
                coeff[l] = (1.0/kk[index_nb[l-6]] + 1.0/kk_old[index_nb[l-6]])/2.0/4.0/(top_h[(l-6)/2]*top_h[(l-6)/2]);
                solver.Set_A(I,I_nb[l],coeff[l]);
                aII += -coeff[l];
            }

            // LOWER
            if (I_nb[4] == -1)
            {
                coeff[10] = (1.0/kk[index] + 1.0/kk_old[index])/2.0/4.0/(top_h[2]*top_h[2]);
                solver.Set_A(I,I_nb[5],coeff[10]);
            }
            else if (I_nb[4] == 1 && I_nb[10] == -1)
            {
                coeff[10] = (1.0/kk[index_nb[4]] + 1.0/kk_old[index_nb[4]])/2.0/4.0/(top_h[2]*top_h[2]);
                solver.Set_A(I,I_nb[4],coeff[10]);
            }
            else
            {
                coeff[10] = (1.0/kk[index_nb[4]] + 1.0/kk_old[index_nb[4]])/2.0/4.0/(top_h[2]*top_h[2]);
                solver.Set_A(I,I_nb[10],coeff[10]);
            }
            aII += -coeff[10];

            // UPPER
            if (I_nb[5] == -1)
            {
                coeff[11] = (1.0/kk[index] + 1.0/kk_old[index])/2.0/4.0/(top_h[2]*top_h[2]);
                solver.Set_A(I,I_nb[4],coeff[11]);
            }
            else if (I_nb[5] == 1 && I_nb[11] == -1)
            {
                coeff[11] = (1.0/kk[index_nb[5]] + 1.0/kk_old[index_nb[5]])/2.0/4.0/(top_h[2]*top_h[2]);
                solver.Set_A(I,I_nb[5],coeff[11]);
            }
            else
            {
                coeff[11] = (1.0/kk[index_nb[5]] + 1.0/kk_old[index_nb[5]])/2.0/4.0/(top_h[2]*top_h[2]);
                solver.Set_A(I,I_nb[11],coeff[11]);
            }
            aII += -coeff[11];

            /*
             * This change reflects the need to treat point with only one
             * interior neighbor (a convex point). Not sure why PETSc cannot
             * handle such case. If we have better understanding, this should
             * be changed back.
             */
            solver.Set_A(I,I,aII);
            solver.Set_b(I,rhs);
        }
        use_neumann_solver = pp_min_status(use_neumann_solver);

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

        start_clock("Before Petsc Solver");
        if (use_neumann_solver)
        {
            printf("\nUsing Neumann Solver for phi!\n");
            solver.Solve_withPureNeumann();

/*
            if (debugging("PETSc"))
            {
                double max, min;
                solver.GetExtremeSingularValues(&max, &min); //max=min=-1, only can be used for GMRES
                (void) printf("The max singular value of A = %lf in poisson_solver3d_Expand_vd\n", max);
                (void) printf("The min singular value of A = %lf in poisson_solver3d_Expand_vd\n", min);
                if ( min != 0.0)
                    (void) printf("The Cond Num of A = %lf in poisson_solver3d_vd\n", max/min);
            }
*/

            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);
            if(rel_residual > 1e-10)
            {
                printf("\n The solution diverges for phi! The residual "
                       "is %le. Solve again using GMRES!\n",rel_residual);
                solver.Reset_x();
                solver.Solve_withPureNeumann_GMRES();

                if (debugging("PETSc"))
                {
                    double max, min;
                    solver.GetExtremeSingularValues(&max, &min);
                    (void) printf("The max singular value of A = %lf in poisson_solver3d_Expand_vd\n", max);
                    (void) printf("The min singular value of A = %lf in poisson_solver3d_Expand_vd\n", min);
                    if ( min != 0.0)
                        (void) printf("The Cond Num of A = %lf in poisson_solver3d_vd\n", max/min);
                }

                solver.GetNumIterations(&num_iter);
                solver.GetFinalRelativeResidualNorm(&rel_residual);
            }

        }
        else
        {
            printf("\nUsing non-Neumann Solver for phi!\n");
            solver.Solve();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

            if(rel_residual > 1e-10)
            {
                printf("\n The solution diverges for phi! The residual "
                       "is %le. Solve again using GMRES!\n",rel_residual);
                solver.Reset_x();
                solver.Solve_GMRES();
                solver.GetNumIterations(&num_iter);
                solver.GetFinalRelativeResidualNorm(&rel_residual);
            }

        }
        stop_clock("After Petsc Solver");

        double *x;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("In poisson_solver3d_Expand_vd(): "
                        "num_iter = %d, rel_residual = %le \n",
                        num_iter, rel_residual);

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            I = ijk_to_I[i][j][k];
            soln[index] = x[I-ilower];
            if (*max_soln < soln[index]) *max_soln = soln[index];
            if (*min_soln > soln[index]) *min_soln = soln[index];
        }
        reflect[0] = reflect[1] = reflect[2] = YES;
        FT_ParallelExchGridArrayBuffer(soln,front,reflect);
        pp_global_max(max_soln,1);
        pp_global_min(min_soln,1);

        if(debugging("step_size"))
        {
            printf("\nThe max value of phi is %.16g\n",*max_soln);
            printf("\nThe min value of phi is %.16g\n",*min_soln);
        }

        FT_FreeThese(1,x);
} /* end poisson_solver3d_Expand_vd */
