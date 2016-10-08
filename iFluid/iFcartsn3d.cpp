/*******************************************************************
  			iFcartsn3d.cpp
 *******************************************************************/
#include "iFluid.h"
#include "solver.h"
#define SLIP 1

//--------------------------------------------------------------------------
// 		   Incompress_Solver_Smooth_3D_Cartesian
//--------------------------------------------------------------------------

void Incompress_Solver_Smooth_3D_Cartesian::computeNewDensity_vd(int flag)
{
        int i, j, k, index;
        COMPONENT comp;
        double density, max_tmp, min_tmp, nu, max_rho, min_rho;
        int indmax[3],indmin[3];

        max_density = -1;
        min_density = HUGE;
        max_rho = (m_rho[0] > m_rho[1]) ? m_rho[0] : m_rho[1];
        min_rho = (m_rho[0] + m_rho[1]) - max_rho;

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            comp = top_comp[index];
            if (!ifluid_comp(comp))
            {
                cell_center[index].m_state.m_rho = 0;
                cell_center[index].m_state.m_rho_old = 0;
                continue;
            }

            if (!flag) //flag==0
            {
                cell_center[index].m_state.m_rho_old = cell_center[index].m_state.m_rho;
                cell_center[index].m_state.m_rho -= m_dt*cell_center[index].m_state.m_rho_adv;
                if (cell_center[index].m_state.m_rho > max_rho)	cell_center[index].m_state.m_rho = max_rho;
                if (cell_center[index].m_state.m_rho < min_rho) cell_center[index].m_state.m_rho = min_rho;
            }
            else //flag==1
            {
                cell_center[index].m_state.m_rho = cell_center[index].m_state.m_rho_old -
                                                   m_dt*cell_center[index].m_state.m_rho_adv;
                if (cell_center[index].m_state.m_rho > max_rho) cell_center[index].m_state.m_rho = max_rho;
                if (cell_center[index].m_state.m_rho < min_rho) cell_center[index].m_state.m_rho = min_rho;
                //update mu
                nu = (m_mu[0]+m_mu[1])/(m_rho[0]+m_rho[1]);
                cell_center[index].m_state.m_mu = cell_center[index].m_state.m_rho*nu;
                cell_center[index].m_state.m_mu_old = cell_center[index].m_state.m_rho_old*nu;
            }

            density = fabs(cell_center[index].m_state.m_rho);
            if (density > max_density)
            {
                max_density = density;
                indmax[0] = i;
                indmax[1] = j;
                indmax[2] = k;
            }
            if (density < min_density)
            {
                min_density = density;
                indmin[0] = i;
                indmin[1] = j;
                indmin[2] = k;
            }
        }
        max_tmp = max_density;
        min_tmp = min_density;

        pp_global_max(&max_density,1);
        pp_global_min(&min_density,1);

	if (debugging("step_size"))
	{
            printf("local max_density after computeNewDensity_vd(%d) "
                   "in cell(%d, %d, %d) of node #%d is: %lf\n", flag,indmax[0],indmax[1],indmax[2],pp_mynode(),max_tmp);
            printf("local min_density after computeNewDensity_vd(%d) "
                   "in cell(%d, %d, %d) of node #%d is: %lf\n", flag,indmin[0],indmin[1],indmin[2],pp_mynode(),min_tmp);
            if (max_tmp == max_density)
                printf("max_density (locates in node #%d) is: %lf\n", pp_mynode(),max_density);
            if (min_tmp == min_density)
        	printf("min_density (locates in node #%d) is: %lf\n", pp_mynode(),min_density);
	}

        //scatter states
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_rho_old;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_rho_old = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_rho;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_rho = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_mu_old;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_mu_old = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_mu;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_mu = array[index];
        }
} /* end computeNewDensity_vd */


void Incompress_Solver_Smooth_3D_Cartesian::computeNewDensity_MAC_constRho_vd(void)
{
        int i, j, k, index;
        COMPONENT comp;

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            comp = top_comp[index];
            if (!ifluid_comp(comp))
            {
                cell_center[index].m_state.m_rho = 0.0;
                cell_center[index].m_state.m_rho_old = 0.0;
                continue;
            }

            if (m_rho[0] == m_rho[1])
                cell_center[index].m_state.m_rho = cell_center[index].m_state.m_rho_old = m_rho[0];
            else
            {
                (void) printf("Density is not constant in computeNewDensity_MAC_constRho_vd()!");
                assert(false);
            }
        }
} /* end computeNewDensity_MAC_constRho_vd */


void Incompress_Solver_Smooth_3D_Cartesian::computeNewConcentration_vd(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
        int indmax[3],indmin[3];
        double coords[MAXD],crx_coords[MAXD];
        double coeff0[6],coeff1[6];
        double Dcoef,rhs;
        double rho_edge[6], rho_old_edge[6], rho_mid, rho0, rho0_old;
        double c_nb[6],c_center;
        int flag[6];  //denote whether there is boundary
        L_STATE state;
        int i,j,k,nb,icoords[MAXD];
        INTERFACE *intfc = front->interf;
        double concentration,max_tmp,min_tmp;
        double *x;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
        POINTER intfc_state;
        HYPER_SURF *hs;
        PetscInt num_iter;
        double rel_residual;
        boolean use_neumann_solver = YES;

        max_concentration = -1;
        min_concentration = HUGE;
        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

        PETSc solver;
        solver.Create(ilower, iupper-1, 7, 7);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index = d_index3d(i,j,k,top_gmax);
            comp = top_comp[index];

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

            Dcoef    = cell_center[index].m_state.m_Dcoef;
            rho0     = cell_center[index].m_state.m_rho;
            rho0_old = cell_center[index].m_state.m_rho_old;
            rho_mid  = 0.5*(rho0 + rho0_old);
            c_center = cell_center[index].m_state.m_c;

            for (nb = 0; nb < 6; nb++)
            {
                if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                {
                    if (wave_type(hs) == DIRICHLET_BOUNDARY ||
                        wave_type(hs) == NEUMANN_BOUNDARY)
                    {
                        flag[nb] = 1;
                        c_nb[nb] = c_center;
                        rho_edge[nb] = rho0;
                        rho_old_edge[nb] = rho0_old;
                    }
                    else
                    {
                        flag[nb]         = 0;
                        c_nb[nb]         = cell_center[index_nb[nb]].m_state.m_c;
                        rho_edge[nb]     = 2.0/(1.0/rho0 + 1.0/cell_center[index_nb[nb]].m_state.m_rho);
                        rho_old_edge[nb] = 2.0/(1.0/rho0_old + 1.0/cell_center[index_nb[nb]].m_state.m_rho_old);
                    }
                }
                else
                {
                    flag[nb]         = 0;
                    c_nb[nb]         = cell_center[index_nb[nb]].m_state.m_c;
                    rho_edge[nb]     = 2.0/(1.0/rho0 + 1.0/cell_center[index_nb[nb]].m_state.m_rho);
                    rho_old_edge[nb] = 2.0/(1.0/rho0_old + 1.0/cell_center[index_nb[nb]].m_state.m_rho_old);
                }
            }

            //Setting the coefficients and matrix for c in concentration eqn.
            coeff0[0] = m_dt/2.0/rho_mid * rho_edge[0] * Dcoef/(top_h[0]*top_h[0]); //west
            coeff0[1] = m_dt/2.0/rho_mid * rho_edge[1] * Dcoef/(top_h[0]*top_h[0]); //east
            coeff0[2] = m_dt/2.0/rho_mid * rho_edge[2] * Dcoef/(top_h[1]*top_h[1]); //south
            coeff0[3] = m_dt/2.0/rho_mid * rho_edge[3] * Dcoef/(top_h[1]*top_h[1]); //north
            coeff0[4] = m_dt/2.0/rho_mid * rho_edge[4] * Dcoef/(top_h[2]*top_h[2]); //lower
            coeff0[5] = m_dt/2.0/rho_mid * rho_edge[5] * Dcoef/(top_h[2]*top_h[2]); //upper

            coeff1[0] = m_dt/2.0/rho_mid * rho_old_edge[0] * Dcoef/(top_h[0]*top_h[0]); //west
            coeff1[1] = m_dt/2.0/rho_mid * rho_old_edge[1] * Dcoef/(top_h[0]*top_h[0]); //east
            coeff1[2] = m_dt/2.0/rho_mid * rho_old_edge[2] * Dcoef/(top_h[1]*top_h[1]); //south
            coeff1[3] = m_dt/2.0/rho_mid * rho_old_edge[3] * Dcoef/(top_h[1]*top_h[1]); //north
            coeff1[4] = m_dt/2.0/rho_mid * rho_old_edge[4] * Dcoef/(top_h[2]*top_h[2]); //lower
            coeff1[5] = m_dt/2.0/rho_mid * rho_old_edge[5] * Dcoef/(top_h[2]*top_h[2]); //upper

            solver.Add_A(I,I,1.0);
            rhs = c_center;

            for (nb = 0; nb < 6; nb++)
            {
                if (flag[nb] == 0) //interior
                {
                    solver.Add_A(I,I_nb[nb],-coeff0[nb]);
                    rhs += coeff1[nb]*c_nb[nb];
                }
                else //homogeneous Neumann BC for c
                    coeff0[nb] = coeff1[nb] = 0.0;

                solver.Add_A(I,I,coeff0[nb]);
                rhs -= coeff1[nb]*c_center;
            }

            rhs -= m_dt*cell_center[index].m_state.m_c_adv; //advection term
            solver.Add_b(I, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

        start_clock("Before Petsc solver in computing concentration");

        solver.Solve_GMRES();

        if (debugging("PETSc"))
        {
            double max, min;
            solver.GetExtremeSingularValues(&max, &min);
            (void) printf("The max singular value of A = %lf in computeNewConcentration_vd\n", max);
            (void) printf("The min singular value of A = %lf in computeNewConcentration_vd\n", min);
            if ( min != 0.0)
                (void) printf("The Cond Num of A = %lf in computeNewConcentration_vd\n", max/min);
        }

        solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);

        //if (rel_residual > 1)
        //{
        //    solver.Reset_x();
        //    solver.Solve_GMRES();
        //    solver.GetNumIterations(&num_iter);
        //    solver.GetFinalRelativeResidualNorm(&rel_residual);
        //}

        stop_clock("After Petsc Solver in computing concentration");

        //get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
        {
            printf("\nIncompress_Solver_Smooth_3D_Cartesian::"
                        "computeNewConcentration_vd: "
                        "num_iter = %d, rel_residual = %le. \n",
                        num_iter,rel_residual);
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_c = x[I-ilower];
                concentration = fabs(cell_center[index].m_state.m_c);
                if (concentration > max_concentration)
                {
                    max_concentration = concentration;
                    indmax[0] = i;
                    indmax[1] = j;
                    indmax[2] = k;
                }
                if (concentration < min_concentration)
                {
                    min_concentration = concentration;
                    indmin[0] = i;
                    indmin[1] = j;
                    indmin[2] = k;
                }
            }
            else
            {
                cell_center[index].m_state.m_c = 0.0;
                concentration = fabs(cell_center[index].m_state.m_c);
                if (concentration > max_concentration)
                {
                    max_concentration = concentration;
                    indmax[0] = -1;
                    indmax[1] = -1;
                    indmax[2] = -1;
                }
                if (concentration < min_concentration)
                {
                    min_concentration = concentration;
                    indmin[0] = -1;
                    indmin[1] = -1;
                    indmin[2] = -1;
                }
            }
        }
        max_tmp = max_concentration;
        min_tmp = min_concentration;

        pp_global_max(&max_concentration,1);
        pp_global_min(&min_concentration,1);

	if (debugging("step_size"))
	{
            (void) printf("local max_concentration after computeNewConcentration_vd "
                          "in cell(%d, %d, %d) of node #%d is: %lf\n",indmax[0],indmax[1],indmax[2],pp_mynode(),max_tmp);
            (void) printf("local min_concentration after computeNewConcentration_vd "
                          "in cell(%d, %d, %d) of node #%d is: %lf\n",indmin[0],indmin[1],indmin[2],pp_mynode(),min_tmp);
            if (max_tmp == max_concentration)
                (void) printf("max_concentration (locates in node #%d) is: %lf\n",pp_mynode(),max_concentration);
            if (min_tmp == min_concentration)
                (void) printf("min_concentration (locates in node #%d) is: %lf\n",pp_mynode(),min_concentration);
	}

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            if (cell_center[index].m_state.m_c < 0)	cell_center[index].m_state.m_c = 0;
            if (cell_center[index].m_state.m_c > 1.0)	cell_center[index].m_state.m_c = 1.0;
            array[index] = cell_center[index].m_state.m_c;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_c = array[index];
        }

        FT_FreeThese(1,x);
} /* end computeNewConcentration_vd */


void Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_velocity_vd(void)
{
        COMPONENT comp;
        int index,index_nb[18],indmax[3],size;
        int I,I_nb[18];
        double coords[MAXD],crx_coords[MAXD];
        double coeff0[6],coeff1[6],coeff2[6],coeff_temp1,coeff_temp2;
        double mu[6],mu_edge[6],mu0,rho,rhs;
        double U0_nb[18],U1_nb[18],U2_nb[18],U0_center,U1_center,U2_center;
        double U0_nb_new[6],U1_nb_new[6],U2_nb_new[6];
        int flag[6]; //denote whether this is dirichlet or neumann boundary
        L_STATE state;
        int i,j,k,l,nb,icoords[MAXD];
        INTERFACE *intfc = front->interf;
        double speed,u,v,w,max_u,max_v,max_w,max_tmp,sum_div,value;
        double *x;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
        double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};
        POINTER intfc_state;
        HYPER_SURF *hs;
        PetscInt num_iter;
        double rel_residual;
        double **vel = iFparams->field->vel;

        max_speed = 0.0;
        max_u = 0.0;
        max_v = 0.0;
        max_w = 0.0;
        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 30, 30);
        // 7u + 9v + 9w for the first equation
        // 7v + 9u + 9w for the second equation
        // 7w + 9u + 9v for the third equation

        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index = d_index3d(i,j,k,top_gmax);
            comp = top_comp[index];

            //6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);

            //xy cut neighbours
            index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
            index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
            index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
            index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);

            //yz cut neighbours
            index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
            index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
            index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
            index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);

            //xz cut neighbours
            index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
            index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
            index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
            index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);


            //6 neighbours of the center cell
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];

            //xy cut neighbours
            I_nb[6] = ijk_to_I[i-1][j-1][k];
            I_nb[7] = ijk_to_I[i+1][j-1][k];
            I_nb[8] = ijk_to_I[i+1][j+1][k];
            I_nb[9] = ijk_to_I[i-1][j+1][k];

            //yz cut neighbours
            I_nb[10] = ijk_to_I[i][j-1][k-1];
            I_nb[11] = ijk_to_I[i][j+1][k-1];
            I_nb[12] = ijk_to_I[i][j+1][k+1];
            I_nb[13] = ijk_to_I[i][j-1][k+1];

            //xz cut neighbours
            I_nb[14] = ijk_to_I[i-1][j][k-1];
            I_nb[15] = ijk_to_I[i+1][j][k-1];
            I_nb[16] = ijk_to_I[i+1][j][k+1];
            I_nb[17] = ijk_to_I[i-1][j][k+1];

            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;

            mu0 = cell_center[index].m_state.m_mu;
            rho = 0.5*(cell_center[index].m_state.m_rho
                       + cell_center[index].m_state.m_rho_old); // rho_mid
            U0_center = cell_center[index].m_state.m_U[0];
            U1_center = cell_center[index].m_state.m_U[1];
            U2_center = cell_center[index].m_state.m_U[2];

            for (nb = 0; nb < 6; nb++)
            {
                if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                {
                    if (wave_type(hs) == DIRICHLET_BOUNDARY ||
                        wave_type(hs) == NEUMANN_BOUNDARY)
                    {
                        flag[nb] = 1;
                        //(void) printf("flag[%d] = 1 in cell(%d,%d,%d)!!\n",nb,i,j,k);
                        U0_nb[nb] = U0_nb_new[nb] = 0.0;
                        U1_nb[nb] = U1_nb_new[nb] = 0.0;
                        U2_nb[nb] = U2_nb_new[nb] = 0.0;
                        mu_edge[nb] = mu0;
                    }
                    else
                    {
                        flag[nb] = 0;
                        U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
                        U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
                        U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
                        mu_edge[nb] = 0.5*(mu0 + cell_center[index_nb[nb]].m_state.m_mu);
                    }
                }
                else
                {
                    flag[nb] = 0;
                    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
                    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
                    U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
                    mu_edge[nb] = 0.5*(mu0 + cell_center[index_nb[nb]].m_state.m_mu);
                }
            }

            for (nb = 6; nb < 18; nb++) //cut corner values, interior
            {
                if (I_nb[nb] != -1)
                {
                    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
                    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
                    U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
                }
                //else
                    //(void) printf("I_nb[%d] = -1 in cell(%d,%d,%d)!!\n",nb,i,j,k);
            }

            // source term
            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);

            /************************************************************************/
            //         Setting the coeffecients for U0 in the first equation
            /************************************************************************/

            coeff0[0] = (2.0/3.0)*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
            coeff0[1] = (2.0/3.0)*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
            coeff0[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
            coeff0[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
            coeff0[4] = 0.5*m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
            coeff0[5] = 0.5*m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

            solver.Add_A(I*3,I*3,1.0);
            rhs = U0_center;

            for (nb = 0; nb < 6; nb++)
            {
                if (flag[nb] == 0) //interior faces
                {
                    solver.Add_A(I*3,I_nb[nb]*3,-coeff0[nb]);
                    rhs += coeff0[nb]*U0_nb[nb];
                }
                else //flag[nb] == 1, use homogeneous Dirichlet B.C. for U
                {
                    coeff0[nb] = 2.0*coeff0[nb];
                    rhs += coeff0[nb]*(U0_nb_new[nb] + U0_nb[nb]);
                }
                solver.Add_A(I*3,I*3,coeff0[nb]);
                rhs -= coeff0[nb]*U0_center;
            }

            //set the coefficients for U1 in the first equation (mu*v_x)_y
            //traverse the four corners
            //(i-1/2,j-1/2,k),(i+1/2,j-1/2,k),(i+1/2,j+1/2,k),(i-1/2,j+1/2,k)

            //corner (i-1/2,j-1/2,k)

            coeff_temp1 = mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho;
            coeff_temp2 = (-2.0/3.0)*mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho;
            if (flag[0] == 0 && flag[2] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U1_nb[0]+U1_nb[2]+U1_nb[6]+U1_center)/8.0;
                solver.Add_A(I*3,I*3+1,       -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[0]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[2]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[6]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else assert(false);

            //corner (i+1/2,j-1/2,k)

            coeff_temp1 = -mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho;
            coeff_temp2 = -(-2.0/3.0)*mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho;
            if (flag[2] == 0 && flag[1] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U1_nb[2]+U1_nb[1]+U1_nb[7]+U1_center)/8.0;
                solver.Add_A(I*3,I*3+1,       -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[2]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[1]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[7]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else assert(false);

            //corner (i+1/2,j+1/2,k)

            coeff_temp1 = mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho;
            coeff_temp2 = (-2.0/3.0)*mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho;
            if (flag[1] == 0 && flag[3] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U1_nb[1]+U1_nb[3]+U1_nb[8]+U1_center)/8.0;
                solver.Add_A(I*3,I*3+1,       -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[1]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[3]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[8]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else assert(false);

            //corner (i-1/2,j+1/2,k)

            coeff_temp1 = -mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho;
            coeff_temp2 = -(-2.0/3.0)*mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho;
            if (flag[0] == 0 && flag[3] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U1_nb[3]+U1_nb[0]+U1_nb[9]+U1_center)/8.0;
                solver.Add_A(I*3,I*3+1,       -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[3]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[0]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[9]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else assert(false);


            //set the coefficients for U2 in the first equation (mu*w_x)_z
            //traverse the four corners
            //(i-1/2,j,k-1/2),(i+1/2,j,k-1/2),(i+1/2,j,k+1/2),(i-1/2,j,k+1/2)

            //corner (i-1/2,j,k-1/2)

            coeff_temp1 = mu_edge[4]*m_dt/(top_h[0]*top_h[2])/rho;
            coeff_temp2 = (-2.0/3.0)*mu_edge[0]*m_dt/(top_h[0]*top_h[2])/rho;
            if (flag[0] == 0 && flag[4] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U2_nb[0]+U2_nb[4]+U2_nb[14]+U2_center)/8.0;
                solver.Add_A(I*3,I*3+2,       -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[0]*3+2, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[4]*3+2, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[14]*3+2, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else //do nothing
            {
            }

            //corner (i+1/2,j,k-1/2)

            coeff_temp1 = -mu_edge[4]*m_dt/(top_h[0]*top_h[2])/rho;
            coeff_temp2 = -(-2.0/3.0)*mu_edge[1]*m_dt/(top_h[0]*top_h[2])/rho;
            if (flag[1] == 0 && flag[4] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U2_nb[1]+U2_nb[4]+U2_nb[15]+U2_center)/8.0;
                solver.Add_A(I*3,I*3+2,       -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[1]*3+2, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[4]*3+2, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[15]*3+2, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else //do nothing
            {
            }

            //corner (i+1/2,j,k+1/2)

            coeff_temp1 = mu_edge[5]*m_dt/(top_h[0]*top_h[2])/rho;
            coeff_temp2 = (-2.0/3.0)*mu_edge[1]*m_dt/(top_h[0]*top_h[2])/rho;
            if (flag[1] == 0 && flag[5] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U2_nb[1]+U2_nb[5]+U2_nb[16]+U2_center)/8.0;
                solver.Add_A(I*3,I*3+2,       -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[1]*3+2, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[5]*3+2, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[16]*3+2, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else //do nothing
            {
            }

            //corner (i-1/2,j,k+1/2)

            coeff_temp1 = -mu_edge[5]*m_dt/(top_h[0]*top_h[2])/rho;
            coeff_temp2 = -(-2.0/3.0)*mu_edge[0]*m_dt/(top_h[0]*top_h[2])/rho;
            if (flag[5] == 0 && flag[0] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U2_nb[5]+U2_nb[0]+U2_nb[17]+U2_center)/8.0;
                solver.Add_A(I*3,I*3+2,       -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[5]*3+2, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[0]*3+2, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3,I_nb[17]*3+2, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else //do nothing
            {
            }

            rhs += m_dt*state.m_U[0];
            rhs += m_dt*cell_center[index].m_state.f_surf[0];
            rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;
            rhs -= m_dt*cell_center[index].m_state.m_adv[0];

            solver.Add_b(I*3, rhs);

            /************************************************************************/
            //       Setting the coeffecients for U1 in the second equation
            /************************************************************************/

            coeff1[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
            coeff1[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
            coeff1[2] = (2.0/3.0)*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
            coeff1[3] = (2.0/3.0)*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
            coeff1[4] = 0.5*m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
            coeff1[5] = 0.5*m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

            solver.Add_A(I*3+1,I*3+1,1.0);
            rhs = U1_center;

            for (nb = 0; nb < 6; nb++)
            {
                if (flag[nb] == 0)
                {
                    solver.Add_A(I*3+1,I_nb[nb]*3+1,-coeff1[nb]);
                    rhs += coeff1[nb]*U1_nb[nb];
                }
                else //flag[nb] == 1
                {
                    coeff1[nb] = 2.0*coeff1[nb];
                    rhs += coeff1[nb]*(U1_nb[nb] + U1_nb_new[nb]);
                }
                solver.Add_A(I*3+1,I*3+1,coeff1[nb]);
                rhs -= coeff1[nb]*U1_center;
            }

            //set the coefficients for U0 in the second equation (mu*u_y)_x
            //traverse the four corners
            //(i-1/2,j-1/2,k),(i+1/2,j-1/2,k),(i+1/2,j+1/2,k),(i-1/2,j+1/2,k)

            //corner (i-1/2,j-1/2,k)

            coeff_temp1 = mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho;
            coeff_temp2 = (-2.0/3.0)*mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho;
            if (flag[0] == 0 && flag[2] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U0_nb[0]+U0_nb[2]+U0_nb[6]+U0_center)/8.0;
                solver.Add_A(I*3+1,I*3,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[0]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[2]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[6]*3,  -(coeff_temp1+coeff_temp2)/8.0);
            }
            else assert(false);

            //corner (i+1/2,j-1/2,k)

            coeff_temp1 = -mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho;
            coeff_temp2 = -(-2.0/3.0)*mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho;
            if (flag[1] == 0 && flag[2] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U0_nb[1]+U0_nb[2]+U0_nb[7]+U0_center)/8.0;
                solver.Add_A(I*3+1,I*3,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[1]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[2]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[7]*3,  -(coeff_temp1+coeff_temp2)/8.0);
            }
            else assert(false);

            //corner (i+1/2,j+1/2,k)

            coeff_temp1 = mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho;
            coeff_temp2 = (-2.0/3.0)*mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho;
            if (flag[1] == 0 && flag[3] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U0_nb[1]+U0_nb[3]+U0_nb[8]+U0_center)/8.0;
                solver.Add_A(I*3+1,I*3,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[1]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[3]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[8]*3,  -(coeff_temp1+coeff_temp2)/8.0);
            }
            else assert(false);

            //corner (i-1/2,j+1/2,k)

            coeff_temp1 = -mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho;
            coeff_temp2 = -(-2.0/3.0)*mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho;
            if (flag[0] == 0 && flag[3] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U0_nb[0]+U0_nb[3]+U0_nb[9]+U0_center)/8.0;
                solver.Add_A(I*3+1,I*3,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[0]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[3]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[9]*3,  -(coeff_temp1+coeff_temp2)/8.0);
            }
            else assert(false);


            //set the coefficients for U2 in the second equation (mu*w_y)_z
            //traverse the four corners
            //(i,j-1/2,k-1/2),(i,j+1/2,k-1/2),(i,j+1/2,k+1/2),(i,j-1/2,k+1/2)

            //corner (i,j-1/2,k-1/2)

            coeff_temp1 = mu_edge[4]*m_dt/(top_h[1]*top_h[2])/rho;
            coeff_temp2 = (-2.0/3.0)*mu_edge[2]*m_dt/(top_h[1]*top_h[2])/rho;
            if (flag[2] == 0 && flag[4] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U2_nb[2]+U2_nb[4]+U2_nb[10]+U2_center)/8.0;
                solver.Add_A(I*3+1,I*3+2,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[2]*3+2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[4]*3+2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[10]*3+2,  -(coeff_temp1+coeff_temp2)/8.0);
            }
            else //do nothing
            {
            }

            //corner (i,j+1/2,k-1/2)

            coeff_temp1 = -mu_edge[4]*m_dt/(top_h[1]*top_h[2])/rho;
            coeff_temp2 = -(-2.0/3.0)*mu_edge[3]*m_dt/(top_h[1]*top_h[2])/rho;
            if (flag[3] == 0 && flag[4] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U2_nb[3]+U2_nb[4]+U2_nb[11]+U2_center)/8.0;
                solver.Add_A(I*3+1,I*3+2,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[3]*3+2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[4]*3+2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[11]*3+2,  -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
            }

            //corner (i,j+1/2,k+1/2)

            coeff_temp1 = mu_edge[5]*m_dt/(top_h[1]*top_h[2])/rho;
            coeff_temp2 = (-2.0/3.0)*mu_edge[3]*m_dt/(top_h[1]*top_h[2])/rho;
            if (flag[3] == 0 && flag[5] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U2_nb[3]+U2_nb[5]+U2_nb[12]+U2_center)/8.0;
                solver.Add_A(I*3+1,I*3+2,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[3]*3+2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[5]*3+2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[12]*3+2,  -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
            }

            //corner (i,j-1/2,k+1/2)

            coeff_temp1 = -mu_edge[5]*m_dt/(top_h[1]*top_h[2])/rho;
            coeff_temp2 = -(-2.0/3.0)*mu_edge[2]*m_dt/(top_h[1]*top_h[2])/rho;
            if (flag[2] == 0 && flag[5] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U2_nb[2]+U2_nb[5]+U2_nb[13]+U2_center)/8.0;

                solver.Add_A(I*3+1,I*3+2,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[2]*3+2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[5]*3+2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+1,I_nb[13]*3+2,  -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
            }

            rhs += m_dt*state.m_U[1];
            rhs += m_dt*cell_center[index].m_state.f_surf[1];
            rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;
            rhs -= m_dt*cell_center[index].m_state.m_adv[1];

            solver.Add_b(I*3+1, rhs);


            /************************************************************************/
            //         Setting the coeffecients of U2 for the third equation
            /************************************************************************/

            coeff2[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
            coeff2[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
            coeff2[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
            coeff2[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
            coeff2[4] = (2.0/3.0)*m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
            coeff2[5] = (2.0/3.0)*m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

            solver.Add_A(I*3+2,I*3+2,1.0);
            rhs = U2_center;

            for (nb = 0; nb < 6; nb++)
            {
                if (flag[nb] == 0)
                {
                    solver.Add_A(I*3+2,I_nb[nb]*3+2,-coeff2[nb]);
                    rhs += coeff2[nb]*U2_nb[nb];
                }
                else //flag[nb] == 1
                {
                    coeff2[nb] = 2.0*coeff2[nb];
                    rhs += coeff2[nb]*(U2_nb[nb] + U2_nb_new[nb]);
                }
                solver.Add_A(I*3+2,I*3+2,coeff2[nb]);
                rhs -= coeff2[nb]*U2_center;
            }

            //set the coefficients for U0 in the thrid equation (mu*u_z)_x
            //traverse the four corners
            //(i-1/2,j,k-1/2),(i+1/2,j,k-1/2),(i+1/2,j,k+1/2),(i-1/2,j,k+1/2)

            //corner (i-1/2,j,k-1/2)

            coeff_temp1 = mu_edge[0]*m_dt/(top_h[0]*top_h[2])/rho;
            coeff_temp2 = (-2.0/3.0)*mu_edge[4]*m_dt/(top_h[0]*top_h[2])/rho;
            if (flag[0] == 0 && flag[4] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U0_nb[0]+U0_nb[4]+U0_nb[14]+U0_center)/8.0;
                solver.Add_A(I*3+2,I*3,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[0]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[4]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[14]*3, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
            }

            //corner (i+1/2,j,k-1/2)

            coeff_temp1 = -mu_edge[1]*m_dt/(top_h[0]*top_h[2])/rho;
            coeff_temp2 = -(-2.0/3.0)*mu_edge[4]*m_dt/(top_h[0]*top_h[2])/rho;
            if (flag[1] == 0 && flag[4] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U0_nb[1]+U0_nb[4]+U0_nb[15]+U0_center)/8.0;
                solver.Add_A(I*3+2,I*3,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[1]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[4]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[15]*3, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
            }

            //corner (i+1/2,j,k+1/2)

            coeff_temp1 = mu_edge[1]*m_dt/(top_h[0]*top_h[2])/rho;
            coeff_temp2 = (-2.0/3.0)*mu_edge[5]*m_dt/(top_h[0]*top_h[2])/rho;
            if (flag[1] == 0 && flag[5] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U0_nb[1]+U0_nb[5]+U0_nb[16]+U0_center)/8.0;
                solver.Add_A(I*3+2,I*3,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[1]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[5]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[16]*3, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
            }

            //corner (i-1/2,j,k+1/2)

            coeff_temp1 = -mu_edge[0]*m_dt/(top_h[0]*top_h[2])/rho;
            coeff_temp2 = -(-2.0/3.0)*mu_edge[5]*m_dt/(top_h[0]*top_h[2])/rho;
            if (flag[0] == 0 && flag[5] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U0_nb[0]+U0_nb[5]+U0_nb[17]+U0_center)/8.0;

                solver.Add_A(I*3+2,I*3,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[0]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[5]*3,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[17]*3, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
            }


            //set the coefficients for U1 in the thrid equation (mu*v_z)_y
            //traverse the four corners
            //(i,j-1/2,k-1/2),(i,j+1/2,k-1/2),(i,j+1/2,k+1/2),(i,j-1/2,k+1/2)

            //corner (i,j-1/2,k-1/2)

            coeff_temp1 = mu_edge[2]*m_dt/(top_h[1]*top_h[2])/rho;
            coeff_temp2 = (-2.0/3.0)*mu_edge[4]*m_dt/(top_h[1]*top_h[2])/rho;
            if (flag[2] == 0 && flag[4] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U1_nb[2]+U1_nb[4]+U1_nb[10]+U1_center)/8.0;

                solver.Add_A(I*3+2,I*3+1,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[2]*3+1,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[4]*3+1,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[10]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
            }

            //corner (i,j+1/2,k-1/2)

            coeff_temp1 = -mu_edge[3]*m_dt/(top_h[1]*top_h[2])/rho;
            coeff_temp2 = -(-2.0/3.0)*mu_edge[4]*m_dt/(top_h[1]*top_h[2])/rho;
            if (flag[3] == 0 && flag[4] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U1_nb[3]+U1_nb[4]+U1_nb[11]+U1_center)/8.0;

                solver.Add_A(I*3+2,I*3+1,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[3]*3+1,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[4]*3+1,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[11]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
            }

            //corner (i,j+1/2,k+1/2)

            coeff_temp1 = mu_edge[3]*m_dt/(top_h[1]*top_h[2])/rho;
            coeff_temp2 = (-2.0/3.0)*mu_edge[5]*m_dt/(top_h[1]*top_h[2])/rho;
            if (flag[3] == 0 && flag[5] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U1_nb[3]+U1_nb[5]+U1_nb[12]+U1_center)/8.0;

                solver.Add_A(I*3+2,I*3+1,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[3]*3+1,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[5]*3+1,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[12]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
            }

            //corner (i,j-1/2,k+1/2)

            coeff_temp1 = -mu_edge[2]*m_dt/(top_h[1]*top_h[2])/rho;
            coeff_temp2 = -(-2.0/3.0)*mu_edge[5]*m_dt/(top_h[1]*top_h[2])/rho;
            if (flag[2] == 0 && flag[5] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U1_nb[2]+U1_nb[5]+U1_nb[13]+U1_center)/8.0;

                solver.Add_A(I*3+2,I*3+1,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[2]*3+1,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[5]*3+1,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*3+2,I_nb[13]*3+1, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
            }

            rhs += m_dt*state.m_U[2];
            rhs += m_dt*cell_center[index].m_state.f_surf[2];
            rhs -= m_dt*cell_center[index].m_state.grad_q[2]/rho;
            rhs -= m_dt*cell_center[index].m_state.m_adv[2];

            solver.Add_b(I*3+2, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

        start_clock("Before Petsc Solve");

        solver.Solve_GMRES();

        if (debugging("PETSc"))
        {
            double max, min;
            solver.GetExtremeSingularValues(&max, &min);
            (void) printf("The max singular value of A = %lf in compDiffWithSmoothProperty_velocity_vd\n", max);
            (void) printf("The min singular value of A = %lf in compDiffWithSmoothProperty_velocity_vd\n", min);
            if ( min != 0.0)
                (void) printf("The Cond Num of A = %lf in compDiffWithSmoothProperty_velocity_vd\n", max/min);
        }

        solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);

        //if (rel_residual > 1)
        //{
        //    printf("\nThe solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
        //    solver.Reset_x();
        //    solver.Solve_GMRES();
        //    solver.GetNumIterations(&num_iter);
        //    solver.GetFinalRelativeResidualNorm(&rel_residual);
        //}

        stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_velocity_vd: "
                        "num_iter = %d, rel_residual = %le. \n",
                        num_iter,rel_residual);

        // store the value of U^n in U_tmp
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_tmp[0] = cell_center[index].m_state.m_U[0];
            cell_center[index].m_state.m_U_tmp[1] = cell_center[index].m_state.m_U[1];
            cell_center[index].m_state.m_U_tmp[2] = cell_center[index].m_state.m_U[2];
        }

        // store the value of U* in U
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
                cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
                cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
                u = fabs(cell_center[index].m_state.m_U[0]);
                v = fabs(cell_center[index].m_state.m_U[1]);
                w = fabs(cell_center[index].m_state.m_U[2]);
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
                        fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                {
                    max_speed = speed;
                    indmax[0] = i;
                    indmax[1] = j;
                    indmax[2] = k;
                }
                if (u > max_u)
                    max_u = u;
                if (v > max_v)
                    max_v = v;
                if (w > max_w)
                    max_w = w;
            }
            else
            {
                (void) printf("I[%d][%d][%d] = -1 !!\n",i,j,k);
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
                cell_center[index].m_state.m_U[2] = 0.0;
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
                        fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                    max_speed = speed;
                if (u > max_u)
                    max_u = u;
                if (v > max_v)
                    max_v = v;
                if (w > max_w)
                    max_w = w;
            }
        }

        //print max and min of speed, u, v, and w
        max_tmp = max_speed;
        pp_global_max(&max_speed,1);
        (void) printf("local max_speed after compDiffWithSmoothProperty_velocity_vd "
                      "in cell(%d, %d, %d) of node #%d is: %lf\n",indmax[0],indmax[1],indmax[2],pp_mynode(),max_tmp);
        if (max_tmp == max_speed)
        (void) printf("max_speed (locates in node #%d) is: %lf\n",pp_mynode(),max_speed);

        max_tmp = max_u;
        pp_global_max(&max_u,1);
        (void) printf("local max_u after compDiffWithSmoothProperty_velocity_vd "
                      "of node #%d is: %lf\n",pp_mynode(),max_tmp);
        if (max_tmp == max_u)
        (void) printf("max_u (locates in node #%d) is: %lf\n",pp_mynode(),max_u);

        max_tmp = max_v;
        pp_global_max(&max_v,1);
        (void) printf("local max_v after compDiffWithSmoothProperty_velocity_vd "
                      "of node #%d is: %lf\n",pp_mynode(),max_tmp);
        if (max_tmp == max_v)
        (void) printf("max_v (locates in node #%d) is: %lf\n",pp_mynode(),max_v);

        max_tmp = max_w;
        pp_global_max(&max_w,1);
        (void) printf("local max_w after compDiffWithSmoothProperty_velocity_vd "
                      "of node #%d is: %lf\n",pp_mynode(),max_tmp);
        if (max_tmp == max_w)
        (void) printf("max_w (locates in node #%d) is: %lf\n",pp_mynode(),max_w);

        for (l = 0; l < 3; ++l)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U_tmp[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U_tmp[l] = array[index];
            }
        }
        for (l = 0; l < 3; ++l)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }

        if(debugging("step_size"))
        {
            value = sum_div = 0.0;

            for (l = 0; l < dim; l++)
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                vel[l][index] = cell_center[index].m_state.m_U[l];
            }

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;
                index = d_index3d(i,j,k,top_gmax);
                source[index] = computeFieldPointDiv(icoords,vel);
            }
            FT_ParallelExchGridArrayBuffer(source,front);

            // store the div(U*) in div_U
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.div_U = source[index];
            }

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(source[index]);
                sum_div += source[index];
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U* is %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U* is %.16g\n",max_value);
            max_value = 0.0;
        }

        FT_FreeThese(1,x);
} /* end compDiffWithSmoothProperty_velocity_vd */


// just for divergence free case, need to modify???????
void Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_velocity_decoupled_vd(void)
{
        COMPONENT comp;
        int index,index_nb[6],indmax[3],size;
        int I,I_nb[6];
        double coords[MAXD],crx_coords[MAXD];
        double coeff0[6],coeff1[6],coeff2[6],coeff_temp1,coeff_temp2;
        double mu[6],mu_edge[6],mu0,rho,rhs;
        double U_nb[6],U_center;
        double U_nb_new[6];
        int flag[6]; //denote whether this is dirichlet or neumann boundary
        L_STATE state;
        int i,j,k,l,nb,icoords[MAXD];
        INTERFACE *intfc = front->interf;
        double speed,u,v,w,max_u,max_v,max_w,max_tmp,sum_div,value;
        double *x;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
        double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};
        POINTER intfc_state;
        HYPER_SURF *hs;
        PetscInt num_iter;
        double rel_residual;
        double **vel = iFparams->field->vel;

        max_speed = 0.0;
        max_u = 0.0;
        max_v = 0.0;
        max_w = 0.0;
        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

    for (l = 0; l < dim; ++l)
    {
        PETSc solver;

        solver.Create(ilower, iupper-1, 7, 7);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();

        for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
                for (i = imin; i <= imax; i++)
                {
                    I = ijk_to_I[i][j][k];
                    if (I == -1) continue;

                    index = d_index3d(i,j,k,top_gmax);
                    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
                    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
                    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
                    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
                    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
                    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

                    icoords[0] = i;
                    icoords[1] = j;
                    icoords[2] = k;
                    comp = top_comp[index];

                    I_nb[0] = ijk_to_I[i-1][j][k]; //west
                    I_nb[1] = ijk_to_I[i+1][j][k]; //east
                    I_nb[2] = ijk_to_I[i][j-1][k]; //south
                    I_nb[3] = ijk_to_I[i][j+1][k]; //north
                    I_nb[4] = ijk_to_I[i][j][k-1]; //lower
                    I_nb[5] = ijk_to_I[i][j][k+1]; //upper


                    mu0 = cell_center[index].m_state.m_mu;
                    rho = 0.5*(cell_center[index].m_state.m_rho +
                               cell_center[index].m_state.m_rho_old);
                    U_center = cell_center[index].m_state.m_U[l];

                    for (nb = 0; nb < 6; nb++)
                    {
                        if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                                comp,&intfc_state,&hs,crx_coords,m_t_old) &&
                                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                        {
                            //old & new boundary condition
                            U_nb[nb] = 0.0;
                            U_nb_new[nb] = 0.0;

                            if (wave_type(hs) == DIRICHLET_BOUNDARY ||
                                    wave_type(hs) == NEUMANN_BOUNDARY)
                                mu[nb] = mu0;
                            else
                                mu[nb] = 1.0/2*(mu0 +
                                        cell_center[index_nb[nb]].m_state.m_mu);
                        }
                        else
                        {
                            U_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[l];
                            mu[nb] = 1.0/2*(mu0 +
                                    cell_center[index_nb[nb]].m_state.m_mu);
                        }
                    }


                    getRectangleCenter(index, coords);
                    computeSourceTerm(coords, state);


                    rhs = 0;
                    double dh[3] = {top_h[0], top_h[1], top_h[2]};
                    for(nb = 0; nb<6; nb++)
                    {
                        // use dh[1]*dh[2] as the face area
                        if(nb<2)
                        {
                            dh[0] = top_h[0];
                            dh[1] = top_h[1];
                            dh[2] = top_h[2];
                        }
                        else if(nb<4 && nb>=2)
                        {
                            dh[0] = top_h[1];
                            dh[1] = top_h[2];
                            dh[2] = top_h[0];
                        }
                        else if(nb<6 && nb>=4)
                        {
                            dh[0] = top_h[2];
                            dh[1] = top_h[0];
                            dh[2] = top_h[1];
                        }

                        if(I_nb[nb]>=0)  // interior
                        {
                            // u^{*}
                            solver.Add_A(I, I,        (-1) * -0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]));
                            solver.Add_A(I, I_nb[nb], (-1) *  0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]));
                            // u^{n}
                            rhs += 0.5*m_dt/rho*mu[nb] *
                                    (U_nb[nb]-U_center) /(dh[0]*dh[0]);
                        }
                        else            // boundary
                        {
                            // u^{*}
                            solver.Add_A(I, I,       (-1) * -0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]) * 2);
                            rhs += 0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]) * 2 * U_nb_new[nb];
                            // u^{n}
                            rhs += 0.5*m_dt/rho*mu[nb] *
                                    (U_nb[nb]-U_center) /(dh[0]*dh[0]) * 2;
                        }
                    }

                    rhs += m_dt*state.m_U[l];
                    rhs += m_dt*cell_center[index].m_state.f_surf[l];
                    rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho;
                    //printf("gradq = %.16g\n",cell_center[index].m_state.grad_q[l]);

                    rhs -= m_dt * cell_center[index].m_state.m_adv[l]; //advection source term

                    solver.Add_A(I, I, 1.0);
                    rhs += U_center;

                    solver.Add_b(I, rhs);
                }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

        start_clock("Before Petsc Solve");

        solver.Solve_GMRES();

        if (debugging("PETSc"))
        {
            double max, min;
            solver.GetExtremeSingularValues(&max, &min);
            (void) printf("The max singular value of A = %lf in compDiffWithSmoothProperty_velocity_decoupled_vd\n", max);
            (void) printf("The min singular value of A = %lf in compDiffWithSmoothProperty_velocity_decoupled_vd\n", min);
            if ( min != 0.0)
                (void) printf("The Cond Num of A = %lf in compDiffWithSmoothProperty_velocity_decoupled_vd\n", max/min);
        }

        PetscInt num_iter;
        double rel_residual;
        solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);

        stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cartesian::"
                          "compDiffWithSmoothProperty_velocity_decoupled_vd: "
                          "num_iter = %d, rel_residual = %le. \n",
                          num_iter,rel_residual);

        // store the value of U^n in U_tmp
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_tmp[l] = cell_center[index].m_state.m_U[l];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[l] = x[I-ilower];
            }
            else
            {
                cell_center[index].m_state.m_U[l] = 0.0;
            }
        }

        //scatter states
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_U_tmp[l];
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_tmp[l] = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_U[l];
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U[l] = array[index];
        }
    }

    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
        index = d_index3d(i,j,k,top_gmax);
        speed = fabs(cell_center[index].m_state.m_U[0]) +
                fabs(cell_center[index].m_state.m_U[1]) +
                fabs(cell_center[index].m_state.m_U[2]);
        if (speed > max_speed)
            max_speed = speed;
    }
    pp_global_max(&max_speed,1);

        if(debugging("step_size"))
        {
            value = sum_div = 0.0;

            for (l = 0; l < dim; l++)
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                vel[l][index] = cell_center[index].m_state.m_U[l];
            }

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;
                index = d_index3d(i,j,k,top_gmax);
                source[index] = computeFieldPointDiv(icoords,vel);
            }
            FT_ParallelExchGridArrayBuffer(source,front);

            // store the div(U*) in div_U
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.div_U = source[index];
            }

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(source[index]);
                sum_div += source[index];
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U* is %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U* is %.16g\n",max_value);
            max_value = 0.0;
        }

        FT_FreeThese(1,x);
} /* end compDiffWithSmoothProperty_velocity_decoupled_vd */


//Viscous terms = mu*laplace(U)
//valid for constant mu and divergence-free case (no cross terms)
void Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_velocity_MAC_decoupled_vd(void)
{
    bool bNoBoundary[6];
    COMPONENT comp;
    int index,index_nb[6],indmax[3],size;
    int I,I_nb[6];
    double coords[MAXD],crx_coords[MAXD],dh[MAXD];
    double mu,rho,rhs;
    double U_nb[6],U_center;
    L_STATE state;
    int i,j,k,l,nb,icoords[MAXD];
    INTERFACE *intfc = front->interf;
    double speed,u,v,w,max_u,max_v,max_w,max_tmp,sum_div,value;
    double *x;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    POINTER intfc_state;
    HYPER_SURF *hs;
    PetscInt num_iter;
    double rel_residual;
    double **vel = iFparams->field->vel;

    max_speed = 0.0;
    max_u = 0.0;
    max_v = 0.0;
    max_w = 0.0;
    setIndexMap();

    size = iupper - ilower;
    FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

    //solve for u, v, and w separately
    for (l = 0; l < dim; ++l)
    {
        PETSc solver;

        solver.Create(ilower, iupper-1, 7, 7);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();

        for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
                for (i = imin; i <= imax; i++)
                {
                    I = ijk_to_I[i][j][k];
                    if (I == -1) continue;

                    index = d_index3d(i,j,k,top_gmax);
                    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
                    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
                    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
                    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
                    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
                    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

                    icoords[0] = i;
                    icoords[1] = j;
                    icoords[2] = k;
                    comp = top_comp[index];

                    I_nb[0] = ijk_to_I[i-1][j][k]; //west
                    I_nb[1] = ijk_to_I[i+1][j][k]; //east
                    I_nb[2] = ijk_to_I[i][j-1][k]; //south
                    I_nb[3] = ijk_to_I[i][j+1][k]; //north
                    I_nb[4] = ijk_to_I[i][j][k-1]; //lower
                    I_nb[5] = ijk_to_I[i][j][k+1]; //upper

                    U_center = cell_center[index].m_state.m_U[l];

                    // 4 directions
                    for (nb = 0; nb < 4; nb++)
                    {
                        bNoBoundary[nb] = YES;
                        U_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[l];
                    }

                    // LOWER
                    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[4],
                            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                    {
                        bNoBoundary[4] = NO;
                        if (l < 2) //for u & v
                            U_nb[4] = -U_center;
                        else //for w
                            U_nb[4] = 0.0;
                    }
                    else
                    {
                        bNoBoundary[4] = YES;
                        U_nb[4] = cell_center[index_nb[4]].m_state.m_U[l];
                    }

                    // UPPER
                    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[5],
                            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                    {
                        bNoBoundary[5] = NO;
                        if (l < 2) //for u & v
                            U_nb[5] = -U_center;
                        else //for w
                            U_nb[5] = -cell_center[index_nb[4]].m_state.m_U[l];
                    }
                    else
                    {
                        bNoBoundary[5] = YES;
                        U_nb[5] = cell_center[index_nb[5]].m_state.m_U[l];
                    }

                    //interpolate mu & rho on the l-th cell face
                    if (!bNoBoundary[2*l+1]) //cells on dir[2*l+1] bdry
                    {
                        mu = cell_center[index].m_state.m_mu;
                        rho = 2.0/(1.0/cell_center[index].m_state.m_rho +
                                   1.0/cell_center[index].m_state.m_rho_old);
                    }
                    else
                    {
                        mu = 0.5*(cell_center[index].m_state.m_mu +
                                  cell_center[index_nb[2*l+1]].m_state.m_mu);
                        rho = 4.0/(1.0/cell_center[index].m_state.m_rho +
                                   1.0/cell_center[index].m_state.m_rho_old +
                                   1.0/cell_center[index_nb[2*l+1]].m_state.m_rho +
                                   1.0/cell_center[index_nb[2*l+1]].m_state.m_rho_old);
                    }

                    //get gravity force terms (constant)
                    getRectangleCenter(index, coords);
                    computeSourceTerm(coords, state);

                    rhs = 0.0;
                    for (nb = 0; nb < 6; nb++)
                    {
                        switch(nb)
                        {
                        case 0:
                        case 1:
                            dh[0] = top_h[0];
                            break;
                        case 2:
                        case 3:
                            dh[0] = top_h[1];
                            break;
                        case 4:
                        case 5:
                            dh[0] = top_h[2];
                            break;
                        default:
                            assert(false);
                        }

                        //if (!bNoBoundary[nb]) //cells on LOWER or UPPER bdry
                        if (I_nb[nb] < 0)  // boundary
                        {
                            if (nb == 4) //cells on LOWER bdry
                            {
                                if (l < 2) //for u & v
                                {
                                    // U^{*}
                                    solver.Add_A(I, I, (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0]) * 2);
                                    // U^{n}
                                    rhs += 0.5*m_dt/rho*mu*(U_nb[nb]-U_center) /(dh[0]*dh[0]);
                                }
                                else //for w
                                {
                                    // U^{*}
                                    solver.Add_A(I, I, (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0]));
                                    // U^{n}
                                    rhs += 0.5*m_dt/rho*mu*(U_nb[nb]-U_center) /(dh[0]*dh[0]);
                                }
                            } //nb == 4
                            else if (nb == 5) //cells on UPPER bdry
                            {
                                if (l < 2) //for u & v
                                {
                                    // U^{*}
                                    solver.Add_A(I, I, (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0]) * 2);
                                    // U^{n}
                                    rhs += 0.5*m_dt/rho*mu*(U_nb[nb]-U_center) /(dh[0]*dh[0]);
                                }
                                else //for w
                                {
                                    // U^{*}
                                    solver.Add_A(I, I,       (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0]));
                                    solver.Add_A(I, I_nb[4], (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0])); //w*[5] = -w*[4]
                                    // U^{n}
                                    rhs += 0.5*m_dt/rho*mu*(U_nb[nb]-U_center) /(dh[0]*dh[0]);
                                }
                            } //nb == 5
                        }
                        //else //bNoBoundary[nb] = YES
                        else            // interior
                        {
                            // U^{*}
                            solver.Add_A(I, I,        (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0]));
                            solver.Add_A(I, I_nb[nb], (-1) *  0.5*m_dt/rho*mu /(dh[0]*dh[0]));
                            // U^{n}
                            rhs += 0.5*m_dt/rho*mu*(U_nb[nb]-U_center) /(dh[0]*dh[0]);
                        }
                    } //loop for nb ends

                    rhs += m_dt*state.m_U[l];
//                  rhs += m_dt*cell_center[index].m_state.f_surf[l];
                    rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho;
                    //printf("gradq = %.16g\n",cell_center[index].m_state.grad_q[l]);

                    rhs -= m_dt * cell_center[index].m_state.m_adv[l]; //advection source term

                    solver.Add_A(I, I, 1.0);
                    rhs += U_center;
                    if (!bNoBoundary[5] && l==2) //for w, rhs = 0 on UPPER bdry
                        rhs = 0.0;

                    solver.Add_b(I, rhs);
                } //loop for (i,j,k) ends

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

        start_clock("Before Petsc Solve");

        solver.Solve_GMRES();

//        solver.Print_A(NULL);
//        solver.Print_b(NULL);

        if (debugging("PETSc"))
        {
            double max, min;
            solver.GetExtremeSingularValues(&max, &min);
            (void) printf("The max singular value of A = %lf in "
                          "compDiffWithSmoothProperty_velocity_MAC_decoupled_vd for U[%d]\n", max,l);
            (void) printf("The min singular value of A = %lf in "
                          "compDiffWithSmoothProperty_velocity_MAC_decoupled_vd for U[%d]\n", min,l);
            if ( min != 0.0)
                (void) printf("The Cond Num of A = %lf in "
                              "compDiffWithSmoothProperty_velocity_MAC_decoupled_vd for U[%d]\n",max/min,l);
        }

        PetscInt num_iter;
        double rel_residual;
        solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);

        stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cartesian::"
                          "compDiffWithSmoothProperty_velocity_MAC_decoupled_vd for U[%d]: "
                          "num_iter = %d, rel_residual = %le. \n",
                          l,num_iter,rel_residual);

        // store the value of U^n in U_tmp
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_tmp[l] = cell_center[index].m_state.m_U[l];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[l] = x[I-ilower];
            }
            else
            {
                cell_center[index].m_state.m_U[l] = 0.0;
            }
        }

        //scatter states
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_U_tmp[l];
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_tmp[l] = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_U[l];
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U[l] = array[index];
        }
    } //the loop for l ends

    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
        index = d_index3d(i,j,k,top_gmax);
        speed = fabs(cell_center[index].m_state.m_U[0]) +
                fabs(cell_center[index].m_state.m_U[1]) +
                fabs(cell_center[index].m_state.m_U[2]);
        if (speed > max_speed)
            max_speed = speed;
    }
    pp_global_max(&max_speed,1);

    if(debugging("step_size"))
    {
        value = sum_div = 0.0;

        for (l = 0; l < dim; l++)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);
        }
        FT_ParallelExchGridArrayBuffer(source,front);

        // store the div(U*) in div_U
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            value = fabs(source[index]);
            sum_div += source[index];
            if(value > max_value)
                max_value = value;
        }
        pp_global_sum(&sum_div,1);
        (void) printf("\nThe summation of divergence of U* is %.16g\n",sum_div);
        pp_global_max(&max_value,1);
        (void) printf("\nThe max value of divergence of U* is %.16g\n",max_value);
        max_value = 0.0;
    }

    FT_FreeThese(1,x);
} /* end compDiffWithSmoothProperty_velocity_MAC_decoupled_vd */


//set w = 0
//Viscous terms = mu*laplace(U)
//valid for constant mu and divergence-free case (no cross terms)
void Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_velocity_MAC_decoupled_zeroW_vd(void)
{
    bool bNoBoundary[6];
    COMPONENT comp;
    int index,index_nb[6],indmax[3],size;
    int I,I_nb[6];
    double coords[MAXD],crx_coords[MAXD],dh[MAXD];
    double mu,rho,rhs;
    double U_nb[6],U_center;
    L_STATE state;
    int i,j,k,l,nb,icoords[MAXD];
    INTERFACE *intfc = front->interf;
    double speed,u,v,w,max_u,max_v,max_w,max_tmp,sum_div,value;
    double *x;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    POINTER intfc_state;
    HYPER_SURF *hs;
    PetscInt num_iter;
    double rel_residual;
    double **vel = iFparams->field->vel;

    max_speed = 0.0;
    max_u = 0.0;
    max_v = 0.0;
    max_w = 0.0;
    setIndexMap();

    size = iupper - ilower;
    FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

    //solve for u and v separately, set w = 0
    for (l = 0; l < (dim-1); ++l)
    {
        PETSc solver;

        solver.Create(ilower, iupper-1, 7, 7);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();

        for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
                for (i = imin; i <= imax; i++)
                {
                    I = ijk_to_I[i][j][k];
                    if (I == -1) continue;

                    index = d_index3d(i,j,k,top_gmax);
                    comp = top_comp[index];

                    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
                    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
                    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
                    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
                    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
                    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

                    icoords[0] = i;
                    icoords[1] = j;
                    icoords[2] = k;

                    I_nb[0] = ijk_to_I[i-1][j][k]; //west
                    I_nb[1] = ijk_to_I[i+1][j][k]; //east
                    I_nb[2] = ijk_to_I[i][j-1][k]; //south
                    I_nb[3] = ijk_to_I[i][j+1][k]; //north
                    I_nb[4] = ijk_to_I[i][j][k-1]; //lower
                    I_nb[5] = ijk_to_I[i][j][k+1]; //upper

                    U_center = cell_center[index].m_state.m_U[l];

                    // 4 directions
                    for (nb = 0; nb < 4; nb++)
                    {
                        bNoBoundary[nb] = YES;
                        U_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[l];
                    }

                    // LOWER
                    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[4],
                            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                    {
                        bNoBoundary[4] = NO;
                        if (l < 2) //for u & v
                            U_nb[4] = -U_center;
                        else //for w
                            U_nb[4] = 0.0;
                    }
                    else
                    {
                        bNoBoundary[4] = YES;
                        U_nb[4] = cell_center[index_nb[4]].m_state.m_U[l];
                    }

                    // UPPER
                    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[5],
                            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                    {
                        bNoBoundary[5] = NO;
                        if (l < 2) //for u & v
                            U_nb[5] = -U_center;
                        else //for w
                            U_nb[5] = -cell_center[index_nb[4]].m_state.m_U[l];
                    }
                    else
                    {
                        bNoBoundary[5] = YES;
                        U_nb[5] = cell_center[index_nb[5]].m_state.m_U[l];
                    }

                    //interpolate mu & rho on the l-th cell face
                    if (!bNoBoundary[2*l+1]) //cells on dir[2*l+1] bdry
                    {
                        mu = cell_center[index].m_state.m_mu;
                        rho = 2.0/(1.0/cell_center[index].m_state.m_rho +
                                   1.0/cell_center[index].m_state.m_rho_old);
                    }
                    else
                    {
                        mu = 0.5*(cell_center[index].m_state.m_mu +
                                  cell_center[index_nb[2*l+1]].m_state.m_mu);
                        rho = 4.0/(1.0/cell_center[index].m_state.m_rho +
                                   1.0/cell_center[index].m_state.m_rho_old +
                                   1.0/cell_center[index_nb[2*l+1]].m_state.m_rho +
                                   1.0/cell_center[index_nb[2*l+1]].m_state.m_rho_old);
                    }

/*
                    //get gravity force terms (constant)
                    getRectangleCenter(index, coords);
                    computeSourceTerm(coords, state);
*/
                    rhs = 0.0;
                    for (nb = 0; nb < 6; nb++)
                    {
                        switch(nb)
                        {
                        case 0:
                        case 1:
                            dh[0] = top_h[0];
                            break;
                        case 2:
                        case 3:
                            dh[0] = top_h[1];
                            break;
                        case 4:
                        case 5:
                            dh[0] = top_h[2];
                            break;
                        default:
                            assert(false);
                        }

                        //if (!bNoBoundary[nb]) //cells on LOWER or UPPER bdry
                        if (I_nb[nb] < 0)  // boundary
                        {
                            if (nb == 4) //cells on LOWER bdry
                            {
                                if (l < 2) //for u & v
                                {
                                    // U^{*}
                                    solver.Add_A(I, I, (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0]) * 2);
                                    // U^{n}
                                    rhs += 0.5*m_dt/rho*mu*(U_nb[nb]-U_center) /(dh[0]*dh[0]);
                                }
                                else //for w
                                {
                                    // U^{*}
                                    solver.Add_A(I, I, (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0]));
                                    // U^{n}
                                    rhs += 0.5*m_dt/rho*mu*(U_nb[nb]-U_center) /(dh[0]*dh[0]);
                                }
                            } //nb == 4
                            else if (nb == 5) //cells on UPPER bdry
                            {
                                if (l < 2) //for u & v
                                {
                                    // U^{*}
                                    solver.Add_A(I, I, (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0]) * 2);
                                    // U^{n}
                                    rhs += 0.5*m_dt/rho*mu*(U_nb[nb]-U_center) /(dh[0]*dh[0]);
                                }
                                else //for w
                                {
                                    // U^{*}
                                    solver.Add_A(I, I,       (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0]));
                                    solver.Add_A(I, I_nb[4], (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0])); //w*[5] = -w*[4]
                                    // U^{n}
                                    rhs += 0.5*m_dt/rho*mu*(U_nb[nb]-U_center) /(dh[0]*dh[0]);
                                }
                            } //nb == 5
                        }
                        //else //bNoBoundary[nb] = YES
                        else            // interior
                        {
                            // U^{*}
                            solver.Add_A(I, I,        (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0]));
                            solver.Add_A(I, I_nb[nb], (-1) *  0.5*m_dt/rho*mu /(dh[0]*dh[0]));
                            // U^{n}
                            rhs += 0.5*m_dt/rho*mu*(U_nb[nb]-U_center) /(dh[0]*dh[0]);
                        }
                    } //loop for nb ends

                    //rhs += m_dt*state.m_U[l];
//                  rhs += m_dt*cell_center[index].m_state.f_surf[l];
                    rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho;
                    //printf("gradq = %.16g\n",cell_center[index].m_state.grad_q[l]);

                    rhs -= m_dt * cell_center[index].m_state.m_adv[l]; //advection source term

                    solver.Add_A(I, I, 1.0);
                    rhs += U_center;
                    if (!bNoBoundary[5] && l==2) //for w, rhs = 0 on UPPER bdry
                        rhs = 0.0;

                    solver.Add_b(I, rhs);
                } //loop for (i,j,k) ends

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

        start_clock("Before Petsc Solve");

        solver.Solve_GMRES();

//        solver.Print_A(NULL);
//        solver.Print_b(NULL);

        if (debugging("PETSc"))
        {
            double max, min;
            solver.GetExtremeSingularValues(&max, &min);
            (void) printf("The max singular value of A = %lf in "
                          "compDiffWithSmoothProperty_velocity_MAC_decoupled_vd for U[%d]\n", max,l);
            (void) printf("The min singular value of A = %lf in "
                          "compDiffWithSmoothProperty_velocity_MAC_decoupled_vd for U[%d]\n", min,l);
            if ( min != 0.0)
                (void) printf("The Cond Num of A = %lf in "
                              "compDiffWithSmoothProperty_velocity_MAC_decoupled_vd for U[%d]\n",max/min,l);
        }

        PetscInt num_iter;
        double rel_residual;
        solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);

        stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cartesian::"
                          "compDiffWithSmoothProperty_velocity_MAC_decoupled_vd for U[%d]: "
                          "num_iter = %d, rel_residual = %le. \n",
                          l,num_iter,rel_residual);

        // store the value of U^n in U_tmp
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_tmp[l] = cell_center[index].m_state.m_U[l];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[l] = x[I-ilower];
            }
            else
            {
                cell_center[index].m_state.m_U[l] = 0.0;
            }
        }

        //scatter states
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_U_tmp[l];
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_tmp[l] = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_U[l];
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U[l] = array[index];
        }
    } //the loop for l ends

    //set w = 0
    for (k = 0; k <= top_gmax[2]; k++)
    for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
    {
        index = d_index3d(i,j,k,top_gmax);
        cell_center[index].m_state.m_U[2] = 0.0;
        cell_center[index].m_state.m_U_tmp[2] = 0.0;
    }

    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
        index = d_index3d(i,j,k,top_gmax);
        speed = fabs(cell_center[index].m_state.m_U[0]) +
                fabs(cell_center[index].m_state.m_U[1]) +
                fabs(cell_center[index].m_state.m_U[2]);
        if (speed > max_speed)
            max_speed = speed;
    }
    pp_global_max(&max_speed,1);

    if(debugging("step_size"))
    {
        value = sum_div = 0.0;

        for (l = 0; l < dim; l++)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);
        }
        FT_ParallelExchGridArrayBuffer(source,front);

        // store the div(U*) in div_U
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            value = fabs(source[index]);
            sum_div += source[index];
            if(value > max_value)
                max_value = value;
        }
        pp_global_sum(&sum_div,1);
        (void) printf("\nThe summation of divergence of U* is %.16g\n",sum_div);
        pp_global_max(&max_value,1);
        (void) printf("\nThe max value of divergence of U* is %.16g\n",max_value);
        max_value = 0.0;
    }

    FT_FreeThese(1,x);
} /* end compDiffWithSmoothProperty_velocity_MAC_decoupled_zeroW_vd */


//set v = 0
//Viscous terms = mu*laplace(U)
//valid for constant mu and divergence-free case (no cross terms)
void Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_velocity_MAC_decoupled_zeroV_vd(void)
{
    bool bNoBoundary[6];
    COMPONENT comp;
    int index,index_nb[6],indmax[3],size;
    int I,I_nb[6];
    double coords[MAXD],crx_coords[MAXD],dh[MAXD];
    double mu,rho,rhs;
    double U_nb[6],U_center;
    L_STATE state;
    int i,j,k,l,nb,icoords[MAXD];
    INTERFACE *intfc = front->interf;
    double speed,u,v,w,max_u,max_v,max_w,max_tmp,sum_div,value;
    double *x;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    POINTER intfc_state;
    HYPER_SURF *hs;
    PetscInt num_iter;
    double rel_residual;
    double **vel = iFparams->field->vel;

    max_speed = 0.0;
    max_u = 0.0;
    max_v = 0.0;
    max_w = 0.0;
    setIndexMap();

    size = iupper - ilower;
    FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

    //solve for u and w separately, set v = 0
    for (l = 0; l < dim; ++l)
    {
        if (l == 1) continue;

        PETSc solver;

        solver.Create(ilower, iupper-1, 7, 7);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();

        for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
                for (i = imin; i <= imax; i++)
                {
                    I = ijk_to_I[i][j][k];
                    if (I == -1) continue;

                    index = d_index3d(i,j,k,top_gmax);
                    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
                    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
                    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
                    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
                    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
                    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

                    icoords[0] = i;
                    icoords[1] = j;
                    icoords[2] = k;
                    comp = top_comp[index];

                    I_nb[0] = ijk_to_I[i-1][j][k]; //west
                    I_nb[1] = ijk_to_I[i+1][j][k]; //east
                    I_nb[2] = ijk_to_I[i][j-1][k]; //south
                    I_nb[3] = ijk_to_I[i][j+1][k]; //north
                    I_nb[4] = ijk_to_I[i][j][k-1]; //lower
                    I_nb[5] = ijk_to_I[i][j][k+1]; //upper

                    U_center = cell_center[index].m_state.m_U[l];

                    // 4 directions
                    for (nb = 0; nb < 4; nb++)
                    {
                        bNoBoundary[nb] = YES;
                        U_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[l];
                    }

                    // LOWER
                    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[4],
                            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                    {
                        bNoBoundary[4] = NO;
                        if (l < 2) //for u & v
                            U_nb[4] = -U_center;
                        else //for w
                            U_nb[4] = 0.0;
                    }
                    else
                    {
                        bNoBoundary[4] = YES;
                        U_nb[4] = cell_center[index_nb[4]].m_state.m_U[l];
                    }

                    // UPPER
                    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[5],
                            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                    {
                        bNoBoundary[5] = NO;
                        if (l < 2) //for u & v
                            U_nb[5] = -U_center;
                        else //for w
                            U_nb[5] = -cell_center[index_nb[4]].m_state.m_U[l];
                    }
                    else
                    {
                        bNoBoundary[5] = YES;
                        U_nb[5] = cell_center[index_nb[5]].m_state.m_U[l];
                    }

                    //interpolate mu & rho on the l-th cell face
                    if (!bNoBoundary[2*l+1]) //cells on dir[2*l+1] bdry
                    {
                        mu = cell_center[index].m_state.m_mu;
                        rho = 2.0/(1.0/cell_center[index].m_state.m_rho +
                                   1.0/cell_center[index].m_state.m_rho_old);
                    }
                    else
                    {
                        mu = 0.5*(cell_center[index].m_state.m_mu +
                                  cell_center[index_nb[2*l+1]].m_state.m_mu);
                        rho = 4.0/(1.0/cell_center[index].m_state.m_rho +
                                   1.0/cell_center[index].m_state.m_rho_old +
                                   1.0/cell_center[index_nb[2*l+1]].m_state.m_rho +
                                   1.0/cell_center[index_nb[2*l+1]].m_state.m_rho_old);
                    }

/*
                    //get gravity force terms (constant)
                    getRectangleCenter(index, coords);
                    computeSourceTerm(coords, state);
*/
                    rhs = 0.0;
                    for (nb = 0; nb < 6; nb++)
                    {
                        switch(nb)
                        {
                        case 0:
                        case 1:
                            dh[0] = top_h[0];
                            break;
                        case 2:
                        case 3:
                            dh[0] = top_h[1];
                            break;
                        case 4:
                        case 5:
                            dh[0] = top_h[2];
                            break;
                        default:
                            assert(false);
                        }

                        //if (!bNoBoundary[nb]) //cells on LOWER or UPPER bdry
                        if (I_nb[nb] < 0)  // boundary
                        {
                            if (nb == 4) //cells on LOWER bdry
                            {
                                if (l < 2) //for u & v
                                {
                                    // U^{*}
                                    solver.Add_A(I, I, (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0]) * 2);
                                    // U^{n}
                                    rhs += 0.5*m_dt/rho*mu*(U_nb[nb]-U_center) /(dh[0]*dh[0]);
                                }
                                else //for w
                                {
                                    // U^{*}
                                    solver.Add_A(I, I, (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0]));
                                    // U^{n}
                                    rhs += 0.5*m_dt/rho*mu*(U_nb[nb]-U_center) /(dh[0]*dh[0]);
                                }
                            } //nb == 4
                            else if (nb == 5) //cells on UPPER bdry
                            {
                                if (l < 2) //for u & v
                                {
                                    // U^{*}
                                    solver.Add_A(I, I, (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0]) * 2);
                                    // U^{n}
                                    rhs += 0.5*m_dt/rho*mu*(U_nb[nb]-U_center) /(dh[0]*dh[0]);
                                }
                                else //for w
                                {
                                    // U^{*}
                                    solver.Add_A(I, I,       (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0]));
                                    solver.Add_A(I, I_nb[4], (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0])); //w*[5] = -w*[4]
                                    // U^{n}
                                    rhs += 0.5*m_dt/rho*mu*(U_nb[nb]-U_center) /(dh[0]*dh[0]);
                                }
                            } //nb == 5
                        }
                        //else //bNoBoundary[nb] = YES
                        else            // interior
                        {
                            // U^{*}
                            solver.Add_A(I, I,        (-1) * -0.5*m_dt/rho*mu /(dh[0]*dh[0]));
                            solver.Add_A(I, I_nb[nb], (-1) *  0.5*m_dt/rho*mu /(dh[0]*dh[0]));
                            // U^{n}
                            rhs += 0.5*m_dt/rho*mu*(U_nb[nb]-U_center) /(dh[0]*dh[0]);
                        }
                    } //loop for nb ends

                    //rhs += m_dt*state.m_U[l];
//                  rhs += m_dt*cell_center[index].m_state.f_surf[l];
                    rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho;
                    //printf("gradq = %.16g\n",cell_center[index].m_state.grad_q[l]);

                    rhs -= m_dt * cell_center[index].m_state.m_adv[l]; //advection source term

                    solver.Add_A(I, I, 1.0);
                    rhs += U_center;
                    if (!bNoBoundary[5] && l==2) //for w, rhs = 0 on UPPER bdry
                        rhs = 0.0;

                    solver.Add_b(I, rhs);
                } //loop for (i,j,k) ends

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

        start_clock("Before Petsc Solve");

        solver.Solve_GMRES();

//        solver.Print_A(NULL);
//        solver.Print_b(NULL);

        if (debugging("PETSc"))
        {
            double max, min;
            solver.GetExtremeSingularValues(&max, &min);
            (void) printf("The max singular value of A = %lf in "
                          "compDiffWithSmoothProperty_velocity_MAC_decoupled_vd for U[%d]\n", max,l);
            (void) printf("The min singular value of A = %lf in "
                          "compDiffWithSmoothProperty_velocity_MAC_decoupled_vd for U[%d]\n", min,l);
            if ( min != 0.0)
                (void) printf("The Cond Num of A = %lf in "
                              "compDiffWithSmoothProperty_velocity_MAC_decoupled_vd for U[%d]\n",max/min,l);
        }

        PetscInt num_iter;
        double rel_residual;
        solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);

        stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cartesian::"
                          "compDiffWithSmoothProperty_velocity_MAC_decoupled_vd for U[%d]: "
                          "num_iter = %d, rel_residual = %le. \n",
                          l,num_iter,rel_residual);

        // store the value of U^n in U_tmp
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_tmp[l] = cell_center[index].m_state.m_U[l];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[l] = x[I-ilower];
            }
            else
            {
                cell_center[index].m_state.m_U[l] = 0.0;
            }
        }

        //scatter states
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_U_tmp[l];
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_tmp[l] = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_U[l];
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U[l] = array[index];
        }
    } //the loop for l ends

    //set v = 0
    for (k = 0; k <= top_gmax[2]; k++)
    for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
    {
        index = d_index3d(i,j,k,top_gmax);
        cell_center[index].m_state.m_U[1] = 0.0;
        cell_center[index].m_state.m_U_tmp[1] = 0.0;
    }

    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
        index = d_index3d(i,j,k,top_gmax);
        speed = fabs(cell_center[index].m_state.m_U[0]) +
                fabs(cell_center[index].m_state.m_U[1]) +
                fabs(cell_center[index].m_state.m_U[2]);
        if (speed > max_speed)
            max_speed = speed;
    }
    pp_global_max(&max_speed,1);

    if(debugging("step_size"))
    {
        value = sum_div = 0.0;

        for (l = 0; l < dim; l++)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);
        }
        FT_ParallelExchGridArrayBuffer(source,front);

        // store the div(U*) in div_U
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            value = fabs(source[index]);
            sum_div += source[index];
            if(value > max_value)
                max_value = value;
        }
        pp_global_sum(&sum_div,1);
        (void) printf("\nThe summation of divergence of U* is %.16g\n",sum_div);
        pp_global_max(&max_value,1);
        (void) printf("\nThe max value of divergence of U* is %.16g\n",max_value);
        max_value = 0.0;
    }

    FT_FreeThese(1,x);
} /* end compDiffWithSmoothProperty_velocity_MAC_decoupled_zeroV_vd */


/*               _                                                                                 _
                |  (4/3*mu*u_x)_x+(mu*u_y)_y+(mu*u_z)_z+(mu*v_x)_y+(mu*w_x)_z-2/3*[mu*(v_y+w_z)]_x  |
viscous terms = |  (mu*v_x)_x+(4/3*mu*v_y)_y+(mu*v_z)_z+(mu*u_y)_x+(mu*w_y)_z-2/3*[mu*(u_x+w_z)]_y  |
                |_ (mu*w_x)_x+(mu*w_y)_y+(4/3*mu*w_z)_z+(mu*u_z)_x+(mu*v_z)_y-2/3*[mu*(u_x+v_y)]_z _|
*/
//valid for variable dynamic viscosity
//NOTE: we use mu_t_old for both U^n and U*, which is a first-order approx.
void Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_velocity_MAC_coupled_vd(void)
{
        int index,index_nb[18],indmax[3],size;
        int I,I_nb[18];
        double coords[MAXD];
        double coeff0[6],coeff1[6],coeff2[6],coeff_temp,coeff_temp_old;
        double coeff0_old[6],coeff1_old[6],coeff2_old[6];
        double mu[6],mu_old[6],rho,rhs;
        double U0_nb[18],U1_nb[18],U2_nb[18],U0_center,U1_center,U2_center;
        double mu_nb[18],mu_t_nb[18],mu_center,mu_t_center;
        double mu_nb_old[18],mu_center_old;
        L_STATE state;
        int i,j,k,l,nb,icoords[MAXD];
        double speed,u,v,w,max_u,max_v,max_w,max_tmp,sum_div,value;
        double *x;
        PetscInt num_iter;
        double rel_residual;
        double **vel = iFparams->field->vel;
        boolean useSGSCellCenter = YES; //use mu_t on the cell center. o.w. use mu_t on 3 cell-faces
        int bNoBoundary[6];//type change
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
        COMPONENT comp;

        max_u = max_v = max_w = max_speed = 0;
        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 20, 20);
        // 7u + 4v + 4w for the 1st momentum equation
        // 4u + 7v + 4w for the 2nd momentum equation
        // 4u + 4v + 7w for the 3rd momentum equation

        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();
        // removal tag: HAOZ REFLECTION BC
        /*
        for (l = 0; l < 3; ++l)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }

        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l];
        }
        // removal tag: HAOZ
        enforceReflectionState(vel);//m_U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U[l] = vel[l][index];
        }// end of Reflection Treatment
        */
        copyMeshStates_vd();

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index = d_index3d(i,j,k,top_gmax);
            comp = cell_center[index].comp;

            //6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);

            //xy cut neighbours
            index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
            index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
            index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
            index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);

            //yz cut neighbours
            index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
            index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
            index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
            index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);

            //xz cut neighbours
            index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
            index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
            index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
            index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);

            //6 neighbours of the center cell
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];

            //xy cut neighbours
            I_nb[6] = ijk_to_I[i-1][j-1][k];
            I_nb[7] = ijk_to_I[i+1][j-1][k];
            I_nb[8] = ijk_to_I[i+1][j+1][k];
            I_nb[9] = ijk_to_I[i-1][j+1][k];

            //yz cut neighbours
            I_nb[10] = ijk_to_I[i][j-1][k-1];
            I_nb[11] = ijk_to_I[i][j+1][k-1];
            I_nb[12] = ijk_to_I[i][j+1][k+1];
            I_nb[13] = ijk_to_I[i][j-1][k+1];

            //xz cut neighbours
            I_nb[14] = ijk_to_I[i-1][j][k-1];
            I_nb[15] = ijk_to_I[i+1][j][k-1];
            I_nb[16] = ijk_to_I[i+1][j][k+1];
            I_nb[17] = ijk_to_I[i-1][j][k+1];

            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;

            U0_center = cell_center[index].m_state.m_U[0];
            U1_center = cell_center[index].m_state.m_U[1];
            U2_center = cell_center[index].m_state.m_U[2];
            mu_center = cell_center[index].m_state.m_mu;
            mu_center_old = cell_center[index].m_state.m_mu_old;
            if (useSGSCellCenter) {
                mu_t_center = cell_center[index].m_state.m_mu_turbulent[3];
            }
            else {
                printf("codes needed for mu_turbulent on cell-faces.\n");
                assert(false);
            }

            //periodic B.C. for WEST,EAST,SOUTH,NORTH
            for (nb = 0; nb < 6; nb++)
            {
                checkBoundaryCondition(dir[nb],icoords,&bNoBoundary[nb],m_t_int,comp);
                if (bNoBoundary[nb] == 3 || !bNoBoundary[nb]) // reflect or periodic or interior
                {
                U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
                U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
                U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
                mu_nb[nb] = cell_center[index_nb[nb]].m_state.m_mu;
                mu_nb_old[nb] = cell_center[index_nb[nb]].m_state.m_mu_old;
                }
                if (bNoBoundary[nb] >= 2) // reflect or neumann
                    I_nb[nb] = -1;
                if (useSGSCellCenter) {
                    mu_t_nb[nb] = cell_center[index_nb[nb]].m_state.m_mu_turbulent[3];
                }
                else {
                    printf("codes needed for mu_turbulent on cell-faces.\n");
                    assert(false);
                }
            }

            //physical B.C. for LOWER & UPPER
            for (nb=4; nb<6; ++nb) {
                if (bNoBoundary[nb] >= 2) { //cells on LOWER/UPPER bdry NEUMANN or REFLECT
                    U0_nb[nb] = -U0_center;
                    U1_nb[nb] = -U1_center;
                    //w[4] = 0
                    if (nb==4)		U2_nb[nb] = 0;
                    //w[5] = -w[4]
                    else if (nb==5)	U2_nb[nb] = -cell_center[index_nb[4]].m_state.m_U[2];
                    else		assert(false);

                    mu_nb[nb] = mu_center;
                    mu_nb_old[nb] = mu_center_old;
                    if (useSGSCellCenter) {
                        mu_t_nb[nb] = mu_t_center;
                    }
                    else {
                        printf("codes needed for mu_turbulent on cell-faces.\n");
                        assert(false);
                    }
                }
                else { //interior cells
                    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
                    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
                    U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
                    mu_nb[nb] = cell_center[index_nb[nb]].m_state.m_mu;
                    mu_nb_old[nb] = cell_center[index_nb[nb]].m_state.m_mu_old;
                    if (useSGSCellCenter) {
                        mu_t_nb[nb] = cell_center[index_nb[nb]].m_state.m_mu_turbulent[3];
                    }
                    else {
                        printf("codes needed for mu_turbulent on cell-faces.\n");
                        assert(false);
                    }
                }
            }

            //gravity force
            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);


            /************************************************************************/
            //         Setting the coeffecients for U in the 1st equation
            /************************************************************************/

            //set rho on u-face
            rho = 4.0/(1.0/cell_center[index].m_state.m_rho +
                       1.0/cell_center[index].m_state.m_rho_old +
                       1.0/cell_center[index_nb[1]].m_state.m_rho +
                       1.0/cell_center[index_nb[1]].m_state.m_rho_old);

            //set mu[0] ~ mu[5]
            if (I_nb[4]>=0 && I_nb[5]>=0) {
                mu_nb[15] = cell_center[index_nb[15]].m_state.m_mu;
                mu_nb[16] = cell_center[index_nb[16]].m_state.m_mu;
                mu_nb_old[15] = cell_center[index_nb[15]].m_state.m_mu_old;
                mu_nb_old[16] = cell_center[index_nb[16]].m_state.m_mu_old;
                if (useSGSCellCenter) {
                    mu_t_nb[15] = cell_center[index_nb[15]].m_state.m_mu_turbulent[3];
                    mu_t_nb[16] = cell_center[index_nb[16]].m_state.m_mu_turbulent[3];
                }
                else {
                    printf("codes needed for mu_turbulent on cell-faces.\n");
                    assert(false);
                }
            }
            else if (I_nb[4] < 0) { //cells on LOWER bdry
                mu_nb[15] = mu_nb[1];
                mu_nb[16] = cell_center[index_nb[16]].m_state.m_mu;
                mu_nb_old[15] = mu_nb_old[1];
                mu_nb_old[16] = cell_center[index_nb[16]].m_state.m_mu_old;
                if (useSGSCellCenter) {
                    mu_t_nb[15] = mu_t_nb[1];
                    mu_t_nb[16] = cell_center[index_nb[16]].m_state.m_mu_turbulent[3];
                }
                else {
                    printf("codes needed for mu_turbulent on cell-faces.\n");
                    assert(false);
                }
            }
            else if (I_nb[5] < 0) { //cells on UPPER bdry
                mu_nb[15] = cell_center[index_nb[15]].m_state.m_mu;
                mu_nb[16] = mu_nb[1];
                mu_nb_old[15] = cell_center[index_nb[15]].m_state.m_mu_old;
                mu_nb_old[16] = mu_nb_old[1];
                if (useSGSCellCenter) {
                    mu_t_nb[15] = cell_center[index_nb[15]].m_state.m_mu_turbulent[3];
                    mu_t_nb[16] = mu_t_nb[1];
                }
                else {
                    printf("codes needed for mu_turbulent on cell-faces.\n");
                    assert(false);
                }
            }

            mu_nb[7] = cell_center[index_nb[7]].m_state.m_mu;
            mu_nb[8] = cell_center[index_nb[8]].m_state.m_mu;
            mu_nb_old[7] = cell_center[index_nb[7]].m_state.m_mu_old;
            mu_nb_old[8] = cell_center[index_nb[8]].m_state.m_mu_old;
            if (useSGSCellCenter) {
                mu_t_nb[7] = cell_center[index_nb[7]].m_state.m_mu_turbulent[3];
                mu_t_nb[8] = cell_center[index_nb[8]].m_state.m_mu_turbulent[3];
            }
            else {
                printf("codes needed for mu_turbulent on cell-faces.\n");
                assert(false);
            }

            mu[0] = mu_center + mu_t_center;
            mu[1] = mu_nb[1] + mu_t_nb[1];
            mu[2] = (mu_center + mu_nb[1] + mu_nb[7] + mu_nb[2] +
                     mu_t_center + mu_t_nb[1] + mu_t_nb[7] + mu_t_nb[2])/4;
            mu[3] = (mu_center + mu_nb[1] + mu_nb[8] + mu_nb[3] +
                     mu_t_center + mu_t_nb[1] + mu_t_nb[8] + mu_t_nb[3])/4;
            mu[4] = (mu_center + mu_nb[1] + mu_nb[15] + mu_nb[4] +
                     mu_t_center + mu_t_nb[1] + mu_t_nb[15] + mu_t_nb[4])/4;
            mu[5] = (mu_center + mu_nb[1] + mu_nb[16] + mu_nb[5] +
                     mu_t_center + mu_t_nb[1] + mu_t_nb[16] + mu_t_nb[5])/4;
            mu_old[0] = mu_center_old + mu_t_center;
            mu_old[1] = mu_nb_old[1] + mu_t_nb[1];
            mu_old[2] = (mu_center_old + mu_nb_old[1] + mu_nb_old[7] + mu_nb_old[2] +
                         mu_t_center + mu_t_nb[1] + mu_t_nb[7] + mu_t_nb[2])/4;
            mu_old[3] = (mu_center_old + mu_nb_old[1] + mu_nb_old[8] + mu_nb_old[3] +
                         mu_t_center + mu_t_nb[1] + mu_t_nb[8] + mu_t_nb[3])/4;
            mu_old[4] = (mu_center_old + mu_nb_old[1] + mu_nb_old[15] + mu_nb_old[4] +
                         mu_t_center + mu_t_nb[1] + mu_t_nb[15] + mu_t_nb[4])/4;
            mu_old[5] = (mu_center_old + mu_nb_old[1] + mu_nb_old[16] + mu_nb_old[5] +
                         mu_t_center + mu_t_nb[1] + mu_t_nb[16] + mu_t_nb[5])/4;

            //set the coeffecients for U0 in the first equation, (4/3*mu*u_x)_x, (mu*u_y)_y and (mu*u_z)_z

            coeff0[0] = 2.0/3*m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            coeff0[1] = 2.0/3*m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            coeff0[2] = 0.5*m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            coeff0[3] = 0.5*m_dt/rho*mu[3]/(top_h[1]*top_h[1]);
            coeff0[4] = 0.5*m_dt/rho*mu[4]/(top_h[2]*top_h[2]);
            coeff0[5] = 0.5*m_dt/rho*mu[5]/(top_h[2]*top_h[2]);
            coeff0_old[0] = 2.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[0]);
            coeff0_old[1] = 2.0/3*m_dt/rho*mu_old[1]/(top_h[0]*top_h[0]);
            coeff0_old[2] = 0.5*m_dt/rho*mu_old[2]/(top_h[1]*top_h[1]);
            coeff0_old[3] = 0.5*m_dt/rho*mu_old[3]/(top_h[1]*top_h[1]);
            coeff0_old[4] = 0.5*m_dt/rho*mu_old[4]/(top_h[2]*top_h[2]);
            coeff0_old[5] = 0.5*m_dt/rho*mu_old[5]/(top_h[2]*top_h[2]);

            solver.Add_A(I*3,I*3,1.0);
            rhs = U0_center;

            // X dir [WEST, EAST] Wall normal to the u* velocity
            if (bNoBoundary[0] == 3) // WEST FACE is on REFLECT BC
            {
                //u[0] = 0
                //no contribution to u[0]
                solver.Add_A(I*3,I*3, coeff0[0]);
                rhs += coeff0_old[0]*(U0_nb[0] - U0_center);
            }
            else if (!bNoBoundary[0]) // WEST FACE is PERIODIC or INTERIOR
            {
                solver.Add_A(I*3,I*3, coeff0[0]);
                solver.Add_A(I*3,I_nb[0]*3, -coeff0[0]);
                rhs += coeff0_old[0]*(U0_nb[0] - U0_center);
            }
            if (bNoBoundary[1] == 3) // EAST FACE is on REFLECT BC
            {
                //u[c] = 0
                //u[1] = - u[0]
                //no contribution to u[c] or u[1]
                solver.Add_A(I*3,I_nb[0]*3,coeff0[1]);
                rhs += coeff0_old[1]*(U0_nb[1] - U0_center);
            }
            else if (!bNoBoundary[1]) // EAST FACE is PERIODIC or INTERIOR
            {
                solver.Add_A(I*3,I*3, coeff0[1]);
                solver.Add_A(I*3,I_nb[1]*3, -coeff0[1]);
                rhs += coeff0_old[1]*(U0_nb[1] - U0_center);
            }
            // Y and Z dir, u* tangential to Wall
            // SLIP flag is applied to tangential velocity on Wall
            if (bNoBoundary[2] == 3) // SOUTH FACE is on REFLECT BC
            {
                if (SLIP)
                {
                    // onWallvel = u*_c (REFLECT, NEUMANN)
                    // u*_nb[2], u*_c no contribution
                    // do nothing for coeff matrix
                    rhs += coeff0_old[2]*(U0_nb[2] - U0_center);
                }
                else
                {
                    // onWallvel = 0
                    printf("NO SLIP CONDITION is not implemented yet in %s!\n", __func__);
                    clean_up(ERROR);
                }
            }
            else if (!bNoBoundary[2]) // SOUTH FACE is PERIODIC or INTERIOR
            {
                solver.Add_A(I*3,I*3, coeff0[2]);
                solver.Add_A(I*3,I_nb[2]*3, -coeff0[2]);
                rhs += coeff0_old[2]*(U0_nb[2] - U0_center);
            }
            if (bNoBoundary[3] == 3) // NORTH FACE is on REFLECT BC
            {
                if (SLIP)
                    rhs += coeff0_old[3]*(U0_nb[3] - U0_center);
                else
                {
                    // onWallvel = 0
                    printf("NO SLIP CONDITION is not implemented yet in %s!\n", __func__);
                    clean_up(ERROR);
                }
            }
            else if (!bNoBoundary[3]) // NORTH FACE is PERIODIC or INTERIOR
            {
                solver.Add_A(I*3,I*3, coeff0[3]);
                solver.Add_A(I*3,I_nb[3]*3, -coeff0[3]);
                rhs += coeff0_old[3]*(U0_nb[3] - U0_center);
            }
            if (bNoBoundary[4] == 3) // LOWER FACE is on REFLECT BC
            {
                if (SLIP)
                    rhs += coeff0_old[4]*(U0_nb[4] - U0_center);
                else
                {
                    printf("NO SLIP CONDITION is not implemented yet in %s!\n", __func__);
                    clean_up(ERROR);
                }
            }
            else if (!bNoBoundary[4]) // LOWER FACE is INTERIOR
            {
                solver.Add_A(I*3,I*3, coeff0[4]);
                solver.Add_A(I*3,I_nb[4]*3, -coeff0[4]);
                rhs += coeff0_old[4]*(U0_nb[4] - U0_center);
            }
            if (bNoBoundary[5] == 3) // UPPER FACE is on REFLECT BC
            {
                if (SLIP)
                    rhs += coeff0_old[5]*(U0_nb[5] - U0_center);
                else
                {
                    printf("NO SLIP CONDITION is not implemented yet in %s!\n", __func__);
                    clean_up(ERROR);
                }
            }
            else if (!bNoBoundary[5]) // UPPER FACE is INTERIOR
            {
                solver.Add_A(I*3,I*3, coeff0[5]);
                solver.Add_A(I*3,I_nb[5]*3, -coeff0[5]);
                rhs += coeff0_old[5]*(U0_nb[5] - U0_center);
            }
            // for NEUMANN BC
            if (bNoBoundary[4] == 2) // LOWER
            {
                solver.Add_A(I*3,I*3, 2.0*coeff0[4]);
                rhs += coeff0_old[4]*(U0_nb[4] - U0_center);
            }
            if (bNoBoundary[5] == 2) // UPPER
            {
                solver.Add_A(I*3,I*3, 2.0*coeff0[5]);
                rhs += coeff0_old[5]*(U0_nb[5] - U0_center);
            }

            //set the coefficients for U1 in the first equation, (mu*v_x)_y and (-2/3)*(mu*v_y)_x
            //traverse the four cells for v, i.e. v_c, v_1, v_2, and v_7

            //concern for boundary condition
            //EAST,SOUTH,NORTH FACEs
            if (!bNoBoundary[1] && !bNoBoundary[2] && !bNoBoundary[3])//PERIODIC or INTERIOR
            {
            //v[1]
            coeff_temp  = 1.0/2*m_dt/rho*mu[3]/(top_h[0]*top_h[1]);
            coeff_temp -= 1.0/3*m_dt/rho*mu[1]/(top_h[0]*top_h[1]);
            solver.Add_A(I*3,I_nb[1]*3+1, -coeff_temp);
            coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[3]/(top_h[0]*top_h[1]);
            coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[1]/(top_h[0]*top_h[1]);
            rhs += coeff_temp_old*U1_nb[1];

            //v[2]
            coeff_temp  = 1.0/2*m_dt/rho*mu[2]/(top_h[0]*top_h[1]);
            coeff_temp -= 1.0/3*m_dt/rho*mu[0]/(top_h[0]*top_h[1]);
            solver.Add_A(I*3,I_nb[2]*3+1, -coeff_temp);
            coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[2]/(top_h[0]*top_h[1]);
            coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[1]);
            rhs += coeff_temp_old*U1_nb[2];

            U1_nb[7] = cell_center[index_nb[7]].m_state.m_U[1];
            //v[7]
            coeff_temp  = -1.0/2*m_dt/rho*mu[2]/(top_h[0]*top_h[1]);
            coeff_temp -= -1.0/3*m_dt/rho*mu[1]/(top_h[0]*top_h[1]);
            solver.Add_A(I*3,I_nb[7]*3+1, -coeff_temp);
            coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[2]/(top_h[0]*top_h[1]);
            coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[1]/(top_h[0]*top_h[1]);
            rhs += coeff_temp_old*U1_nb[7];

            //v[c]
            coeff_temp  = -1.0/2*m_dt/rho*mu[3]/(top_h[0]*top_h[1]);
            coeff_temp -= -1.0/3*m_dt/rho*mu[0]/(top_h[0]*top_h[1]);
            solver.Add_A(I*3,I*3+1, -coeff_temp);
            coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[3]/(top_h[0]*top_h[1]);
            coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[1]);
            rhs += coeff_temp_old*U1_center;
            }
            else if ((bNoBoundary[2] == 3) && !bNoBoundary[1]) // ONLY SOUTH is on BC
            {
                // v[2] = 0 and v[7] = 0

                //v[1]
                coeff_temp  = 1.0/2*m_dt/rho*mu[3]/(top_h[0]*top_h[1]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[1]/(top_h[0]*top_h[1]);
                solver.Add_A(I*3,I_nb[1]*3+1, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[3]/(top_h[0]*top_h[1]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[1]/(top_h[0]*top_h[1]);
                rhs += coeff_temp_old*U1_nb[1];

                //v[c]
                coeff_temp  = -1.0/2*m_dt/rho*mu[3]/(top_h[0]*top_h[1]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[0]/(top_h[0]*top_h[1]);
                solver.Add_A(I*3,I*3+1, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[3]/(top_h[0]*top_h[1]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[1]);
                rhs += coeff_temp_old*U1_center;
            }
            else if ((bNoBoundary[3] == 3) && !bNoBoundary[1]) // ONLY NORTH is on BC
            {
                // v[c] = 0 and v[1] = 0

                //v[7]
                U1_nb[7] = cell_center[index_nb[7]].m_state.m_U[1];
                coeff_temp  = -1.0/2*m_dt/rho*mu[2]/(top_h[0]*top_h[1]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[1]/(top_h[0]*top_h[1]);
                solver.Add_A(I*3,I_nb[7]*3+1, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[2]/(top_h[0]*top_h[1]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[1]/(top_h[0]*top_h[1]);
                rhs += coeff_temp_old*U1_nb[7];

                //v[2]
                coeff_temp  = 1.0/2*m_dt/rho*mu[2]/(top_h[0]*top_h[1]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[0]/(top_h[0]*top_h[1]);
                solver.Add_A(I*3,I_nb[2]*3+1, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[2]/(top_h[0]*top_h[1]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[1]);
                rhs += coeff_temp_old*U1_nb[2];

            }
            else if (bNoBoundary[1] == 3) // EAST is on BC
            {
                //hardwired for slip condition
                double temp1 = -1.0/3*m_dt/rho*mu[1]/(top_h[0]*top_h[1]);
                temp1 += 1.0/3*m_dt/rho*mu[0]/(top_h[0]*top_h[1]);
                double temp2 = -temp1;
                if (bNoBoundary[2] == 3) // SOUTH Added
                {
                    //v[c]
                    solver.Add_A(I*3,I*3+1, -temp1);
                    coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[3]/(top_h[0]*top_h[1]);
                    coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[1]);
                    rhs += coeff_temp_old*U1_center;
                }
                else if (bNoBoundary[3] == 3) // NORTH Added
                {
                    //v[2]
                    solver.Add_A(I*3,I_nb[2]*3+1, -temp2);
                    coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[2]/(top_h[0]*top_h[1]);
                    coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[1]);
                    rhs += coeff_temp_old*U1_nb[2];
                }
                else
                {
                    //v[c]
                    solver.Add_A(I*3,I*3+1, -temp1);
                    coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[3]/(top_h[0]*top_h[1]);
                    coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[1]);
                    rhs += coeff_temp_old*U1_center;
                    //v[2]
                    solver.Add_A(I*3,I_nb[2]*3+1, -temp2);
                    coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[2]/(top_h[0]*top_h[1]);
                    coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[1]);
                    rhs += coeff_temp_old*U1_nb[2];
                }
            }

            //set the coefficients for U2 in the first equation, (mu*w_x)_z and (-2/3)*(mu*w_z)_x
            //traverse the four cells for w, i.e. w_c, w_1, w_4, and w_15
            //consider boundary condition: EAST, LOWER, UPPER
            if (!bNoBoundary[1] && !bNoBoundary[4] && !bNoBoundary[5]) // PERIODIC or INTERIOR
            {
                //w[1]
                coeff_temp  = 1.0/2*m_dt/rho*mu[5]/(top_h[0]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[1]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3,I_nb[1]*3+2, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[5]/(top_h[0]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[1]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U2_nb[1];

                U2_nb[4] = cell_center[index_nb[4]].m_state.m_U[2];
                //w[4]
                coeff_temp  = 1.0/2*m_dt/rho*mu[4]/(top_h[0]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[0]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3,I_nb[4]*3+2, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U2_nb[4];

                U2_nb[15] = cell_center[index_nb[15]].m_state.m_U[2];
                //w[15]
                coeff_temp  = -1.0/2*m_dt/rho*mu[4]/(top_h[0]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[1]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3,I_nb[15]*3+2, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[1]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U2_nb[15];

                //w[c]
                coeff_temp  = -1.0/2*m_dt/rho*mu[5]/(top_h[0]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[0]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3,I*3+2, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[5]/(top_h[0]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U2_center;
            }
            else if ((bNoBoundary[4] == 3) && !bNoBoundary[1]) // ONLY LOWER is on BC
            {
                //w[4] = 0 w[15] = 0
                //w[c]
                coeff_temp  = -1.0/2*m_dt/rho*mu[5]/(top_h[0]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[0]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3,I*3+2, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[5]/(top_h[0]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U2_center;

                //w[1]
                coeff_temp  = 1.0/2*m_dt/rho*mu[5]/(top_h[0]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[1]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3,I_nb[1]*3+2, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[5]/(top_h[0]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[1]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U2_nb[1];
            }
            else if ((bNoBoundary[5] == 3) && !bNoBoundary[1]) // ONLY UPPER is on BC
            {
                //w[c] 0 w[1] = 0

                U2_nb[4] = cell_center[index_nb[4]].m_state.m_U[2];
                //w[4]
                coeff_temp  = 1.0/2*m_dt/rho*mu[4]/(top_h[0]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[0]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3,I_nb[4]*3+2, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U2_nb[4];

                U2_nb[15] = cell_center[index_nb[15]].m_state.m_U[2];
                //w[15]
                coeff_temp  = -1.0/2*m_dt/rho*mu[4]/(top_h[0]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[1]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3,I_nb[15]*3+2, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[1]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U2_nb[15];
            }
            else if (bNoBoundary[1] == 3) // EAST is on BC
            {
                //hardwired for slip condition
                double temp1 = -1.0/3*m_dt/rho*mu[1]/(top_h[0]*top_h[2]);
                temp1 += 1.0/3*m_dt/rho*mu[0]/(top_h[0]*top_h[2]);
                double temp2 = -temp1;
                if (bNoBoundary[4] == 3) // LOWER added
                {
                    //w[c]
                    solver.Add_A(I*3,I*3+2, -temp1);
                    coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[5]/(top_h[0]*top_h[2]);
                    coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                    rhs += coeff_temp_old*U2_center;
                }
                else if (bNoBoundary[5] == 3) // UPPER added
                {
                    //w[4]
                    solver.Add_A(I*3,I_nb[4]*3+2, -temp2);
                    coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                    coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                    rhs += coeff_temp_old*U2_nb[4];
                }
                else
                {
                    //w[c]
                    solver.Add_A(I*3,I*3+2, -temp1);
                    coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[5]/(top_h[0]*top_h[2]);
                    coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                    rhs += coeff_temp_old*U2_center;

                    //w[4]
                    solver.Add_A(I*3,I_nb[4]*3+2, -temp2);
                    coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                    coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                    rhs += coeff_temp_old*U2_nb[4];
                }
            }
            if (bNoBoundary[4] == 2) //cells on LOWER bdry
            {
                //w[4] = w[15] = 0

                //w[1]
                coeff_temp  = 1.0/2*m_dt/rho*mu[5]/(top_h[0]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[1]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3,I_nb[1]*3+2, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[5]/(top_h[0]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[1]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U2_nb[1];

                //w[c]
                coeff_temp  = -1.0/2*m_dt/rho*mu[5]/(top_h[0]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[0]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3,I*3+2, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[5]/(top_h[0]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U2_center;
            }

            //add other terms to rhs
            rhs += m_dt*state.m_U[0];
//            rhs += m_dt*cell_center[index].m_state.f_surf[0];
            rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;
            rhs -= m_dt*cell_center[index].m_state.m_adv[0];

            solver.Add_b(I*3, rhs);


            /************************************************************************/
            //       Setting the coeffecients for U in the 2nd equation
            /************************************************************************/


            //set rho on v-face
            rho = 4.0/(1.0/cell_center[index].m_state.m_rho +
                       1.0/cell_center[index].m_state.m_rho_old +
                       1.0/cell_center[index_nb[3]].m_state.m_rho +
                       1.0/cell_center[index_nb[3]].m_state.m_rho_old);

            //set mu[0] ~ mu[5]
            if (I_nb[4]>=0 && I_nb[5]>=0) {
                mu_nb[11] = cell_center[index_nb[11]].m_state.m_mu;
                mu_nb[12] = cell_center[index_nb[12]].m_state.m_mu;
                mu_nb_old[11] = cell_center[index_nb[11]].m_state.m_mu_old;
                mu_nb_old[12] = cell_center[index_nb[12]].m_state.m_mu_old;
                if (useSGSCellCenter) {
                    mu_t_nb[11] = cell_center[index_nb[11]].m_state.m_mu_turbulent[3];
                    mu_t_nb[12] = cell_center[index_nb[12]].m_state.m_mu_turbulent[3];
                }
                else {
                    printf("codes needed for mu_turbulent on cell-faces.\n");
                    assert(false);
                }
            }
            else if (I_nb[4] < 0) { //cells on LOWER bdry
                mu_nb[11] = mu_nb[3];
                mu_nb[12] = cell_center[index_nb[12]].m_state.m_mu;
                mu_nb_old[11] = mu_nb_old[3];
                mu_nb_old[12] = cell_center[index_nb[12]].m_state.m_mu_old;
                if (useSGSCellCenter) {
                    mu_t_nb[11] = mu_t_nb[3];
                    mu_t_nb[12] = cell_center[index_nb[12]].m_state.m_mu_turbulent[3];
                }
                else {
                    printf("codes needed for mu_turbulent on cell-faces.\n");
                    assert(false);
                }
            }
            else if (I_nb[5] < 0) { //cells on UPPER bdry
                mu_nb[11] = cell_center[index_nb[11]].m_state.m_mu;
                mu_nb[12] = mu_nb[3];
                mu_nb_old[11] = cell_center[index_nb[11]].m_state.m_mu_old;
                mu_nb_old[12] = mu_nb_old[3];
                if (useSGSCellCenter) {
                    mu_t_nb[11] = cell_center[index_nb[11]].m_state.m_mu_turbulent[3];
                    mu_t_nb[12] = mu_t_nb[3];
                }
                else {
                    printf("codes needed for mu_turbulent on cell-faces.\n");
                    assert(false);
                }
            }

            mu_nb[8] = cell_center[index_nb[8]].m_state.m_mu;
            mu_nb[9] = cell_center[index_nb[9]].m_state.m_mu;
            mu_nb_old[8] = cell_center[index_nb[8]].m_state.m_mu_old;
            mu_nb_old[9] = cell_center[index_nb[9]].m_state.m_mu_old;
            if (useSGSCellCenter) {
                mu_t_nb[8] = cell_center[index_nb[8]].m_state.m_mu_turbulent[3];
                mu_t_nb[9] = cell_center[index_nb[9]].m_state.m_mu_turbulent[3];
            }
            else {
                printf("codes needed for mu_turbulent on cell-faces.\n");
                assert(false);
            }

            mu[0] = (mu_center + mu_nb[0] + mu_nb[9] + mu_nb[3] +
                     mu_t_center + mu_t_nb[0] + mu_t_nb[9] + mu_t_nb[3])/4;
            mu[1] = (mu_center + mu_nb[1] + mu_nb[8] + mu_nb[3] +
                     mu_t_center + mu_t_nb[1] + mu_t_nb[8] + mu_t_nb[3])/4;
            mu[2] = mu_center + mu_t_center;
            mu[3] = mu_nb[3] + mu_t_nb[3];
            mu[4] = (mu_center + mu_nb[3] + mu_nb[11] + mu_nb[4] +
                     mu_t_center + mu_t_nb[3] + mu_t_nb[11] + mu_t_nb[4])/4;
            mu[5] = (mu_center + mu_nb[3] + mu_nb[12] + mu_nb[5] +
                     mu_t_center + mu_t_nb[3] + mu_t_nb[12] + mu_t_nb[5])/4;
            mu_old[0] = (mu_center_old + mu_nb_old[0] + mu_nb_old[9] + mu_nb_old[3] +
                         mu_t_center + mu_t_nb[0] + mu_t_nb[9] + mu_t_nb[3])/4;
            mu_old[1] = (mu_center_old + mu_nb_old[1] + mu_nb_old[8] + mu_nb_old[3] +
                         mu_t_center + mu_t_nb[1] + mu_t_nb[8] + mu_t_nb[3])/4;
            mu_old[2] = mu_center_old + mu_t_center;
            mu_old[3] = mu_nb_old[3] + mu_t_nb[3];
            mu_old[4] = (mu_center_old + mu_nb_old[3] + mu_nb_old[11] + mu_nb_old[4] +
                         mu_t_center + mu_t_nb[3] + mu_t_nb[11] + mu_t_nb[4])/4;
            mu_old[5] = (mu_center_old + mu_nb_old[3] + mu_nb_old[12] + mu_nb_old[5] +
                         mu_t_center + mu_t_nb[3] + mu_t_nb[12] + mu_t_nb[5])/4;

            //set the coeffecients for U1 in the second equation, (mu*v_x)_x, (4/3*mu*v_y)_y and (mu*v_z)_z

            coeff1[0] = 0.5*m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            coeff1[1] = 0.5*m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            coeff1[2] = 2.0/3*m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            coeff1[3] = 2.0/3*m_dt/rho*mu[3]/(top_h[1]*top_h[1]);
            coeff1[4] = 0.5*m_dt/rho*mu[4]/(top_h[2]*top_h[2]);
            coeff1[5] = 0.5*m_dt/rho*mu[5]/(top_h[2]*top_h[2]);
            coeff1_old[0] = 0.5*m_dt/rho*mu_old[0]/(top_h[0]*top_h[0]);
            coeff1_old[1] = 0.5*m_dt/rho*mu_old[1]/(top_h[0]*top_h[0]);
            coeff1_old[2] = 2.0/3*m_dt/rho*mu_old[2]/(top_h[1]*top_h[1]);
            coeff1_old[3] = 2.0/3*m_dt/rho*mu_old[3]/(top_h[1]*top_h[1]);
            coeff1_old[4] = 0.5*m_dt/rho*mu_old[4]/(top_h[2]*top_h[2]);
            coeff1_old[5] = 0.5*m_dt/rho*mu_old[5]/(top_h[2]*top_h[2]);

            solver.Add_A(I*3+1,I*3+1,1.0);
            rhs = U1_center;
            /*
            for (nb = 0; nb < 6; ++nb)
            {
                if (I_nb[nb] >= 0) //interior cells
                {
                    solver.Add_A(I*3+1,I*3+1, coeff1[nb]);
                    solver.Add_A(I*3+1,I_nb[nb]*3+1, -coeff1[nb]);
                    rhs += coeff1_old[nb]*(U1_nb[nb] - U1_center);
                }
                else //cells on boundary
                {
                    solver.Add_A(I*3+1,I*3+1, 2.0*coeff1[nb]);
                    rhs += coeff1_old[nb]*(U1_nb[nb] - U1_center);
                }
            }*/

            //Y dir [SOUTH, NORTH] Wall normal to the v* velocity
            if (bNoBoundary[0] == 3) // WEST FACE is on REFLECT BC
            {
                if (SLIP)
                {
                    // do nothing for coeff matrix
                    rhs += coeff1_old[0]*(U1_nb[0] - U1_center);
                }
                else
                {
                    printf("NO SLIP CONDITION is not implmeneted yet in %s\n", __func__);
                    clean_up(ERROR);
                }
            }
            else if (!bNoBoundary[0]) // WEST FACE PERIODIC or INTERIOR
            {
                solver.Add_A(I*3+1,I*3+1, coeff1[0]);
                solver.Add_A(I*3+1,I_nb[0]*3+1, -coeff1[0]);
                rhs += coeff1_old[0]*(U1_nb[0] - U1_center);
            }
            if (bNoBoundary[1] == 3) // EAST FACE is on REFLECT BC
            {
                if (SLIP)
                {
                    // do nothing for coeff matrix
                    rhs += coeff1_old[1]*(U1_nb[1] - U1_center);
                }
                else
                {
                    printf("NO SLIP CONDITION is not implmeneted yet in %s\n", __func__);
                    clean_up(ERROR);
                }
            }
            else if (!bNoBoundary[1]) // EAST FACE PERIODIC or INTERIOR
            {
                solver.Add_A(I*3+1,I*3+1, coeff1[1]);
                solver.Add_A(I*3+1,I_nb[1]*3+1, -coeff1[1]);
                rhs += coeff1_old[1]*(U1_nb[1] - U1_center);
            }
            if (bNoBoundary[2] == 3) // SOUTH FACE is on REFLECT BC
            {
                //v[2] = 0
                solver.Add_A(I*3+1,I*3+1, coeff1[2]);
                rhs += coeff1_old[2]*(U1_nb[2] - U1_center);
            }
            else if (!bNoBoundary[2]) // SOUTH FACE PERIODIC or INTERIOR
            {
                solver.Add_A(I*3+1,I*3+1, coeff1[2]);
                solver.Add_A(I*3+1,I_nb[2]*3+1, -coeff1[2]);
                rhs += coeff1_old[2]*(U1_nb[2] - U1_center);
            }
            if (bNoBoundary[3] == 3) // NORTH FACE is on REFLECT BC
            {
                //v[c] = 0
                //v[3] = -v[2]
                solver.Add_A(I*3+1,I_nb[2]*3+1, coeff1[3]);
                rhs += coeff1_old[3]*(U1_nb[3] - U1_center);
            }
            else if (!bNoBoundary[3]) // NORTH FACE PERIODIC or INTERIOR
            {
                solver.Add_A(I*3+1,I*3+1, coeff1[3]);
                solver.Add_A(I*3+1,I_nb[3]*3+1, -coeff1[3]);
                rhs += coeff1_old[3]*(U1_nb[3] - U1_center);
            }
            if (bNoBoundary[4] == 3) // LOWER FACE is on REFLECT BC
            {
                if (SLIP)
                {
                    // do nothing for coeff matrix
                    rhs += coeff1_old[4]*(U1_nb[4] - U1_center);
                }
                else
                {
                    printf("NO SLIP CONDITION is not implmeneted yet in %s\n", __func__);
                    clean_up(ERROR);
                }
            }
            else if (!bNoBoundary[4]) // LOWER FACE PERIODIC or INTERIOR
            {
                solver.Add_A(I*3+1,I*3+1, coeff1[4]);
                solver.Add_A(I*3+1,I_nb[4]*3+1, -coeff1[4]);
                rhs += coeff1_old[4]*(U1_nb[4] - U1_center);
            }
            if (bNoBoundary[5] == 3) // UPPER FACE is on REFLECT BC
            {
                if (SLIP)
                {
                    // do nothing for coeff matrix
                    rhs += coeff1_old[5]*(U1_nb[5] - U1_center);
                }
                else
                {
                    printf("NO SLIP CONDITION is not implmeneted yet in %s\n", __func__);
                    clean_up(ERROR);
                }
            }
            else if (!bNoBoundary[5]) // UPPER FACE PERIODIC or INTERIOR
            {
                solver.Add_A(I*3+1,I*3+1, coeff1[5]);
                solver.Add_A(I*3+1,I_nb[5]*3+1, -coeff1[5]);
                rhs += coeff1_old[5]*(U1_nb[5] - U1_center);
            }
            // NEUMANN BC
            if(bNoBoundary[4] == 2)
            {
                solver.Add_A(I*3+1,I*3+1, 2.0*coeff1[nb]);
                rhs += coeff1_old[nb]*(U1_nb[nb] - U1_center);
            }
            if (bNoBoundary[5] == 2)
            {
                solver.Add_A(I*3+1,I*3+1, 2.0*coeff1[nb]);
                rhs += coeff1_old[nb]*(U1_nb[nb] - U1_center);
            }

            //set the coefficients for U0 in the second equation, (mu*u_y)_x and (-2/3)*(mu*u_x)_y
            //traverse the four cells for u, i.e. u_c, u_0, u_3, and u_9
            //consider NORTH, WEST, EAST

            if (!bNoBoundary[3] && !bNoBoundary[0] && !bNoBoundary[1])//PERIODIC or INTERIOR
            {
            //u[0]
            coeff_temp  = 1.0/2*m_dt/rho*mu[0]/(top_h[0]*top_h[1]);
            coeff_temp -= 1.0/3*m_dt/rho*mu[2]/(top_h[0]*top_h[1]);
            solver.Add_A(I*3+1,I_nb[0]*3, -coeff_temp);
            coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[0]/(top_h[0]*top_h[1]);
            coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[2]/(top_h[0]*top_h[1]);
            rhs += coeff_temp_old*U0_nb[0];

            //u[3]
            coeff_temp  = 1.0/2*m_dt/rho*mu[1]/(top_h[0]*top_h[1]);
            coeff_temp -= 1.0/3*m_dt/rho*mu[3]/(top_h[0]*top_h[1]);
            solver.Add_A(I*3+1,I_nb[3]*3, -coeff_temp);
            coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[1]/(top_h[0]*top_h[1]);
            coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[3]/(top_h[0]*top_h[1]);
            rhs += coeff_temp_old*U0_nb[3];

            U0_nb[9] = cell_center[index_nb[9]].m_state.m_U[0];
            //u[9]
            coeff_temp  = -1.0/2*m_dt/rho*mu[0]/(top_h[0]*top_h[1]);
            coeff_temp -= -1.0/3*m_dt/rho*mu[3]/(top_h[0]*top_h[1]);
            solver.Add_A(I*3+1,I_nb[9]*3, -coeff_temp);
            coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[0]/(top_h[0]*top_h[1]);
            coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[3]/(top_h[0]*top_h[1]);
            rhs += coeff_temp_old*U0_nb[9];

            //u[c]
            coeff_temp  = -1.0/2*m_dt/rho*mu[1]/(top_h[0]*top_h[1]);
            coeff_temp -= -1.0/3*m_dt/rho*mu[2]/(top_h[0]*top_h[1]);
            solver.Add_A(I*3+1,I*3, -coeff_temp);
            coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[1]/(top_h[0]*top_h[1]);
            coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[2]/(top_h[0]*top_h[1]);
            rhs += coeff_temp_old*U0_center;
            }
            else if ((bNoBoundary[0] == 3) && !bNoBoundary[3])//WEST FACE is on BC
            {
                // u[9] = 0 u[0] = 0

                //u[c]
                coeff_temp  = -1.0/2*m_dt/rho*mu[1]/(top_h[0]*top_h[1]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[2]/(top_h[0]*top_h[1]);
                solver.Add_A(I*3+1,I*3, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[1]/(top_h[0]*top_h[1]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[2]/(top_h[0]*top_h[1]);
                rhs += coeff_temp_old*U0_center;

                //u[3]
                coeff_temp  = 1.0/2*m_dt/rho*mu[1]/(top_h[0]*top_h[1]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[3]/(top_h[0]*top_h[1]);
                solver.Add_A(I*3+1,I_nb[3]*3, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[1]/(top_h[0]*top_h[1]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[3]/(top_h[0]*top_h[1]);
                rhs += coeff_temp_old*U0_nb[3];
            }
            else if ((bNoBoundary[1] == 3) && !bNoBoundary[3])//EAST FACE is on BC
            {
                //u[3] = 0 u[c] = 0

                U0_nb[9] = cell_center[index_nb[9]].m_state.m_U[0];
                //u[9]
                coeff_temp  = -1.0/2*m_dt/rho*mu[0]/(top_h[0]*top_h[1]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[3]/(top_h[0]*top_h[1]);
                solver.Add_A(I*3+1,I_nb[9]*3, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[0]/(top_h[0]*top_h[1]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[3]/(top_h[0]*top_h[1]);
                rhs += coeff_temp_old*U0_nb[9];

                //u[0]
                coeff_temp  = 1.0/2*m_dt/rho*mu[0]/(top_h[0]*top_h[1]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[2]/(top_h[0]*top_h[1]);
                solver.Add_A(I*3+1,I_nb[0]*3, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[0]/(top_h[0]*top_h[1]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[2]/(top_h[0]*top_h[1]);
                rhs += coeff_temp_old*U0_nb[0];
            }
            else if (bNoBoundary[3] == 3)// NORTH FACE is on BC
            {
                // hardwire for slip condition
                double temp1 = 1.0/3*m_dt/rho*mu[3]/(top_h[0]*top_h[1]);
                temp1 -= -1.0/3*m_dt/rho*mu[2]/(top_h[0]*top_h[1]);
                double temp2 = - temp1;
                if (bNoBoundary[0] == 3) // WEST added
                {
                    // u[c]
                    solver.Add_A(I*3+1,I*3, -temp2);
                    coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[1]/(top_h[0]*top_h[1]);
                    coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[2]/(top_h[0]*top_h[1]);
                    rhs += coeff_temp_old*U0_center;
                }
                else if (bNoBoundary[1] == 3) // EAST added
                {
                    // u[0]
                    solver.Add_A(I*3+1,I_nb[0]*3, -temp1);
                    coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[0]/(top_h[0]*top_h[1]);
                    coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[2]/(top_h[0]*top_h[1]);
                    rhs += coeff_temp_old*U0_nb[0];
                }
                else
                {
                    // u[c]
                    solver.Add_A(I*3+1,I*3, -temp2);
                    coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[1]/(top_h[0]*top_h[1]);
                    coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[2]/(top_h[0]*top_h[1]);
                    rhs += coeff_temp_old*U0_center;

                    // u[0]
                    solver.Add_A(I*3+1,I_nb[0]*3, -temp1);
                    coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[0]/(top_h[0]*top_h[1]);
                    coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[2]/(top_h[0]*top_h[1]);
                    rhs += coeff_temp_old*U0_nb[0];

                }
            }

            //set the coefficients for U2 in the second equation, (mu*w_y)_z and (-2/3)*(mu*w_z)_y
            //traverse the four cells for w, i.e. w_c, w_3, w_4, and w_11
            //consider LOWER, UPPER, EAST
            if (!bNoBoundary[4] && !bNoBoundary[5] && !bNoBoundary[1]) // PERIODIC or INTERIOR
            {
                //w[3]
                coeff_temp  = 1.0/2*m_dt/rho*mu[5]/(top_h[1]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[3]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+1,I_nb[3]*3+2, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[5]/(top_h[1]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[3]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U2_nb[3];

                U2_nb[4] = cell_center[index_nb[4]].m_state.m_U[2];
                //w[4]
                coeff_temp  = 1.0/2*m_dt/rho*mu[4]/(top_h[1]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[2]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+1,I_nb[4]*3+2, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U2_nb[4];

                U2_nb[11] = cell_center[index_nb[11]].m_state.m_U[2];
                //w[11]
                coeff_temp  = -1.0/2*m_dt/rho*mu[4]/(top_h[1]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[3]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+1,I_nb[11]*3+2, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[3]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U2_nb[11];

                //w[c]
                coeff_temp  = -1.0/2*m_dt/rho*mu[5]/(top_h[1]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[2]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+1,I*3+2, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[5]/(top_h[1]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U2_center;
            }
            else if ((bNoBoundary[4] >= 2) && !bNoBoundary[1]) // LOWER FACE is on BC
            {
                // w[4] = 0 w[11] = 0

                //w[c]
                coeff_temp  = -1.0/2*m_dt/rho*mu[5]/(top_h[1]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[2]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+1,I*3+2, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[5]/(top_h[1]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U2_center;

                //w[3]
                coeff_temp  = 1.0/2*m_dt/rho*mu[5]/(top_h[1]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[3]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+1,I_nb[3]*3+2, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[5]/(top_h[1]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[3]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U2_nb[3];
            }
            else if ((bNoBoundary[5] == 3) && !bNoBoundary[1]) // UPPER FACE is on BC
            {
                // w[c] = 0 w[3] = 0

                U2_nb[4] = cell_center[index_nb[4]].m_state.m_U[2];
                //w[4]
                coeff_temp  = 1.0/2*m_dt/rho*mu[4]/(top_h[1]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[2]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+1,I_nb[4]*3+2, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U2_nb[4];

                U2_nb[11] = cell_center[index_nb[11]].m_state.m_U[2];
                //w[11]
                coeff_temp  = -1.0/2*m_dt/rho*mu[4]/(top_h[1]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[3]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+1,I_nb[11]*3+2, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[3]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U2_nb[11];
            }
            else if (bNoBoundary[1] == 3) // EAST FACE
            {
                //hardwired for slip condition
                double temp1= 1.0/3*m_dt/rho*mu[3]/(top_h[1]*top_h[2]);
                temp1 -= 1.0/3*m_dt/rho*mu[2]/(top_h[1]*top_h[2]);
                double temp2 = - temp1;
                if (bNoBoundary[4] == 3) // LOWER added
                {
                    // w[c]
                    solver.Add_A(I*3+1,I*3+2, -temp2);
                    coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[5]/(top_h[1]*top_h[2]);
                    coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                    rhs += coeff_temp_old*U2_center;
                }
                else if (bNoBoundary[5] == 3) // UPPER added
                {
                    //w[4]
                    solver.Add_A(I*3+1,I_nb[4]*3+2, -temp1);
                    coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                    coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                    rhs += coeff_temp_old*U2_nb[4];
                }
                else
                {
                    // w[c]
                    solver.Add_A(I*3+1,I*3+2, -temp2);
                    coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[5]/(top_h[1]*top_h[2]);
                    coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                    rhs += coeff_temp_old*U2_center;

                    //w[4]
                    solver.Add_A(I*3+1,I_nb[4]*3+2, -temp1);
                    coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                    coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                    rhs += coeff_temp_old*U2_nb[4];
                }
            }

            //add other terms to rhs
            rhs += m_dt*state.m_U[1];
//            rhs += m_dt*cell_center[index].m_state.f_surf[1];
            rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;
            rhs -= m_dt*cell_center[index].m_state.m_adv[1];

            solver.Add_b(I*3+1, rhs);


            /************************************************************************/
            //         Setting the coeffecients of U for the 3rd equation
            /************************************************************************/

            //set rho on w-face
            //set mu[0] ~ mu[5]
            if (I_nb[5] < 0) //cells on UPPER bdry
            {
                rho = 2.0/(1.0/cell_center[index].m_state.m_rho +
                           1.0/cell_center[index].m_state.m_rho_old);

                mu_nb[12] = mu_nb[3];
                mu_nb[13] = mu_nb[2];
                mu_nb[16] = mu_nb[1];
                mu_nb[17] = mu_nb[0];
                mu_nb_old[12] = mu_nb_old[3];
                mu_nb_old[13] = mu_nb_old[2];
                mu_nb_old[16] = mu_nb_old[1];
                mu_nb_old[17] = mu_nb_old[0];
                if (useSGSCellCenter) {
                    mu_t_nb[12] = mu_t_nb[3];
                    mu_t_nb[13] = mu_t_nb[2];
                    mu_t_nb[16] = mu_t_nb[1];
                    mu_t_nb[17] = mu_t_nb[0];
                }
                else {
                    printf("codes needed for mu_turbulent on cell-faces.\n");
                    assert(false);
                }
            }
            else //other cells
            {
                rho = 4.0/(1.0/cell_center[index].m_state.m_rho +
                           1.0/cell_center[index].m_state.m_rho_old +
                           1.0/cell_center[index_nb[5]].m_state.m_rho +
                           1.0/cell_center[index_nb[5]].m_state.m_rho_old);

                mu_nb[12] = cell_center[index_nb[12]].m_state.m_mu;
                mu_nb[13] = cell_center[index_nb[13]].m_state.m_mu;
                mu_nb[16] = cell_center[index_nb[16]].m_state.m_mu;
                mu_nb[17] = cell_center[index_nb[17]].m_state.m_mu;
                mu_nb_old[12] = cell_center[index_nb[12]].m_state.m_mu_old;
                mu_nb_old[13] = cell_center[index_nb[13]].m_state.m_mu_old;
                mu_nb_old[16] = cell_center[index_nb[16]].m_state.m_mu_old;
                mu_nb_old[17] = cell_center[index_nb[17]].m_state.m_mu_old;
                if (useSGSCellCenter) {
                    mu_t_nb[12] = cell_center[index_nb[12]].m_state.m_mu_turbulent[3];
                    mu_t_nb[13] = cell_center[index_nb[13]].m_state.m_mu_turbulent[3];
                    mu_t_nb[16] = cell_center[index_nb[16]].m_state.m_mu_turbulent[3];
                    mu_t_nb[17] = cell_center[index_nb[17]].m_state.m_mu_turbulent[3];
                }
                else {
                    printf("codes needed for mu_turbulent on cell-faces.\n");
                    assert(false);
                }
            }

            mu[0] = (mu_center + mu_nb[0] + mu_nb[17] + mu_nb[5] +
                     mu_t_center + mu_t_nb[0] + mu_t_nb[17] + mu_t_nb[5])/4;
            mu[1] = (mu_center + mu_nb[1] + mu_nb[16] + mu_nb[5] +
                     mu_t_center + mu_t_nb[1] + mu_t_nb[16] + mu_t_nb[5])/4;
            mu[2] = (mu_center + mu_nb[2] + mu_nb[13] + mu_nb[5] +
                     mu_t_center + mu_t_nb[2] + mu_t_nb[13] + mu_t_nb[5])/4;
            mu[3] = (mu_center + mu_nb[3] + mu_nb[12] + mu_nb[5] +
                     mu_t_center + mu_t_nb[3] + mu_t_nb[12] + mu_t_nb[5])/4;
            mu[4] = mu_center + mu_t_center;
            mu[5] = mu_nb[5] + mu_t_nb[5];
            mu_old[0] = (mu_center_old + mu_nb_old[0] + mu_nb_old[17] + mu_nb_old[5] +
                         mu_t_center + mu_t_nb[0] + mu_t_nb[17] + mu_t_nb[5])/4;
            mu_old[1] = (mu_center_old + mu_nb_old[1] + mu_nb_old[16] + mu_nb_old[5] +
                         mu_t_center + mu_t_nb[1] + mu_t_nb[16] + mu_t_nb[5])/4;
            mu_old[2] = (mu_center_old + mu_nb_old[2] + mu_nb_old[13] + mu_nb_old[5] +
                         mu_t_center + mu_t_nb[2] + mu_t_nb[13] + mu_t_nb[5])/4;
            mu_old[3] = (mu_center_old + mu_nb_old[3] + mu_nb_old[12] + mu_nb_old[5] +
                         mu_t_center + mu_t_nb[3] + mu_t_nb[12] + mu_t_nb[5])/4;
            mu_old[4] = mu_center_old + mu_t_center;
            mu_old[5] = mu_nb_old[5] + mu_t_nb[5];

            //set the coeffecients for U2 in the third equation, (mu*w_x)_x, (mu*w_y)_y and (4/3*mu*w_z)_z

            coeff2[0] = 0.5*m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            coeff2[1] = 0.5*m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            coeff2[2] = 0.5*m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            coeff2[3] = 0.5*m_dt/rho*mu[3]/(top_h[1]*top_h[1]);
            coeff2[4] = 2.0/3*m_dt/rho*mu[4]/(top_h[2]*top_h[2]);
            coeff2[5] = 2.0/3*m_dt/rho*mu[5]/(top_h[2]*top_h[2]);
            coeff2_old[0] = 0.5*m_dt/rho*mu_old[0]/(top_h[0]*top_h[0]);
            coeff2_old[1] = 0.5*m_dt/rho*mu_old[1]/(top_h[0]*top_h[0]);
            coeff2_old[2] = 0.5*m_dt/rho*mu_old[2]/(top_h[1]*top_h[1]);
            coeff2_old[3] = 0.5*m_dt/rho*mu_old[3]/(top_h[1]*top_h[1]);
            coeff2_old[4] = 2.0/3*m_dt/rho*mu_old[4]/(top_h[2]*top_h[2]);
            coeff2_old[5] = 2.0/3*m_dt/rho*mu_old[5]/(top_h[2]*top_h[2]);

            solver.Add_A(I*3+2,I*3+2,1.0);
            rhs = U2_center;
            /*
            for (nb = 0; nb < 6; ++nb)
            {
                if (I_nb[nb] >= 0) //interior cells
                {
                    solver.Add_A(I*3+2,I*3+2, coeff2[nb]);
                    solver.Add_A(I*3+2,I_nb[nb]*3+2, -coeff2[nb]);
                    rhs += coeff2_old[nb]*(U2_nb[nb] - U2_center);
                }
                else //cells on boundary
                {
                    if (nb == 4) //LOWER bdry
                    {
                        //w[4] = 0

                        solver.Add_A(I*3+2,I*3+2, coeff2[nb]);
                        rhs += coeff2_old[nb]*(U2_nb[nb] - U2_center);
                    }
                    else if (nb == 5) //UPPER bdry
                    {
                        //w[5] = -w[4]

                        solver.Add_A(I*3+2,I*3+2, coeff2[nb]);
                        solver.Add_A(I*3+2,I_nb[4]*3+2, coeff2[nb]);
                        rhs += coeff2_old[nb]*(U2_nb[nb] - U2_center);
                    }
                    else assert(false);
                }
            }
            */
            // Z dir [LOWER, UPPER] Wall normal to the w* velocity
            if (bNoBoundary[0] == 3) // WEST FACE is on REFLECT BC
            {
                if (SLIP)
                {
                    //onWallvel = w[c]
                    //do nothing for coeff matrix
                    rhs += coeff2_old[0]*(U2_nb[0] - U2_center);
                }
                else
                {
                    printf("NO SLIP CONDITION was implemented.\n");
                    clean_up(ERROR);
                }
            }
            else // PERIODIC or INTERIOR
            {
                solver.Add_A(I*3+2,I*3+2, coeff2[0]);
                solver.Add_A(I*3+2,I_nb[0]*3+2, -coeff2[0]);
                rhs += coeff2_old[0]*(U2_nb[0] - U2_center);
            }
            if (bNoBoundary[1] == 3) // EAST FACE is on REFLECT BC
            {
                if (SLIP)
                {
                    //onWallvel = w[c]
                    //do nothing for coeff matrix
                    rhs += coeff2_old[1]*(U2_nb[1] - U2_center);
                }
                else
                {
                    printf("NO SLIP CONDITION was implemented.\n");
                    clean_up(ERROR);
                }
            }
            else // PERIODIC or INTERIOR
            {
                solver.Add_A(I*3+2,I*3+2, coeff2[1]);
                solver.Add_A(I*3+2,I_nb[1]*3+2, -coeff2[1]);
                rhs += coeff2_old[1]*(U2_nb[1] - U2_center);
            }
            if (bNoBoundary[2] == 3) // SOUTH FACE is on REFLECT BC
            {
                if (SLIP)
                {
                    //onWallvel = w[c]
                    //do nothing for coeff matrix
                    rhs += coeff2_old[2]*(U2_nb[2] - U2_center);
                }
                else
                {
                    printf("NO SLIP CONDITION was implemented.\n");
                    clean_up(ERROR);
                }
            }
            else // PERIODIC or INTERIOR
            {
                solver.Add_A(I*3+2,I*3+2, coeff2[1]);
                solver.Add_A(I*3+2,I_nb[2]*3+2, -coeff2[2]);
                rhs += coeff2_old[1]*(U2_nb[2] - U2_center);
            }
            if (bNoBoundary[3] == 3) // NORTH FACE is on REFLECT BC
            {
                if (SLIP)
                {
                    //onWallvel = w[c]
                    //do nothing for coeff matrix
                    rhs += coeff2_old[3]*(U2_nb[3] - U2_center);
                }
                else
                {
                    printf("NO SLIP CONDITION was implemented.\n");
                    clean_up(ERROR);
                }
            }
            else // PERIODIC or INTERIOR
            {
                solver.Add_A(I*3+2,I*3+2, coeff2[3]);
                solver.Add_A(I*3+2,I_nb[3]*3+2, -coeff2[3]);
                rhs += coeff2_old[3]*(U2_nb[3] - U2_center);
            }
            if (bNoBoundary[4] == 3) // LOWER FACE is on REFLECT BC
            {
                //w[4] = 0
                solver.Add_A(I*3+2,I*3+2, coeff2[4]);
                rhs += coeff2_old[4]*(U2_nb[4] - U2_center);
            }
            else
            {
                solver.Add_A(I*3+2,I*3+2, coeff2[3]);
                solver.Add_A(I*3+2,I_nb[3]*3+2, -coeff2[3]);
                rhs += coeff2_old[3]*(U2_nb[3] - U2_center);
            }
            if (bNoBoundary[5] == 3) // UPPER FACE is on REFLECT BC
            {
                //w[5] = -w[4]
                //w[c] = 0
                solver.Add_A(I*3+2,I_nb[4]*3+2, coeff2[5]);
                rhs += coeff2_old[5]*(U2_nb[5] - U2_center);
            }
            else
            {
                solver.Add_A(I*3+2,I*3+2, coeff2[5]);
                solver.Add_A(I*3+2,I_nb[5]*3+2, -coeff2[5]);
                rhs += coeff2_old[3]*(U2_nb[5] - U2_center);
            }
            //set the coefficients for U0 in the third equation, (mu*u_z)_x and (-2/3)*(mu*u_x)_z
            //traverse the four cells for u, i.e. u_c, u_0, u_5, and u_17
            //consider WEST, EAST and UPPER
            if(!bNoBoundary[0] && !bNoBoundary[1] && !bNoBoundary[5])// PERIODIC or INTERIOR
            {
                //u[0]
                coeff_temp  = 1.0/2*m_dt/rho*mu[0]/(top_h[0]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[4]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[0]*3, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U0_nb[0];

                U0_nb[5] = cell_center[index_nb[5]].m_state.m_U[0];
                //u[5]
                coeff_temp  = 1.0/2*m_dt/rho*mu[1]/(top_h[0]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[5]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[5]*3, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[1]/(top_h[0]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[5]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U0_nb[5];

                U0_nb[17] = cell_center[index_nb[17]].m_state.m_U[0];
                //u[17]
                coeff_temp  = -1.0/2*m_dt/rho*mu[0]/(top_h[0]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[5]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[17]*3, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[5]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U0_nb[17];

                //u[c]
                coeff_temp  = -1.0/2*m_dt/rho*mu[1]/(top_h[0]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[4]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3+2,I*3, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[1]/(top_h[0]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U0_center;
            }
            else if ((bNoBoundary[0] == 3) && !bNoBoundary[5])// ONLY WEST FACE is on BC
            {
                // u[17] = 0, u[0] = 0

                U0_nb[5] = cell_center[index_nb[5]].m_state.m_U[0];
                //u[5]
                coeff_temp  = 1.0/2*m_dt/rho*mu[1]/(top_h[0]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[5]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[5]*3, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[1]/(top_h[0]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[5]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U0_nb[5];

                //u[c]
                coeff_temp  = -1.0/2*m_dt/rho*mu[1]/(top_h[0]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[4]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3+2,I*3, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[1]/(top_h[0]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U0_center;
            }
            else if ((bNoBoundary[1] == 3) && !bNoBoundary[5])// ONLY EAST FACE is on BC
            {
                //u[5] = 0, u[c] = 0

                //u[0]
                coeff_temp  = 1.0/2*m_dt/rho*mu[0]/(top_h[0]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[4]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[0]*3, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U0_nb[0];

                U0_nb[17] = cell_center[index_nb[17]].m_state.m_U[0];
                //u[17]
                coeff_temp  = -1.0/2*m_dt/rho*mu[0]/(top_h[0]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[5]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[17]*3, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[5]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U0_nb[17];

            }
            else if (bNoBoundary[5] == 3) // UPPER FACE is on BC
            {
                //hardwired for slip condition
                double temp1 = -1.0/3*m_dt/rho*mu[5]/(top_h[0]*top_h[2]);
                temp1 += 1.0/3*m_dt/rho*mu[4]/(top_h[0]*top_h[2]);
                double temp2 = -temp1;
                if (bNoBoundary[0] == 3) // WEST added
                {
                    //u[c]
                    solver.Add_A(I*3+2,I*3, -temp1);
                    coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[1]/(top_h[0]*top_h[2]);
                    coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                    rhs += coeff_temp_old*U0_center;
                }
                else if (bNoBoundary[1] == 3) // EAST added
                {
                    //u[0]
                    solver.Add_A(I*3+2,I_nb[0]*3, -temp2);
                    coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                    coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                    rhs += coeff_temp_old*U0_nb[0];
                }
                else
                {
                    //u[c]
                    solver.Add_A(I*3+2,I*3, -temp1);
                    coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[1]/(top_h[0]*top_h[2]);
                    coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                    rhs += coeff_temp_old*U0_center;

                    //u[0]
                    solver.Add_A(I*3+2,I_nb[0]*3, -temp2);
                    coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                    coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                    rhs += coeff_temp_old*U0_nb[0];
                }
            }

            if (bNoBoundary[5] == 2) //cells on UPPER bdry
            {
                //u[5]=-u[c] & u[17]=-u[0]

                //u[0]
                coeff_temp  = 1.0/2*m_dt/rho*mu[0]/(top_h[0]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[4]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[0]*3, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U0_nb[0];

                U0_nb[5] = -U0_center;
                //u[5]
                coeff_temp  = 1.0/2*m_dt/rho*mu[1]/(top_h[0]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[5]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3+2,I*3, coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[1]/(top_h[0]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[5]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U0_nb[5];

                U0_nb[17] = -U0_nb[0];
                //u[17]
                coeff_temp  = -1.0/2*m_dt/rho*mu[0]/(top_h[0]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[5]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[0]*3, coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[0]/(top_h[0]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[5]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U0_nb[17];

                //u[c]
                coeff_temp  = -1.0/2*m_dt/rho*mu[1]/(top_h[0]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[4]/(top_h[0]*top_h[2]);
                solver.Add_A(I*3+2,I*3, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[1]/(top_h[0]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[4]/(top_h[0]*top_h[2]);
                rhs += coeff_temp_old*U0_center;
            }

            //set the coefficients for U1 in the third equation (mu*v_z)_y and (-2/3)*(mu*v_y)_z
            //traverse the four cells for v, i.e. v_c, v_2, v_5, and v_13
            //consider UPPER, SOUTH, NORTH
            if (!bNoBoundary[2] && !bNoBoundary[3] && !bNoBoundary[5]) // PERIODIC or INTERIOR
            {
                //v[2]
                coeff_temp  = 1.0/2*m_dt/rho*mu[2]/(top_h[1]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[4]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[2]*3+1, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U1_nb[2];

                U1_nb[5] = cell_center[index_nb[5]].m_state.m_U[1];
                //v[5]
                coeff_temp  = 1.0/2*m_dt/rho*mu[3]/(top_h[1]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[5]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[5]*3+1, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[3]/(top_h[1]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[5]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U1_nb[5];

                U1_nb[13] = cell_center[index_nb[13]].m_state.m_U[1];
                //v[13]
                coeff_temp  = -1.0/2*m_dt/rho*mu[2]/(top_h[1]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[5]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[13]*3+1, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[5]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U1_nb[13];

                //v[c]
                coeff_temp  = -1.0/2*m_dt/rho*mu[3]/(top_h[1]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[4]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+2,I*3+1, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[3]/(top_h[1]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U1_center;
            }
            else if ((bNoBoundary[2] == 3) && !bNoBoundary[5]) // ONLY SOUTH is on BC
            {
                //v[13] = 0, v[2] = 0

                U1_nb[5] = cell_center[index_nb[5]].m_state.m_U[1];
                //v[5]
                coeff_temp  = 1.0/2*m_dt/rho*mu[3]/(top_h[1]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[5]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[5]*3+1, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[3]/(top_h[1]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[5]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U1_nb[5];

                //v[c]
                coeff_temp  = -1.0/2*m_dt/rho*mu[3]/(top_h[1]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[4]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+2,I*3+1, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[3]/(top_h[1]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U1_center;
            }
            else if ((bNoBoundary[3] == 3) && !bNoBoundary[5]) // ONLY NORTH is on BC
            {
                //v[5] = 0, v[c] = 0

                U1_nb[13] = cell_center[index_nb[13]].m_state.m_U[1];
                //v[13]
                coeff_temp  = -1.0/2*m_dt/rho*mu[2]/(top_h[1]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[5]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[13]*3+1, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[5]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U1_nb[13];

                //v[2]
                coeff_temp  = 1.0/2*m_dt/rho*mu[2]/(top_h[1]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[4]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[2]*3+1, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U1_nb[2];
            }
            else if (bNoBoundary[5] == 3) // UPPER
            {
                double temp1 = -1.0/3*m_dt/rho*mu[5]/(top_h[1]*top_h[2]);
                temp1 += 1.0/3*m_dt/rho*mu[4]/(top_h[1]*top_h[2]);
                double temp2 = -temp1;
                if (bNoBoundary[2] == 3) // SOUTH added
                {
                    //v[c]
                    solver.Add_A(I*3+2,I*3+1, -temp1);
                    coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[3]/(top_h[1]*top_h[2]);
                    coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                    rhs += coeff_temp_old*U1_center;
                }
                else if (bNoBoundary[3] == 3) // NORTH added
                {
                    //v[2]
                    solver.Add_A(I*3+2,I_nb[2]*3+1, -temp1);
                    coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                    coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                    rhs += coeff_temp_old*U1_nb[2];
                }
                else
                {
                    //v[c]
                    solver.Add_A(I*3+2,I*3+1, -temp1);
                    coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[3]/(top_h[1]*top_h[2]);
                    coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                    rhs += coeff_temp_old*U1_center;

                    //v[2]
                    solver.Add_A(I*3+2,I_nb[2]*3+1, -temp1);
                    coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                    coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                    rhs += coeff_temp_old*U1_nb[2];
                }
            }
            if (bNoBoundary[5] == 2) //cells on UPPER bdry
            {
                //v[5]=-v[c] & v[13]=-v[2]

                //v[2]
                coeff_temp  = 1.0/2*m_dt/rho*mu[2]/(top_h[1]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[4]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[2]*3+1, -coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U1_nb[2];

                U1_nb[5] = -U1_center;
                //v[5]
                coeff_temp  = 1.0/2*m_dt/rho*mu[3]/(top_h[1]*top_h[2]);
                coeff_temp -= 1.0/3*m_dt/rho*mu[5]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+2,I*3+1, coeff_temp);
                coeff_temp_old  = 1.0/2*m_dt/rho*mu_old[3]/(top_h[1]*top_h[2]);
                coeff_temp_old -= 1.0/3*m_dt/rho*mu_old[5]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U1_nb[5];

                U1_nb[13] = -U1_nb[2];
                //v[13]
                coeff_temp  = -1.0/2*m_dt/rho*mu[2]/(top_h[1]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[5]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+2,I_nb[2]*3+1, coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[2]/(top_h[1]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[5]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U1_nb[13];

                //v[c]
                coeff_temp  = -1.0/2*m_dt/rho*mu[3]/(top_h[1]*top_h[2]);
                coeff_temp -= -1.0/3*m_dt/rho*mu[4]/(top_h[1]*top_h[2]);
                solver.Add_A(I*3+2,I*3+1, -coeff_temp);
                coeff_temp_old  = -1.0/2*m_dt/rho*mu_old[3]/(top_h[1]*top_h[2]);
                coeff_temp_old -= -1.0/3*m_dt/rho*mu_old[4]/(top_h[1]*top_h[2]);
                rhs += coeff_temp_old*U1_center;
            }

            //add other terms to rhs
            rhs += m_dt*state.m_U[2];
//            rhs += m_dt*cell_center[index].m_state.f_surf[2];
            rhs -= m_dt*cell_center[index].m_state.grad_q[2]/rho;
            rhs -= m_dt*cell_center[index].m_state.m_adv[2];

            solver.Add_b(I*3+2, rhs);
        } //loop for (i,j,k) ends

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

        start_clock("Before Petsc Solve");

        solver.Solve_GMRES();

        if (debugging("PETSc"))
        {
            double max, min;
            solver.GetExtremeSingularValues(&max, &min);
            (void) printf("The max singular value of A = %lf in "
                          "compDiffWithSmoothProperty_velocity_MAC_coupled_vd\n", max);
            (void) printf("The min singular value of A = %lf in "
                          "compDiffWithSmoothProperty_velocity_MAC_coupled_vd\n", min);
            if ( min != 0.0)
                (void) printf("The Cond Num of A = %lf in "
                              "compDiffWithSmoothProperty_velocity_MAC_coupled_vd\n", max/min);
        }

        solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);

        //if (rel_residual > 1)
        //{
        //    printf("\nThe solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
        //    solver.Reset_x();
        //    solver.Solve_GMRES();
        //    solver.GetNumIterations(&num_iter);
        //    solver.GetFinalRelativeResidualNorm(&rel_residual);
        //}

        stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cartesian::"
                          "compDiffWithSmoothProperty_velocity_MAC_coupled_vd: "
                          "num_iter = %d, rel_residual = %le. \n",num_iter,rel_residual);

        // store the value of U^n in U_tmp
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_tmp[0] = cell_center[index].m_state.m_U[0];
            cell_center[index].m_state.m_U_tmp[1] = cell_center[index].m_state.m_U[1];
            cell_center[index].m_state.m_U_tmp[2] = cell_center[index].m_state.m_U[2];
        }

        // store the value of U* in U
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
                cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
                cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
                u = fabs(cell_center[index].m_state.m_U[0]);
                v = fabs(cell_center[index].m_state.m_U[1]);
                w = fabs(cell_center[index].m_state.m_U[2]);
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
                        fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                {
                    max_speed = speed;
                    indmax[0] = i;
                    indmax[1] = j;
                    indmax[2] = k;
                }
                if (u > max_u)
                    max_u = u;
                if (v > max_v)
                    max_v = v;
                if (w > max_w)
                    max_w = w;
            }
            else
            {
                (void) printf("I[%d][%d][%d] = -1\n",i,j,k);
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
                cell_center[index].m_state.m_U[2] = 0.0;
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
                        fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                    max_speed = speed;
                if (u > max_u)
                    max_u = u;
                if (v > max_v)
                    max_v = v;
                if (w > max_w)
                    max_w = w;
            }
        }

        //print max and min of speed, u, v, and w
        max_tmp = max_speed;
        pp_global_max(&max_speed,1);
        if (debugging("step_size"))
	{
            (void) printf("local max_speed after compDiffWithSmoothProperty_velocity_MAC_coupled_vd "
                          "in cell(%d, %d, %d) of node #%d is: %lf\n",
                          indmax[0],indmax[1],indmax[2],pp_mynode(),max_tmp);
            if (max_tmp == max_speed)
                (void) printf("max_speed (locates in node #%d) is: %lf\n",pp_mynode(),max_speed);
	}

        max_tmp = max_u;
        pp_global_max(&max_u,1);
        if (debugging("step_size"))
        {
            (void) printf("local max_u after compDiffWithSmoothProperty_velocity_MAC_coupled_vd "
                          "of node #%d is: %lf\n",pp_mynode(),max_tmp);
            if (max_tmp == max_u)
                (void) printf("max_u (locates in node #%d) is: %lf\n",pp_mynode(),max_u);
	}

        max_tmp = max_v;
        pp_global_max(&max_v,1);
        if (debugging("step_size"))
        {
            (void) printf("local max_v after compDiffWithSmoothProperty_velocity_MAC_coupled_vd "
                          "of node #%d is: %lf\n",pp_mynode(),max_tmp);
            if (max_tmp == max_v)
                (void) printf("max_v (locates in node #%d) is: %lf\n",pp_mynode(),max_v);
	}

        max_tmp = max_w;
        pp_global_max(&max_w,1);
        if (debugging("step_size"))
        {
            (void) printf("local max_w after compDiffWithSmoothProperty_velocity_MAC_coupled_vd "
                          "of node #%d is: %lf\n",pp_mynode(),max_tmp);
            if (max_tmp == max_w)
        	(void) printf("max_w (locates in node #%d) is: %lf\n",pp_mynode(),max_w);
	}

        for (l = 0; l < 3; ++l)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U_tmp[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U_tmp[l] = array[index];
                vecarray[l][index] = array[index];
            }
        }
        /*
         * removal tag: HAOZ
         * */
        enforceReflectionState(vecarray); // on m_U_tmp
        for (l = 0; l < 3; l++)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_tmp[l] = vecarray[l][index];
        }
        for (l = 0; l < 3; ++l)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }

        //calc div_U for U*
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l];
        }
        // removal tag: HAOZ
        enforceReflectionState(vel);//m_U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U[l] = vel[l][index];
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel); //removal tag: adjust Reflection B.C. afterwards
        }
        //removal tag: HAOZ. adjust
        FT_ParallelExchGridArrayBuffer(source,front);
        // store the div(U*) in div_U
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        if (debugging("step_size"))
        {
            value = sum_div = 0;
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(source[index]);
                sum_div += source[index];
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U* is %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U* is %.16g\n",max_value);
            max_value = 0.0;
        }

        FT_FreeThese(1,x);
} /* end compDiffWithSmoothProperty_velocity_MAC_coupled_vd */


void Incompress_Solver_Smooth_3D_Cartesian::computeProjection_vd(void)
{
        int index;
        int i,j,k,l,icoords[MAXD];
        double P_max,P_min;
        double **vel = iFparams->field->vel;
        double sum_div;
        double value;
        double diffusion[5],diffusion_old[5];

        sum_div = 0.0;
        max_value = -1;

/*
        // project U*-U^n
        for (l = 0; l < dim; l++)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l] -
                            cell_center[index].m_state.m_U_tmp[l]; //(U*-U^n)
            vel[l][index] /= accum_dt; //(U*-U^n)/dt
        }

        // Compute velocity divergence
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv(icoords,vel); //div[(U*-U^n)/dt]
            getDivU_coupled_vd(icoords,diffusion,2);           //divergence constraint, i.e. S^{n+1}
            getDivU_coupled_vd(icoords,diffusion_old,3);       //divergence constraint, i.e. S^n
            source[index] -= (diffusion[3]-diffusion_old[3])/accum_dt;
                                                               //div[(U*-U^n)/dt]-[S^{n+1}-S^n]/dt
            diff_coeff[index] = cell_center[index].m_state.m_rho;
            diff_coeff_old[index] = cell_center[index].m_state.m_rho_old;
        }
*/

        // project U*, instead of (U*-U^n)
        for (l = 0; l < dim; l++)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l]; //U*
            vel[l][index] /= accum_dt; //(1/dt)*U*
        }

        // compute velocity divergence
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv(icoords,vel); //div((1/dt)*U*)
            getDivU_coupled_vd(icoords,diffusion,2);           //divergence constraint S^{n+1}
            source[index] -= diffusion[3]/accum_dt;            //div((1/dt)*U*)-S^{n+1}/dt
            diff_coeff[index] = cell_center[index].m_state.m_rho;
            diff_coeff_old[index] = cell_center[index].m_state.m_rho_old;
        }

        FT_ParallelExchGridArrayBuffer(source,front);
        FT_ParallelExchGridArrayBuffer(diff_coeff,front);
        FT_ParallelExchGridArrayBuffer(diff_coeff_old,front);

        poisson_solver3d_vd(front,ilower,iupper,ijk_to_I,source,diff_coeff,diff_coeff_old,
                        array,&P_max,&P_min);

        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_phi = array[index];
        }
} /* end computeProjection_vd */


void Incompress_Solver_Smooth_3D_Cartesian::computeProjection_MAC_vd(void)
{
        int index;
        int i,j,k,l,icoords[MAXD];
        double P_max,P_min;
        double **vel = iFparams->field->vel;
        double sum_div;
        double value;
        double diffusion;

        sum_div = 0;
        max_value = -1.0;

        //project U^{*}, instead of (U^{*}-U^{n})
        //get div(U^{*}) and divergence constraint S^{n+1}
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);

            source[index] = cell_center[index].m_state.div_U/accum_dt; //div(U^{*}/dt)
            getDivU_MAC_vd(icoords,&diffusion,2,bGhostCell); //divergence constraint S^{n+1}
            source[index] -= diffusion/accum_dt; //div(U^{*}/dt) - S^{n+1}/dt

            diff_coeff[index] = cell_center[index].m_state.m_rho;
            diff_coeff_old[index] = cell_center[index].m_state.m_rho_old;
        }

        FT_ParallelExchGridArrayBuffer(source,front);
        FT_ParallelExchGridArrayBuffer(diff_coeff,front);
        FT_ParallelExchGridArrayBuffer(diff_coeff_old,front);

        poisson_solver3d_vd(front,ilower,iupper,ijk_to_I,source,diff_coeff,diff_coeff_old,
                        array,&P_max,&P_min);

        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_phi = array[index];
        }
} /* end computeProjection_MAC_vd */


//solve PPE for pressure P^{n+1/2} and set q = P
void Incompress_Solver_Smooth_3D_Cartesian::computeProjection_MAC_PPE_vd(void)
{
} /* end computeProjection_MAC_PPE_vd */


void Incompress_Solver_Smooth_3D_Cartesian::computeNewVelocity_vd(void)
{
} /* end computeNewVelocity_vd */


// computeNewVelocity_fullMAC_vd
void Incompress_Solver_Smooth_3D_Cartesian::computeNewVelocity_fullMAC_vd(void)
{
        bool bNoBoundary;
        int i,j,k,l,index0,index,tol_num;
        double grad_phi[3],rho;
        COMPONENT comp;
        double speed,sum_div,value,diffusion;
        int icoords[MAXD];
        double **vel = iFparams->field->vel;
        int I, index_nb[6];
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        double t;
        max_speed = 0;

        //update U^{n+1} using U^{*}
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_phi;
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);

            comp = top_comp[index];
            if (!ifluid_comp(comp))
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
                cell_center[index].m_state.m_U[2] = 0.0;
                continue;
            }

            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            // UPPER
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
                    comp,&intfc_state,&hs,crx_coords,m_t_int) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                bNoBoundary = NO;
            else
                bNoBoundary = YES;

            if (!bNoBoundary) //cells on UPPER bdry
            {
                for (l = 0; l < 2; ++l)
                {
                    rho = 4.0/(1.0/cell_center[index].m_state.m_rho +
                               1.0/cell_center[index].m_state.m_rho_old +
                               1.0/cell_center[index_nb[2*l+1]].m_state.m_rho +
                               1.0/cell_center[index_nb[2*l+1]].m_state.m_rho_old);
                    grad_phi[l] = (cell_center[index_nb[2*l+1]].m_state.m_phi -
                                   cell_center[index].m_state.m_phi)/top_h[l];
                    cell_center[index].m_state.m_U[l] -= accum_dt/rho*grad_phi[l];
                }
                cell_center[index].m_state.m_U[2] = 0.0;
            }
            else //other cells
            {
                for (l = 0; l < 3; ++l)
                {
                    rho = 4.0/(1.0/cell_center[index].m_state.m_rho +
                               1.0/cell_center[index].m_state.m_rho_old +
                               1.0/cell_center[index_nb[2*l+1]].m_state.m_rho +
                               1.0/cell_center[index_nb[2*l+1]].m_state.m_rho_old);
                    grad_phi[l] = (cell_center[index_nb[2*l+1]].m_state.m_phi -
                                   cell_center[index].m_state.m_phi)/top_h[l];
                    cell_center[index].m_state.m_U[l] -= accum_dt/rho*grad_phi[l];
                }
            }

            speed = fabs(cell_center[index].m_state.m_U[0]) +
                    fabs(cell_center[index].m_state.m_U[1]) +
                    fabs(cell_center[index].m_state.m_U[2]);
            if (speed > max_speed)
                max_speed = speed;
        }
        //scatter states
        for (l = 0; l < 3; ++l)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }
        pp_global_max(&max_speed,1);

        //get div(U^{n+1}) and store the values in div_U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l];
        }
        //removal tag: HAOZ
        enforceReflectionState(vel);//m_U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U[l] = vel[l][index];
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);//removal tag: adjust Reflection B.C. afterwards
        }
        //removal tag: HAOZ. adjust
        FT_ParallelExchGridArrayBuffer(source,front);
        //store the values of div(U^{n+1}) in div_U
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        if (debugging("step_size"))
        {
            value = sum_div = max_value = 0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                sum_div += cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U^{n+1} "
                          "is: %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U^{n+1} "
                          "is: %.16g\n",max_value);
            max_value = 0;
        }

        if (debugging("step_size"))
        {
            value = sum_div = max_value = 0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;

                //get the divergence constraint S^{n+1}
                getDivU_MAC_vd(icoords,&diffusion,2,bGhostCell);
                source[index] = diffusion;
            }
            FT_ParallelExchGridArrayBuffer(source,front);

            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.div_U_tmp = source[index];
            }

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U_tmp);
                sum_div += cell_center[index].m_state.div_U_tmp;
                if(value > max_value)
                    max_value = value;
            }

            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of S^{n+1} is %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of S^{n+1} is %.16g\n",max_value);
            max_value = 0;
        }
} /* end computeNewVelocity_fullMAC_vd */


// computeNewVelocity_fullMAC_vd
void Incompress_Solver_Smooth_3D_Cartesian::computeNewVelocity_fullMAC_zeroW_vd(void)
{
        bool bNoBoundary[6];
        int i,j,k,l,index0,index,tol_num;
        double grad_phi[3],rho;
        COMPONENT comp;
        double speed,sum_div,value,diffusion;
        int icoords[MAXD];
        double **vel = iFparams->field->vel;

        int I, index_nb[6];
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        double t;
        max_speed = 0.0;

        //update U^{n+1} using U^{*}
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_phi;
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);

            comp = top_comp[index];
            if (!ifluid_comp(comp))
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
                cell_center[index].m_state.m_U[2] = 0.0;
                continue;
            }

            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            // 4 directions
            bNoBoundary[0] = YES;
            bNoBoundary[1] = YES;
            bNoBoundary[2] = YES;
            bNoBoundary[3] = YES;
            // LOWER
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
                    comp,&intfc_state,&hs,crx_coords,m_t_int) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                bNoBoundary[4] = NO;
            else
                bNoBoundary[4] = YES;
            // UPPER
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
                    comp,&intfc_state,&hs,crx_coords,m_t_int) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                bNoBoundary[5] = NO;
            else
                bNoBoundary[5] = YES;

            if (!bNoBoundary[5]) //cells on UPPER bdry
            {
                for (l = 0; l < 2; ++l)
                {
                    rho = 4.0/(1.0/cell_center[index].m_state.m_rho +
                               1.0/cell_center[index].m_state.m_rho_old +
                               1.0/cell_center[index_nb[2*l+1]].m_state.m_rho +
                               1.0/cell_center[index_nb[2*l+1]].m_state.m_rho_old);
                    grad_phi[l] = (cell_center[index_nb[2*l+1]].m_state.m_phi -
                                   cell_center[index].m_state.m_phi)/top_h[l];
                    cell_center[index].m_state.m_U[l] -= accum_dt/rho*grad_phi[l];
                }
                cell_center[index].m_state.m_U[2] = 0.0;
            }
            else //other cells
            {
                for (l = 0; l < 3; ++l)
                {
                    rho = 4.0/(1.0/cell_center[index].m_state.m_rho +
                               1.0/cell_center[index].m_state.m_rho_old +
                               1.0/cell_center[index_nb[2*l+1]].m_state.m_rho +
                               1.0/cell_center[index_nb[2*l+1]].m_state.m_rho_old);
                    grad_phi[l] = (cell_center[index_nb[2*l+1]].m_state.m_phi -
                                   cell_center[index].m_state.m_phi)/top_h[l];
                    cell_center[index].m_state.m_U[l] -= accum_dt/rho*grad_phi[l];
                }
            }

            //set w = 0
            cell_center[index].m_state.m_U[2] = 0.0;

            speed = fabs(cell_center[index].m_state.m_U[0]) +
                    fabs(cell_center[index].m_state.m_U[1]) +
                    fabs(cell_center[index].m_state.m_U[2]);
            if (speed > max_speed)
                max_speed = speed;
        }
        //scatter states
        for (l = 0; l < 3; ++l)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }
        pp_global_max(&max_speed,1);

        //get div(U^{n+1}) and store the values in div_U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l];
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);
        }
        FT_ParallelExchGridArrayBuffer(source,front);
        //store the values of div(U^{n+1}) in div_U
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        if(debugging("step_size"))
        {
            value = sum_div = 0.0;
            max_value = 0.0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                sum_div += cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U^{n+1} "
                          "is: %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U^{n+1} "
                          "is: %.16g\n",max_value);
            max_value = 0.0;
        }
        if (debugging("step_size"))
        {
            value = sum_div = 0.0;
            max_value = 0.0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;
                //get the divergence constraint S^{n+1}
                getDivU_MAC_vd(icoords,&diffusion,2,bGhostCell);
                source[index] = diffusion;
            }
            FT_ParallelExchGridArrayBuffer(source,front);

            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.div_U_tmp = source[index];
            }

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U_tmp);
                sum_div += cell_center[index].m_state.div_U_tmp;
                if(value > max_value)
                    max_value = value;
            }

            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of S^{n+1} is %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of S^{n+1} is %.16g\n",max_value);
            max_value = 0.0;
        }
} /* end computeNewVelocity_fullMAC_zeroW_vd */


// computeNewVelocity_fullMAC_vd
void Incompress_Solver_Smooth_3D_Cartesian::computeNewVelocity_fullMAC_zeroV_vd(void)
{
        bool bNoBoundary[6];
        int i,j,k,l,index0,index,tol_num;
        double grad_phi[3],rho;
        COMPONENT comp;
        double speed,sum_div,value,diffusion;
        int icoords[MAXD];
        double **vel = iFparams->field->vel;

        int I, index_nb[6];
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        double t;
        max_speed = 0.0;

        //update U^{n+1} using U^{*}
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_phi;
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);

            comp = top_comp[index];
            if (!ifluid_comp(comp))
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
                cell_center[index].m_state.m_U[2] = 0.0;
                continue;
            }

            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            // 4 directions
            bNoBoundary[0] = YES;
            bNoBoundary[1] = YES;
            bNoBoundary[2] = YES;
            bNoBoundary[3] = YES;
            // LOWER
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
                    comp,&intfc_state,&hs,crx_coords,m_t_int) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                bNoBoundary[4] = NO;
            else
                bNoBoundary[4] = YES;
            // UPPER
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
                    comp,&intfc_state,&hs,crx_coords,m_t_int) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                bNoBoundary[5] = NO;
            else
                bNoBoundary[5] = YES;

            if (!bNoBoundary[5]) //cells on UPPER bdry
            {
                for (l = 0; l < 2; ++l)
                {
                    rho = 4.0/(1.0/cell_center[index].m_state.m_rho +
                               1.0/cell_center[index].m_state.m_rho_old +
                               1.0/cell_center[index_nb[2*l+1]].m_state.m_rho +
                               1.0/cell_center[index_nb[2*l+1]].m_state.m_rho_old);
                    grad_phi[l] = (cell_center[index_nb[2*l+1]].m_state.m_phi -
                                   cell_center[index].m_state.m_phi)/top_h[l];
                    cell_center[index].m_state.m_U[l] -= accum_dt/rho*grad_phi[l];
                }
                cell_center[index].m_state.m_U[2] = 0.0;
            }
            else //other cells
            {
                for (l = 0; l < 3; ++l)
                {
                    rho = 4.0/(1.0/cell_center[index].m_state.m_rho +
                               1.0/cell_center[index].m_state.m_rho_old +
                               1.0/cell_center[index_nb[2*l+1]].m_state.m_rho +
                               1.0/cell_center[index_nb[2*l+1]].m_state.m_rho_old);
                    grad_phi[l] = (cell_center[index_nb[2*l+1]].m_state.m_phi -
                                   cell_center[index].m_state.m_phi)/top_h[l];
                    cell_center[index].m_state.m_U[l] -= accum_dt/rho*grad_phi[l];
                }
            }

            //set v = 0
            cell_center[index].m_state.m_U[1] = 0.0;

            speed = fabs(cell_center[index].m_state.m_U[0]) +
                    fabs(cell_center[index].m_state.m_U[1]) +
                    fabs(cell_center[index].m_state.m_U[2]);
            if (speed > max_speed)
                max_speed = speed;
        }
        //scatter states
        for (l = 0; l < 3; ++l)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }
        pp_global_max(&max_speed,1);

        //get div(U^{n+1}) and store the values in div_U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l];
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);
        }
        FT_ParallelExchGridArrayBuffer(source,front);
        //store the values of div(U^{n+1}) in div_U
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        if(debugging("step_size"))
        {
            value = sum_div = 0.0;
            max_value = 0.0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                sum_div += cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U^{n+1} "
                          "is: %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U^{n+1} "
                          "is: %.16g\n",max_value);
            max_value = 0.0;
        }
        if (debugging("step_size"))
        {
            value = sum_div = 0.0;
            max_value = 0.0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;
                //get the divergence constraint S^{n+1}
                getDivU_MAC_vd(icoords,&diffusion,2,bGhostCell);
                source[index] = diffusion;
            }
            FT_ParallelExchGridArrayBuffer(source,front);

            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.div_U_tmp = source[index];
            }

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U_tmp);
                sum_div += cell_center[index].m_state.div_U_tmp;
                if(value > max_value)
                    max_value = value;
            }

            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of S^{n+1} is %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of S^{n+1} is %.16g\n",max_value);
            max_value = 0.0;
        }
} /* end computeNewVelocity_fullMAC_zeroV_vd */


// ComputeNewVelocity_Symmetry_rhoFilter_vd
void Incompress_Solver_Smooth_3D_Cartesian::computeNewVelocity_Filter_vd(void)
{
} /* end computeNewVelocity_Filter_vd */


// Symmetric rho-weighted filter
void Incompress_Solver_Smooth_3D_Cartesian::compFilter_Ud_vd(void)
{
} /* end compFilter_Ud_vd */


void Incompress_Solver_Smooth_3D_Cartesian::compAdvectionTerm_coupled_vd(int flag)
{
} /* end compAdvectionTerm_coupled_vd */


void Incompress_Solver_Smooth_3D_Cartesian::compAdvectionTerm_decoupled_vd(int flag)
{
} /* end compAdvectionTerm_decoupled_vd */


//valid for variable dynamic viscosity
void Incompress_Solver_Smooth_3D_Cartesian::compAdvectionTerm_MAC_decoupled_vd(int flag, boolean bGhostCell)
{
    int I, i, j, k, l, index;
    int icoords[MAXD], iCoords[MAXD], index_nb[6];
    L_STATE state;
    double U[6], c[6], rho[6];
    double diffusion, value, sum_div, c_bar;
    double **vel = iFparams->field->vel;

    setIndexMap();

    if (!flag)
    {
        // (1) get U_center_bar on cell centers, and U_face_bar on cell faces
        // removal tag: HAOZ
        // REFLECTION B.C.
        // split the calculation of face vals and center vals from the triple loop.
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index = d_index3d(i,j,k,top_gmax);
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;

            getCellCenterVelocityBar_MAC_decoupled_vd(icoords,m_t_int);
            getCellFaceVelocityBar_MAC_decoupled_vd(icoords,m_t_int);
        }
        //scatter states
        for (l=0; l<dim; l++)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U_center_bar[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U_center_bar[l] = array[index];
                vecarray[l][index] = array[index];
            }

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U_face_bar[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U_face_bar[l] = array[index];
            }
        }

        //get div(U_face_bar) and store the values in div_U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U_face_bar[l];
        }
        // removal tag: HAOZ
        // REFLECTION B.C.
        enforceReflectionState(vecarray);//m_U_center_bar
        for (l = 0; l < 3; l++)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_center_bar[l] = vecarray[l][index];
        }
        enforceReflectionState(vel);//m_U_face_bar
        for (l = 0; l < 3; l++)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_face_bar[l] = vel[l][index];
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);//removal tag: adjust Reflection B.C. afterwards
        }
        // removal tag: HAOZ adjust
        FT_ParallelExchGridArrayBuffer(source,front);
        //store the values of div(U_face_bar) in div_U
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        if (debugging("step_size"))
        {
            value = sum_div = max_value = 0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                sum_div += cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U_face_bar "
                          "BEFORE MAC projection is: %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U_face_bar "
                          "BEFORE MAC projection is: %.16g\n",max_value);
            max_value = 0;
        }

        //get phi_mac_0 in MAC-projection step using div(U_face_bar) and S^{n}
        computeMacPhi_MAC_vd(flag);
        //get U_face_bar_0 using phi_mac_0 and rho^{n}
        computeNewUFaceBar_MAC_vd(m_t_int,flag);

        //get div(U_face_bar_0) and store the values in div_U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U_face_bar[l];
        }
        //removal tag: HAOZ
        enforceReflectionState(vel); // on m_U_face_bar
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_face_bar[l] = vel[l][index];
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);//removal tag: adjust Reflection B.C. afterwards
        }
        // removal tag: HAOZ adjust
        FT_ParallelExchGridArrayBuffer(source,front);
        //store div(U_face_bar_0) in div_U
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        if (debugging("step_size"))
        {
            value = sum_div = max_value = 0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                sum_div += cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U_face_bar "
                          "AFTER 1st MAC projection is: %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U_face_bar "
                          "AFTER 1st MAC projection is: %.16g\n",max_value);
            max_value = 0;
        }

        // (2) get U_edge_bar on cell edges
        // (3) get convection terms of NS eqn's
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index = d_index3d(i,j,k,top_gmax);
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            //get U_edge_bar and convection terms
            getCellEdgeVelocityBar_MAC_decoupled_vd(icoords,m_t_int);
        }
        //scatter states
        for (l=0; l<dim; l++)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_adv[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_adv[l] = array[index];
            }
        }

        if (debugging("step_size"))
        {
            value = sum_div = max_value = 0;

            for (l = 0; l < dim; ++l)
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                vel[l][index] = cell_center[index].m_state.m_adv[l];
            }
            /*
             * removal tag: HAOZ
             *
             * m_adv was treated as vel fashion.
             *
             * So, adv term need REFLECTION TREATMENT.
             * */
            enforceReflectionState(vel);// on m_adv terms
            for (l = 0; l < 3; l++)
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_adv[l] = vel[l][index];
            }
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;
                index = d_index3d(i,j,k,top_gmax);
                source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);//removal tag: adjust Reflection B.C. afterwards
            }
            // removal tag: HAOZ adjust
            FT_ParallelExchGridArrayBuffer(source,front);

            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.div_U_tmp = source[index];
            }

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U_tmp);
                sum_div += cell_center[index].m_state.div_U_tmp;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of N(U^{n+1/2}) "
                          "is: %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of N(U^{n+1/2}) "
                          "is: %.16g\n",max_value);
            max_value = 0;
        }

        // (4) get rho_bar and c_bar by U_face_bar_0
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index = d_index3d(i,j,k,top_gmax);
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            getScalarBar_MAC_vd(icoords,m_t_int,bGhostCell);
        }
        //scatter states
        for (l=0; l<6; l++)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_rho_bar[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_rho_bar[l] = array[index];
            }
        }
        for (l=0; l<6; l++)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_c_bar[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_c_bar[l] = array[index];
            }
        }

        // (5) get advection terms of continuity & concentration eqn's
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            cell_center[index].m_state.m_rho_adv = 0;
            cell_center[index].m_state.m_c_adv = 0;

            c_bar = 0;
            for (l = 0; l < 6; l++)
            {
                rho[l] = cell_center[index].m_state.m_rho_bar[l];
                c[l] = cell_center[index].m_state.m_c_bar[l];
                c_bar += c[l];
            }
            c_bar /= 6.0;

            U[0] = cell_center[index_nb[0]].m_state.m_U_face_bar[0];
            U[1] = cell_center[index].m_state.m_U_face_bar[0];
            U[2] = cell_center[index_nb[2]].m_state.m_U_face_bar[1];
            U[3] = cell_center[index].m_state.m_U_face_bar[1];
            U[5] = cell_center[index].m_state.m_U_face_bar[2];
            if (ijk_to_I[i][j][k-1] < 0) //cells on LOWER bdry
                U[4] = 0;
            else //other cells
                U[4] = cell_center[index_nb[4]].m_state.m_U_face_bar[2];

            //m_rho_adv = div(rho*U)
            cell_center[index].m_state.m_rho_adv += (rho[1]*U[1] - rho[0]*U[0])/top_h[0];
            cell_center[index].m_state.m_rho_adv += (rho[3]*U[3] - rho[2]*U[2])/top_h[1];
            cell_center[index].m_state.m_rho_adv += (rho[5]*U[5] - rho[4]*U[4])/top_h[2];

            //m_c_adv = div(c*U) - c*div(U)
            cell_center[index].m_state.m_c_adv += (c[1]*U[1] - c[0]*U[0])/top_h[0];
            cell_center[index].m_state.m_c_adv += (c[3]*U[3] - c[2]*U[2])/top_h[1];
            cell_center[index].m_state.m_c_adv += (c[5]*U[5] - c[4]*U[4])/top_h[2];
            cell_center[index].m_state.m_c_adv -= c_bar*cell_center[index].m_state.div_U;
        }
        //scatter states
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_rho_adv;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_rho_adv = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_c_adv;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_c_adv = array[index];
        }
    }
    else //flag == 1
    {
        //get phi_mac_1 in MAC-projection step using div(U_face_bar_0) and S^{n+1/2}
        computeMacPhi_MAC_vd(flag);
        //get U_face_bar_1 using phi_mac_1 and rho^{n}
        computeNewUFaceBar_MAC_vd(m_t_int,flag);

        //get div(U_face_bar_1) and store the values in div_U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U_face_bar[l];
        }
        // removal tag: HAOZ
        enforceReflectionState(vel); // on m_U_face_bar
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_face_bar[l] = vel[l][index];
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);//removal tag: adjust Reflection B.C. afterwards
        }
        // removal tag: HAOZ adjust
        FT_ParallelExchGridArrayBuffer(source,front);
        //store div(U_face_bar_1) in div_U
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        if (debugging("step_size"))
        {
            value = sum_div = max_value = 0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                sum_div += cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U_face_bar "
                          "AFTER 2nd MAC projection is: %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U_face_bar "
                          "AFTER 2nd MAC projection is: %.16g\n",max_value);
            max_value = 0;
        }

        //get S^{n+1/2}
        if (debugging("step_size"))
        {
            value = sum_div = max_value = 0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;

                getDivU_MAC_vd(icoords,&diffusion,flag,bGhostCell);
                source[index] = diffusion;
            }
            FT_ParallelExchGridArrayBuffer(source,front);

            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.div_U_tmp = source[index];
            }

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U_tmp);
                sum_div += cell_center[index].m_state.div_U_tmp;
                if(value > max_value)
                    max_value = value;
            }

            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of S^{n+1/2} is %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of S^{n+1/2} is %.16g\n",max_value);
            max_value = 0;
        }

        //set rho to be rho^n
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_rho = cell_center[index].m_state.m_rho_old;
        }

        //update rho_bar and c_bar by U_face_bar_1
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index = d_index3d(i,j,k,top_gmax);
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            getScalarBar_MAC_vd(icoords,m_t_int,bGhostCell);
        }
        //scatter states
        for (l=0; l<6; l++)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_rho_bar[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_rho_bar[l] = array[index];
            }
        }
        for (l=0; l<6; l++)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_c_bar[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_c_bar[l] = array[index];
            }
        }

        //get advection terms of continuity & concentration eqn's
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            cell_center[index].m_state.m_rho_adv = 0;
            cell_center[index].m_state.m_c_adv = 0;

            c_bar = 0;
            for (l=0; l<6; l++)
            {
                //TODO: use rho^(n+1/2) instead of updated rho_bar
                rho[l] = cell_center[index].m_state.m_rho_bar[l];
                c[l] = cell_center[index].m_state.m_c_bar[l];
                c_bar += c[l];
            }
            c_bar /= 6.0;

            U[0] = cell_center[index_nb[0]].m_state.m_U_face_bar[0];
            U[1] = cell_center[index].m_state.m_U_face_bar[0];
            U[2] = cell_center[index_nb[2]].m_state.m_U_face_bar[1];
            U[3] = cell_center[index].m_state.m_U_face_bar[1];
            U[5] = cell_center[index].m_state.m_U_face_bar[2];
            if (ijk_to_I[i][j][k-1] < 0) //cells on LOWER bdry
                U[4] = 0;
            else //other cells
                U[4] = cell_center[index_nb[4]].m_state.m_U_face_bar[2];

            cell_center[index].m_state.m_rho_adv += (rho[1]*U[1] - rho[0]*U[0])/top_h[0];
            cell_center[index].m_state.m_rho_adv += (rho[3]*U[3] - rho[2]*U[2])/top_h[1];
            cell_center[index].m_state.m_rho_adv += (rho[5]*U[5] - rho[4]*U[4])/top_h[2];

            //m_c_adv = div(cU) - c*div(U)
            cell_center[index].m_state.m_c_adv += (c[1]*U[1] - c[0]*U[0])/top_h[0];
            cell_center[index].m_state.m_c_adv += (c[3]*U[3] - c[2]*U[2])/top_h[1];
            cell_center[index].m_state.m_c_adv += (c[5]*U[5] - c[4]*U[4])/top_h[2];
            cell_center[index].m_state.m_c_adv -= c_bar*cell_center[index].m_state.div_U;
        }

        //scatter states
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_rho_adv;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_rho_adv = array[index];
        }
        //scatter states
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_c_adv;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_c_adv = array[index];
        }
    }
} /* end compAdvectionTerm_MAC_decoupled_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getVelocityBar_coupled_vd(
        int *icoords,
        double m_t_int)
{
} /* end getVelocityBar_coupled_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getVelocityBar_decoupled_vd(
        int *icoords,
        double m_t_int)
{
} /* end getVelocityBar_decoupled_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getCellCenterVelocityBar_MAC_decoupled_vd(
        int *icoords,
        double m_t_int)
{
    int bNoBoundary[6];//type change
    int ICoords[3];
    int index,index_nb[18];
    int i,j,k;
    L_STATE sl, sr;
    L_STATE state_center_bar;
    L_STATE state_center_hat;
    L_STATE state_center_hat_l, state_center_hat_r;

    COMPONENT comp;
    double w_nb, crx_coords[MAXD];

    int nb;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

    i = icoords[0];
    j = icoords[1];
    k = icoords[2];
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    comp = top_comp[index];

    //6 neighbours of the center cell
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    //xy cut neighbours
    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);

    //yz cut neighbours
    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);

    //xz cut neighbours
    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);

    //TODO && FIXME: Replace this with general formula, remove
    /*
    // 4 directions
    bNoBoundary[0] = YES;
    bNoBoundary[1] = YES;
    bNoBoundary[2] = YES;
    bNoBoundary[3] = YES;
    // LOWER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
            comp,&intfc_state,&hs,crx_coords,m_t_int) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[4] = NO;
    else
        bNoBoundary[4] = YES;
    // UPPER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
            comp,&intfc_state,&hs,crx_coords,m_t_int) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[5] = NO;
    else
        bNoBoundary[5] = YES;
    */
    for (nb = 0; nb < 6; nb++)
    {
        checkBoundaryCondition(dir[nb],icoords,&bNoBoundary[nb],m_t_int,comp);
    }


    //////////////////////// Get the U_hat & U_bar on cell centers //////////////////////////

    // WEST & EAST, get u_hat & u_bar
    ICoords[0] = icoords[0] - 1;
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2];
    if ((!bNoBoundary[0] && !bNoBoundary[1]) || bNoBoundary[0] == 3 || bNoBoundary[1] == 3)// PERIODIC or INTERIOR CELLs, bNoBoundary = 0
    {
        //u_hat
        getCenterVelocity_MAC_middleStep_hat_vd(ICoords,EAST,state_center_hat_l);
        getCenterVelocity_MAC_middleStep_hat_vd(icoords,WEST,state_center_hat_r);
        //u_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_X,EAST,sl,state_center_hat_l);
        getVelocity_MAC_middleStep_bar_decoupled_vd(icoords,COORD_X,WEST,sr,state_center_hat_r);
        getRiemannSolution_MAC_CenterVelocity_vd(COORD_X,sl,sr,state_center_bar,icoords);
    }
    else assert(false);;


    // SOUTH & NORTH, get v_hat & v_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1] - 1;
    ICoords[2] = icoords[2];
    if ((!bNoBoundary[2] && !bNoBoundary[3]) || bNoBoundary[2] == 3 || bNoBoundary[3] == 3)// PERIODIC or INTERIOR CELLs, bNoBoundary = 0
    {
        //v_hat
        getCenterVelocity_MAC_middleStep_hat_vd(ICoords,NORTH,state_center_hat_l);
        getCenterVelocity_MAC_middleStep_hat_vd(icoords,SOUTH,state_center_hat_r);
        //v_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Y,NORTH,sl,state_center_hat_l);
        getVelocity_MAC_middleStep_bar_decoupled_vd(icoords,COORD_Y,SOUTH,sr,state_center_hat_r);
        getRiemannSolution_MAC_CenterVelocity_vd(COORD_Y,sl,sr,state_center_bar,icoords);
    }
    else assert(false);


    // LOWER & UPPER, get w_hat & w_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2] - 1;
    if ((!bNoBoundary[4] && !bNoBoundary[5]) || bNoBoundary[4] == 3 || bNoBoundary[5] == 3)// PERIODIC or INTERIOR CELLs, bNoBoundary = 0
    {
        //w_hat
        getCenterVelocity_MAC_middleStep_hat_vd(ICoords,UPPER,state_center_hat_l);
        getCenterVelocity_MAC_middleStep_hat_vd(icoords,LOWER,state_center_hat_r);
        //w_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Z,UPPER,sl,state_center_hat_l);
        getVelocity_MAC_middleStep_bar_decoupled_vd(icoords,COORD_Z,LOWER,sr,state_center_hat_r);
        getRiemannSolution_MAC_CenterVelocity_vd(COORD_Z,sl,sr,state_center_bar,icoords);
    }
    else if (bNoBoundary[4] == 2) //cells on LOWER NEUMANN bdry no slip
    {
        //w_hat
        getCenterVelocity_MAC_middleStep_hat_vd(icoords,LOWER,state_center_hat_r);
        //w_bar
        //cells on the out side of LOWER bdry (i.e. icoords[2] = -1)
        w_nb = cell_center[index].m_state.m_U[2];
        sl.m_U[2] = 0.0 + top_h[2]/2.0*(w_nb - 0.0)/top_h[2] + 0.0;//no slip ONLY
        getVelocity_MAC_middleStep_bar_decoupled_vd(icoords,COORD_Z,LOWER,sr,state_center_hat_r);
        getRiemannSolution_MAC_CenterVelocity_vd(COORD_Z,sl,sr,state_center_bar,icoords);
    }
    else if (bNoBoundary[5] == 2) //cells on UPPER NEUMANN bdry no slip
    {
        //w_hat
        getCenterVelocity_MAC_middleStep_hat_vd(ICoords,UPPER,state_center_hat_l);
        //w_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Z,UPPER,sl,state_center_hat_l);
        //cells on UPPER bdry
        w_nb = cell_center[index_nb[4]].m_state.m_U[2];
        sr.m_U[2] = 0.0 - top_h[2]/2.0*(0.0 - w_nb)/top_h[2] + 0.0;
        getRiemannSolution_MAC_CenterVelocity_vd(COORD_Z,sl,sr,state_center_bar,icoords);
    }
    else assert(false);

    // store U_bar on cell centers
    cell_center[index].m_state.m_U_center_bar[0] = state_center_bar.m_U[0];
    cell_center[index].m_state.m_U_center_bar[1] = state_center_bar.m_U[1];
    cell_center[index].m_state.m_U_center_bar[2] = state_center_bar.m_U[2];
} /* end getCellCenterVelocityBar_MAC_decoupled_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getCellFaceVelocityBar_MAC_decoupled_vd(
        int *icoords,
        double m_t_int)
{
    int bNoBoundary[6];
    int index;
    int i,j,k,l;
    L_STATE state_face_bar[3];
    L_STATE state_face_hat[3];

    COMPONENT comp;
    double crx_coords[MAXD];
    int nb;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

    i = icoords[0];
    j = icoords[1];
    k = icoords[2];
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    comp = cell_center[index].comp;

    for (nb = 0; nb < 6; nb++)
    {
        checkBoundaryCondition(dir[nb],icoords,&bNoBoundary[nb],m_t_int,comp);
    }
    /*
    // 4 directions
    bNoBoundary[0] = YES;
    bNoBoundary[1] = YES;
    bNoBoundary[2] = YES;
    bNoBoundary[3] = YES;
    // LOWER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
            comp,&intfc_state,&hs,crx_coords,m_t_int) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[4] = NO;
    else
        bNoBoundary[4] = YES;
    // UPPER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
            comp,&intfc_state,&hs,crx_coords,m_t_int) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[5] = NO;
    else
        bNoBoundary[5] = YES;
    */

    //////////////////////// Get the U_hat & U_bar on cell faces //////////////////////////

    //u-faces, get u_hat & u_bar
    //u_hat
    getFaceVelocity_MAC_middleStep_hat_vd(icoords,COORD_X,state_face_hat[0]);
    //u_bar
    getVelocity_MAC_middleStep_bar_decoupled_vd(icoords,COORD_X,CENTER,state_face_bar[0],state_face_hat[0]);

    //v-faces, get v_hat & v_bar
    //v_hat
    getFaceVelocity_MAC_middleStep_hat_vd(icoords,COORD_Y,state_face_hat[1]);
    //v_bar
    getVelocity_MAC_middleStep_bar_decoupled_vd(icoords,COORD_Y,CENTER,state_face_bar[1],state_face_hat[1]);


    //w-faces, get w_hat & w_bar
    if (bNoBoundary[5]==2) //cells on UPPER bdry
        state_face_bar[2].m_U[2] = 0.0;
    else //other cells
    {
        //w_hat
        getFaceVelocity_MAC_middleStep_hat_vd(icoords,COORD_Z,state_face_hat[2]);
        //w_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(icoords,COORD_Z,CENTER,state_face_bar[2],state_face_hat[2]);
    }

    // store U_face_bar on cell faces
    for(l = 0; l < 3; ++l)
        cell_center[index].m_state.m_U_face_bar[l] = state_face_bar[l].m_U[l];
} /* end getCellFaceVelocityBar_MAC_decoupled_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getCellEdgeVelocityBar_MAC_decoupled_vd(
        int *icoords,
        double m_t_int)
{
    int bNoBoundary[6],nb;
    int ICoords[3],ICOORDs[3];
    int index,index_nb[18];
    int i, j, k;
    double dx, dy, dz;
    double tmp, tmp_nb;
    L_STATE sl, sr;
    L_STATE SL, SR;
    L_STATE ul, ur;
    L_STATE vl, vr;
    L_STATE wl, wr;
    L_STATE state_edge_bar[18];
    L_STATE state_edge_hat_l[18], state_edge_hat_r[18];
    COMPONENT comp;
    double crx_coords[MAXD];
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

    i = icoords[0];
    j = icoords[1];
    k = icoords[2];
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    comp = top_comp[index];

    //6 neighbours of the center cell
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    //xy cut neighbours
    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);

    //yz cut neighbours
    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);

    //xz cut neighbours
    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);

    for (nb = 0; nb < 6; nb++) // general way to determine bc type
    {
        checkBoundaryCondition(dir[nb],icoords,&bNoBoundary[nb],m_t_int,comp);
    }

    /*
    // 4 directions
    bNoBoundary[0] = YES;
    bNoBoundary[1] = YES;
    bNoBoundary[2] = YES;
    bNoBoundary[3] = YES;

    // LOWER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
            comp,&intfc_state,&hs,crx_coords,m_t_int) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[4] = NO;
    else
        bNoBoundary[4] = YES;

    // UPPER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
            comp,&intfc_state,&hs,crx_coords,m_t_int) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[5] = NO;
    else
        bNoBoundary[5] = YES;
    */


    //***********************************************************************************//
    //////////////////////// Get the U_hat & U_bar on cell edges //////////////////////////
    //***********************************************************************************//

    ///////////////////// Get the (u,v)_hat & (u,v)_bar on cell edge[7] ///////////////////
    // SOUTH & NORTH, get u_hat & u_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1] - 1;
    ICoords[2] = icoords[2];
    ICOORDs[0] = icoords[0];
    ICOORDs[1] = icoords[1];
    ICOORDs[2] = icoords[2];
    //u_hat
    getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_X,NORTH,state_edge_hat_l[7]);
    getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_X,SOUTH,state_edge_hat_r[7]);
    //u_bar
    getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_X,NORTH,ul,state_edge_hat_l[7]);
    getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_X,SOUTH,ur,state_edge_hat_r[7]);

    // WEST & EAST, get v_hat & v_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1] - 1;
    ICoords[2] = icoords[2];
    ICOORDs[0] = icoords[0] + 1;
    ICOORDs[1] = icoords[1] - 1;
    ICOORDs[2] = icoords[2];
    //v_hat
    getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_Y,EAST,state_edge_hat_l[7]);
    getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_Y,WEST,state_edge_hat_r[7]);
    //v_bar
    getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Y,EAST,vl,state_edge_hat_l[7]);
    getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_Y,WEST,vr,state_edge_hat_r[7]);

    //get (u,v)_bar on cell edges by solving an approx. RP, ambiguities????
    getRiemannSolution_MAC_EdgeVelocity_vd(COORD_X,COORD_Y,ul,ur,vl,vr,state_edge_bar[7]);


    ///////////////////// Get the (u,v)_hat & (u,v)_bar on cell edge[8] ///////////////////
    // SOUTH & NORTH, get u_hat & u_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2];
    ICOORDs[0] = icoords[0];
    ICOORDs[1] = icoords[1] + 1;
    ICOORDs[2] = icoords[2];
    //u_hat
    getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_X,NORTH,state_edge_hat_l[8]);
    getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_X,SOUTH,state_edge_hat_r[8]);
    //u_bar
    getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_X,NORTH,ul,state_edge_hat_l[8]);
    getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_X,SOUTH,ur,state_edge_hat_r[8]);

    // WEST & EAST, get v_hat & v_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2];
    ICOORDs[0] = icoords[0] + 1;
    ICOORDs[1] = icoords[1];
    ICOORDs[2] = icoords[2];
    //v_hat
    getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_Y,EAST,state_edge_hat_l[8]);
    getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_Y,WEST,state_edge_hat_r[8]);
    //v_bar
    getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Y,EAST,vl,state_edge_hat_l[8]);
    getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_Y,WEST,vr,state_edge_hat_r[8]);

    //get (u,v)_bar on cell edge by solving an approx. RP
    getRiemannSolution_MAC_EdgeVelocity_vd(COORD_X,COORD_Y,ul,ur,vl,vr,state_edge_bar[8]);


    ///////////////////// Get the (u,v)_hat & (u,v)_bar on cell edge[9] ///////////////////
    // SOUTH & NORTH, get u_hat & u_bar
    ICoords[0] = icoords[0] - 1;
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2];
    ICOORDs[0] = icoords[0] - 1;
    ICOORDs[1] = icoords[1] + 1;
    ICOORDs[2] = icoords[2];
    //u_hat
    getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_X,NORTH,state_edge_hat_l[9]);
    getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_X,SOUTH,state_edge_hat_r[9]);
    //u_bar
    getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_X,NORTH,ul,state_edge_hat_l[9]);
    getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_X,SOUTH,ur,state_edge_hat_r[9]);

    // WEST & EAST, get v_hat & v_bar
    ICoords[0] = icoords[0] - 1;
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2];
    ICOORDs[0] = icoords[0];
    ICOORDs[1] = icoords[1];
    ICOORDs[2] = icoords[2];
    //v_hat
    getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_Y,EAST,state_edge_hat_l[9]);
    getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_Y,WEST,state_edge_hat_r[9]);
    //v_bar
    getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Y,EAST,vl,state_edge_hat_l[9]);
    getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_Y,WEST,vr,state_edge_hat_r[9]);

    //get (u,v)_bar on cell edge by solving an approx. RP
    getRiemannSolution_MAC_EdgeVelocity_vd(COORD_X,COORD_Y,ul,ur,vl,vr,state_edge_bar[9]);


    ///////////////////// Get the (v,w)_hat & (v,w)_bar on cell edge[11] ///////////////////
    // LOWER & UPPER, get v_hat & v_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2] - 1;
    ICOORDs[0] = icoords[0];
    ICOORDs[1] = icoords[1];
    ICOORDs[2] = icoords[2];
    if (bNoBoundary[4] == 2) //cells on LOWER bdry
    {
        vl.m_U[1] = 0.0;
        vr.m_U[1] = 0.0;
    }
    else //other cells
    {
        //v_hat
        getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_Y,UPPER,state_edge_hat_l[11]);
        getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_Y,LOWER,state_edge_hat_r[11]);
        //v_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Y,UPPER,vl,state_edge_hat_l[11]);
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_Y,LOWER,vr,state_edge_hat_r[11]);
    }

    // SOUTH & NORTH, get w_hat & w_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2] - 1;
    ICOORDs[0] = icoords[0];
    ICOORDs[1] = icoords[1] + 1;
    ICOORDs[2] = icoords[2] - 1;
    if (bNoBoundary[4] == 2) //cells on LOWER bdry
    {
        wl.m_U[2] = 0.0;
        wr.m_U[2] = 0.0;
    }
    else //other cells
    {
        //w_hat
        getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_Z,NORTH,state_edge_hat_l[11]);
        getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_Z,SOUTH,state_edge_hat_r[11]);
        //w_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Z,NORTH,wl,state_edge_hat_l[11]);
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_Z,SOUTH,wr,state_edge_hat_r[11]);
    }

    //get (v,w)_bar on cell edge by solving an approx. RP
    getRiemannSolution_MAC_EdgeVelocity_vd(COORD_Y,COORD_Z,vl,vr,wl,wr,state_edge_bar[11]);


    ///////////////////// Get the (v,w)_hat & (v,w)_bar on cell edge[12] ///////////////////
    // LOWER & UPPER, get v_hat & v_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2];
    ICOORDs[0] = icoords[0];
    ICOORDs[1] = icoords[1];
    ICOORDs[2] = icoords[2] + 1;
    if (bNoBoundary[5] == 2) //cells on UPPER bdry
    {
        vl.m_U[1] = 0.0;
        vr.m_U[1] = 0.0;
    }
    else //other cells
    {
        //v_hat
        getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_Y,UPPER,state_edge_hat_l[12]);
        getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_Y,LOWER,state_edge_hat_r[12]);
        //v_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Y,UPPER,vl,state_edge_hat_l[12]);
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_Y,LOWER,vr,state_edge_hat_r[12]);
    }

    // SOUTH & NORTH, get w_hat & w_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2];
    ICOORDs[0] = icoords[0];
    ICOORDs[1] = icoords[1] + 1;
    ICOORDs[2] = icoords[2];
    if (bNoBoundary[5] == 2) //cells on UPPER bdry
    {
        wl.m_U[2] = 0.0;
        wr.m_U[2] = 0.0;
    }
    else //other cells
    {
        //w_hat
        getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_Z,NORTH,state_edge_hat_l[12]);
        getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_Z,SOUTH,state_edge_hat_r[12]);
        //w_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Z,NORTH,wl,state_edge_hat_l[12]);
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_Z,SOUTH,wr,state_edge_hat_r[12]);
    }

    //get (v,w)_bar on cell edge by solving an approx. RP
    getRiemannSolution_MAC_EdgeVelocity_vd(COORD_Y,COORD_Z,vl,vr,wl,wr,state_edge_bar[12]);


    ///////////////////// Get the (v,w)_hat & (v,w)_bar on cell edge[13] ///////////////////
    // LOWER & UPPER, get v_hat & v_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1] - 1;
    ICoords[2] = icoords[2];
    ICOORDs[0] = icoords[0];
    ICOORDs[1] = icoords[1] - 1;
    ICOORDs[2] = icoords[2] + 1;
    if (bNoBoundary[5] == 2) //cells on UPPER bdry
    {
        vl.m_U[1] = 0.0;
        vr.m_U[1] = 0.0;
    }
    else //other cells
    {
        //v_hat
        getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_Y,UPPER,state_edge_hat_l[13]);
        getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_Y,LOWER,state_edge_hat_r[13]);
        //v_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Y,UPPER,vl,state_edge_hat_l[13]);
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_Y,LOWER,vr,state_edge_hat_r[13]);
    }

    // SOUTH & NORTH, get w_hat & w_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1] - 1;
    ICoords[2] = icoords[2];
    ICOORDs[0] = icoords[0];
    ICOORDs[1] = icoords[1];
    ICOORDs[2] = icoords[2];
    if (bNoBoundary[5] == 2) //cells on UPPER bdry
    {
        wl.m_U[2] = 0.0;
        wr.m_U[2] = 0.0;
    }
    else //other cells
    {
        //w_hat
        getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_Z,NORTH,state_edge_hat_l[13]);
        getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_Z,SOUTH,state_edge_hat_r[13]);
        //w_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Z,NORTH,wl,state_edge_hat_l[13]);
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_Z,SOUTH,wr,state_edge_hat_r[13]);
    }

    //get (v,w)_bar on cell edge by solving an approx. RP
    getRiemannSolution_MAC_EdgeVelocity_vd(COORD_Y,COORD_Z,vl,vr,wl,wr,state_edge_bar[13]);


    ///////////////////// Get the (w,u)_hat & (w,u)_bar on cell edge[15] ///////////////////
    // WEST & EAST, get w_hat & w_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2] - 1;
    ICOORDs[0] = icoords[0] + 1;
    ICOORDs[1] = icoords[1];
    ICOORDs[2] = icoords[2] - 1;
    if (bNoBoundary[4] == 2) //cells on LOWER bdry
    {
        wl.m_U[2] = 0.0;
        wr.m_U[2] = 0.0;
    }
    else //other cells
    {
        //w_hat
        getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_Z,EAST,state_edge_hat_l[15]);
        getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_Z,WEST,state_edge_hat_r[15]);
        //w_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Z,EAST,wl,state_edge_hat_l[15]);
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_Z,WEST,wr,state_edge_hat_r[15]);
    }

    // LOWER & UPPER, get u_hat & u_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2] - 1;
    ICOORDs[0] = icoords[0];
    ICOORDs[1] = icoords[1];
    ICOORDs[2] = icoords[2];
    if (bNoBoundary[4] == 2) //cells on LOWER bdry
    {
        ul.m_U[0] = 0.0;
        ur.m_U[0] = 0.0;
    }
    else //other cells
    {
        //u_hat
        getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_X,UPPER,state_edge_hat_l[15]);
        getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_X,LOWER,state_edge_hat_r[15]);
        //u_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_X,UPPER,ul,state_edge_hat_l[15]);
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_X,LOWER,ur,state_edge_hat_r[15]);
    }

    //get (w,u)_bar on cell edge by solving an approx. RP
    getRiemannSolution_MAC_EdgeVelocity_vd(COORD_Z,COORD_X,wl,wr,ul,ur,state_edge_bar[15]);


    ///////////////////// Get the (w,u)_hat & (w,u)_bar on cell edge[16] ///////////////////
    // WEST & EAST, get w_hat & w_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2];
    ICOORDs[0] = icoords[0] + 1;
    ICOORDs[1] = icoords[1];
    ICOORDs[2] = icoords[2];
    if (bNoBoundary[5] == 2) //cells on UPPER bdry
    {
        wl.m_U[2] = 0.0;
        wr.m_U[2] = 0.0;
    }
    else //other cells
    {
        //w_hat
        getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_Z,EAST,state_edge_hat_l[16]);
        getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_Z,WEST,state_edge_hat_r[16]);
        //w_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Z,EAST,wl,state_edge_hat_l[16]);
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_Z,WEST,wr,state_edge_hat_r[16]);
    }

    // LOWER & UPPER, get u_hat & u_bar
    ICoords[0] = icoords[0];
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2];
    ICOORDs[0] = icoords[0];
    ICOORDs[1] = icoords[1];
    ICOORDs[2] = icoords[2] + 1;
    if (bNoBoundary[5] == 2) //cells on UPPER bdry
    {
        ul.m_U[0] = 0.0;
        ur.m_U[0] = 0.0;
    }
    else //other cells
    {
        //u_hat
        getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_X,UPPER,state_edge_hat_l[16]);
        getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_X,LOWER,state_edge_hat_r[16]);
        //u_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_X,UPPER,ul,state_edge_hat_l[16]);
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_X,LOWER,ur,state_edge_hat_r[16]);
    }

    //get (w,u)_bar on cell edge by solving an approx. RP
    getRiemannSolution_MAC_EdgeVelocity_vd(COORD_Z,COORD_X,wl,wr,ul,ur,state_edge_bar[16]);


    ///////////////////// Get the (w,u)_hat & (w,u)_bar on cell edge[17] ///////////////////
    // WEST & EAST, get w_hat & w_bar
    ICoords[0] = icoords[0] - 1;
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2];
    ICOORDs[0] = icoords[0];
    ICOORDs[1] = icoords[1];
    ICOORDs[2] = icoords[2];
    if (bNoBoundary[5] == 2) //cells on UPPER bdry
    {
        wl.m_U[2] = 0.0;
        wr.m_U[2] = 0.0;
    }
    else //other cells
    {
        //w_hat
        getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_Z,EAST,state_edge_hat_l[17]);
         getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_Z,WEST,state_edge_hat_r[17]);
        //w_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Z,EAST,wl,state_edge_hat_l[17]);
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_Z,WEST,wr,state_edge_hat_r[17]);
    }

    // LOWER & UPPER, get u_hat & u_bar
    ICoords[0] = icoords[0] - 1;
    ICoords[1] = icoords[1];
    ICoords[2] = icoords[2];
    ICOORDs[0] = icoords[0] - 1;
    ICOORDs[1] = icoords[1];
    ICOORDs[2] = icoords[2] + 1;
    if (bNoBoundary[5] == 2) //cells on UPPER bdry
    {
        ul.m_U[0] = 0.0;
        ur.m_U[0] = 0.0;
    }
    else //other cells
    {
        //u_hat
        getEdgeVelocity_MAC_middleStep_hat_vd(ICoords,COORD_X,UPPER,state_edge_hat_l[17]);
        getEdgeVelocity_MAC_middleStep_hat_vd(ICOORDs,COORD_X,LOWER,state_edge_hat_r[17]);
        //u_bar
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_X,UPPER,ul,state_edge_hat_l[17]);
        getVelocity_MAC_middleStep_bar_decoupled_vd(ICOORDs,COORD_X,LOWER,ur,state_edge_hat_r[17]);
    }

    //get (w,u)_bar on cell edge by solving an approx. RP
    getRiemannSolution_MAC_EdgeVelocity_vd(COORD_Z,COORD_X,wl,wr,ul,ur,state_edge_bar[17]);


    //***********************************************************************************//
    //                           get convection terms of NS eqn's
    //***********************************************************************************//
    dx = top_h[0];
    dy = top_h[1];
    dz = top_h[2];
    cell_center[index].m_state.m_adv[0] = 0.0;
    cell_center[index].m_state.m_adv[1] = 0.0;
    cell_center[index].m_state.m_adv[2] = 0.0;

    //u*u_x + v*u_y + w*u_z on u-faces
    tmp = cell_center[index].m_state.m_U_center_bar[0];
    tmp_nb = cell_center[index_nb[1]].m_state.m_U_center_bar[0];
    cell_center[index].m_state.m_adv[0] += 0.5*(tmp_nb + tmp)*(tmp_nb - tmp)/dx;
    cell_center[index].m_state.m_adv[0] += 0.5*(state_edge_bar[8].m_U[1] + state_edge_bar[7].m_U[1])*
                                               (state_edge_bar[8].m_U[0] - state_edge_bar[7].m_U[0])/dy;
    cell_center[index].m_state.m_adv[0] += 0.5*(state_edge_bar[16].m_U[2] + state_edge_bar[15].m_U[2])*
                                               (state_edge_bar[16].m_U[0] - state_edge_bar[15].m_U[0])/dz;

    //u*v_x + v*v_y + w*v_z on v-faces
    tmp = cell_center[index].m_state.m_U_center_bar[1];
    tmp_nb = cell_center[index_nb[3]].m_state.m_U_center_bar[1];
    cell_center[index].m_state.m_adv[1] += 0.5*(state_edge_bar[8].m_U[0] + state_edge_bar[9].m_U[0])*
                                               (state_edge_bar[8].m_U[1] - state_edge_bar[9].m_U[1])/dx;
    cell_center[index].m_state.m_adv[1] += 0.5*(tmp_nb + tmp)*(tmp_nb - tmp)/dy;
    cell_center[index].m_state.m_adv[1] += 0.5*(state_edge_bar[12].m_U[2] + state_edge_bar[11].m_U[2])*
                                               (state_edge_bar[12].m_U[1] - state_edge_bar[11].m_U[1])/dz;

    //u*w_x + v*w_y + w*w_z on w-faces
    tmp = cell_center[index].m_state.m_U_center_bar[2];
    tmp_nb = cell_center[index_nb[5]].m_state.m_U_center_bar[2];
    if (bNoBoundary[5] == 2) //cells on UPPER bdry
        cell_center[index].m_state.m_adv[2] += 0.0;
    else //other cells
    {
        cell_center[index].m_state.m_adv[2] += 0.5*(state_edge_bar[16].m_U[0] + state_edge_bar[17].m_U[0])*
                                                   (state_edge_bar[16].m_U[2] - state_edge_bar[17].m_U[2])/dx;
        cell_center[index].m_state.m_adv[2] += 0.5*(state_edge_bar[12].m_U[1] + state_edge_bar[13].m_U[1])*
                                                   (state_edge_bar[12].m_U[2] - state_edge_bar[13].m_U[2])/dy;
        cell_center[index].m_state.m_adv[2] += 0.5*(tmp_nb + tmp)*(tmp_nb - tmp)/dz;
    }
} /* end getCellEdgeVelocityBar_MAC_decoupled_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getScalarBar_coupled_vd(
        int *icoords,
        double m_t_int)
{
} /* end getScalarBar_coupled_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getScalarBar_MAC_vd(
        int *icoords,
        double m_t_int,
        boolean bGhostCell)
{
    int bNoBoundary[6],nb;
    int ICoords[3];
    int index;
    L_STATE sl, sr, state_west_bar, state_east_bar, state_south_bar, state_north_bar, state_lower_bar, state_upper_bar;
    L_STATE state_west_hat, state_east_hat, state_south_hat, state_north_hat, state_lower_hat, state_upper_hat;
    L_STATE state_west_hat_l, state_west_hat_r;
    L_STATE state_east_hat_l, state_east_hat_r;
    L_STATE state_south_hat_l, state_south_hat_r;
    L_STATE state_north_hat_l, state_north_hat_r;
    L_STATE state_lower_hat_l, state_lower_hat_r;
    L_STATE state_upper_hat_l, state_upper_hat_r;
    COMPONENT comp;
    double crx_coords[MAXD];
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    comp = top_comp[index];

    for (nb = 0; nb < 6; nb++)
    {
        checkBoundaryCondition(dir[nb],icoords,&bNoBoundary[nb],m_t_int,comp);
        if (bNoBoundary[nb] == 1) // no DIRICHLET BC
            clean_up(ERROR);
    }

    /*
    // 4 directions
    bNoBoundary[0] = YES;
    bNoBoundary[1] = YES;
    bNoBoundary[2] = YES;
    bNoBoundary[3] = YES;

    // LOWER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
            comp,&intfc_state,&hs,crx_coords,m_t_int) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[4] = NO;
    else
        bNoBoundary[4] = YES;

    // UPPER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
            comp,&intfc_state,&hs,crx_coords,m_t_int) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[5] = NO;
    else
        bNoBoundary[5] = YES;
    */

    //////////////////////// Get the scalar_hat on six surfaces first //////////////////////////

    // WEST
    if (bNoBoundary[0]>=2)
    {
        getFaceScalar_MAC_middleStep_hat_vd(icoords,COORD_X,WEST,state_west_hat,bGhostCell);
    }
    else
    {
        ICoords[0] = icoords[0] - 1;
        ICoords[1] = icoords[1];
        ICoords[2] = icoords[2];
        getFaceScalar_MAC_middleStep_hat_vd(ICoords,COORD_X,EAST,state_west_hat_l,bGhostCell);
        getFaceScalar_MAC_middleStep_hat_vd(icoords,COORD_X,WEST,state_west_hat_r,bGhostCell);
    }

    // EAST
    if (bNoBoundary[1]>=2)
    {
        getFaceScalar_MAC_middleStep_hat_vd(icoords,COORD_X,EAST,state_east_hat,bGhostCell);
    }
    else
    {
        ICoords[0] = icoords[0] + 1;
        ICoords[1] = icoords[1];
        ICoords[2] = icoords[2];
        getFaceScalar_MAC_middleStep_hat_vd(ICoords,COORD_X,WEST,state_east_hat_r,bGhostCell);
        getFaceScalar_MAC_middleStep_hat_vd(icoords,COORD_X,EAST,state_east_hat_l,bGhostCell);
    }

    // SOUTH
    if (bNoBoundary[2]>=2)
    {
        getFaceScalar_MAC_middleStep_hat_vd(icoords,COORD_Y,SOUTH,state_east_hat,bGhostCell);
    }
    else
    {
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1] - 1;
        ICoords[2] = icoords[2];
        getFaceScalar_MAC_middleStep_hat_vd(ICoords,COORD_Y,NORTH,state_south_hat_l,bGhostCell);
        getFaceScalar_MAC_middleStep_hat_vd(icoords,COORD_Y,SOUTH,state_south_hat_r,bGhostCell);
    }

    // NORTH
    if (bNoBoundary[2]>=2)
    {
        getFaceScalar_MAC_middleStep_hat_vd(icoords,COORD_Y,NORTH,state_east_hat,bGhostCell);
    }
    else
    {
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1] + 1;
        ICoords[2] = icoords[2];
        getFaceScalar_MAC_middleStep_hat_vd(ICoords,COORD_Y,SOUTH,state_north_hat_r,bGhostCell);
        getFaceScalar_MAC_middleStep_hat_vd(icoords,COORD_Y,NORTH,state_north_hat_l,bGhostCell);
    }

    // LOWER
    if(bNoBoundary[4]>=2) //Neumann BC for scalar_hat
    {
        getFaceScalar_MAC_middleStep_hat_vd(icoords,COORD_Z,LOWER,state_lower_hat,bGhostCell);
    }
    else
    {
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1];
        ICoords[2] = icoords[2] - 1;
        getFaceScalar_MAC_middleStep_hat_vd(ICoords,COORD_Z,UPPER,state_lower_hat_l,bGhostCell);
        getFaceScalar_MAC_middleStep_hat_vd(icoords,COORD_Z,LOWER,state_lower_hat_r,bGhostCell);
    }

    // UPPER
    if(bNoBoundary[5]>=2) //Neumann BC for scalar_hat
    {
        getFaceScalar_MAC_middleStep_hat_vd(icoords,COORD_Z,UPPER,state_upper_hat,bGhostCell);
    }
    else
    {
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1];
        ICoords[2] = icoords[2] + 1;
        getFaceScalar_MAC_middleStep_hat_vd(ICoords,COORD_Z,LOWER,state_upper_hat_r,bGhostCell);
        getFaceScalar_MAC_middleStep_hat_vd(icoords,COORD_Z,UPPER,state_upper_hat_l,bGhostCell);
    }

    ///////////////////////////// get the scalar_bar on six surfaces ////////////////////////////////

    // WEST
    if(bNoBoundary[0]>=2) //homogeneous Neumann BC for scalar
    {
        getScalar_MAC_middleStep_bar_decoupled_vd(icoords,COORD_X,WEST,state_west_bar,state_west_hat,bGhostCell);
    }
    else
    {
        ICoords[0] = icoords[0] - 1;
        ICoords[1] = icoords[1];
        ICoords[2] = icoords[2];
        getScalar_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_X,EAST,sl,state_west_hat_l,bGhostCell); //stability?????
        getScalar_MAC_middleStep_bar_decoupled_vd(icoords,COORD_X,WEST,sr,state_west_hat_r,bGhostCell);
        getRiemannSolution_MAC_Scalar_vd(sl,sr,state_west_bar,icoords,WEST);
    }

    // EAST
    if(bNoBoundary[1]>=2) //homogeneous Neumann BC for scalar
    {
        getScalar_MAC_middleStep_bar_decoupled_vd(icoords,COORD_X,EAST,state_east_bar,state_east_hat,bGhostCell);
    }
    else
    {
        ICoords[0] = icoords[0] + 1;
        ICoords[1] = icoords[1];
        ICoords[2] = icoords[2];
        getScalar_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_X,WEST,sr,state_east_hat_r,bGhostCell);
        getScalar_MAC_middleStep_bar_decoupled_vd(icoords,COORD_X,EAST,sl,state_east_hat_l,bGhostCell);
        getRiemannSolution_MAC_Scalar_vd(sl,sr,state_east_bar,icoords,EAST);
    }

    // SOUTH
    if(bNoBoundary[2]>=2) //homogeneous Neumann BC for scalar
    {
        getScalar_MAC_middleStep_bar_decoupled_vd(icoords,COORD_Y,SOUTH,state_south_bar,state_south_hat,bGhostCell);
    }
    else
    {
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1] - 1;
        ICoords[2] = icoords[2];
        getScalar_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Y,NORTH,sl,state_south_hat_l,bGhostCell);
        getScalar_MAC_middleStep_bar_decoupled_vd(icoords,COORD_Y,SOUTH,sr,state_south_hat_r,bGhostCell);
        getRiemannSolution_MAC_Scalar_vd(sl,sr,state_south_bar,icoords,SOUTH);
    }

    // NORTH
    if(bNoBoundary[3]>=2) //homogeneous Neumann BC for scalar
    {
        getScalar_MAC_middleStep_bar_decoupled_vd(icoords,COORD_Y,NORTH,state_north_bar,state_north_hat,bGhostCell);
    }
    else
    {
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1] + 1;
        ICoords[2] = icoords[2];
        getScalar_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Y,SOUTH,sr,state_north_hat_r,bGhostCell);
        getScalar_MAC_middleStep_bar_decoupled_vd(icoords,COORD_Y,NORTH,sl,state_north_hat_l,bGhostCell);
        getRiemannSolution_MAC_Scalar_vd(sl,sr,state_north_bar,icoords,NORTH);
    }

    // LOWER
    if(bNoBoundary[4]>=2) //homogeneous Neumann BC for scalar
    {
        getScalar_MAC_middleStep_bar_decoupled_vd(icoords,COORD_Z,LOWER,state_lower_bar,state_lower_hat,bGhostCell);
    }
    else
    {
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1];
        ICoords[2] = icoords[2] - 1;
        getScalar_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Z,UPPER,sl,state_lower_hat_l,bGhostCell);
        getScalar_MAC_middleStep_bar_decoupled_vd(icoords,COORD_Z,LOWER,sr,state_lower_hat_r,bGhostCell);
        getRiemannSolution_MAC_Scalar_vd(sl,sr,state_lower_bar,icoords,LOWER);
    }

    // UPPER
    if(!bNoBoundary[5]>=2) //homogeneous Neumann BC for scalar
    {
        getScalar_MAC_middleStep_bar_decoupled_vd(icoords,COORD_Z,UPPER,state_upper_bar,state_upper_hat,bGhostCell);
    }
    else
    {
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1];
        ICoords[2] = icoords[2] + 1;
        getScalar_MAC_middleStep_bar_decoupled_vd(ICoords,COORD_Z,LOWER,sr,state_upper_hat_r,bGhostCell);
        getScalar_MAC_middleStep_bar_decoupled_vd(icoords,COORD_Z,UPPER,sl,state_upper_hat_l,bGhostCell);
        getRiemannSolution_MAC_Scalar_vd(sl,sr,state_upper_bar,icoords,UPPER);
    }

    //store states bar
    cell_center[index].m_state.m_rho_bar[0] = state_west_bar.m_rho;
    cell_center[index].m_state.m_rho_bar[1] = state_east_bar.m_rho;
    cell_center[index].m_state.m_rho_bar[2] = state_south_bar.m_rho;
    cell_center[index].m_state.m_rho_bar[3] = state_north_bar.m_rho;
    cell_center[index].m_state.m_rho_bar[4] = state_lower_bar.m_rho;
    cell_center[index].m_state.m_rho_bar[5] = state_upper_bar.m_rho;

    cell_center[index].m_state.m_c_bar[0] = state_west_bar.m_c;
    cell_center[index].m_state.m_c_bar[1] = state_east_bar.m_c;
    cell_center[index].m_state.m_c_bar[2] = state_south_bar.m_c;
    cell_center[index].m_state.m_c_bar[3] = state_north_bar.m_c;
    cell_center[index].m_state.m_c_bar[4] = state_lower_bar.m_c;
    cell_center[index].m_state.m_c_bar[5] = state_upper_bar.m_c;
} /* end getScalarBar_MAC_vd */


void Incompress_Solver_Smooth_3D_Cartesian::computeMacPhi_vd(int flag)
{
        int index;
        int i,j,k,icoords[MAXD];
        double P_max,P_min;
        double diffusion[5];

        /* Compute velocity divergence */
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);

            //get the divergence constraint: S^n or S^{n+1/2}
            getDivU_coupled_vd(icoords,diffusion,flag);
            source[index] = (cell_center[index].m_state.div_U - diffusion[3])/(m_dt/2.0);

            if(!flag)
                diff_coeff[index] = cell_center[index].m_state.m_rho; //rho_n
            else
                diff_coeff[index] = cell_center[index].m_state.m_rho_old; //rho_n
        }
        FT_ParallelExchGridArrayBuffer(source,front);
        FT_ParallelExchGridArrayBuffer(diff_coeff,front);

        poisson_solver3d_MacPhi_vd(front,ilower,iupper,ijk_to_I,source,diff_coeff,
                        array,&P_max,&P_min);

        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_phi = array[index];
        }
} /* end computeMacPhi_vd */


void Incompress_Solver_Smooth_3D_Cartesian::computeMacPhi_MAC_vd(int flag)
{
    int index;
    int i,j,k,icoords[MAXD];
    double P_max,P_min;
    double diffusion;

    if (m_dt < min_dt) {
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_phi = 0;
        }
        return;
    }

    //get source terms for MAC-projection
    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
        icoords[0] = i;
        icoords[1] = j;
        icoords[2] = k;
        index = d_index3d(i,j,k,top_gmax);

        //get the divergence constraint: S^{n} for flag=0 or S^{n+1/2} for flag=1
        //div_U = div(U_face_bar) for flag=0 or div_U = div(U_face_bar_0) for flag=1
        getDivU_MAC_vd(icoords,&diffusion,flag,bGhostCell);
        source[index] = (cell_center[index].m_state.div_U - diffusion)/(m_dt/2.0);

        if (!flag)
            diff_coeff[index] = cell_center[index].m_state.m_rho;     //rho_n
        else
            diff_coeff[index] = cell_center[index].m_state.m_rho_old; //rho_n
    }
    FT_ParallelExchGridArrayBuffer(source,front);
    FT_ParallelExchGridArrayBuffer(diff_coeff,front);

    if (debugging("trace"))
        printf("enter poisson_solver3d_MacPhi_vd in computeMacPhi_MAC_vd()\n");
    /*
     * removal tag: HAOZ
     * REFLECTION B.C. needs HOMOGENOUS NEUMANN B.C. for Poisson solver
     * TODO && FIXME
     * */
    poisson_solver3d_MacPhi_vd(front,ilower,iupper,ijk_to_I,source,diff_coeff,
                    array,&P_max,&P_min);
    if (debugging("trace"))
        printf("leave poisson_solver3d_MacPhi_vd in computeMacPhi_MAC_vd()\n");

    for (k = 0; k <= top_gmax[2]; k++)
    for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
    {
        index = d_index3d(i,j,k,top_gmax);
        cell_center[index].m_state.m_phi = array[index];
    }
} /* end computeMacPhi_MAC_vd */


void Incompress_Solver_Smooth_3D_Cartesian::computeNewU_MAC_vd(double t)
{
} /* end computeNewU_MAC_vd */


void Incompress_Solver_Smooth_3D_Cartesian::computeNewUFaceBar_MAC_vd(double t, int flag)
{
        int i, j, k, l, nb, index;
        double tmp, phi, rho, rho_nb[6];
        COMPONENT comp;
        int icoords[MAXD];
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        int ind[6],index_nb[6];

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            //6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);

            comp = cell_center[index].comp;
            phi = cell_center[index].m_state.m_phi;

            if (!ifluid_comp(comp))
            {
                for (l=0; l<3; l++)
                    cell_center[index].m_state.m_U_face_bar[l] = 0.0;
                continue;
            }

            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;

            // 4 directions
            ind[0] = 1;
            ind[1] = 1;
            ind[2] = 1;
            ind[3] = 1;
            // LOWER
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
                    comp,&intfc_state,&hs,crx_coords,t) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                ind[4] = 0;
            else
                ind[4] = 1;
            // UPPER
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
                    comp,&intfc_state,&hs,crx_coords,t) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                ind[5] = 0;
            else
                ind[5] = 1;

            if (!flag) //flag == 0
            {
                rho = cell_center[index].m_state.m_rho;
                for (nb = 0; nb < 6; nb++)
                {
                    if (ind[nb] == 0)
                        rho_nb[nb] = rho;
                    else
                        rho_nb[nb] = cell_center[index_nb[nb]].m_state.m_rho;
                }
            }
            else //flag == 1
            {
                rho = cell_center[index].m_state.m_rho_old;
                for (nb = 0; nb < 6; nb++)
                {
                    if (ind[nb] == 0)
                        rho_nb[nb] = rho;
                    else
                        rho_nb[nb] = cell_center[index_nb[nb]].m_state.m_rho_old;
                }
            }

            if (ind[5] == 0) //cells on UPPER bdry
            {
                //u-faces
                tmp = m_dt/4.0*(1.0/rho_nb[1] + 1.0/rho);
                cell_center[index].m_state.m_U_face_bar[0] -= (cell_center[index_nb[1]].m_state.m_phi - phi)/top_h[0]*tmp;

                //v-faces
                tmp = m_dt/4.0*(1.0/rho_nb[3] + 1.0/rho);
                cell_center[index].m_state.m_U_face_bar[1] -= (cell_center[index_nb[3]].m_state.m_phi - phi)/top_h[1]*tmp;

                //w-faces
               cell_center[index].m_state.m_U_face_bar[2] = 0.0;
            }
            else //other cells
            {
                //u-faces
                tmp = m_dt/4.0*(1.0/rho_nb[1] + 1.0/rho);
                cell_center[index].m_state.m_U_face_bar[0] -= (cell_center[index_nb[1]].m_state.m_phi - phi)/top_h[0]*tmp;

                //v-faces
                tmp = m_dt/4.0*(1.0/rho_nb[3] + 1.0/rho);
                cell_center[index].m_state.m_U_face_bar[1] -= (cell_center[index_nb[3]].m_state.m_phi - phi)/top_h[1]*tmp;

                //w-faces
                tmp = m_dt/4.0*(1.0/rho_nb[5] + 1.0/rho);
                cell_center[index].m_state.m_U_face_bar[2] -= (cell_center[index_nb[5]].m_state.m_phi - phi)/top_h[2]*tmp;
            }
        }

        // scatter states
        for (l = 0; l < 3; l++)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U_face_bar[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U_face_bar[l] = array[index];
            }
        }
} /* end computeNewUFaceBar_MAC_vd */


void Incompress_Solver_Smooth_3D_Cartesian::computeNewUbar_vd(double t, int flag)
{
} /* end computeNewUbar_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getFaceVelocity_middleStep_hat_vd(
        int *icoords,
        GRID_DIRECTION dir,
        L_STATE &state_hat)
{
    L_STATE  state_orig;
    int index;
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    state_orig = cell_center[index].m_state;
    double sL;

    double dx = 0, dy = 0, dz = 0;
    double slope_x_limited[5] = {0,0,0,0,0}, slope_y_limited[5] = {0,0,0,0,0}, slope_z_limited[5] = {0,0,0,0,0};

    switch(dir)
    {
    case WEST:
        getLimitedSlope_vd(icoords,COORD_X,slope_x_limited);
        dx = top_h[0];
        sL = 1.0;
        state_orig.m_U[0] = std::min(state_orig.m_U[0], 0.0);
        state_hat.m_U[0] = state_orig.m_U[0] + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[0];
        state_hat.m_U[1] = state_orig.m_U[1] + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[1];
        state_hat.m_U[2] = state_orig.m_U[2] + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[2];
        break;


    case EAST:
        getLimitedSlope_vd(icoords,COORD_X,slope_x_limited);
        dx = top_h[0];
        sL = 1.0;
        state_orig.m_U[0] = std::max(state_orig.m_U[0], 0.0);
        state_hat.m_U[0] = state_orig.m_U[0] + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[0];
        state_hat.m_U[1] = state_orig.m_U[1] + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[1];
        state_hat.m_U[2] = state_orig.m_U[2] + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[2];
        break;


    case SOUTH:
        getLimitedSlope_vd(icoords,COORD_Y,slope_y_limited);
        dy = top_h[1];
        sL = 1.0;
        state_orig.m_U[1] = std::min(state_orig.m_U[1], 0.0);
        state_hat.m_U[0] = state_orig.m_U[0] + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[0];
        state_hat.m_U[1] = state_orig.m_U[1] + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[1];
        state_hat.m_U[2] = state_orig.m_U[2] + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[2];
        break;

    case NORTH:
        getLimitedSlope_vd(icoords,COORD_Y,slope_y_limited);
        dy = top_h[1];
        sL = 1.0;
        state_orig.m_U[1] = std::max(state_orig.m_U[1], 0.0);
        state_hat.m_U[0] = state_orig.m_U[0] + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[0];
        state_hat.m_U[1] = state_orig.m_U[1] + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[1];
        state_hat.m_U[2] = state_orig.m_U[2] + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[2];
        break;

    case LOWER:
        getLimitedSlope_vd(icoords,COORD_Z,slope_z_limited);
        dz = top_h[2];
        sL = 1.0;
        state_orig.m_U[2] = std::min(state_orig.m_U[2], 0.0);
        state_hat.m_U[0] = state_orig.m_U[0] + sL*(-dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[0];
        state_hat.m_U[1] = state_orig.m_U[1] + sL*(-dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[1];
        state_hat.m_U[2] = state_orig.m_U[2] + sL*(-dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[2];
        break;

    case UPPER:
        getLimitedSlope_vd(icoords,COORD_Z,slope_z_limited);
        dz = top_h[2];
        sL = 1.0;
        state_orig.m_U[2] = std::max(state_orig.m_U[2], 0.0);
        state_hat.m_U[0] = state_orig.m_U[0] + sL*(dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[0];
        state_hat.m_U[1] = state_orig.m_U[1] + sL*(dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[1];
        state_hat.m_U[2] = state_orig.m_U[2] + sL*(dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[2];
        break;
    default:
        assert(false);
    }
} /* end getFaceVelocity_middleStep_hat_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getCenterVelocity_MAC_middleStep_hat_vd(
        int *icoords,
        GRID_DIRECTION dir,
        L_STATE &state_hat)
{
    L_STATE  state_orig;
    int index;
    double sL = 1.0;
    double dx, dy, dz;
    double slope_limited[3] = {0.0,0.0,0.0};

    state_hat.m_U[0] = 0.0;
    state_hat.m_U[1] = 0.0;
    state_hat.m_U[2] = 0.0;

    dx = top_h[0];
    dy = top_h[1];
    dz = top_h[2];

    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    state_orig = cell_center[index].m_state;

    switch(dir)
    {
    case WEST:
        getLimitedSlope_Velocity_MAC_vd(icoords,COORD_X,slope_limited); //u_x, u_y, and u_z
        state_orig.m_U[0] = std::min(state_orig.m_U[0], 0.0);
        state_hat.m_U[0] = state_orig.m_U[0] + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_limited[0];
        break;

    case EAST:
        getLimitedSlope_Velocity_MAC_vd(icoords,COORD_X,slope_limited); //u_x, u_y, and u_z
        state_orig.m_U[0] = std::max(state_orig.m_U[0], 0.0);
        state_hat.m_U[0] = state_orig.m_U[0] + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_limited[0];
        break;

    case SOUTH:
        getLimitedSlope_Velocity_MAC_vd(icoords,COORD_Y,slope_limited); //v_x, v_y, and v_z
        state_orig.m_U[1] = std::min(state_orig.m_U[1], 0.0);
        state_hat.m_U[1] = state_orig.m_U[1] + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_limited[1];
        break;

    case NORTH:
        getLimitedSlope_Velocity_MAC_vd(icoords,COORD_Y,slope_limited); //v_x, v_y, and v_z
        state_orig.m_U[1] = std::max(state_orig.m_U[1], 0.0);
        state_hat.m_U[1] = state_orig.m_U[1] + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_limited[1];
        break;

    case LOWER:
        getLimitedSlope_Velocity_MAC_vd(icoords,COORD_Z,slope_limited); //w_x, w_y, and w_z
        state_orig.m_U[2] = std::min(state_orig.m_U[2], 0.0);
        state_hat.m_U[2] = state_orig.m_U[2] + sL*(-dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_limited[2];
        break;

    case UPPER:
        getLimitedSlope_Velocity_MAC_vd(icoords,COORD_Z,slope_limited); //w_x, w_y, and w_z
        state_orig.m_U[2] = std::max(state_orig.m_U[2], 0.0);
        state_hat.m_U[2] = state_orig.m_U[2] + sL*(dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_limited[2];
        break;
    default:
        assert(false);
    }
} /* end getCenterVelocity_MAC_middleStep_hat_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getFaceVelocity_MAC_middleStep_hat_vd(
        int *icoords,
        EBM_COORD xyz,
        L_STATE &state_hat)
{
    L_STATE  state_orig;
    int index;
    double slope_limited[3] = {0.0, 0.0, 0.0};

    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    state_orig = cell_center[index].m_state;
    state_hat.m_U[0] = 0.0;
    state_hat.m_U[1] = 0.0;
    state_hat.m_U[2] = 0.0;

    getLimitedSlope_Velocity_MAC_vd(icoords,xyz,slope_limited);
    //state_orig.m_U[xyz] = std::min(state_orig.m_U[xyz], 0.0);
    state_hat.m_U[xyz] = state_orig.m_U[xyz] - m_dt/2.0*state_orig.m_U[xyz]*slope_limited[xyz];
} /* end getFaceVelocity_MAC_middleStep_hat_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getEdgeVelocity_MAC_middleStep_hat_vd(
        int *icoords,
        EBM_COORD xyz,
        GRID_DIRECTION dir,
        L_STATE &state_hat)
{
    L_STATE  state_orig;
    int index, l;
    double sL;
    double slope_limited[3] = {0.0, 0.0, 0.0};

    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    state_orig = cell_center[index].m_state;
    state_hat.m_U[0] = 0.0;
    state_hat.m_U[1] = 0.0;
    state_hat.m_U[2] = 0.0;

    if ((dir == EAST) || (dir == NORTH) || (dir == UPPER))
        sL = 1.0;
    else
        sL = -1.0;

    if ((dir == WEST) || (dir == EAST))
        l = 0;
    else if ((dir == SOUTH) || (dir == NORTH))
        l = 1;
    else //dir == LOWER or dir == UPPER
        l = 2;

    // removal tag: HAOZ
    // slope limiter function need further REFLECTION BOUNDARY CONDITION treatment
    getLimitedSlope_Velocity_MAC_vd(icoords,xyz,slope_limited);
    //state_orig.m_U[xyz] = std::min(state_orig.m_U[xyz], 0.0);
    state_hat.m_U[xyz] = state_orig.m_U[xyz] + 0.5*(sL*top_h[l]*slope_limited[l] -
                         m_dt*state_orig.m_U[xyz]*slope_limited[xyz]);
} /* end getEdgeVelocity_MAC_middleStep_hat_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getCenterVelocity_MAC_zBdry_middleStep_hat_vd(
        int *icoords,
        GRID_DIRECTION dir,
        L_STATE &state_hat)
{
    L_STATE  state_orig, state_nb;
    int sL, index_nb;
    double dz = top_h[2];
    double slope_limited;

    if (dir == UPPER)
        sL = 1;
    else if (dir == LOWER)
        sL = -1;
    else
        assert(false);

    index_nb = d_index3d(icoords[0],icoords[1],icoords[2]+sL,top_gmax);
    state_nb = cell_center[index_nb].m_state;
    state_orig.m_U[2] = 0.0;
    slope_limited = sL*(state_nb.m_U[2] - state_orig.m_U[2])/dz; //w_z on LOWER or UPPER bdry (one-sided)
    state_hat.m_U[2] = state_orig.m_U[2] + (sL*dz/2.0 - m_dt/2.0*state_orig.m_U[2])*slope_limited;
} /* end getCenterVelocity_MAC_zBdry_middleStep_hat_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getFaceScalar_middleStep_hat_vd(
        int *icoords,
        GRID_DIRECTION dir,
        L_STATE &state_hat)
{
    L_STATE  state_orig;
    int index;
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    state_orig = cell_center[index].m_state;
    double sL;

    double dx = 0, dy = 0, dz = 0;
    double slope_x_limited[5] = {0,0,0,0,0}, slope_y_limited[5] = {0,0,0,0,0}, slope_z_limited[5] = {0,0,0,0,0};

    switch(dir)
    {
    case WEST:
        getLimitedSlope_vd(icoords,COORD_X,slope_x_limited);
        dx = top_h[0];
        sL = 1.0;
        state_orig.m_U[0] = std::min(state_orig.m_U[0], 0.0);
        state_hat.m_rho  = state_orig.m_rho  + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[3];
        state_hat.m_c    = state_orig.m_c    + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[4];
        break;


    case EAST:
        getLimitedSlope_vd(icoords,COORD_X,slope_x_limited);
        dx = top_h[0];
        sL = 1.0;
        state_orig.m_U[0] = std::max(state_orig.m_U[0], 0.0);
        state_hat.m_rho  = state_orig.m_rho  + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[3];
        state_hat.m_c    = state_orig.m_c    + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[4];
        break;


    case SOUTH:
        getLimitedSlope_vd(icoords,COORD_Y,slope_y_limited);
        dy = top_h[1];
        sL = 1.0;
        state_orig.m_U[1] = std::min(state_orig.m_U[1], 0.0);
        state_hat.m_rho  = state_orig.m_rho  + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[3];
        state_hat.m_c    = state_orig.m_c    + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[4];
        break;

    case NORTH:
        getLimitedSlope_vd(icoords,COORD_Y,slope_y_limited);
        dy = top_h[1];
        sL = 1.0;
        state_orig.m_U[1] = std::max(state_orig.m_U[1], 0.0);
        state_hat.m_rho  = state_orig.m_rho  + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[3];
        state_hat.m_c    = state_orig.m_c    + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[4];
        break;

    case LOWER:
        getLimitedSlope_vd(icoords,COORD_Z,slope_z_limited);
        dz = top_h[2];
        sL = 1.0;
        state_orig.m_U[2] = std::min(state_orig.m_U[2], 0.0);
        state_hat.m_rho  = state_orig.m_rho  + sL*(-dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[3];
        state_hat.m_c    = state_orig.m_c    + sL*(-dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[4];
        break;

    case UPPER:
        getLimitedSlope_vd(icoords,COORD_Z,slope_z_limited);
        dz = top_h[2];
        sL = 1.0;
        state_orig.m_U[2] = std::max(state_orig.m_U[2], 0.0);
        state_hat.m_rho  = state_orig.m_rho  + sL*(dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[3];
        state_hat.m_c    = state_orig.m_c    + sL*(dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[4];
        break;
    default:
        assert(false);
    }
} /* end getFaceScalar_middleStep_hat_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getFaceScalar_MAC_middleStep_hat_vd(
        int *icoords,
        EBM_COORD xyz,
        GRID_DIRECTION dir,
        L_STATE &state_hat,
        boolean bGhostCell)
{
    L_STATE  state_orig;
    int index, index_nb[6];
    double sL, U_center;
    double slope_limited[2] = {0, 0};
    int bNoBoundary[6],nb;
    COMPONENT comp;
    double crx_coords[MAXD];
    GRID_DIRECTION dir_n[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    comp = top_comp[index];

    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

    for (nb = 0; nb < 6; nb++)
    {
        checkBoundaryCondition(dir_n[nb],icoords,&bNoBoundary[nb],m_t_old,comp);
    }

    /*
    // 4 directions
    bNoBoundary[0] = YES;
    bNoBoundary[1] = YES;
    bNoBoundary[2] = YES;
    bNoBoundary[3] = YES;
    // LOWER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[4] = NO;
    else
        bNoBoundary[4] = YES;
    // UPPER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[5] = NO;
    else
        bNoBoundary[5] = YES;
    */

    state_orig = cell_center[index].m_state;
    state_hat.m_rho = 0;
    state_hat.m_c = 0;

    if ((dir == EAST) || (dir == NORTH) || (dir == UPPER))
        sL = 1.0;
    else
        sL = -1.0;

    if (xyz == COORD_X) //average u at cell centers
        U_center = 0.5*(state_orig.m_U[0] + cell_center[index_nb[0]].m_state.m_U[0]);
    else if (xyz == COORD_Y) //average v at cell centers
        U_center = 0.5*(state_orig.m_U[1] + cell_center[index_nb[2]].m_state.m_U[1]);
    else //average w at cell centers
    {
        if (bNoBoundary[4]==2) //cells on LOWER boundary
            U_center = 0.5*(state_orig.m_U[2] + 0);
        else if (bNoBoundary[5]==2) //cells on UPPER boundary
            U_center = 0.5*(0 + cell_center[index_nb[4]].m_state.m_U[2]);
        else //other cells
            U_center = 0.5*(state_orig.m_U[2] + cell_center[index_nb[4]].m_state.m_U[2]);
    }

    // removal tag: HAOZ scalar version of slope limiter need REFLECTION BOUNDARY CONDITION treatment
    getLimitedSlope_Scalar_MAC_vd(icoords,xyz,slope_limited,bGhostCell);
    state_hat.m_rho = state_orig.m_rho + 0.5*(sL*top_h[xyz] - m_dt*U_center)*slope_limited[0];
    state_hat.m_c = state_orig.m_c + 0.5*(sL*top_h[xyz] - m_dt*U_center)*slope_limited[1];
} /* end getFaceScalar_MAC_middleStep_hat_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getFaceVelocity_middleStep_bar_vd(
        int *icoords,
        GRID_DIRECTION dir,
        L_STATE &state_bar,
        double transverseD[5],
        L_STATE state_hat)
{
    int index;
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double rho = cell_center[index].m_state.m_rho;

    double diffusion[5];

    getViscousTerm_coupled_vd(icoords,diffusion); //get viscous terms of NS eqn's

    state_bar.m_U[0] = state_hat.m_U[0] + m_dt/2.0 * (-transverseD[0] + diffusion[0]/rho);
    state_bar.m_U[1] = state_hat.m_U[1] + m_dt/2.0 * (-transverseD[1] + diffusion[1]/rho);
    state_bar.m_U[2] = state_hat.m_U[2] + m_dt/2.0 * (-transverseD[2] + diffusion[2]/rho);

    double coords[3];
    L_STATE source_term;

    //add gravity term
    getRectangleCenter(index, coords);
    computeSourceTerm_Adv(coords, source_term);
    state_bar.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_bar.m_U[1] += m_dt/2 * source_term.m_U[1];
    state_bar.m_U[2] += m_dt/2 * source_term.m_U[2];

    //add pressure gradient
    state_bar.m_U[0] -= m_dt/2 * cell_center[index].m_state.grad_q[0]/rho;
    state_bar.m_U[1] -= m_dt/2 * cell_center[index].m_state.grad_q[1]/rho;
    state_bar.m_U[2] -= m_dt/2 * cell_center[index].m_state.grad_q[2]/rho;
} /* end getFaceVelocity_middleStep_bar_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getFaceVelocity_middleStep_bar_decoupled_vd(
        int *icoords,
        GRID_DIRECTION dir,
        L_STATE &state_bar,
        double transverseD[5],
        L_STATE state_hat)
{
    int index;
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double rho = cell_center[index].m_state.m_rho;

    double diffusion[5];

    getViscousTerm_decoupled_vd(icoords,diffusion); //get viscous terms of NS eqn's

    state_bar.m_U[0] = state_hat.m_U[0] + m_dt/2.0 * (-transverseD[0] + diffusion[0]/rho);
    state_bar.m_U[1] = state_hat.m_U[1] + m_dt/2.0 * (-transverseD[1] + diffusion[1]/rho);
    state_bar.m_U[2] = state_hat.m_U[2] + m_dt/2.0 * (-transverseD[2] + diffusion[2]/rho);

    double coords[3];
    L_STATE source_term;

    //add gravity term
    getRectangleCenter(index, coords);
    computeSourceTerm_Adv(coords, source_term);
    state_bar.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_bar.m_U[1] += m_dt/2 * source_term.m_U[1];
    state_bar.m_U[2] += m_dt/2 * source_term.m_U[2];

    //add pressure gradient
    state_bar.m_U[0] -= m_dt/2 * cell_center[index].m_state.grad_q[0]/rho;
    state_bar.m_U[1] -= m_dt/2 * cell_center[index].m_state.grad_q[1]/rho;
    state_bar.m_U[2] -= m_dt/2 * cell_center[index].m_state.grad_q[2]/rho;
} /* end getFaceVelocity_middleStep_bar_decoupled_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getVelocity_MAC_middleStep_bar_decoupled_vd(
        int *icoords,
        EBM_COORD xyz,
        GRID_DIRECTION dir,
        L_STATE &state_bar,
        L_STATE state_hat)
{
    int index, index_nb;
    bool bNoBoundary;
    double rho;
    double coords[3], diffusion[3], transverseD[3];
    L_STATE source_term;
    COMPONENT comp;
    double crx_coords[MAXD];
    POINTER intfc_state;
    HYPER_SURF *hs;

    state_bar.m_U[xyz] = 0.0;

    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    comp = top_comp[index];

    // UPPER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary = NO;
    else
        bNoBoundary = YES;

    //get transverse derivatives terms of NS eqn's
    getTransverseDTerm_Velocity_MAC_vd(icoords, transverseD);
    //get viscous terms of NS eqn's
    getViscousTerm_MAC_decoupled_vd(icoords, xyz, diffusion);
    //add gravity term (constant)
    getRectangleCenter(index, coords);
    computeSourceTerm_Adv(coords, source_term);

    switch(xyz)
    {
    case COORD_X:
        index_nb = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
        break;
    case COORD_Y:
        index_nb = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
        break;
    case COORD_Z:
        index_nb = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);
        break;
    default:
        assert(false);
    }

    if (!bNoBoundary && xyz==2) //for w on UPPER boundary cells
        rho = cell_center[index].m_state.m_rho;
    else
        rho = 2.0/(1.0/cell_center[index].m_state.m_rho +
                   1.0/cell_center[index_nb].m_state.m_rho);

    state_bar.m_U[xyz] += state_hat.m_U[xyz] + m_dt/2.0*(-transverseD[xyz] +
                          diffusion[xyz]/rho);
    state_bar.m_U[xyz] += m_dt/2.0*source_term.m_U[xyz];
    state_bar.m_U[xyz] -= m_dt/2.0*cell_center[index].m_state.grad_q[xyz]/rho;
} /* end getVelocity_MAC_middleStep_bar_decoupled_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getScalar_MAC_middleStep_bar_decoupled_vd(
        int *icoords,
        EBM_COORD xyz,
        GRID_DIRECTION dir,
        L_STATE &state_bar,
        L_STATE state_hat,
        boolean bGhostCell)
{
    int index;
    double rho, diff, diffC, transverseD[2];

    state_bar.m_rho = state_bar.m_c = 0;

    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    rho = cell_center[index].m_state.m_rho;

    //get transverse derivatives terms of continuity & concentration eqn's
    getTransverseDTerm_Scalar_MAC_vd(icoords,xyz,transverseD,bGhostCell); //stability???
    //get the divergence constraint S^n
    getDivU_MAC_vd(icoords,&diff,0,NO);
    //get the diffusion term of concentration eqn
    getDiffusionC_MAC_vd(icoords,&diffC,NO);

    state_bar.m_rho += state_hat.m_rho + m_dt/2.0*(-transverseD[0]-diff*rho);
    state_bar.m_c += state_hat.m_c + m_dt/2.0*(-transverseD[1]+diffC/rho);
} /* end getScalar_MAC_middleStep_bar_decoupled_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getFaceScalar_middleStep_bar_vd(
        int *icoords,
        GRID_DIRECTION dir,
        L_STATE &state_bar,
        double transverseD[5],
        L_STATE state_hat)
{
    int index;
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double rho = cell_center[index].m_state.m_rho;

    double diffusion[5];

    getDivU_coupled_vd(icoords,diffusion,0);      //get the divergence constraint S^n
    getDiffusionC_coupled_vd(icoords,diffusion);  //get the diffusion term of concentration eqn

    state_bar.m_rho  = state_hat.m_rho  + m_dt/2.0 * (-transverseD[3] - diffusion[3]*rho);
    state_bar.m_c    = state_hat.m_c    + m_dt/2.0 * (-transverseD[4] + diffusion[4]/rho);
} /* end getFaceScalar_middleStep_bar_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getViscousTerm_coupled_tmp_vd(
        int *icoords,
        double diffusion[5])
{
} /* end getViscousTerm_coupled_tmp_vd */


/*               _                                                                                 _
                |  (4/3*mu*u_x)_x+(mu*u_y)_y+(mu*u_z)_z+(mu*v_x)_y+(mu*w_x)_z-2/3*[mu*(v_y+w_z)]_x  |
viscous terms = |  (mu*v_x)_x+(4/3*mu*v_y)_y+(mu*v_z)_z+(mu*u_y)_x+(mu*w_y)_z-2/3*[mu*(u_x+w_z)]_y  |
                |_ (mu*w_x)_x+(mu*w_y)_y+(4/3*mu*w_z)_z+(mu*u_z)_x+(mu*v_z)_y-2/3*[mu*(u_x+v_y)]_z _|
*/
//valid for variable mu case
void Incompress_Solver_Smooth_3D_Cartesian::getViscousTerm_coupled_vd(
        int *icoords,
        double diffusion[5])
{
    int index,index_nb[18];
    double mu_edge[6],mu0;
    L_STATE Unb;
    double U0_nb[18],U1_nb[18],U2_nb[18],U0_center,U1_center,U2_center;
    // U0[6] -- U0[17] U1[6] -- U1[17]  U2[6] -- U2[17] are corner values
    int nb;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    bool bNoBoundary[6];
    double dh[3],dh0[3],dh1[3];
    double coords[MAXD];

    diffusion[0] = 0.0;
    diffusion[1] = 0.0;
    diffusion[2] = 0.0;

    int i = icoords[0];
    int j = icoords[1];
    int k = icoords[2];

    index = d_index3d(i,j,k,top_gmax);

    //6 neighbours of the center cell
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    //xy cut neighbours
    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);

    //yz cut neighbours
    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);

    //xz cut neighbours
    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);

    mu0 = cell_center[index].m_state.m_mu;
    U0_center = cell_center[index].m_state.m_U[0];
    U1_center = cell_center[index].m_state.m_U[1];
    U2_center = cell_center[index].m_state.m_U[2];

    for (nb = 0; nb < 6; nb++)
    {
        bNoBoundary[nb] = getNeighborOrBoundaryState_vd(icoords,dir[nb],Unb,m_t_old);
        U0_nb[nb] = Unb.m_U[0];
        U1_nb[nb] = Unb.m_U[1];
        U2_nb[nb] = Unb.m_U[2];
        if(!bNoBoundary[nb])
            mu_edge[nb] = mu0;
        else
            mu_edge[nb] = 0.5*(mu0 + cell_center[index_nb[nb]].m_state.m_mu);
    }

    // non-cross derivative terms
    dh[0] = top_h[0];
    dh[1] = top_h[1];
    dh[2] = top_h[2];

    if (bNoBoundary[0])
        dh0[0] = top_h[0];
    else
        dh0[0] = top_h[0]/2.0;
    if (bNoBoundary[1])
        dh1[0] = top_h[0];
    else
        dh1[0] = top_h[0]/2.0;

    if (bNoBoundary[2])
        dh0[1] = top_h[1];
    else
        dh0[1] = top_h[1]/2.0;
    if (bNoBoundary[3])
        dh1[1] = top_h[1];
    else
        dh1[1] = top_h[1]/2.0;

    if (bNoBoundary[4])
        dh0[2] = top_h[2];
    else
        dh0[2] = top_h[2]/2.0;
    if (bNoBoundary[5])
        dh1[2] = top_h[2];
    else
        dh1[2] = top_h[2]/2.0;

    diffusion[0] += 4.0/3.0*(mu_edge[1]*(U0_nb[1]-U0_center)/dh1[0] - mu_edge[0]*(U0_center-U0_nb[0])/dh0[0])/dh[0];//(4/3*mu*u_x)_x
    diffusion[1] +=         (mu_edge[1]*(U1_nb[1]-U1_center)/dh1[0] - mu_edge[0]*(U1_center-U1_nb[0])/dh0[0])/dh[0];//(mu*v_x)_x
    diffusion[2] +=         (mu_edge[1]*(U2_nb[1]-U2_center)/dh1[0] - mu_edge[0]*(U2_center-U2_nb[0])/dh0[0])/dh[0];//(mu*w_x)_x

    diffusion[0] +=         (mu_edge[3]*(U0_nb[3]-U0_center)/dh1[1] - mu_edge[2]*(U0_center-U0_nb[2])/dh0[1])/dh[1];//(mu*u_y)_y
    diffusion[1] += 4.0/3.0*(mu_edge[3]*(U1_nb[3]-U1_center)/dh1[1] - mu_edge[2]*(U1_center-U1_nb[2])/dh0[1])/dh[1];//(4/3*mu*v_y)_y
    diffusion[2] +=         (mu_edge[3]*(U2_nb[3]-U2_center)/dh1[1] - mu_edge[2]*(U2_center-U2_nb[2])/dh0[1])/dh[1];//(mu*w_y)_y

    diffusion[0] +=         (mu_edge[5]*(U0_nb[5]-U0_center)/dh1[2] - mu_edge[4]*(U0_center-U0_nb[4])/dh0[2])/dh[2];//(mu*u_z)_z
    diffusion[1] +=         (mu_edge[5]*(U1_nb[5]-U1_center)/dh1[2] - mu_edge[4]*(U1_center-U1_nb[4])/dh0[2])/dh[2];//(mu*v_z)_z
    diffusion[2] += 4.0/3.0*(mu_edge[5]*(U2_nb[5]-U2_center)/dh1[2] - mu_edge[4]*(U2_center-U2_nb[4])/dh0[2])/dh[2];//(4/3*mu*w_z)_z


    //cross derivative terms

    //traverse the corners on 3 cut planes

    //corner (i-1/2,j-1/2,k)
    if (bNoBoundary[0] && bNoBoundary[2])
    {
        U0_nb[6] = (U0_nb[0]+U0_nb[2]+cell_center[index_nb[6]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[6] = (U1_nb[0]+U1_nb[2]+cell_center[index_nb[6]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[6] = (U2_nb[0]+U2_nb[2]+cell_center[index_nb[6]].m_state.m_U[2]+U2_center)/4.0;
    }
    else assert(false);

    //corner (i+1/2,j-1/2,k)
    if (bNoBoundary[1] && bNoBoundary[2])
    {
        U0_nb[7] = (U0_nb[1]+U0_nb[2]+cell_center[index_nb[7]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[7] = (U1_nb[1]+U1_nb[2]+cell_center[index_nb[7]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[7] = (U2_nb[1]+U2_nb[2]+cell_center[index_nb[7]].m_state.m_U[2]+U2_center)/4.0;
    }
    else assert(false);

    //corner (i+1/2,j+1/2,k)
    if (bNoBoundary[1] && bNoBoundary[3])
    {
        U0_nb[8] = (U0_nb[1]+U0_nb[3]+cell_center[index_nb[8]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[8] = (U1_nb[1]+U1_nb[3]+cell_center[index_nb[8]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[8] = (U2_nb[1]+U2_nb[3]+cell_center[index_nb[8]].m_state.m_U[2]+U2_center)/4.0;
    }
    else assert(false);

    //corner (i-1/2,j+1/2,k)
    if (bNoBoundary[0] && bNoBoundary[3])
    {
        U0_nb[9] = (U0_nb[0]+U0_nb[3]+cell_center[index_nb[9]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[9] = (U1_nb[0]+U1_nb[3]+cell_center[index_nb[9]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[9] = (U2_nb[0]+U2_nb[3]+cell_center[index_nb[9]].m_state.m_U[2]+U2_center)/4.0;
    }
    else assert(false);

    //corner (i,j-1/2,k-1/2)
    if (bNoBoundary[4])
    {
        U0_nb[10] = (U0_nb[2]+U0_nb[4]+cell_center[index_nb[10]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[10] = (U1_nb[2]+U1_nb[4]+cell_center[index_nb[10]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[10] = (U2_nb[2]+U2_nb[4]+cell_center[index_nb[10]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        U0_nb[10] = 0.0;
        U1_nb[10] = 0.0;
        U2_nb[10] = 0.0;
    }

    //corner (i,j+1/2,k-1/2)
    if (bNoBoundary[4])
    {
        U0_nb[11] = (U0_nb[3]+U0_nb[4]+cell_center[index_nb[11]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[11] = (U1_nb[3]+U1_nb[4]+cell_center[index_nb[11]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[11] = (U2_nb[3]+U2_nb[4]+cell_center[index_nb[11]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        U0_nb[11] = 0.0;
        U1_nb[11] = 0.0;
        U2_nb[11] = 0.0;
    }

    //corner (i,j+1/2,k+1/2)
    if (bNoBoundary[5])
    {
        U0_nb[12] = (U0_nb[3]+U0_nb[5]+cell_center[index_nb[12]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[12] = (U1_nb[3]+U1_nb[5]+cell_center[index_nb[12]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[12] = (U2_nb[3]+U2_nb[5]+cell_center[index_nb[12]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        U0_nb[12] = 0.0;
        U1_nb[12] = 0.0;
        U2_nb[12] = 0.0;
    }


    //corner (i,j-1/2,k+1/2)
    if (bNoBoundary[5])
    {
        U0_nb[13] = (U0_nb[2]+U0_nb[5]+cell_center[index_nb[13]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[13] = (U1_nb[2]+U1_nb[5]+cell_center[index_nb[13]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[13] = (U2_nb[2]+U2_nb[5]+cell_center[index_nb[13]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        U0_nb[13] = 0.0;
        U1_nb[13] = 0.0;
        U2_nb[13] = 0.0;
    }

    //corner (i-1/2,j,k-1/2)
    if (bNoBoundary[4])
    {
        U0_nb[14] = (U0_nb[0]+U0_nb[4]+cell_center[index_nb[14]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[14] = (U1_nb[0]+U1_nb[4]+cell_center[index_nb[14]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[14] = (U2_nb[0]+U2_nb[4]+cell_center[index_nb[14]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        U0_nb[14] = 0.0;
        U1_nb[14] = 0.0;
        U2_nb[14] = 0.0;
    }

    //corner (i+1/2,j,k-1/2)
    if (bNoBoundary[4])
    {
        U0_nb[15] = (U0_nb[1]+U0_nb[4]+cell_center[index_nb[15]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[15] = (U1_nb[1]+U1_nb[4]+cell_center[index_nb[15]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[15] = (U2_nb[1]+U2_nb[4]+cell_center[index_nb[15]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        U0_nb[15] = 0.0;
        U1_nb[15] = 0.0;
        U2_nb[15] = 0.0;
    }

    //corner (i+1/2,j,k+1/2)
    if (bNoBoundary[5])
    {
        U0_nb[16] = (U0_nb[1]+U0_nb[5]+cell_center[index_nb[16]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[16] = (U1_nb[1]+U1_nb[5]+cell_center[index_nb[16]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[16] = (U2_nb[1]+U2_nb[5]+cell_center[index_nb[16]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        U0_nb[16] = 0.0;
        U1_nb[16] = 0.0;
        U2_nb[16] = 0.0;
    }

    //corner (i-1/2,j,k+1/2)
    if (bNoBoundary[5])
    {
        U0_nb[17] = (U0_nb[0]+U0_nb[5]+cell_center[index_nb[17]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[17] = (U1_nb[0]+U1_nb[5]+cell_center[index_nb[17]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[17] = (U2_nb[0]+U2_nb[5]+cell_center[index_nb[17]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        U0_nb[17] = 0.0;
        U1_nb[17] = 0.0;
        U2_nb[17] = 0.0;
    }

    diffusion[0] += (mu_edge[2]*U1_nb[6]-mu_edge[2]*U1_nb[7]+mu_edge[3]*U1_nb[8]-mu_edge[3]*U1_nb[9])/(top_h[0]*top_h[1]);//(mu*v_x)_y
    diffusion[0] += (mu_edge[4]*U2_nb[14]-mu_edge[4]*U2_nb[15]+mu_edge[5]*U2_nb[16]-mu_edge[5]*U2_nb[17])/(top_h[0]*top_h[2]);//(mu*w_x)_z
    diffusion[0] += -(2.0/3.0)*(mu_edge[0]*U1_nb[6]-mu_edge[1]*U1_nb[7]+mu_edge[1]*U1_nb[8]-mu_edge[0]*U1_nb[9])/(top_h[0]*top_h[1]);
                         //(-2/3)*(mu*v_y)_x
    diffusion[0] += -(2.0/3.0)*(mu_edge[0]*U2_nb[14]-mu_edge[1]*U2_nb[15]+mu_edge[1]*U2_nb[16]-mu_edge[0]*U2_nb[17])/(top_h[0]*top_h[2]);
                         //(-2/3)*(mu*w_z)_x

    diffusion[1] += (mu_edge[0]*U0_nb[6]-mu_edge[1]*U0_nb[7]+mu_edge[1]*U0_nb[8]-mu_edge[0]*U0_nb[9])/(top_h[0]*top_h[1]);//(mu*u_y)_x
    diffusion[1] += (mu_edge[4]*U2_nb[10]-mu_edge[4]*U2_nb[11]+mu_edge[5]*U2_nb[12]-mu_edge[5]*U2_nb[13])/(top_h[1]*top_h[2]);//(mu*w_y)_z
    diffusion[1] += -(2.0/3.0)*(mu_edge[2]*U0_nb[6]-mu_edge[2]*U0_nb[7]+mu_edge[3]*U0_nb[8]-mu_edge[3]*U0_nb[9])/(top_h[0]*top_h[1]);
                         //(-2/3)*(mu*u_x)_y
    diffusion[1] += -(2.0/3.0)*(mu_edge[2]*U2_nb[10]-mu_edge[3]*U2_nb[11]+mu_edge[3]*U2_nb[12]-mu_edge[2]*U2_nb[13])/(top_h[1]*top_h[2]);
                         //(-2/3)*(mu*w_z)_y

    diffusion[2] += (mu_edge[0]*U0_nb[14]-mu_edge[1]*U0_nb[15]+mu_edge[1]*U0_nb[16]-mu_edge[0]*U0_nb[17])/(top_h[0]*top_h[2]);//(mu*u_z)_x
    diffusion[2] += (mu_edge[2]*U1_nb[10]-mu_edge[3]*U1_nb[11]+mu_edge[3]*U1_nb[12]-mu_edge[2]*U1_nb[13])/(top_h[1]*top_h[2]);//(mu*v_z)_y
    diffusion[2] += -(2.0/3.0)*(mu_edge[4]*U0_nb[14]-mu_edge[4]*U0_nb[15]+mu_edge[5]*U0_nb[16]-mu_edge[5]*U0_nb[17])/(top_h[0]*top_h[2]);
                         //(-2/3)*(mu*u_x)_z
    diffusion[2] += -(2.0/3.0)*(mu_edge[4]*U1_nb[10]-mu_edge[4]*U1_nb[11]+mu_edge[5]*U1_nb[12]-mu_edge[5]*U1_nb[13])/(top_h[1]*top_h[2]);
                         //(-2/3)*(mu*v_y)_z
} /* end getViscousTerm_coupled_vd */


// just for divergence free case, need to modify???????
void Incompress_Solver_Smooth_3D_Cartesian::getViscousTerm_decoupled_vd(
        int *icoords,
        double diffusion[5])
{
    int index,index_nb[6];
    double mu0;
    L_STATE Unb;
    double U0_nb[6],U1_nb[6],U2_nb[6],U0_center,U1_center,U2_center;
    int nb, ICoords[3];
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    bool bNoBoundary[6];
    double dh[3],dh0[3],dh1[3];
    double coords[MAXD], diff0[5],diff1[5];

    diffusion[0] = 0.0;
    diffusion[1] = 0.0;
    diffusion[2] = 0.0;

    int i = icoords[0];
    int j = icoords[1];
    int k = icoords[2];

    index = d_index3d(i,j,k,top_gmax);

    //6 neighbours of the center cell
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    mu0 = cell_center[index].m_state.m_mu;
    U0_center = cell_center[index].m_state.m_U[0];
    U1_center = cell_center[index].m_state.m_U[1];
    U2_center = cell_center[index].m_state.m_U[2];

    for (nb = 0; nb < 6; nb++)
    {
        bNoBoundary[nb] = getNeighborOrBoundaryState_vd(icoords,dir[nb],Unb,m_t_old);
        U0_nb[nb] = Unb.m_U[0];
        U1_nb[nb] = Unb.m_U[1];
        U2_nb[nb] = Unb.m_U[2];
    }

    // non-cross derivative terms
    dh[0] = top_h[0];
    dh[1] = top_h[1];
    dh[2] = top_h[2];

    if (bNoBoundary[0])
        dh0[0] = top_h[0];
    else
        dh0[0] = top_h[0]/2.0;
    if (bNoBoundary[1])
        dh1[0] = top_h[0];
    else
        dh1[0] = top_h[0]/2.0;

    if (bNoBoundary[2])
        dh0[1] = top_h[1];
    else
        dh0[1] = top_h[1]/2.0;
    if (bNoBoundary[3])
        dh1[1] = top_h[1];
    else
        dh1[1] = top_h[1]/2.0;

    if (bNoBoundary[4])
        dh0[2] = top_h[2];
    else
        dh0[2] = top_h[2]/2.0;
    if (bNoBoundary[5])
        dh1[2] = top_h[2];
    else
        dh1[2] = top_h[2]/2.0;

    diffusion[0] += mu0*((U0_nb[1]-U0_center)/dh1[0] - (U0_center-U0_nb[0])/dh0[0])/dh[0];//mu*u_xx
    diffusion[1] += mu0*((U1_nb[1]-U1_center)/dh1[0] - (U1_center-U1_nb[0])/dh0[0])/dh[0];//mu*v_xx
    diffusion[2] += mu0*((U2_nb[1]-U2_center)/dh1[0] - (U2_center-U2_nb[0])/dh0[0])/dh[0];//mu*w_xx

    diffusion[0] += mu0*((U0_nb[3]-U0_center)/dh1[1] - (U0_center-U0_nb[2])/dh0[1])/dh[1];//mu*u_yy
    diffusion[1] += mu0*((U1_nb[3]-U1_center)/dh1[1] - (U1_center-U1_nb[2])/dh0[1])/dh[1];//mu*v_yy
    diffusion[2] += mu0*((U2_nb[3]-U2_center)/dh1[1] - (U2_center-U2_nb[2])/dh0[1])/dh[1];//mu*w_yy

    diffusion[0] += mu0*((U0_nb[5]-U0_center)/dh1[2] - (U0_center-U0_nb[4])/dh0[2])/dh[2];//mu*u_zz
    diffusion[1] += mu0*((U1_nb[5]-U1_center)/dh1[2] - (U1_center-U1_nb[4])/dh0[2])/dh[2];//mu*v_zz
    diffusion[2] += mu0*((U2_nb[5]-U2_center)/dh1[2] - (U2_center-U2_nb[4])/dh0[2])/dh[2];//mu*w_zz
} /* end getViscousTerm_decoupled_vd */


/*
//valid for constant mu and divergence-free case (no cross terms)
void Incompress_Solver_Smooth_3D_Cartesian::getViscousTerm_MAC_decoupled_vd(
        int *icoords,
        EBM_COORD xyz,
        double diffusion[3])
{
    int nb,index,index_nb[6];
    double mu0;
    double U0_nb[6],U1_nb[6],U2_nb[6],U0_center,U1_center,U2_center;
    bool bNoBoundary[6];
    double dh[3];
    COMPONENT comp;
    double crx_coords[MAXD];
    POINTER intfc_state;
    HYPER_SURF *hs;

    diffusion[0] = 0.0;
    diffusion[1] = 0.0;
    diffusion[2] = 0.0;

    int i = icoords[0];
    int j = icoords[1];
    int k = icoords[2];

    index = d_index3d(i,j,k,top_gmax);
    comp = cell_center[index].comp;

    //6 neighbours of the center cell
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    U0_center = cell_center[index].m_state.m_U[0];
    U1_center = cell_center[index].m_state.m_U[1];
    U2_center = cell_center[index].m_state.m_U[2];

    // 4 directions
    bNoBoundary[0] = YES;
    bNoBoundary[1] = YES;
    bNoBoundary[2] = YES;
    bNoBoundary[3] = YES;

    // LOWER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[4] = NO; //cells on LOWER bdry
    else
        bNoBoundary[4] = YES;

    // UPPER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[5] = NO; //cells on UPPER bdry
    else
        bNoBoundary[5] = YES;

    for (nb = 0; nb < 4; nb++) //4 directions only
    {
        U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
        U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
        U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
    }

    // non-cross derivative terms
    dh[0] = top_h[0];
    dh[1] = top_h[1];
    dh[2] = top_h[2];

    switch(xyz)
    {
    case COORD_X:
        if (bNoBoundary[4] && bNoBoundary[5])
        {
            U0_nb[4] = cell_center[index_nb[4]].m_state.m_U[0];
            U0_nb[5] = cell_center[index_nb[5]].m_state.m_U[0];
        }
        else if (!bNoBoundary[4]) //cells on LOWER bdry
        {
            U0_nb[4] = -U0_center;
            U0_nb[5] = cell_center[index_nb[5]].m_state.m_U[0];
        }
        else if (!bNoBoundary[5]) //cells on UPPER bdry
        {
            U0_nb[4] = cell_center[index_nb[4]].m_state.m_U[0];
            U0_nb[5] = -U0_center;
        }
        else assert(false);
        mu0 = 0.5*(cell_center[index].m_state.m_mu + cell_center[index_nb[1]].m_state.m_mu);
        diffusion[0] += mu0*((U0_nb[1]-U0_center) - (U0_center-U0_nb[0]))/dh[0]/dh[0];//mu*u_xx
        diffusion[0] += mu0*((U0_nb[3]-U0_center) - (U0_center-U0_nb[2]))/dh[1]/dh[1];//mu*u_yy
        diffusion[0] += mu0*((U0_nb[5]-U0_center) - (U0_center-U0_nb[4]))/dh[2]/dh[2];//mu*u_zz
        break;

    case COORD_Y:
        if (bNoBoundary[4] && bNoBoundary[5])
        {
            U1_nb[4] = cell_center[index_nb[4]].m_state.m_U[1];
            U1_nb[5] = cell_center[index_nb[5]].m_state.m_U[1];
        }
        else if (!bNoBoundary[4]) //cells on LOWER bdry
        {
            U1_nb[4] = -U1_center;
            U1_nb[5] = cell_center[index_nb[5]].m_state.m_U[1];
        }
        else if (!bNoBoundary[5]) //cells on UPPER bdry
        {
            U1_nb[4] = cell_center[index_nb[4]].m_state.m_U[1];
            U1_nb[5] = -U1_center;
        }
        else assert(false);
        mu0 = 0.5*(cell_center[index].m_state.m_mu + cell_center[index_nb[3]].m_state.m_mu);
        diffusion[1] += mu0*((U1_nb[1]-U1_center) - (U1_center-U1_nb[0]))/dh[0]/dh[0];//mu*v_xx
        diffusion[1] += mu0*((U1_nb[3]-U1_center) - (U1_center-U1_nb[2]))/dh[1]/dh[1];//mu*v_yy
        diffusion[1] += mu0*((U1_nb[5]-U1_center) - (U1_center-U1_nb[4]))/dh[2]/dh[2];//mu*v_zz
        break;

    case COORD_Z:
        if (bNoBoundary[4] && bNoBoundary[5])
        {
            U2_nb[4] = cell_center[index_nb[4]].m_state.m_U[2];
            U2_nb[5] = cell_center[index_nb[5]].m_state.m_U[2];
            mu0 = 0.5*(cell_center[index].m_state.m_mu + cell_center[index_nb[5]].m_state.m_mu);
            diffusion[2] += mu0*((U2_nb[1]-U2_center) - (U2_center-U2_nb[0]))/dh[0]/dh[0];//mu*w_xx
            diffusion[2] += mu0*((U2_nb[3]-U2_center) - (U2_center-U2_nb[2]))/dh[1]/dh[1];//mu*w_yy
            diffusion[2] += mu0*((U2_nb[5]-U2_center) - (U2_center-U2_nb[4]))/dh[2]/dh[2];//mu*w_zz
        }
        else if (!bNoBoundary[4]) //cells on LOWER bdry
        {
            U2_nb[4] = 0.0;
            U2_nb[5] = cell_center[index_nb[5]].m_state.m_U[2];
            mu0 = 0.5*(cell_center[index].m_state.m_mu + cell_center[index_nb[5]].m_state.m_mu);
            diffusion[2] += mu0*((U2_nb[1]-U2_center) - (U2_center-U2_nb[0]))/dh[0]/dh[0];//mu*w_xx
            diffusion[2] += mu0*((U2_nb[3]-U2_center) - (U2_center-U2_nb[2]))/dh[1]/dh[1];//mu*w_yy
            diffusion[2] += mu0*((U2_nb[5]-U2_center) - (U2_center-U2_nb[4]))/dh[2]/dh[2];//mu*w_zz
        }
        else if (!bNoBoundary[5]) //cells on UPPER bdry
        {
            diffusion[2] += 0.0;//mu*w_xx
            diffusion[2] += 0.0;//mu*w_yy
            diffusion[2] += 0.0;//mu*w_zz
        }
        else assert(false);
        break;

    default:
        assert(false);
    }
} //end getViscousTerm_MAC_decoupled_vd
*/


/*
//viscous terms = mu*[laplace(U) + grad(div(U))/3]
//valid for constant mu
void Incompress_Solver_Smooth_3D_Cartesian::getViscousTerm_MAC_decoupled_vd(
        int *icoords,
        EBM_COORD xyz,
        double diffusion[3])
{
    int nb,index,index_nb[18];
    double mu0;
    double U0_nb[18],U1_nb[18],U2_nb[18],U0_center,U1_center,U2_center;
    bool bNoBoundary[6];
    double dh[3];
    COMPONENT comp;
    double crx_coords[MAXD];
    POINTER intfc_state;
    HYPER_SURF *hs;


    diffusion[0] = 0.0;
    diffusion[1] = 0.0;
    diffusion[2] = 0.0;

    int i = icoords[0];
    int j = icoords[1];
    int k = icoords[2];
    index = d_index3d(i,j,k,top_gmax);
    comp = top_comp[index];

    //6 neighbours of the center cell
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    //xy cut neighbours
    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);

    //yz cut neighbours
    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);

    //xz cut neighbours
    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);

    U0_center = cell_center[index].m_state.m_U[0];
    U1_center = cell_center[index].m_state.m_U[1];
    U2_center = cell_center[index].m_state.m_U[2];

    // 4 directions
    bNoBoundary[0] = YES;
    bNoBoundary[1] = YES;
    bNoBoundary[2] = YES;
    bNoBoundary[3] = YES;

    // LOWER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[4] = NO; //cells on LOWER bdry
    else
        bNoBoundary[4] = YES;

    // UPPER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[5] = NO; //cells on UPPER bdry
    else
        bNoBoundary[5] = YES;

    for (nb = 0; nb < 4; nb++) //4 directions only
    {
        U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
        U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
        U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
    }

    // non-cross derivative terms
    dh[0] = top_h[0];
    dh[1] = top_h[1];
    dh[2] = top_h[2];

    switch(xyz)
    {
    case COORD_X:
        if (bNoBoundary[4] && bNoBoundary[5])
        {
            U0_nb[4] = cell_center[index_nb[4]].m_state.m_U[0];
            U0_nb[5] = cell_center[index_nb[5]].m_state.m_U[0];

            U2_nb[1] = cell_center[index_nb[1]].m_state.m_U[2];
            U2_nb[4] = cell_center[index_nb[4]].m_state.m_U[2];
            U2_nb[15] = cell_center[index_nb[15]].m_state.m_U[2];
        }
        else if (!bNoBoundary[4]) //cells on LOWER bdry
        {
            U0_nb[4] = -U0_center;
            U0_nb[5] = cell_center[index_nb[5]].m_state.m_U[0];

            U2_nb[1] = cell_center[index_nb[1]].m_state.m_U[2];
            U2_nb[4] = 0.0;
            U2_nb[15] = 0.0;
        }
        else if (!bNoBoundary[5]) //cells on UPPER bdry
        {
            U0_nb[4] = cell_center[index_nb[4]].m_state.m_U[0];
            U0_nb[5] = -U0_center;

            U2_nb[1] = cell_center[index_nb[1]].m_state.m_U[2];
            U2_nb[4] = cell_center[index_nb[4]].m_state.m_U[2];
            U2_nb[15] = cell_center[index_nb[15]].m_state.m_U[2];
        }
        else assert(false);
        mu0 = 0.5*(cell_center[index].m_state.m_mu + cell_center[index_nb[1]].m_state.m_mu);
        U1_nb[1] = cell_center[index_nb[1]].m_state.m_U[1];
        U1_nb[2] = cell_center[index_nb[2]].m_state.m_U[1];
        U1_nb[7] = cell_center[index_nb[7]].m_state.m_U[1];
        diffusion[0] += mu0*4.0/3.0*((U0_nb[1]-U0_center) - (U0_center-U0_nb[0]))/dh[0]/dh[0];//mu*(4/3)*u_xx
        diffusion[0] += mu0*((U0_nb[3]-U0_center) - (U0_center-U0_nb[2]))/dh[1]/dh[1];//mu*u_yy
        diffusion[0] += mu0*((U0_nb[5]-U0_center) - (U0_center-U0_nb[4]))/dh[2]/dh[2];//mu*u_zz
        diffusion[0] += mu0/3.0*((U1_nb[1]-U1_nb[7]) - (U1_center-U1_nb[2]))/dh[0]/dh[1];//mu*(1/3)*v_yx
        diffusion[0] += mu0/3.0*((U2_nb[1]-U2_nb[15]) - (U2_center-U2_nb[4]))/dh[0]/dh[2];//mu*(1/3)*w_zx
        break;

    case COORD_Y:
        if (bNoBoundary[4] && bNoBoundary[5])
        {
            U1_nb[4] = cell_center[index_nb[4]].m_state.m_U[1];
            U1_nb[5] = cell_center[index_nb[5]].m_state.m_U[1];

            U2_nb[3] = cell_center[index_nb[3]].m_state.m_U[2];
            U2_nb[4] = cell_center[index_nb[4]].m_state.m_U[2];
            U2_nb[11] = cell_center[index_nb[11]].m_state.m_U[2];
        }
        else if (!bNoBoundary[4]) //cells on LOWER bdry
        {
            U1_nb[4] = -U1_center;
            U1_nb[5] = cell_center[index_nb[5]].m_state.m_U[1];

            U2_nb[3] = cell_center[index_nb[3]].m_state.m_U[2];
            U2_nb[4] = 0.0;
            U2_nb[11] = 0.0;
        }
        else if (!bNoBoundary[5]) //cells on UPPER bdry
        {
            U1_nb[4] = cell_center[index_nb[4]].m_state.m_U[1];
            U1_nb[5] = -U1_center;

            U2_nb[3] = cell_center[index_nb[3]].m_state.m_U[2];
            U2_nb[4] = cell_center[index_nb[4]].m_state.m_U[2];
            U2_nb[11] = cell_center[index_nb[11]].m_state.m_U[2];
        }
        else assert(false);
        mu0 = 0.5*(cell_center[index].m_state.m_mu + cell_center[index_nb[3]].m_state.m_mu);
        U0_nb[0] = cell_center[index_nb[0]].m_state.m_U[0];
        U0_nb[3] = cell_center[index_nb[3]].m_state.m_U[0];
        U0_nb[9] = cell_center[index_nb[9]].m_state.m_U[0];
        diffusion[1] += mu0*((U1_nb[1]-U1_center) - (U1_center-U1_nb[0]))/dh[0]/dh[0];//mu*v_xx
        diffusion[1] += mu0*4.0/3.0*((U1_nb[3]-U1_center) - (U1_center-U1_nb[2]))/dh[1]/dh[1];//mu*(4/3)*v_yy
        diffusion[1] += mu0*((U1_nb[5]-U1_center) - (U1_center-U1_nb[4]))/dh[2]/dh[2];//mu*v_zz
        diffusion[1] += mu0/3.0*((U0_nb[3]-U0_nb[9]) - (U0_center-U0_nb[0]))/dh[0]/dh[1];//mu*(1/3)*u_xy
        diffusion[1] += mu0/3.0*((U2_nb[3]-U2_nb[11]) - (U2_center-U2_nb[4]))/dh[1]/dh[2];//mu*(1/3)*w_zy
        break;

    case COORD_Z:
        if (bNoBoundary[4] && bNoBoundary[5])
        {
            U2_nb[4] = cell_center[index_nb[4]].m_state.m_U[2];
            U2_nb[5] = cell_center[index_nb[5]].m_state.m_U[2];

            U0_nb[0] = cell_center[index_nb[0]].m_state.m_U[0];
            U0_nb[5] = cell_center[index_nb[5]].m_state.m_U[0];
            U0_nb[17] = cell_center[index_nb[17]].m_state.m_U[0];

            U1_nb[2] = cell_center[index_nb[2]].m_state.m_U[1];
            U1_nb[5] = cell_center[index_nb[5]].m_state.m_U[1];
            U1_nb[13] = cell_center[index_nb[13]].m_state.m_U[1];

            mu0 = 0.5*(cell_center[index].m_state.m_mu + cell_center[index_nb[5]].m_state.m_mu);
            diffusion[2] += mu0*((U2_nb[1]-U2_center) - (U2_center-U2_nb[0]))/dh[0]/dh[0];//mu*w_xx
            diffusion[2] += mu0*((U2_nb[3]-U2_center) - (U2_center-U2_nb[2]))/dh[1]/dh[1];//mu*w_yy
            diffusion[2] += mu0*4.0/3.0*((U2_nb[5]-U2_center) - (U2_center-U2_nb[4]))/dh[2]/dh[2];//mu*(4/3)*w_zz
            diffusion[2] += mu0/3.0*((U0_nb[5]-U0_nb[17]) - (U0_center-U0_nb[0]))/dh[0]/dh[2];//mu*(1/3)*u_xz
            diffusion[2] += mu0/3.0*((U1_nb[5]-U1_nb[13]) - (U1_center-U1_nb[2]))/dh[1]/dh[2];//mu*(1/3)*v_yz
        }
        else if (!bNoBoundary[4]) //cells on LOWER bdry
        {
            U2_nb[4] = 0.0;
            U2_nb[5] = cell_center[index_nb[5]].m_state.m_U[2];

            U0_nb[0] = cell_center[index_nb[0]].m_state.m_U[0];
            U0_nb[5] = cell_center[index_nb[5]].m_state.m_U[0];
            U0_nb[17] = cell_center[index_nb[17]].m_state.m_U[0];

            U1_nb[2] = cell_center[index_nb[2]].m_state.m_U[1];
            U1_nb[5] = cell_center[index_nb[5]].m_state.m_U[1];
            U1_nb[13] = cell_center[index_nb[13]].m_state.m_U[1];

            mu0 = 0.5*(cell_center[index].m_state.m_mu + cell_center[index_nb[5]].m_state.m_mu);
            diffusion[2] += mu0*((U2_nb[1]-U2_center) - (U2_center-U2_nb[0]))/dh[0]/dh[0];//mu*w_xx
            diffusion[2] += mu0*((U2_nb[3]-U2_center) - (U2_center-U2_nb[2]))/dh[1]/dh[1];//mu*w_yy
            diffusion[2] += mu0*4.0/3.0*((U2_nb[5]-U2_center) - (U2_center-U2_nb[4]))/dh[2]/dh[2];//mu*(4/3)*w_zz
            diffusion[2] += mu0/3.0*((U0_nb[5]-U0_nb[17]) - (U0_center-U0_nb[0]))/dh[0]/dh[2];//mu*(1/3)*u_xz
            diffusion[2] += mu0/3.0*((U1_nb[5]-U1_nb[13]) - (U1_center-U1_nb[2]))/dh[1]/dh[2];//mu*(1/3)*v_yz
        }
        else if (!bNoBoundary[5]) //cells on UPPER bdry
        {
            U0_nb[0] = cell_center[index_nb[0]].m_state.m_U[0];
            U0_nb[5] = -U0_center;
            U0_nb[17] = -U0_nb[0];

            U1_nb[2] = cell_center[index_nb[2]].m_state.m_U[1];
            U1_nb[5] = -U1_center;
            U1_nb[13] = -U1_nb[2];

            mu0 = cell_center[index].m_state.m_mu;
            diffusion[2] += 0.0;//mu*w_xx
            diffusion[2] += 0.0;//mu*w_yy
            diffusion[2] += 0.0;//mu*(4/3)*w_zz
            diffusion[2] += mu0/3.0*((U0_nb[5]-U0_nb[17]) - (U0_center-U0_nb[0]))/dh[0]/dh[2];//mu*(1/3)*u_xz
            diffusion[2] += mu0/3.0*((U1_nb[5]-U1_nb[13]) - (U1_center-U1_nb[2]))/dh[1]/dh[2];//mu*(1/3)*v_yz
        }
        else assert(false);
        break;

    default:
        assert(false);
    }
} //end getViscousTerm_MAC_decoupled_vd
*/


/*               _                                                                                 _
                |  (4/3*mu*u_x)_x+(mu*u_y)_y+(mu*u_z)_z+(mu*v_x)_y+(mu*w_x)_z-2/3*[mu*(v_y+w_z)]_x  |
viscous terms = |  (mu*v_x)_x+(4/3*mu*v_y)_y+(mu*v_z)_z+(mu*u_y)_x+(mu*w_y)_z-2/3*[mu*(u_x+w_z)]_y  |
                |_ (mu*w_x)_x+(mu*w_y)_y+(4/3*mu*w_z)_z+(mu*u_z)_x+(mu*v_z)_y-2/3*[mu*(u_x+v_y)]_z _|
*/
//valid for variable dynamic viscosity
void Incompress_Solver_Smooth_3D_Cartesian::getViscousTerm_MAC_decoupled_vd(
        int *icoords,
        EBM_COORD xyz,
        double diffusion[MAXD])
{
    int nb,index,index_nb[18];
    double mu[2*MAXD];
    double U0_nb[18],U1_nb[18],U2_nb[18],U0_center,U1_center,U2_center;
    double mu_nb[18], mu_t_nb[18], mu_center, mu_t_center;
    int bNoBoundary[2*MAXD];
    double dh[MAXD], crx_coords[MAXD];
    COMPONENT comp;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    int i=icoords[0];
    int j=icoords[1];
    int k=icoords[2];
    bool useSGSCellCenter = YES; //use mu_t on the cell center. o.w. use mu_t on 3 cell-faces

    index = d_index3d(i,j,k,top_gmax);
    comp = top_comp[index];

    //6 neighbours of the center cell
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    //xy cut neighbours
    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);

    //yz cut neighbours
    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);

    //xz cut neighbours
    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);

    U0_center = cell_center[index].m_state.m_U[0];
    U1_center = cell_center[index].m_state.m_U[1];
    U2_center = cell_center[index].m_state.m_U[2];
    mu_center = cell_center[index].m_state.m_mu;
    if (useSGSCellCenter) {
        mu_t_center = cell_center[index].m_state.m_mu_turbulent[3];
    }
    else {
        printf("codes needed for mu_turbulent on cell-faces.\n");
        assert(false);
    }

    for (nb = 0; nb < 6; nb++)
    {
        checkBoundaryCondition(dir[nb],icoords,&bNoBoundary[nb],m_t_old,comp);
    }

    for (nb=0; nb<6; ++nb) {
        if(!bNoBoundary[nb] || bNoBoundary[nb] == 3) // reflect
        {
        U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
        U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
        U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
        mu_nb[nb] = cell_center[index_nb[nb]].m_state.m_mu;
        }
        if (useSGSCellCenter) {
            mu_t_nb[nb] = cell_center[index_nb[nb]].m_state.m_mu_turbulent[3];
        }
        else {
            printf("codes needed for mu_turbulent on cell-faces.\n");
            assert(false);
        }
    }

    //physical B.C. for LOWER & UPPER
    for (nb=4; nb<6; ++nb) {

        if (bNoBoundary[nb] == 2) { //cells on LOWER/UPPER bdry
            U0_nb[nb] = -U0_center;
            U1_nb[nb] = -U1_center;
            //w[4] = 0
            if (nb==4)          U2_nb[nb] = 0;
            //w[5] = -w[4]
            else if (nb==5)     U2_nb[nb] = -cell_center[index_nb[4]].m_state.m_U[2];
            else                assert(false);

            mu_nb[nb] = mu_center;
            if (useSGSCellCenter) {
                mu_t_nb[nb] = mu_t_center;
            }
            else {
                printf("codes needed for mu_turbulent on cell-faces.\n");
                assert(false);
            }
        }
    }

    for (i=0; i<MAXD; ++i) {
        diffusion[i] = 0;
        dh[i] = top_h[i];
    }

    //TODO: implement mu[] for mu_t on 3 cell-faces
    switch(xyz)
    {
    //u-face
    case COORD_X:
        if (bNoBoundary[2] == 3) // SOUTH DIR is Reflect
        {
            //printf("if no-slip, need implement\n");
        }
        if (bNoBoundary[3] == 3) // NORTH DIR is Reflect
        {
            //printf("if no-slip, need implement\n");
        }
        if ((!bNoBoundary[4] && !bNoBoundary[5]) || bNoBoundary[4] == 3 || bNoBoundary[5] == 3) // reflect
        {
            U2_nb[15] = cell_center[index_nb[15]].m_state.m_U[2];
            mu_nb[15] = cell_center[index_nb[15]].m_state.m_mu;
            mu_nb[16] = cell_center[index_nb[16]].m_state.m_mu;
            if (useSGSCellCenter) {
                mu_t_nb[15] = cell_center[index_nb[15]].m_state.m_mu_turbulent[3];
                mu_t_nb[16] = cell_center[index_nb[16]].m_state.m_mu_turbulent[3];
            }
            else {
                printf("codes needed for mu_turbulent on cell-faces.\n");
                assert(false);
            }
        }
        else if (bNoBoundary[4] == 2) //cells on LOWER bdry
        {
            U2_nb[15] = 0;
            mu_nb[15] = mu_nb[1];
            mu_nb[16] = cell_center[index_nb[16]].m_state.m_mu;
            if (useSGSCellCenter) {
                mu_t_nb[15] = mu_t_nb[1];
                mu_t_nb[16] = cell_center[index_nb[16]].m_state.m_mu_turbulent[3];
            }
            else {
                printf("codes needed for mu_turbulent on cell-faces.\n");
                assert(false);
            }
        }
        else if (bNoBoundary[5] == 2) //cells on UPPER bdry
        {
            U2_nb[15] = cell_center[index_nb[15]].m_state.m_U[2];
            mu_nb[15] = cell_center[index_nb[15]].m_state.m_mu;
            mu_nb[16] = mu_nb[1];
            if (useSGSCellCenter) {
                mu_t_nb[15] = cell_center[index_nb[15]].m_state.m_mu_turbulent[3];
                mu_t_nb[16] = mu_t_nb[1];
            }
            else {
                printf("codes needed for mu_turbulent on cell-faces.\n");
                assert(false);
            }
        }

        U1_nb[7] = cell_center[index_nb[7]].m_state.m_U[1];
        mu_nb[7] = cell_center[index_nb[7]].m_state.m_mu;
        mu_nb[8] = cell_center[index_nb[8]].m_state.m_mu;
        if (useSGSCellCenter) {
            mu_t_nb[7] = cell_center[index_nb[7]].m_state.m_mu_turbulent[3];
            mu_t_nb[8] = cell_center[index_nb[8]].m_state.m_mu_turbulent[3];
        }
        else {
            printf("codes needed for mu_turbulent on cell-faces.\n");
            assert(false);
        }
        mu[0] = mu_center + mu_t_center;
        mu[1] = mu_nb[1] + mu_t_nb[1];
        mu[2] = (mu_center + mu_nb[1] + mu_nb[7] + mu_nb[2] +
                 mu_t_center + mu_t_nb[1] + mu_t_nb[7] + mu_t_nb[2])/4;
        mu[3] = (mu_center + mu_nb[1] + mu_nb[8] + mu_nb[3] +
                 mu_t_center + mu_t_nb[1] + mu_t_nb[8] + mu_t_nb[3])/4;
        mu[4] = (mu_center + mu_nb[1] + mu_nb[15] + mu_nb[4] +
                 mu_t_center + mu_t_nb[1] + mu_t_nb[15] + mu_t_nb[4])/4;
        mu[5] = (mu_center + mu_nb[1] + mu_nb[16] + mu_nb[5] +
                 mu_t_center + mu_t_nb[1] + mu_t_nb[16] + mu_t_nb[5])/4;
        diffusion[0] += 4.0/3*(mu[1]*(U0_nb[1]-U0_center)-mu[0]*(U0_center-U0_nb[0]))/dh[0]/dh[0];//(4/3*mu*u_x)_x
        diffusion[0] += (mu[3]*(U0_nb[3]-U0_center)-mu[2]*(U0_center-U0_nb[2]))/dh[1]/dh[1];//(mu*u_y)_y
        diffusion[0] += (mu[5]*(U0_nb[5]-U0_center)-mu[4]*(U0_center-U0_nb[4]))/dh[2]/dh[2];//(mu*u_z)_z
        diffusion[0] += (mu[3]*(U1_nb[1]-U1_center)-mu[2]*(U1_nb[7]-U1_nb[2]))/dh[0]/dh[1];//(mu*v_x)_y
        diffusion[0] += (mu[5]*(U2_nb[1]-U2_center)-mu[4]*(U2_nb[15]-U2_nb[4]))/dh[0]/dh[2];//(mu*w_x)_z
        diffusion[0] += -2.0/3*(mu[1]*(U1_nb[1]-U1_nb[7])-mu[0]*(U1_center-U1_nb[2]))/dh[0]/dh[1];//(-2/3)*(mu*v_y)_x
        diffusion[0] += -2.0/3*(mu[1]*(U2_nb[1]-U2_nb[15])-mu[0]*(U2_center-U2_nb[4]))/dh[0]/dh[2];//(-2/3)*(mu*w_z)_x
        break;

    //v-face
    case COORD_Y:
        if ((!bNoBoundary[4] && !bNoBoundary[5]) || bNoBoundary[4] == 3 || bNoBoundary[5] == 3)
        {
            U2_nb[11] = cell_center[index_nb[11]].m_state.m_U[2];
            mu_nb[11] = cell_center[index_nb[11]].m_state.m_mu;
            mu_nb[12] = cell_center[index_nb[12]].m_state.m_mu;
            if (useSGSCellCenter) {
                mu_t_nb[11] = cell_center[index_nb[11]].m_state.m_mu_turbulent[3];
                mu_t_nb[12] = cell_center[index_nb[12]].m_state.m_mu_turbulent[3];
            }
            else {
                printf("codes needed for mu_turbulent on cell-faces.\n");
                assert(false);
            }
        }
        else if (bNoBoundary[4] == 2) //cells on LOWER bdry
        {
            U2_nb[11] = 0;
            mu_nb[11] = mu_nb[3];
            mu_nb[12] = cell_center[index_nb[12]].m_state.m_mu;
            if (useSGSCellCenter) {
                mu_t_nb[11] = mu_t_nb[3];
                mu_t_nb[12] = cell_center[index_nb[12]].m_state.m_mu_turbulent[3];
            }
            else {
                printf("codes needed for mu_turbulent on cell-faces.\n");
                assert(false);
            }
        }
        else if (bNoBoundary[5] == 2) //cells on UPPER bdry
        {
            U2_nb[11] = cell_center[index_nb[11]].m_state.m_U[2];
            mu_nb[11] = cell_center[index_nb[11]].m_state.m_mu;
            mu_nb[12] = mu_nb[3];
            if (useSGSCellCenter) {
                mu_t_nb[11] = cell_center[index_nb[11]].m_state.m_mu_turbulent[3];
                mu_t_nb[12] = mu_t_nb[3];
            }
            else {
                printf("codes needed for mu_turbulent on cell-faces.\n");
                assert(false);
            }
        }

        U0_nb[9] = cell_center[index_nb[9]].m_state.m_U[0];
        mu_nb[8] = cell_center[index_nb[8]].m_state.m_mu;
        mu_nb[9] = cell_center[index_nb[9]].m_state.m_mu;
        if (useSGSCellCenter) {
            mu_t_nb[8] = cell_center[index_nb[8]].m_state.m_mu_turbulent[3];
            mu_t_nb[9] = cell_center[index_nb[9]].m_state.m_mu_turbulent[3];
        }
        else {
            printf("codes needed for mu_turbulent on cell-faces.\n");
            assert(false);
        }
        mu[0] = (mu_center + mu_nb[0] + mu_nb[9] + mu_nb[3] +
                 mu_t_center + mu_t_nb[0] + mu_t_nb[9] + mu_t_nb[3])/4;
        mu[1] = (mu_center + mu_nb[1] + mu_nb[8] + mu_nb[3] +
                 mu_t_center + mu_t_nb[1] + mu_t_nb[8] + mu_t_nb[3])/4;
        mu[2] = mu_center + mu_t_center;
        mu[3] = mu_nb[3] + mu_t_nb[3];
        mu[4] = (mu_center + mu_nb[3] + mu_nb[11] + mu_nb[4] +
                 mu_t_center + mu_t_nb[3] + mu_t_nb[11] + mu_t_nb[4])/4;
        mu[5] = (mu_center + mu_nb[3] + mu_nb[12] + mu_nb[5] +
                 mu_t_center + mu_t_nb[3] + mu_t_nb[12] + mu_t_nb[5])/4;
        diffusion[1] += (mu[1]*(U1_nb[1]-U1_center)-mu[0]*(U1_center-U1_nb[0]))/dh[0]/dh[0];//(mu*v_x)_x
        diffusion[1] += 4.0/3*(mu[3]*(U1_nb[3]-U1_center)-mu[2]*(U1_center-U1_nb[2]))/dh[1]/dh[1];//(4/3*mu*v_y)_y
        diffusion[1] += (mu[5]*(U1_nb[5]-U1_center)-mu[4]*(U1_center-U1_nb[4]))/dh[2]/dh[2];//(mu*v_z)_z
        diffusion[1] += (mu[1]*(U0_nb[3]-U0_center)-mu[0]*(U0_nb[9]-U0_nb[0]))/dh[0]/dh[1];//(mu*u_y)_x
        diffusion[1] += (mu[5]*(U2_nb[3]-U2_center)-mu[4]*(U2_nb[11]-U2_nb[4]))/dh[1]/dh[2];//(mu*w_y)_z
        diffusion[1] += -2.0/3*(mu[3]*(U0_nb[3]-U0_nb[9])-mu[2]*(U0_center-U0_nb[0]))/dh[0]/dh[1];//(-2/3)*(mu*u_x)_y
        diffusion[1] += -2.0/3*(mu[3]*(U2_nb[3]-U2_nb[11])-mu[2]*(U2_center-U2_nb[4]))/dh[1]/dh[2];//(-2/3)*(mu*w_z)_y
        break;

    //w-face
    case COORD_Z:
        if (bNoBoundary[5] == 2) //cells on UPPER bdry
        {
            U0_nb[17] = -U0_nb[0];
            U1_nb[13] = -U1_nb[2];
            mu_nb[12] = mu_nb[3];
            mu_nb[13] = mu_nb[2];
            mu_nb[16] = mu_nb[1];
            mu_nb[17] = mu_nb[0];
            if (useSGSCellCenter) {
                mu_t_nb[12] = mu_t_nb[3];
                mu_t_nb[13] = mu_t_nb[2];
                mu_t_nb[16] = mu_t_nb[1];
                mu_t_nb[17] = mu_t_nb[0];
            }
            else {
                printf("codes needed for mu_turbulent on cell-faces.\n");
                assert(false);
            }
        }
        else //other cells
        {
            U0_nb[17] = cell_center[index_nb[17]].m_state.m_U[0];
            U1_nb[13] = cell_center[index_nb[13]].m_state.m_U[1];
            mu_nb[12] = cell_center[index_nb[12]].m_state.m_mu;
            mu_nb[13] = cell_center[index_nb[13]].m_state.m_mu;
            mu_nb[16] = cell_center[index_nb[16]].m_state.m_mu;
            mu_nb[17] = cell_center[index_nb[17]].m_state.m_mu;
            if (useSGSCellCenter) {
                mu_t_nb[12] = cell_center[index_nb[12]].m_state.m_mu_turbulent[3];
                mu_t_nb[13] = cell_center[index_nb[13]].m_state.m_mu_turbulent[3];
                mu_t_nb[16] = cell_center[index_nb[16]].m_state.m_mu_turbulent[3];
                mu_t_nb[17] = cell_center[index_nb[17]].m_state.m_mu_turbulent[3];
            }
            else {
                printf("codes needed for mu_turbulent on cell-faces.\n");
                assert(false);
            }
        }

        mu[0] = (mu_center + mu_nb[0] + mu_nb[17] + mu_nb[5] +
                 mu_t_center + mu_t_nb[0] + mu_t_nb[17] + mu_t_nb[5])/4;
        mu[1] = (mu_center + mu_nb[1] + mu_nb[16] + mu_nb[5] +
                 mu_t_center + mu_t_nb[1] + mu_t_nb[16] + mu_t_nb[5])/4;
        mu[2] = (mu_center + mu_nb[2] + mu_nb[13] + mu_nb[5] +
                 mu_t_center + mu_t_nb[2] + mu_t_nb[13] + mu_t_nb[5])/4;
        mu[3] = (mu_center + mu_nb[3] + mu_nb[12] + mu_nb[5] +
                 mu_t_center + mu_t_nb[3] + mu_t_nb[12] + mu_t_nb[5])/4;
        mu[4] = mu_center + mu_t_center;
        mu[5] = mu_nb[5] + mu_t_nb[5];
        diffusion[2] += (mu[1]*(U2_nb[1]-U2_center)-mu[0]*(U2_center-U2_nb[0]))/dh[0]/dh[0];//(mu*w_x)_x
        diffusion[2] += (mu[3]*(U2_nb[3]-U2_center)-mu[2]*(U2_center-U2_nb[2]))/dh[1]/dh[1];//(mu*w_y)_y
        diffusion[2] += 4.0/3*(mu[5]*(U2_nb[5]-U2_center)-mu[4]*(U2_center-U2_nb[4]))/dh[2]/dh[2];//(4/3*mu*w_z)_z
        diffusion[2] += (mu[1]*(U0_nb[5]-U0_center)-mu[0]*(U0_nb[17]-U0_nb[0]))/dh[0]/dh[2];//(mu*u_z)_x
        diffusion[2] += (mu[3]*(U1_nb[5]-U1_center)-mu[2]*(U1_nb[13]-U1_nb[2]))/dh[1]/dh[2];//(mu*v_z)_y
        diffusion[2] += -2.0/3*(mu[5]*(U0_nb[5]-U0_nb[17])-mu[4]*(U0_center-U0_nb[0]))/dh[0]/dh[2];//(-2/3)*(mu*u_x)_z
        diffusion[2] += -2.0/3*(mu[5]*(U1_nb[5]-U1_nb[13])-mu[4]*(U1_center-U1_nb[2]))/dh[1]/dh[2];//(-2/3)*(mu*v_y)_z
        break;

    default:
        assert(false);
    }
} //end getViscousTerm_MAC_decoupled_vd


//Stability issue proposed by M. Monion, refer to JCP 123 pp.435-449(1996)
void Incompress_Solver_Smooth_3D_Cartesian::getTransverseDTerm_Velocity_MAC_vd(
        int *icoords,
        double transverseD[3])
{
    int i, j, k;
    int bNoBoundary[6];
    int index,index_nb[18],ICoords[MAXD];
    double tmp, dx, dy, dz, rho, rho_nb;
    COMPONENT comp;
    double crx_coords[MAXD], diff[3], diff_nb[3];
    int nb;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

    transverseD[0] = 0.0;
    transverseD[1] = 0.0;
    transverseD[2] = 0.0;

    dx = top_h[0];
    dy = top_h[1];
    dz = top_h[2];

    i = icoords[0];
    j = icoords[1];
    k = icoords[2];
    index = d_index3d(i,j,k,top_gmax);
    comp = top_comp[index];

    //6 neighbours of the center cell
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    //xy cut neighbours
    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);

    //yz cut neighbours
    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);

    //xz cut neighbours
    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);

    for (nb = 0; nb < 6; nb++)
    {
        checkBoundaryCondition(dir[nb],icoords,&bNoBoundary[nb],m_t_old,comp);
    }
    /*
    // 4 directions
    bNoBoundary[0] = YES;
    bNoBoundary[1] = YES;
    bNoBoundary[2] = YES;
    bNoBoundary[3] = YES;

    // LOWER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[4] = NO; //cells on LOWER bdry
    else
        bNoBoundary[4] = YES;

    // UPPER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[5] = NO; //cells on UPPER bdry
    else
        bNoBoundary[5] = YES;
    */


///////////////////////////////////////////////////
//  transverseD terms on u_faces = v*u_y + w*u_z
///////////////////////////////////////////////////

    //v*u_y
    tmp = (cell_center[index].m_state.m_U[1] +
           cell_center[index_nb[1]].m_state.m_U[1] +
           cell_center[index_nb[7]].m_state.m_U[1] +
           cell_center[index_nb[2]].m_state.m_U[1])/4.0;
    tmp /= dy;
    if (tmp > 0.0)
    {
        transverseD[0] += tmp*(cell_center[index].m_state.m_U[0] - cell_center[index_nb[2]].m_state.m_U[0]);

        //Stability issue proposed by M. Monion, refer to JCP 123 pp.435-449(1996)
        getViscousTerm_MAC_decoupled_vd(icoords,COORD_X,diff);
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1] - 1;
        ICoords[2] = icoords[2];
        getViscousTerm_MAC_decoupled_vd(ICoords,COORD_X,diff_nb);
        rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[1]].m_state.m_rho);
        rho_nb = 2.0/(1.0/cell_center[index_nb[2]].m_state.m_rho + 1.0/cell_center[index_nb[7]].m_state.m_rho);
        transverseD[0] += tmp*m_dt/2.0*(diff[0]/rho - diff_nb[0]/rho_nb);
    }
    else //tmp <= 0.0
    {
        transverseD[0] += tmp*(cell_center[index_nb[3]].m_state.m_U[0] - cell_center[index].m_state.m_U[0]);

        getViscousTerm_MAC_decoupled_vd(icoords,COORD_X,diff);
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1] + 1;
        ICoords[2] = icoords[2];
        getViscousTerm_MAC_decoupled_vd(ICoords,COORD_X,diff_nb);
        rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[1]].m_state.m_rho);
        rho_nb = 2.0/(1.0/cell_center[index_nb[3]].m_state.m_rho + 1.0/cell_center[index_nb[8]].m_state.m_rho);
        transverseD[0] += tmp*m_dt/2.0*(diff_nb[0]/rho_nb - diff[0]/rho);
    }

    //w*u_z
    if ((!bNoBoundary[4] && !bNoBoundary[5]) || bNoBoundary[4] == 3 || bNoBoundary[5] == 3)
    {
        tmp = (cell_center[index].m_state.m_U[2] +
               cell_center[index_nb[1]].m_state.m_U[2] +
               cell_center[index_nb[15]].m_state.m_U[2] +
               cell_center[index_nb[4]].m_state.m_U[2])/4.0;
        tmp /= dz;
        if (tmp > 0.0)
        {
            transverseD[0] += tmp*(cell_center[index].m_state.m_U[0] - cell_center[index_nb[4]].m_state.m_U[0]);

            getViscousTerm_MAC_decoupled_vd(icoords,COORD_X,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] - 1;
            getViscousTerm_MAC_decoupled_vd(ICoords,COORD_X,diff_nb);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[1]].m_state.m_rho);
            rho_nb = 2.0/(1.0/cell_center[index_nb[4]].m_state.m_rho + 1.0/cell_center[index_nb[15]].m_state.m_rho);
            transverseD[0] += tmp*m_dt/2.0*(diff[0]/rho - diff_nb[0]/rho_nb);
        }
        else //tmp <= 0.0
        {
            transverseD[0] += tmp*(cell_center[index_nb[5]].m_state.m_U[0] - cell_center[index].m_state.m_U[0]);

            getViscousTerm_MAC_decoupled_vd(icoords,COORD_X,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] + 1;
            getViscousTerm_MAC_decoupled_vd(ICoords,COORD_X,diff_nb);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[1]].m_state.m_rho);
            rho_nb = 2.0/(1.0/cell_center[index_nb[5]].m_state.m_rho + 1.0/cell_center[index_nb[16]].m_state.m_rho);
            transverseD[0] += tmp*m_dt/2.0*(diff_nb[0]/rho_nb - diff[0]/rho);
        }
    }
    else if (bNoBoundary[4] == 2) //cells on LOWER boundary
    {
        tmp = (cell_center[index].m_state.m_U[2] +
               cell_center[index_nb[1]].m_state.m_U[2] + 0.0 + 0.0)/4.0;
        tmp /= dz;
        if (tmp > 0.0)
        {
            transverseD[0] += tmp*(2.0*cell_center[index].m_state.m_U[0]);

            //laplace of u or v at (i,j,-1) equals to minus of that at (i,j,0)
            getViscousTerm_MAC_decoupled_vd(icoords,COORD_X,diff);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[1]].m_state.m_rho);
            transverseD[0] += tmp*m_dt/2.0*(2.0*diff[0]/rho);
        }
        else //tmp <= 0.0
        {
            transverseD[0] += tmp*(cell_center[index_nb[5]].m_state.m_U[0] - cell_center[index].m_state.m_U[0]);

            getViscousTerm_MAC_decoupled_vd(icoords,COORD_X,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] + 1;
            getViscousTerm_MAC_decoupled_vd(ICoords,COORD_X,diff_nb);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[1]].m_state.m_rho);
            rho_nb = 2.0/(1.0/cell_center[index_nb[5]].m_state.m_rho + 1.0/cell_center[index_nb[16]].m_state.m_rho);
            transverseD[0] += tmp*m_dt/2.0*(diff_nb[0]/rho_nb - diff[0]/rho);
        }
    }
    else if (bNoBoundary[5] == 2) //cells on UPPER boundary
    {
        tmp = (cell_center[index].m_state.m_U[2] +
               cell_center[index_nb[1]].m_state.m_U[2] +
               cell_center[index_nb[15]].m_state.m_U[2] +
               cell_center[index_nb[4]].m_state.m_U[2])/4.0;
        tmp /= dz;
        if (tmp > 0.0)
        {
            transverseD[0] += tmp*(cell_center[index].m_state.m_U[0] - cell_center[index_nb[4]].m_state.m_U[0]);

            getViscousTerm_MAC_decoupled_vd(icoords,COORD_X,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] - 1;
            getViscousTerm_MAC_decoupled_vd(ICoords,COORD_X,diff_nb);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[1]].m_state.m_rho);
            rho_nb = 2.0/(1.0/cell_center[index_nb[4]].m_state.m_rho + 1.0/cell_center[index_nb[15]].m_state.m_rho);
            transverseD[0] += tmp*m_dt/2.0*(diff[0]/rho - diff_nb[0]/rho_nb);
        }
        else //tmp <= 0.0
        {
            transverseD[0] += tmp*(-2.0*cell_center[index].m_state.m_U[0]);

            getViscousTerm_MAC_decoupled_vd(icoords,COORD_X,diff);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[1]].m_state.m_rho);
            transverseD[0] += tmp*m_dt/2.0*(-2.0*diff[0]/rho);
        }
    }
    else assert(false);


///////////////////////////////////////////////////
//  transverseD terms on v_faces = u*v_x + w*v_z
///////////////////////////////////////////////////

    //u*v_x
    tmp = (cell_center[index].m_state.m_U[0] +
           cell_center[index_nb[0]].m_state.m_U[0] +
           cell_center[index_nb[9]].m_state.m_U[0] +
           cell_center[index_nb[3]].m_state.m_U[0])/4.0;
    tmp /= dx;
    if (tmp > 0.0)
    {
        transverseD[1] += tmp*(cell_center[index].m_state.m_U[1] - cell_center[index_nb[0]].m_state.m_U[1]);

        getViscousTerm_MAC_decoupled_vd(icoords,COORD_Y,diff);
        ICoords[0] = icoords[0] - 1;
        ICoords[1] = icoords[1];
        ICoords[2] = icoords[2];
        getViscousTerm_MAC_decoupled_vd(ICoords,COORD_Y,diff_nb);
        rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[3]].m_state.m_rho);
        rho_nb = 2.0/(1.0/cell_center[index_nb[0]].m_state.m_rho + 1.0/cell_center[index_nb[9]].m_state.m_rho);
        transverseD[1] += tmp*m_dt/2.0*(diff[1]/rho - diff_nb[1]/rho_nb);
    }
    else //tmp <= 0.0
    {
        transverseD[1] += tmp*(cell_center[index_nb[1]].m_state.m_U[1] - cell_center[index].m_state.m_U[1]);

        getViscousTerm_MAC_decoupled_vd(icoords,COORD_Y,diff);
        ICoords[0] = icoords[0] + 1;
        ICoords[1] = icoords[1];
        ICoords[2] = icoords[2];
        getViscousTerm_MAC_decoupled_vd(ICoords,COORD_Y,diff_nb);
        rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[3]].m_state.m_rho);
        rho_nb = 2.0/(1.0/cell_center[index_nb[1]].m_state.m_rho + 1.0/cell_center[index_nb[8]].m_state.m_rho);
        transverseD[1] += tmp*m_dt/2.0*(diff_nb[1]/rho_nb - diff[1]/rho);
    }

    //w*v_z
    if ((!bNoBoundary[4] && !bNoBoundary[5]) || bNoBoundary[4] == 3 || bNoBoundary[5] == 3)
    {
        tmp = (cell_center[index].m_state.m_U[2] +
               cell_center[index_nb[3]].m_state.m_U[2] +
               cell_center[index_nb[11]].m_state.m_U[2] +
               cell_center[index_nb[4]].m_state.m_U[2])/4.0;
        tmp /= dz;
        if (tmp > 0.0)
        {
            transverseD[1] += tmp*(cell_center[index].m_state.m_U[1] - cell_center[index_nb[4]].m_state.m_U[1]);

            getViscousTerm_MAC_decoupled_vd(icoords,COORD_Y,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] - 1;
            getViscousTerm_MAC_decoupled_vd(ICoords,COORD_Y,diff_nb);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[3]].m_state.m_rho);
            rho_nb = 2.0/(1.0/cell_center[index_nb[4]].m_state.m_rho + 1.0/cell_center[index_nb[11]].m_state.m_rho);
            transverseD[1] += tmp*m_dt/2.0*(diff[1]/rho - diff_nb[1]/rho_nb);
        }
        else //tmp <= 0.0
        {
            transverseD[1] += tmp*(cell_center[index_nb[5]].m_state.m_U[1] - cell_center[index].m_state.m_U[1]);

            getViscousTerm_MAC_decoupled_vd(icoords,COORD_Y,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] + 1;
            getViscousTerm_MAC_decoupled_vd(ICoords,COORD_Y,diff_nb);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[3]].m_state.m_rho);
            rho_nb = 2.0/(1.0/cell_center[index_nb[5]].m_state.m_rho + 1.0/cell_center[index_nb[12]].m_state.m_rho);
            transverseD[1] += tmp*m_dt/2.0*(diff_nb[1]/rho_nb - diff[1]/rho);
        }
    }
    else if (bNoBoundary[4] == 2) //cells on LOWER boundary
    {
        tmp = (cell_center[index].m_state.m_U[2] +
               cell_center[index_nb[3]].m_state.m_U[2] + 0.0 + 0.0)/4.0;
        tmp /= dz;
        if (tmp > 0.0)
        {
            transverseD[1] += tmp*(2.0*cell_center[index].m_state.m_U[1]);

            //laplace of u or v at (i,j,-1) equals to minus of that at (i,j,0)
            getViscousTerm_MAC_decoupled_vd(icoords,COORD_Y,diff);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[3]].m_state.m_rho);
            transverseD[1] += tmp*m_dt/2.0*(2.0*diff[1]/rho);
        }
        else //tmp <= 0.0
        {
            transverseD[1] += tmp*(cell_center[index_nb[5]].m_state.m_U[1] - cell_center[index].m_state.m_U[1]);

            getViscousTerm_MAC_decoupled_vd(icoords,COORD_Y,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] + 1;
            getViscousTerm_MAC_decoupled_vd(ICoords,COORD_Y,diff_nb);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[3]].m_state.m_rho);
            rho_nb = 2.0/(1.0/cell_center[index_nb[5]].m_state.m_rho + 1.0/cell_center[index_nb[12]].m_state.m_rho);
            transverseD[1] += tmp*m_dt/2.0*(diff_nb[1]/rho_nb - diff[1]/rho);
        }
    }
    else if (bNoBoundary[5] == 2) //cells on UPPER boundary
    {
        tmp = (cell_center[index].m_state.m_U[2] +
               cell_center[index_nb[3]].m_state.m_U[2] +
               cell_center[index_nb[11]].m_state.m_U[2] +
               cell_center[index_nb[4]].m_state.m_U[2])/4.0;
        tmp /= dz;
        if (tmp > 0.0)
        {
            transverseD[1] += tmp*(cell_center[index].m_state.m_U[1] - cell_center[index_nb[4]].m_state.m_U[1]);

            getViscousTerm_MAC_decoupled_vd(icoords,COORD_Y,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] - 1;
            getViscousTerm_MAC_decoupled_vd(ICoords,COORD_Y,diff_nb);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[3]].m_state.m_rho);
            rho_nb = 2.0/(1.0/cell_center[index_nb[4]].m_state.m_rho + 1.0/cell_center[index_nb[11]].m_state.m_rho);
            transverseD[1] += tmp*m_dt/2.0*(diff[1]/rho - diff_nb[1]/rho_nb);
        }
        else //tmp <= 0.0
        {
            transverseD[1] += tmp*(-2.0*cell_center[index].m_state.m_U[1]);

            getViscousTerm_MAC_decoupled_vd(icoords,COORD_Y,diff);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[3]].m_state.m_rho);
            transverseD[1] += tmp*m_dt/2.0*(-2.0*diff[1]/rho);
        }
    }
    else assert(false);


///////////////////////////////////////////////////
//  transverseD terms on w_faces = u*w_x + v*w_y
///////////////////////////////////////////////////

    //u*w_x
    if (bNoBoundary[5] == 2) //cells on UPPER boundary
    {
        tmp = 0.0; //u = 0 on UPPER bdry
        transverseD[2] += 0.0;
    }
    else //other cells
    {
        tmp = (cell_center[index].m_state.m_U[0] +
               cell_center[index_nb[0]].m_state.m_U[0] +
               cell_center[index_nb[17]].m_state.m_U[0] +
               cell_center[index_nb[5]].m_state.m_U[0])/4.0;
        tmp /= dx;
        if (tmp > 0.0)
        {
            transverseD[2] += tmp*(cell_center[index].m_state.m_U[2] - cell_center[index_nb[0]].m_state.m_U[2]);

            getViscousTerm_MAC_decoupled_vd(icoords,COORD_Z,diff);
            ICoords[0] = icoords[0] - 1;
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2];
            getViscousTerm_MAC_decoupled_vd(ICoords,COORD_Z,diff_nb);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[5]].m_state.m_rho);
            rho_nb = 2.0/(1.0/cell_center[index_nb[0]].m_state.m_rho + 1.0/cell_center[index_nb[17]].m_state.m_rho);
            transverseD[2] += tmp*m_dt/2.0*(diff[2]/rho - diff_nb[2]/rho_nb);
        }
        else //tmp <= 0.0
        {
            transverseD[2] += tmp*(cell_center[index_nb[1]].m_state.m_U[2] - cell_center[index].m_state.m_U[2]);

            getViscousTerm_MAC_decoupled_vd(icoords,COORD_Z,diff);
            ICoords[0] = icoords[0] + 1;
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2];
            getViscousTerm_MAC_decoupled_vd(ICoords,COORD_Z,diff_nb);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[5]].m_state.m_rho);
            rho_nb = 2.0/(1.0/cell_center[index_nb[1]].m_state.m_rho + 1.0/cell_center[index_nb[16]].m_state.m_rho);
            transverseD[2] += tmp*m_dt/2.0*(diff_nb[2]/rho_nb - diff[2]/rho);
        }
    }

    //v*w_y
    if (bNoBoundary[5] == 2) //UPPER boundary
    {
        tmp = 0.0; //v = 0 on UPPER bdry
        transverseD[2] += 0.0;
    }
    else //other cells
    {
        tmp = (cell_center[index].m_state.m_U[1] +
               cell_center[index_nb[2]].m_state.m_U[1] +
               cell_center[index_nb[13]].m_state.m_U[1] +
               cell_center[index_nb[5]].m_state.m_U[1])/4.0;
        tmp /= dy;
        if (tmp > 0.0)
        {
            transverseD[2] += tmp*(cell_center[index].m_state.m_U[2] - cell_center[index_nb[2]].m_state.m_U[2]);

            getViscousTerm_MAC_decoupled_vd(icoords,COORD_Z,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1] - 1;
            ICoords[2] = icoords[2];
            getViscousTerm_MAC_decoupled_vd(ICoords,COORD_Z,diff_nb);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[5]].m_state.m_rho);
            rho_nb = 2.0/(1.0/cell_center[index_nb[2]].m_state.m_rho + 1.0/cell_center[index_nb[13]].m_state.m_rho);
            transverseD[2] += tmp*m_dt/2.0*(diff[2]/rho - diff_nb[2]/rho_nb);
        }
        else //tmp <= 0.0
        {
            transverseD[2] += tmp*(cell_center[index_nb[3]].m_state.m_U[2] - cell_center[index].m_state.m_U[2]);

            getViscousTerm_MAC_decoupled_vd(icoords,COORD_Z,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1] + 1;
            ICoords[2] = icoords[2];
            getViscousTerm_MAC_decoupled_vd(ICoords,COORD_Z,diff_nb);
            rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[5]].m_state.m_rho);
            rho_nb = 2.0/(1.0/cell_center[index_nb[3]].m_state.m_rho + 1.0/cell_center[index_nb[12]].m_state.m_rho);
            transverseD[2] += tmp*m_dt/2.0*(diff_nb[2]/rho_nb - diff[2]/rho);
        }
    }
} //end getTransverseDTerm_Velocity_MAC_vd


/*
//BCG, refer to JCP 85 pp257-283(1989)
void Incompress_Solver_Smooth_3D_Cartesian::getTransverseDTerm_Velocity_MAC_vd(
        int *icoords,
        double transverseD[3])
{
    int i, j, k;
    bool bNoBoundary[6];
    int index,index_nb[18],ICoords[MAXD];
    double tmp, dx, dy, dz, rho, rho_nb;
    COMPONENT comp;
    double crx_coords[MAXD], diff[3], diff_nb[3];
    POINTER intfc_state;
    HYPER_SURF *hs;

    transverseD[0] = 0.0;
    transverseD[1] = 0.0;
    transverseD[2] = 0.0;

    dx = top_h[0];
    dy = top_h[1];
    dz = top_h[2];

    i = icoords[0];
    j = icoords[1];
    k = icoords[2];
    index = d_index3d(i,j,k,top_gmax);
    comp = cell_center[index].comp;

    //6 neighbours of the center cell
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    //xy cut neighbours
    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);

    //yz cut neighbours
    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);

    //xz cut neighbours
    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);

    // 4 directions
    bNoBoundary[0] = YES;
    bNoBoundary[1] = YES;
    bNoBoundary[2] = YES;
    bNoBoundary[3] = YES;

    // LOWER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[4] = NO; //cells on LOWER bdry
    else
        bNoBoundary[4] = YES;

    // UPPER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[5] = NO; //cells on UPPER bdry
    else
        bNoBoundary[5] = YES;


///////////////////////////////////////////////////
//  transverseD terms on u_faces = v*u_y + w*u_z
///////////////////////////////////////////////////

    //v*u_y
    tmp = (cell_center[index].m_state.m_U[1] +
           cell_center[index_nb[1]].m_state.m_U[1] +
           cell_center[index_nb[7]].m_state.m_U[1] +
           cell_center[index_nb[2]].m_state.m_U[1])/4.0;
    tmp /= dy;
    if (tmp > 0.0)
    {
        transverseD[0] += tmp*(cell_center[index].m_state.m_U[0] - cell_center[index_nb[2]].m_state.m_U[0]);

        //Stability issue proposed by M. Monion, refer to JCP 123 pp.435-449(1996)
        getLimitedSlope_Velocity_MAC_vd(icoords,COORD_X,diff);
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1] - 1;
        ICoords[2] = icoords[2];
        getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_X,diff_nb);
        transverseD[0] += tmp*dy*(1.0 - m_dt*tmp)/2.0*(diff[1] - diff_nb[1]);
    }
    else //tmp <= 0.0
    {
        transverseD[0] += tmp*(cell_center[index_nb[3]].m_state.m_U[0] - cell_center[index].m_state.m_U[0]);

        getLimitedSlope_Velocity_MAC_vd(icoords,COORD_X,diff);
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1] + 1;
        ICoords[2] = icoords[2];
        getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_X,diff_nb);
        transverseD[0] += tmp*dy*(1.0 + m_dt*tmp)/2.0*(diff[1] - diff_nb[1]);
    }

    //w*u_z
    if (bNoBoundary[4] && bNoBoundary[5])
    {
        tmp = (cell_center[index].m_state.m_U[2] +
               cell_center[index_nb[1]].m_state.m_U[2] +
               cell_center[index_nb[15]].m_state.m_U[2] +
               cell_center[index_nb[4]].m_state.m_U[2])/4.0;
        tmp /= dz;
        if (tmp > 0.0)
        {
            transverseD[0] += tmp*(cell_center[index].m_state.m_U[0] - cell_center[index_nb[4]].m_state.m_U[0]);

            getLimitedSlope_Velocity_MAC_vd(icoords,COORD_X,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] - 1;
            getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_X,diff_nb);
            transverseD[0] += tmp*dz*(1.0 - m_dt*tmp)/2.0*(diff[2] - diff_nb[2]);
        }
        else //tmp <= 0.0
        {
            transverseD[0] += tmp*(cell_center[index_nb[5]].m_state.m_U[0] - cell_center[index].m_state.m_U[0]);

            getLimitedSlope_Velocity_MAC_vd(icoords,COORD_X,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] + 1;
            getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_X,diff_nb);
            transverseD[0] += tmp*dz*(1.0 + m_dt*tmp)/2.0*(diff[2] - diff_nb[2]);
        }
    }
    else if (!bNoBoundary[4]) //cells on LOWER boundary
    {
        tmp = (cell_center[index].m_state.m_U[2] +
               cell_center[index_nb[1]].m_state.m_U[2] + 0.0 + 0.0)/4.0;
        tmp /= dz;
        if (tmp > 0.0)
        {
            transverseD[0] += tmp*(2.0*cell_center[index].m_state.m_U[0]);

            //u_z or v_z at (i,j,-1) equals to that at (i,j,0)
            transverseD[0] += 0.0;
        }
        else //tmp <= 0.0
        {
            transverseD[0] += tmp*(cell_center[index_nb[5]].m_state.m_U[0] - cell_center[index].m_state.m_U[0]);

            getLimitedSlope_Velocity_MAC_vd(icoords,COORD_X,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] + 1;
            getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_X,diff_nb);
            transverseD[0] += tmp*dz*(1.0 + m_dt*tmp)/2.0*(diff[2] - diff_nb[2]);
        }
    }
    else if (!bNoBoundary[5]) //cells on UPPER boundary
    {
        tmp = (cell_center[index].m_state.m_U[2] +
               cell_center[index_nb[1]].m_state.m_U[2] +
               cell_center[index_nb[15]].m_state.m_U[2] +
               cell_center[index_nb[4]].m_state.m_U[2])/4.0;
        tmp /= dz;
        if (tmp > 0.0)
        {
            transverseD[0] += tmp*(cell_center[index].m_state.m_U[0] - cell_center[index_nb[4]].m_state.m_U[0]);

            getLimitedSlope_Velocity_MAC_vd(icoords,COORD_X,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] - 1;
            getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_X,diff_nb);
            transverseD[0] += tmp*dz*(1.0 - m_dt*tmp)/2.0*(diff[2] - diff_nb[2]);
        }
        else //tmp <= 0.0
        {
            transverseD[0] += tmp*(-2.0*cell_center[index].m_state.m_U[0]);

            //u_z or v_z at (i,j,N_z + 1) equals to that at (i,j,N_z)
            transverseD[0] += 0.0;
        }
    }
    else assert(false);


///////////////////////////////////////////////////
//  transverseD terms on v_faces = u*v_x + w*v_z
///////////////////////////////////////////////////

    //u*v_x
    tmp = (cell_center[index].m_state.m_U[0] +
           cell_center[index_nb[0]].m_state.m_U[0] +
           cell_center[index_nb[9]].m_state.m_U[0] +
           cell_center[index_nb[3]].m_state.m_U[0])/4.0;
    tmp /= dx;
    if (tmp > 0.0)
    {
        transverseD[1] += tmp*(cell_center[index].m_state.m_U[1] - cell_center[index_nb[0]].m_state.m_U[1]);

        getLimitedSlope_Velocity_MAC_vd(icoords,COORD_Y,diff);
        ICoords[0] = icoords[0] - 1;
        ICoords[1] = icoords[1];
        ICoords[2] = icoords[2];
        getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_Y,diff_nb);
        transverseD[1] += tmp*dx*(1.0 - m_dt*tmp)/2.0*(diff[0] - diff_nb[0]);
    }
    else //tmp <= 0.0
    {
        transverseD[1] += tmp*(cell_center[index_nb[1]].m_state.m_U[1] - cell_center[index].m_state.m_U[1]);

        getLimitedSlope_Velocity_MAC_vd(icoords,COORD_Y,diff);
        ICoords[0] = icoords[0] + 1;
        ICoords[1] = icoords[1];
        ICoords[2] = icoords[2];
        getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_Y,diff_nb);
        transverseD[1] += tmp*dx*(1.0 + m_dt*tmp)/2.0*(diff[0] - diff_nb[0]);
    }

    //w*v_z
    if (bNoBoundary[4] && bNoBoundary[5])
    {
        tmp = (cell_center[index].m_state.m_U[2] +
               cell_center[index_nb[3]].m_state.m_U[2] +
               cell_center[index_nb[11]].m_state.m_U[2] +
               cell_center[index_nb[4]].m_state.m_U[2])/4.0;
        tmp /= dz;
        if (tmp > 0.0)
        {
            transverseD[1] += tmp*(cell_center[index].m_state.m_U[1] - cell_center[index_nb[4]].m_state.m_U[1]);

            getLimitedSlope_Velocity_MAC_vd(icoords,COORD_Y,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] - 1;
            getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_Y,diff_nb);
            transverseD[1] += tmp*dz*(1.0 - m_dt*tmp)/2.0*(diff[2] - diff_nb[2]);
        }
        else //tmp <= 0.0
        {
            transverseD[1] += tmp*(cell_center[index_nb[5]].m_state.m_U[1] - cell_center[index].m_state.m_U[1]);

            getLimitedSlope_Velocity_MAC_vd(icoords,COORD_Y,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] + 1;
            getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_Y,diff_nb);
            transverseD[1] += tmp*dz*(1.0 + m_dt*tmp)/2.0*(diff[2] - diff_nb[2]);
        }
    }
    else if (!bNoBoundary[4]) //cells on LOWER boundary
    {
        tmp = (cell_center[index].m_state.m_U[2] +
               cell_center[index_nb[3]].m_state.m_U[2] + 0.0 + 0.0)/4.0;
        tmp /= dz;
        if (tmp > 0.0)
        {
            transverseD[1] += tmp*(2.0*cell_center[index].m_state.m_U[1]);

            //u_z or v_z at (i,j,-1) equals to that at (i,j,0)
            transverseD[1] += 0.0;
        }
        else //tmp <= 0.0
        {
            transverseD[1] += tmp*(cell_center[index_nb[5]].m_state.m_U[1] - cell_center[index].m_state.m_U[1]);

            getLimitedSlope_Velocity_MAC_vd(icoords,COORD_Y,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] + 1;
            getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_Y,diff_nb);
            transverseD[1] += tmp*dz*(1.0 + m_dt*tmp)/2.0*(diff[2] - diff_nb[2]);
        }
    }
    else if (!bNoBoundary[5]) //cells on UPPER boundary
    {
        tmp = (cell_center[index].m_state.m_U[2] +
               cell_center[index_nb[3]].m_state.m_U[2] +
               cell_center[index_nb[11]].m_state.m_U[2] +
               cell_center[index_nb[4]].m_state.m_U[2])/4.0;
        tmp /= dz;
        if (tmp > 0.0)
        {
            transverseD[1] += tmp*(cell_center[index].m_state.m_U[1] - cell_center[index_nb[4]].m_state.m_U[1]);

            getLimitedSlope_Velocity_MAC_vd(icoords,COORD_Y,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2] - 1;
            getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_Y,diff_nb);
            transverseD[1] += tmp*dz*(1.0 - m_dt*tmp)/2.0*(diff[2] - diff_nb[2]);
        }
        else //tmp <= 0.0
        {
            transverseD[1] += tmp*(-2.0*cell_center[index].m_state.m_U[1]);

            //u_z or v_z at (i,j,N_z + 1) equals to that at (i,j,N_z)
            transverseD[1] += 0.0;
        }
    }
    else assert(false);


///////////////////////////////////////////////////
//  transverseD terms on w_faces = u*w_x + v*w_y
///////////////////////////////////////////////////

    //u*w_x
    if (!bNoBoundary[5]) //cells on UPPER boundary
    {
        tmp = 0.0; //u = 0 on UPPER bdry
        transverseD[2] += 0.0;
    }
    else //other cells
    {
        tmp = (cell_center[index].m_state.m_U[0] +
               cell_center[index_nb[0]].m_state.m_U[0] +
               cell_center[index_nb[17]].m_state.m_U[0] +
               cell_center[index_nb[5]].m_state.m_U[0])/4.0;
        tmp /= dx;
        if (tmp > 0.0)
        {
            transverseD[2] += tmp*(cell_center[index].m_state.m_U[2] - cell_center[index_nb[0]].m_state.m_U[2]);

            getLimitedSlope_Velocity_MAC_vd(icoords,COORD_Z,diff);
            ICoords[0] = icoords[0] - 1;
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2];
            getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_Z,diff_nb);
            transverseD[2] += tmp*dx*(1.0 - m_dt*tmp)/2.0*(diff[0] - diff_nb[0]);
        }
        else //tmp <= 0.0
        {
            transverseD[2] += tmp*(cell_center[index_nb[1]].m_state.m_U[2] - cell_center[index].m_state.m_U[2]);

            getLimitedSlope_Velocity_MAC_vd(icoords,COORD_Z,diff);
            ICoords[0] = icoords[0] + 1;
            ICoords[1] = icoords[1];
            ICoords[2] = icoords[2];
            getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_Z,diff_nb);
            transverseD[2] += tmp*dx*(1.0 + m_dt*tmp)/2.0*(diff[0] - diff_nb[0]);
        }
    }

    //v*w_y
    if (!bNoBoundary[5]) //UPPER boundary
    {
        tmp = 0.0; //v = 0 on UPPER bdry
        transverseD[2] += 0.0;
    }
    else //other cells
    {
        tmp = (cell_center[index].m_state.m_U[1] +
               cell_center[index_nb[2]].m_state.m_U[1] +
               cell_center[index_nb[13]].m_state.m_U[1] +
               cell_center[index_nb[5]].m_state.m_U[1])/4.0;
        tmp /= dy;
        if (tmp > 0.0)
        {
            transverseD[2] += tmp*(cell_center[index].m_state.m_U[2] - cell_center[index_nb[2]].m_state.m_U[2]);

            getLimitedSlope_Velocity_MAC_vd(icoords,COORD_Z,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1] - 1;
            ICoords[2] = icoords[2];
            getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_Z,diff_nb);
            transverseD[2] += tmp*dy*(1.0 - m_dt*tmp)/2.0*(diff[1] - diff_nb[1]);
        }
        else //tmp <= 0.0
        {
            transverseD[2] += tmp*(cell_center[index_nb[3]].m_state.m_U[2] - cell_center[index].m_state.m_U[2]);

            getLimitedSlope_Velocity_MAC_vd(icoords,COORD_Z,diff);
            ICoords[0] = icoords[0];
            ICoords[1] = icoords[1] + 1;
            ICoords[2] = icoords[2];
            getLimitedSlope_Velocity_MAC_vd(ICoords,COORD_Z,diff_nb);
            transverseD[2] += tmp*dy*(1.0 + m_dt*tmp)/2.0*(diff[1] - diff_nb[1]);
        }
    }
} //end getTransverseDTerm_Velocity_MAC_vd
*/


//Stability issue proposed by M. Monion, refer to JCP 123 pp.435-449(1996)???
void Incompress_Solver_Smooth_3D_Cartesian::getTransverseDTerm_Scalar_MAC_vd(
        int *icoords,
        EBM_COORD xyz,
        double transverseD[2],
        boolean bGhostCell)
{
    int i, j, k;
    int bNoBoundary[6], nb;
    int index,index_nb[6];
    double tmp, dx, dy, dz;
    COMPONENT comp;
    double crx_coords[MAXD];
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

    transverseD[0] = 0;
    transverseD[1] = 0;

    dx = top_h[0];
    dy = top_h[1];
    dz = top_h[2];

    i = icoords[0];
    j = icoords[1];
    k = icoords[2];
    index = d_index3d(i,j,k,top_gmax);
    comp = top_comp[index];

    //6 neighbours of the center cell
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    for (nb = 0; nb < 6; nb++)
    {
        checkBoundaryCondition(dir[nb],icoords,&bNoBoundary[nb],m_t_old,comp);
    }
    /*
    //4 directions
    bNoBoundary[0] = YES;
    bNoBoundary[1] = YES;
    bNoBoundary[2] = YES;
    bNoBoundary[3] = YES;

    //LOWER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[4] = NO; //cells on LOWER bdry
    else
        bNoBoundary[4] = YES;

    //UPPER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[5] = NO; //cells on UPPER bdry
    else
        bNoBoundary[5] = YES;
    */


    switch(xyz)
    {
    //transverseD terms in x-dirc = v*scalar_y + w*scalar_z
    case COORD_X:
        //v*rho_y & v*c_y
        tmp = (cell_center[index].m_state.m_U[1] +
               cell_center[index_nb[2]].m_state.m_U[1])/2.0;
        if (bNoBoundary[2] == 3 || bNoBoundary[3] == 3)
            dy = 0.5 * dy;
        tmp /= dy;
        if (tmp > 0)
        {
            if (bGhostCell && (top_comp[index] != top_comp[index_nb[2]]))
            {
                transverseD[0] += tmp*0.0;
                transverseD[1] += tmp*0.0;
            }
            else
            {
                transverseD[0] += tmp*(cell_center[index].m_state.m_rho - cell_center[index_nb[2]].m_state.m_rho);
                transverseD[1] += tmp*(cell_center[index].m_state.m_c - cell_center[index_nb[2]].m_state.m_c);
            }
        }
        else //tmp <= 0
        {
            if (bGhostCell && (top_comp[index] != top_comp[index_nb[3]]))
	    {
                transverseD[0] += tmp*0.0;
                transverseD[1] += tmp*0.0;
            }
            else
            {
                transverseD[0] += tmp*(cell_center[index_nb[3]].m_state.m_rho - cell_center[index].m_state.m_rho);
                transverseD[1] += tmp*(cell_center[index_nb[3]].m_state.m_c - cell_center[index].m_state.m_c);
            }
        }

        //w*rho_z & w*c_z
        if ((!bNoBoundary[4] && !bNoBoundary[5]) || bNoBoundary[4] == 3 || bNoBoundary[5] == 3)
        {
            tmp = (cell_center[index].m_state.m_U[2] +
                   cell_center[index_nb[4]].m_state.m_U[2])/2.0;
            if (bNoBoundary[4] == 3 || bNoBoundary[5] == 3)
                dz = 0.5 * dz;
            tmp /= dz;
            if (tmp > 0)
            {
                if (bGhostCell && (top_comp[index] != top_comp[index_nb[4]]))
		{
                    transverseD[0] += tmp*0.0;
                    transverseD[1] += tmp*0.0;
                }
                else
                {
                    transverseD[0] += tmp*(cell_center[index].m_state.m_rho - cell_center[index_nb[4]].m_state.m_rho);
                    transverseD[1] += tmp*(cell_center[index].m_state.m_c - cell_center[index_nb[4]].m_state.m_c);
                }
            }
            else //tmp <= 0
            {
                if (bGhostCell && (top_comp[index] != top_comp[index_nb[5]]))
                {
                    transverseD[0] += tmp*0.0;
                    transverseD[1] += tmp*0.0;
                }
                else
                {
                    transverseD[0] += tmp*(cell_center[index_nb[5]].m_state.m_rho - cell_center[index].m_state.m_rho);
                    transverseD[1] += tmp*(cell_center[index_nb[5]].m_state.m_c - cell_center[index].m_state.m_c);
                }
            }
        }
        else if (bNoBoundary[4]==2) //cells on LOWER boundary
        {
            tmp = (cell_center[index].m_state.m_U[2] + 0.0)/2.0;
            tmp /= dz;
            if (tmp > 0)
            {
                transverseD[0] += 0.0;
                transverseD[1] += 0.0;
            }
            else //tmp <= 0
            {
                if (bGhostCell && (top_comp[index] != top_comp[index_nb[5]]))
                {
                    transverseD[0] += tmp*0.0;
                    transverseD[1] += tmp*0.0;
                }
                else
                {
                    transverseD[0] += tmp*(cell_center[index_nb[5]].m_state.m_rho - cell_center[index].m_state.m_rho);
                    transverseD[1] += tmp*(cell_center[index_nb[5]].m_state.m_c - cell_center[index].m_state.m_c);
                }
            }
        }
        else if (bNoBoundary[5]==2) //cells on UPPER boundary
        {
            tmp = (0.0 + cell_center[index_nb[4]].m_state.m_U[2])/2.0;
            tmp /= dz;
            if (tmp > 0)
            {
                if (bGhostCell && (top_comp[index] != top_comp[index_nb[4]]))
                {
                    transverseD[0] += tmp*0.0;
                    transverseD[1] += tmp*0.0;
                }
                else
                {
                    transverseD[0] += tmp*(cell_center[index].m_state.m_rho - cell_center[index_nb[4]].m_state.m_rho);
                    transverseD[1] += tmp*(cell_center[index].m_state.m_c - cell_center[index_nb[4]].m_state.m_c);
                }
            }
            else //tmp <= 0
            {
                transverseD[0] += 0.0;
                transverseD[1] += 0.0;
            }
        }
        else assert(false);
        break;

    //transverseD terms in y-dirc = u*scalar_x + w*scalar_z
    case COORD_Y:
        //u*rho_x & u*c_x
        tmp = (cell_center[index].m_state.m_U[0] +
               cell_center[index_nb[0]].m_state.m_U[0])/2.0;
        if (bNoBoundary[0] == 3 || bNoBoundary[1] == 3)
            dx = 0.5 * dx;
        tmp /= dx;
        if (tmp > 0.0)
        {
            if (bGhostCell && (top_comp[index] != top_comp[index_nb[0]]))
            {
                transverseD[0] += tmp*0.0;
                transverseD[1] += tmp*0.0;
            }
            else
            {
                transverseD[0] += tmp*(cell_center[index].m_state.m_rho - cell_center[index_nb[0]].m_state.m_rho);
                transverseD[1] += tmp*(cell_center[index].m_state.m_c - cell_center[index_nb[0]].m_state.m_c);
            }
        }
        else //tmp <= 0.0
        {
            if (bGhostCell && (top_comp[index] != top_comp[index_nb[1]]))
            {
                transverseD[0] += tmp*0.0;
                transverseD[1] += tmp*0.0;
            }
            else
            {
                transverseD[0] += tmp*(cell_center[index_nb[1]].m_state.m_rho - cell_center[index].m_state.m_rho);
                transverseD[1] += tmp*(cell_center[index_nb[1]].m_state.m_c - cell_center[index].m_state.m_c);
            }
        }

        //w*rho_z & w*c_z
        if ((!bNoBoundary[4] && !bNoBoundary[5]) || bNoBoundary[4] == 3 || bNoBoundary[5] == 3)
        {
            tmp = (cell_center[index].m_state.m_U[2] +
                   cell_center[index_nb[4]].m_state.m_U[2])/2.0;
            if (bNoBoundary[4] == 3 || bNoBoundary[5] == 3)
                dz = 0.5 * dz;
            tmp /= dz;
            if (tmp > 0.0)
            {
                if (bGhostCell && (top_comp[index] != top_comp[index_nb[4]]))
                {
                    transverseD[0] += tmp*0.0;
                    transverseD[1] += tmp*0.0;
                }
                else
                {
                    transverseD[0] += tmp*(cell_center[index].m_state.m_rho - cell_center[index_nb[4]].m_state.m_rho);
                    transverseD[1] += tmp*(cell_center[index].m_state.m_c - cell_center[index_nb[4]].m_state.m_c);
                }
            }
            else //tmp <= 0.0
            {
                if (bGhostCell && (top_comp[index] != top_comp[index_nb[5]]))
                {
                    transverseD[0] += tmp*0.0;
                    transverseD[1] += tmp*0.0;
                }
                else
                {
                    transverseD[0] += tmp*(cell_center[index_nb[5]].m_state.m_rho - cell_center[index].m_state.m_rho);
                    transverseD[1] += tmp*(cell_center[index_nb[5]].m_state.m_c - cell_center[index].m_state.m_c);
                }
            }
        }
        else if (bNoBoundary[4]==2) //cells on LOWER boundary
        {
            tmp = (cell_center[index].m_state.m_U[2] + 0.0)/2.0;
            tmp /= dz;
            if (tmp > 0.0)
            {
                transverseD[0] += 0.0;
                transverseD[1] += 0.0;
            }
            else //tmp <= 0.0
            {
                if (bGhostCell && (top_comp[index] != top_comp[index_nb[5]]))
                {
                    transverseD[0] += tmp*0.0;
                    transverseD[1] += tmp*0.0;
                }
                else
                {
                    transverseD[0] += tmp*(cell_center[index_nb[5]].m_state.m_rho - cell_center[index].m_state.m_rho);
                    transverseD[1] += tmp*(cell_center[index_nb[5]].m_state.m_c - cell_center[index].m_state.m_c);
                }
            }
        }
        else if (bNoBoundary[5]==2) //cells on UPPER boundary
        {
            tmp = (0.0 + cell_center[index_nb[4]].m_state.m_U[2])/2.0;
            tmp /= dz;
            if (tmp > 0.0)
            {
                if (bGhostCell && (top_comp[index] != top_comp[index_nb[4]]))
                {
                    transverseD[0] += tmp*0.0;
                    transverseD[1] += tmp*0.0;
                }
                else
                {
                    transverseD[0] += tmp*(cell_center[index].m_state.m_rho - cell_center[index_nb[4]].m_state.m_rho);
                    transverseD[1] += tmp*(cell_center[index].m_state.m_c - cell_center[index_nb[4]].m_state.m_c);
                }
            }
            else //tmp <= 0.0
            {
                transverseD[0] += 0.0;
                transverseD[1] += 0.0;
            }
        }
        else assert(false);
        break;

    //transverseD terms on z-dirc = u*scalar_x + v*scalar_y
    case COORD_Z:
        //u*rho_x & u*c_x
        tmp = (cell_center[index].m_state.m_U[0] +
               cell_center[index_nb[0]].m_state.m_U[0])/2.0;
        if (bNoBoundary[0] == 3 || bNoBoundary[1] == 3)
            dx = 0.5 * dx;
        tmp /= dx;
        if (tmp > 0.0)
        {
            if (bGhostCell && (top_comp[index] != top_comp[index_nb[0]]))
            {
                transverseD[0] += tmp*0.0;
                transverseD[1] += tmp*0.0;
            }
            else
            {
                transverseD[0] += tmp*(cell_center[index].m_state.m_rho - cell_center[index_nb[0]].m_state.m_rho);
                transverseD[1] += tmp*(cell_center[index].m_state.m_c - cell_center[index_nb[0]].m_state.m_c);
            }
        }
        else //tmp <= 0.0
        {
            if (bGhostCell && (top_comp[index] != top_comp[index_nb[1]]))
            {
                transverseD[0] += tmp*0.0;
                transverseD[1] += tmp*0.0;
            }
            else
            {
                transverseD[0] += tmp*(cell_center[index_nb[1]].m_state.m_rho - cell_center[index].m_state.m_rho);
                transverseD[1] += tmp*(cell_center[index_nb[1]].m_state.m_c - cell_center[index].m_state.m_c);
            }
        }

        //v*rho_y & v*c_y
        tmp = (cell_center[index].m_state.m_U[1] +
               cell_center[index_nb[2]].m_state.m_U[1])/2.0;
        if (bNoBoundary[2] == 3 || bNoBoundary[3] == 3)
            dy = 0.5 * dy;
        tmp /= dy;
        if (tmp > 0.0)
        {
            if (bGhostCell && (top_comp[index] != top_comp[index_nb[2]]))
            {
                transverseD[0] += tmp*0.0;
                transverseD[1] += tmp*0.0;
            }
            else
            {
                transverseD[0] += tmp*(cell_center[index].m_state.m_rho - cell_center[index_nb[2]].m_state.m_rho);
                transverseD[1] += tmp*(cell_center[index].m_state.m_c - cell_center[index_nb[2]].m_state.m_c);
            }
        }
        else //tmp <= 0.0
        {
            if (bGhostCell && (top_comp[index] != top_comp[index_nb[3]]))
            {
                transverseD[0] += tmp*0.0;
                transverseD[1] += tmp*0.0;
            }
            else
            {
                transverseD[0] += tmp*(cell_center[index_nb[3]].m_state.m_rho - cell_center[index].m_state.m_rho);
                transverseD[1] += tmp*(cell_center[index_nb[3]].m_state.m_c - cell_center[index].m_state.m_c);
            }
        }
        break;

    default:
        assert(false);
    }
} /* end getTransverseDTerm_Scalar_MAC_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getDivU_coupled_vd(
        int *icoords,
        double diffusion[5],
        int flag)
{
    int index,index_nb[6];
    double Dcoef[6],Dcoef_edge[6],Dcoef0;
    double rho_edge[6];
    L_STATE statenb;
    double rho_nb[6],rho_center;
    int nb;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    bool bNoBoundary[6];
    double dh[3],dh0[3],dh1[3];

    diffusion[3] = 0.0;

    int i = icoords[0];
    int j = icoords[1];
    int k = icoords[2];

    index = d_index3d(i,j,k,top_gmax);
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    Dcoef0 = cell_center[index].m_state.m_Dcoef;
    if (flag == 0 || flag == 2)
        rho_center = cell_center[index].m_state.m_rho;
    if (flag == 1)
        rho_center = 0.5*(cell_center[index].m_state.m_rho +
                     cell_center[index].m_state.m_rho_old);
    if (flag == 3)
        rho_center = cell_center[index].m_state.m_rho_old;

    for (nb = 0; nb < 6; nb++)
    {
        if (flag == 0 || flag == 3)
            bNoBoundary[nb] = getNeighborOrBoundaryState_vd(icoords,dir[nb],statenb,m_t_old);
        if (flag == 1)
            bNoBoundary[nb] = getNeighborOrBoundaryState_vd(icoords,dir[nb],statenb,m_t_int);
        if (flag == 2)
            bNoBoundary[nb] = getNeighborOrBoundaryState_vd(icoords,dir[nb],statenb,m_t_new);

        // Homogeneous Neumann B.C. for rho
        if (flag == 0 || flag == 2)
            rho_nb[nb] = statenb.m_rho;
        if (flag == 1)
            rho_nb[nb] = 0.5*(statenb.m_rho + statenb.m_rho_old);
        if (flag == 3)
            rho_nb[nb] = statenb.m_rho_old;

        if(!bNoBoundary[nb])
        {
            Dcoef[nb] = Dcoef0;
            Dcoef_edge[nb] = Dcoef0;
            rho_edge[nb] = rho_nb[nb];
        }
        else
        {
            Dcoef[nb] = cell_center[index_nb[nb]].m_state.m_Dcoef;
            Dcoef_edge[nb] = 0.5*(cell_center[index_nb[nb]].m_state.m_Dcoef + Dcoef0);
            rho_edge[nb] = 2.0/(1.0/rho_nb[nb] + 1.0/rho_center);
        }
    }

    dh[0] = top_h[0];
    dh[1] = top_h[1];
    dh[2] = top_h[2];

    if (bNoBoundary[0])
        dh0[0] = top_h[0];
    else
        dh0[0] = top_h[0]/2.0;
    if (bNoBoundary[1])
        dh1[0] = top_h[0];
    else
        dh1[0] = top_h[0]/2.0;

    if (bNoBoundary[2])
        dh0[1] = top_h[1];
    else
        dh0[1] = top_h[1]/2.0;
    if (bNoBoundary[3])
        dh1[1] = top_h[1];
    else
        dh1[1] = top_h[1]/2.0;

    if (bNoBoundary[4])
        dh0[2] = top_h[2];
    else
        dh0[2] = top_h[2]/2.0;
    if (bNoBoundary[5])
        dh1[2] = top_h[2];
    else
        dh1[2] = top_h[2]/2.0;

    // -div(Dcoef/rho*grad(rho))
    diffusion[3] -= (Dcoef_edge[1]/rho_edge[1]*(rho_nb[1]-rho_center)/dh1[0] - Dcoef_edge[0]/rho_edge[0]*(rho_center-rho_nb[0])/dh0[0])/dh[0];
                                                  // -(Dcoef/rho*rho_x)_x
    diffusion[3] -= (Dcoef_edge[3]/rho_edge[3]*(rho_nb[3]-rho_center)/dh1[1] - Dcoef_edge[2]/rho_edge[2]*(rho_center-rho_nb[2])/dh0[1])/dh[1];
                                                  // -(Dcoef/rho*rho_y)_y
    diffusion[3] -= (Dcoef_edge[5]/rho_edge[5]*(rho_nb[5]-rho_center)/dh1[2] - Dcoef_edge[4]/rho_edge[4]*(rho_center-rho_nb[4])/dh0[2])/dh[2];
                                                  // -(Dcoef/rho*rho_z)_z
} /* end getDivU_coupled_vd */


//valid for variable diffusity
void Incompress_Solver_Smooth_3D_Cartesian::getDivU_MAC_vd(
        int *icoords,
        double *diff,
        int flag,
        boolean bGhostCell)
{
    int nb,index,index_nb[6];
    L_STATE statenb;
    double t;
    double rho,rho_nb[6],rho_face[6];
    double Dcoef,Dcoef_t,Dcoef_face[6];
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    int bNoBoundary[6];
    double dh[3],dh0[3],dh1[3];
    int i = icoords[0];
    int j = icoords[1];
    int k = icoords[2];
    boolean useSGSCellCenter = YES;

    //use no ghost cells for diff calc
    bGhostCell = NO;

    index = d_index3d(i,j,k,top_gmax);
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    //Dcoef_center & rho_center
    Dcoef = cell_center[index].m_state.m_Dcoef;
    if (flag==0 || flag==2)
        rho = cell_center[index].m_state.m_rho;
    else if (flag==1)
        rho = (cell_center[index].m_state.m_rho + cell_center[index].m_state.m_rho_old)/2;
    else if (flag==3)
        rho = cell_center[index].m_state.m_rho_old;
    else assert(false);

    for (nb=0; nb<6; ++nb)
    {
        if (flag==0 || flag==3)		t = m_t_old;
        else if (flag==1)		t = m_t_int;
        else if (flag==2)		t = m_t_new;
        else				assert(false);

        /*
         * removal tag: HAOZ
         *
         * REFLECTION BOUNDARY CONDITION
         *
         *
         * */
        if (!bGhostCell)
            bNoBoundary[nb] = getNeighborOrBoundaryScalar_MAC_vd(icoords,dir[nb],statenb,t);

        //homogeneous Neumann B.C. for rho and Dcoef
        if (flag==0 || flag==2)
            rho_nb[nb] = statenb.m_rho;
        else if (flag==1)
            rho_nb[nb] = (statenb.m_rho + statenb.m_rho_old)/2;
        else if (flag==3)
            rho_nb[nb] = statenb.m_rho_old;
        else assert(false);

        //TODO: implement Dcoef_face[] for Dcoef_t on 3 cell-faces
        //Dcoef_t at m_t_old, which is a first-order approx. in time
        if (bNoBoundary[nb]==2) {
            rho_face[nb] = rho_nb[nb];
            if (useSGSCellCenter) {
                Dcoef_face[nb] = statenb.m_Dcoef + statenb.m_Dcoef_turbulent[3];
            }
            else {
                printf("codes needed for Dcoef_turbulent on cell-faces.\n");
                assert(false);
            }
        }
        else {
            rho_face[nb] = 2.0/(1.0/rho_nb[nb] + 1.0/rho);
            if (useSGSCellCenter) {
                Dcoef_t = cell_center[index].m_state.m_Dcoef_turbulent[3];
                Dcoef_face[nb] = (statenb.m_Dcoef + Dcoef +
                                  statenb.m_Dcoef_turbulent[3] + Dcoef_t)/2;
            }
            else {
                printf("codes needed for Dcoef_turbulent on cell-faces.\n");
                assert(false);
            }
        }
    }

    dh[0] = top_h[0];
    dh[1] = top_h[1];
    dh[2] = top_h[2];

    for (nb = 0; nb < 6; nb++) // no DIRICHLET here
        if (bNoBoundary[nb] == 1)
            clean_up(ERROR);

    if (bNoBoundary[0] < 1)
        dh0[0] = top_h[0];
    else
        dh0[0] = top_h[0]/2;
    if (bNoBoundary[1] < 1)
        dh1[0] = top_h[0];
    else
        dh1[0] = top_h[0]/2;

    if (bNoBoundary[2] < 1)
        dh0[1] = top_h[1];
    else
        dh0[1] = top_h[1]/2;
    if (bNoBoundary[3] < 1)
        dh1[1] = top_h[1];
    else
        dh1[1] = top_h[1]/2;

    if (bNoBoundary[4] < 1)
        dh0[2] = top_h[2];
    else
        dh0[2] = top_h[2]/2;
    if (bNoBoundary[5] < 1)
        dh1[2] = top_h[2];
    else
        dh1[2] = top_h[2]/2;

/***** divergence constraint -div(Dcoef/rho*grad(rho)) on cell centers *****/
    *diff = 0;
    //-(Dcoef/rho*rho_x)_x
    *diff -= (Dcoef_face[1]/rho_face[1]*(rho_nb[1]-rho)/dh1[0] -
              Dcoef_face[0]/rho_face[0]*(rho-rho_nb[0])/dh0[0])/dh[0];
    //-(Dcoef/rho*rho_y)_y
    *diff -= (Dcoef_face[3]/rho_face[3]*(rho_nb[3]-rho)/dh1[1] -
              Dcoef_face[2]/rho_face[2]*(rho-rho_nb[2])/dh0[1])/dh[1];
    //-(Dcoef/rho*rho_z)_z
    *diff -= (Dcoef_face[5]/rho_face[5]*(rho_nb[5]-rho)/dh1[2] -
              Dcoef_face[4]/rho_face[4]*(rho-rho_nb[4])/dh0[2])/dh[2];
} /* end getDivU_MAC_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getDiffusionC_coupled_vd(
        int *icoords,
        double diffusion[5])
{
    int index,index_nb[6];
    double Dcoef[6],Dcoef_edge[6],Dcoef0;
    double rho_edge[6],c_edge[6];
    L_STATE statenb;
    double rho_nb[6],c_nb[6],rho_center,c_center;
    int nb;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    bool bNoBoundary[6];
    double dh[3],dh0[3],dh1[3];

    diffusion[4] = 0;

    int i = icoords[0];
    int j = icoords[1];
    int k = icoords[2];
    index = d_index3d(i,j,k,top_gmax);
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    Dcoef0 = cell_center[index].m_state.m_Dcoef;
    c_center = cell_center[index].m_state.m_c;
    rho_center = cell_center[index].m_state.m_rho;

    for (nb = 0; nb < 6; nb++)
    {
        bNoBoundary[nb] = getNeighborOrBoundaryState_vd(icoords,dir[nb],statenb,m_t_old);
        // homogeneous Neumann B.C. for c and rho
        c_nb[nb] = statenb.m_c;
        rho_nb[nb] = statenb.m_rho;

        if(!bNoBoundary[nb])
        {
            Dcoef[nb] = Dcoef0;
            Dcoef_edge[nb] = Dcoef0;
            rho_edge[nb] = rho_nb[nb];
        }
        else
        {
            Dcoef[nb] = cell_center[index_nb[nb]].m_state.m_Dcoef;
            Dcoef_edge[nb] = 0.5*(cell_center[index_nb[nb]].m_state.m_Dcoef + Dcoef0);
            rho_edge[nb] = 2.0/(1.0/rho_nb[nb] + 1.0/rho_center);
        }
    }

    dh[0] = top_h[0];
    dh[1] = top_h[1];
    dh[2] = top_h[2];

    if (bNoBoundary[0])
        dh0[0] = top_h[0];
    else
        dh0[0] = top_h[0]/2.0;
    if (bNoBoundary[1])
        dh1[0] = top_h[0];
    else
        dh1[0] = top_h[0]/2.0;

    if (bNoBoundary[2])
        dh0[1] = top_h[1];
    else
        dh0[1] = top_h[1]/2.0;
    if (bNoBoundary[3])
        dh1[1] = top_h[1];
    else
        dh1[1] = top_h[1]/2.0;

    if (bNoBoundary[4])
        dh0[2] = top_h[2];
    else
        dh0[2] = top_h[2]/2.0;
    if (bNoBoundary[5])
        dh1[2] = top_h[2];
    else
        dh1[2] = top_h[2]/2.0;

    // div(Dcoef*rho*grad(c))
    diffusion[4] += (Dcoef_edge[1]*rho_edge[1]*(c_nb[1]-c_center)/dh1[0] - Dcoef_edge[0]*rho_edge[0]*(c_center-c_nb[0])/dh0[0])/dh[0];
                                                     // (Dcoef*rho*c_x)_x
    diffusion[4] += (Dcoef_edge[3]*rho_edge[3]*(c_nb[3]-c_center)/dh1[1] - Dcoef_edge[2]*rho_edge[2]*(c_center-c_nb[2])/dh0[1])/dh[1];
                                                     // (Dcoef*rho*c_y)_y
    diffusion[4] += (Dcoef_edge[5]*rho_edge[5]*(c_nb[5]-c_center)/dh1[2] - Dcoef_edge[4]*rho_edge[4]*(c_center-c_nb[4])/dh0[2])/dh[2];
                                                     // (Dcoef*rho*c_z)_z
} /* end getDiffusionC_coupled_vd */


//valid for variable diffusity
void Incompress_Solver_Smooth_3D_Cartesian::getDiffusionC_MAC_vd(
        int *icoords,
        double *diff,
        boolean bGhostCell)
{
    int nb,index,index_nb[6];
    double rho,rho_nb[6],rho_face[6];
    double Dcoef,Dcoef_t,Dcoef_face[6];
    double c,c_nb[6];
    L_STATE statenb;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    int bNoBoundary[6];
    double dh[3],dh0[3],dh1[3];
    int i = icoords[0];
    int j = icoords[1];
    int k = icoords[2];
    boolean useSGSCellCenter = YES;

    //use no ghost cells for diff calc
    bGhostCell = NO;

    index = d_index3d(i,j,k,top_gmax);
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    c = cell_center[index].m_state.m_c;
    rho = cell_center[index].m_state.m_rho;
    Dcoef = cell_center[index].m_state.m_Dcoef;

    for (nb=0; nb<6; ++nb)
    {
        if (!bGhostCell)
            bNoBoundary[nb] = getNeighborOrBoundaryScalar_MAC_vd(icoords,dir[nb],statenb,m_t_old);

        //homogeneous Neumann B.C. for c, rho, and Dcoef
        c_nb[nb] = statenb.m_c;
        rho_nb[nb] = statenb.m_rho;

        //TODO: implement Dcoef_face[] for Dcoef_t on 3 cell-faces
        if (bNoBoundary[nb]==2) { //cells on LOWER/UPPER bdry
            rho_face[nb] = rho_nb[nb];
            if (useSGSCellCenter) {
                Dcoef_face[nb] = statenb.m_Dcoef + statenb.m_Dcoef_turbulent[3];
            }
            else {
                printf("codes needed for Dcoef_turbulent on cell-faces.\n");
                assert(false);
            }
        }
        else { //interior cells
            rho_face[nb] = 2.0/(1.0/rho_nb[nb] + 1.0/rho);
            if (useSGSCellCenter) {
                Dcoef_t = cell_center[index].m_state.m_Dcoef_turbulent[3];
                Dcoef_face[nb] = (statenb.m_Dcoef + Dcoef +
                                  statenb.m_Dcoef_turbulent[3] + Dcoef_t)/2;
            }
            else {
                printf("codes needed for Dcoef_turbulent on cell-faces.\n");
                assert(false);
            }
        }
    }

    dh[0] = top_h[0];
    dh[1] = top_h[1];
    dh[2] = top_h[2];

    for (nb = 0; nb < 6; nb++)
        if (bNoBoundary[nb] == 1) // NO DIRICHLET
            clean_up(ERROR);

    if (bNoBoundary[0]<1)
        dh0[0] = top_h[0];
    else
        dh0[0] = top_h[0]/2.0;
    if (bNoBoundary[1]<1)
        dh1[0] = top_h[0];
    else
        dh1[0] = top_h[0]/2.0;

    if (bNoBoundary[2]<1)
        dh0[1] = top_h[1];
    else
        dh0[1] = top_h[1]/2.0;
    if (bNoBoundary[3]<1)
        dh1[1] = top_h[1];
    else
        dh1[1] = top_h[1]/2.0;

    if (bNoBoundary[4]<1)
        dh0[2] = top_h[2];
    else
        dh0[2] = top_h[2]/2.0;
    if (bNoBoundary[5]<1)
        dh1[2] = top_h[2];
    else
        dh1[2] = top_h[2]/2.0;

//// diffusion term of concentration eqn is: div(Dcoef*rho*grad(c)) ////
    *diff = 0;
    //(Dcoef*rho*c_x)_x
    *diff += (Dcoef_face[1]*rho_face[1]*(c_nb[1]-c)/dh1[0] -
              Dcoef_face[0]*rho_face[0]*(c-c_nb[0])/dh0[0])/dh[0];
    //(Dcoef*rho*c_y)_y
    *diff += (Dcoef_face[3]*rho_face[3]*(c_nb[3]-c)/dh1[1] -
              Dcoef_face[2]*rho_face[2]*(c-c_nb[2])/dh0[1])/dh[1];
    //(Dcoef*rho*c_z)_z
    *diff += (Dcoef_face[5]*rho_face[5]*(c_nb[5]-c)/dh1[2] -
              Dcoef_face[4]*rho_face[4]*(c-c_nb[4])/dh0[2])/dh[2];
} /* end getDiffusionC_MAC_vd */


// Minmod slope limiter
void Incompress_Solver_Smooth_3D_Cartesian::getLimitedSlope_vd(
        int *icoords,
        EBM_COORD xyz,
        double slope[5])
{
    double dh0, dh1;
    L_STATE U0, U1, U2;

    U1 = cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_state;

    bool bNoBoundary[2];

    if(xyz==COORD_X)
    {
        bNoBoundary[0] = getNeighborOrBoundaryState_vd(icoords,WEST,U0,m_t_old);
        if(bNoBoundary[0])
            dh0 = top_h[0];
        else
            dh0 = top_h[0]/2;

        bNoBoundary[1] = getNeighborOrBoundaryState_vd(icoords,EAST,U2,m_t_old);
        if(bNoBoundary[1])
            dh1 = top_h[0];
        else
            dh1 = top_h[0]/2;
    }
    else if(xyz==COORD_Y)
    {
        bNoBoundary[0] = getNeighborOrBoundaryState_vd(icoords,SOUTH,U0,m_t_old);
        if(bNoBoundary[0])
            dh0 = top_h[1];
        else
            dh0 = top_h[1]/2;

        bNoBoundary[1] = getNeighborOrBoundaryState_vd(icoords,NORTH,U2,m_t_old);
        if(bNoBoundary[1])
            dh1 = top_h[1];
        else
            dh1 = top_h[1]/2;
    }
    else //xyz == COORD_Z
    {
        bNoBoundary[0] = getNeighborOrBoundaryState_vd(icoords,LOWER,U0,m_t_old);
        if(bNoBoundary[0])
            dh0 = top_h[2];
        else
            dh0 = top_h[2]/2;

        bNoBoundary[1] = getNeighborOrBoundaryState_vd(icoords,UPPER,U2,m_t_old);
        if(bNoBoundary[1])
            dh1 = top_h[2];
        else
            dh1 = top_h[2]/2;
    }
    slope[0] = EBM_minmod((U1.m_U[0]-U0.m_U[0])/dh0, (U2.m_U[0]-U1.m_U[0])/dh1);
    slope[1] = EBM_minmod((U1.m_U[1]-U0.m_U[1])/dh0, (U2.m_U[1]-U1.m_U[1])/dh1);
    slope[2] = EBM_minmod((U1.m_U[2]-U0.m_U[2])/dh0, (U2.m_U[2]-U1.m_U[2])/dh1);
    slope[3] = EBM_minmod((U1.m_rho -U0.m_rho) /dh0, (U2.m_rho -U1.m_rho) /dh1);
    slope[4] = EBM_minmod((U1.m_c   -U0.m_c)   /dh0, (U2.m_c   -U1.m_c)   /dh1);
} /* end getLimitedSlope_vd */


// Minmod slope limiter
void Incompress_Solver_Smooth_3D_Cartesian::getLimitedSlope_Scalar_MAC_vd(
        int *icoords,
        EBM_COORD xyz,
        double slope[2],
        boolean bGhostCell)
{
    double dh0, dh1;
    L_STATE U0, U1, U2;

    U1 = cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_state;

    int bNoBoundary[2];

    if (xyz==COORD_X)
    {
        if (!bGhostCell)
            bNoBoundary[0] = getNeighborOrBoundaryScalar_MAC_vd(icoords,WEST,U0,m_t_old);
        else
            bNoBoundary[0] = getNeighborOrBoundaryScalar_MAC_GhostCell_vd(icoords,WEST,U0,m_t_old);//need to change
        if(bNoBoundary[0]<1)
            dh0 = top_h[0];
        else
            dh0 = top_h[0]/2;

        if (!bGhostCell)
            bNoBoundary[1] = getNeighborOrBoundaryScalar_MAC_vd(icoords,EAST,U2,m_t_old);
        else
            bNoBoundary[1] = getNeighborOrBoundaryScalar_MAC_GhostCell_vd(icoords,EAST,U2,m_t_old);
        if(bNoBoundary[1]<1)
            dh1 = top_h[0];
        else
            dh1 = top_h[0]/2;
    }
    else if (xyz==COORD_Y)
    {
        if (!bGhostCell)
            bNoBoundary[0] = getNeighborOrBoundaryScalar_MAC_vd(icoords,SOUTH,U0,m_t_old);
        else
            bNoBoundary[0] = getNeighborOrBoundaryScalar_MAC_GhostCell_vd(icoords,SOUTH,U0,m_t_old);
        if(bNoBoundary[0]<1)
            dh0 = top_h[1];
        else
            dh0 = top_h[1]/2;

        if (!bGhostCell)
            bNoBoundary[1] = getNeighborOrBoundaryScalar_MAC_vd(icoords,NORTH,U2,m_t_old);
        else
            bNoBoundary[1] = getNeighborOrBoundaryScalar_MAC_GhostCell_vd(icoords,NORTH,U2,m_t_old);
        if(bNoBoundary[1]<1)
            dh1 = top_h[1];
        else
            dh1 = top_h[1]/2;
    }
    else if (xyz==COORD_Z)
    {
        if (!bGhostCell)
            bNoBoundary[0] = getNeighborOrBoundaryScalar_MAC_vd(icoords,LOWER,U0,m_t_old);
        else
            bNoBoundary[0] = getNeighborOrBoundaryScalar_MAC_GhostCell_vd(icoords,LOWER,U0,m_t_old);
        if(bNoBoundary[0]<1)
            dh0 = top_h[2];
        else
            dh0 = top_h[2]/2;

        if (!bGhostCell)
            bNoBoundary[1] = getNeighborOrBoundaryScalar_MAC_vd(icoords,UPPER,U2,m_t_old);
        else
            bNoBoundary[1] = getNeighborOrBoundaryScalar_MAC_GhostCell_vd(icoords,UPPER,U2,m_t_old);
        if(bNoBoundary[1]<1)
            dh1 = top_h[2];
        else
            dh1 = top_h[2]/2;
    }
    else assert(false);

    slope[0] = EBM_minmod((U1.m_rho-U0.m_rho)/dh0, (U2.m_rho-U1.m_rho)/dh1);
    slope[1] = EBM_minmod((U1.m_c-U0.m_c)/dh0, (U2.m_c-U1.m_c)/dh1);
} /* end getLimitedSlope_Scalar_MAC_vd */


//Minmod slope limiter
void Incompress_Solver_Smooth_3D_Cartesian::getLimitedSlope_Velocity_MAC_vd(
        int *icoords,
        EBM_COORD xyz,
        double slope[3])
{
    int i, j, k;
    int bNoBoundary[6];//type change
    int index,index_nb[6];
    COMPONENT comp;
    double dx, dy, dz;
    L_STATE U0, U1, U2;
    double crx_coords[MAXD];
    double onWallvel; // if on wall. average interior and exterior

    int nb;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

    dx = top_h[0];
    dy = top_h[1];
    dz = top_h[2];

    i = icoords[0];
    j = icoords[1];
    k = icoords[2];
    index = d_index3d(i,j,k,top_gmax);
    comp = top_comp[index];

    //6 neighbours of the center cell
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    //TODO && FIXME: remove
    /*
    // 4 directions
    bNoBoundary[0] = YES;
    bNoBoundary[1] = YES;
    bNoBoundary[2] = YES;
    bNoBoundary[3] = YES;
    // removal tag: HAOZ
    // wrt REFLECTION
    // LOWER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[4] = NO;
    else
        bNoBoundary[4] = YES;
    // UPPER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[5] = NO;
    else
        bNoBoundary[5] = YES;
    */
    for (nb = 0; nb < 6; nb++)
    {
        checkBoundaryCondition(dir[nb],icoords,&bNoBoundary[nb],m_t_int,comp);
    }

    //Limited slopes on u-faces, v-faces, or w-faces
    U1 = cell_center[index].m_state;

    //u-faces
    if (xyz==COORD_X)
    {
        //u_x
        if ((!bNoBoundary[0] && !bNoBoundary[1]) || bNoBoundary[0] == 3 || bNoBoundary[1] == 3) // PERIODIC, INTERIOR, REFLECT
        {
            U0 = cell_center[index_nb[0]].m_state;
            U2 = cell_center[index_nb[1]].m_state;
            slope[0] = EBM_minmod((U1.m_U[0]-U0.m_U[0])/dx, (U2.m_U[0]-U1.m_U[0])/dx);
            //TODO && FIXME: check the velocity normal to the wall is zero
            /*
            if (bNoBoundary[0] == 3)
                printf("U0.m_U[0] = %24.24f\n", U0.m_U[0]);
            if (bNoBoundary[1] == 3)
                printf("U1.m_U[0] = %24.24f\n", U1.m_U[0]);
            */
        }
        else assert(false);

        //u_y
        if (!bNoBoundary[2] && !bNoBoundary[3])
        {
            U0 = cell_center[index_nb[2]].m_state;
            U2 = cell_center[index_nb[3]].m_state;
            slope[1] = EBM_minmod((U1.m_U[0]-U0.m_U[0])/dy, (U2.m_U[0]-U1.m_U[0])/dy);
        }
        else if (bNoBoundary[2] == 3)
        {
            //printf("WEST REFLECTION\n");
            U0 = cell_center[index_nb[2]].m_state;
            U2 = cell_center[index_nb[3]].m_state;
            onWallvel = 0.5 * (U0.m_U[0] + U1.m_U[0]);
            slope[1] = EBM_minmod(2.0*(U1.m_U[0]-onWallvel)/dy, (U2.m_U[0]-U1.m_U[0])/dy);
        }
        else if (bNoBoundary[3] == 3)
        {
            //printf("EAST REFLECTION\n");
            U0 = cell_center[index_nb[2]].m_state;
            U2 = cell_center[index_nb[3]].m_state;
            onWallvel = 0.5 * (U2.m_U[0] + U1.m_U[0]);
            slope[1] = EBM_minmod(2.0*(onWallvel-U1.m_U[0])/dy, (U1.m_U[0]-U0.m_U[0])/dy);
        }
        else assert(false);

        //u_z
        if (!bNoBoundary[4] && !bNoBoundary[5])
        {
            U0 = cell_center[index_nb[4]].m_state;
            U2 = cell_center[index_nb[5]].m_state;
            slope[2] = EBM_minmod((U1.m_U[0]-U0.m_U[0])/dz, (U2.m_U[0]-U1.m_U[0])/dz);
        }
        else if (bNoBoundary[4] == 2) //cells on LOWER boundary
        {
            U2 = cell_center[index_nb[5]].m_state;
            slope[2] = EBM_minmod((2.0*U1.m_U[0])/dz, (U2.m_U[0]-U1.m_U[0])/dz);
        }
        else if (bNoBoundary[5] == 2) //cells on UPPER boundary
        {
            U0 = cell_center[index_nb[4]].m_state;
            slope[2] = EBM_minmod((U1.m_U[0]-U0.m_U[0])/dz, (-2.0*U1.m_U[0])/dz);
        }
        else if (bNoBoundary[4] == 3) // reflect
        {
            U0 = cell_center[index_nb[4]].m_state;
            U2 = cell_center[index_nb[5]].m_state;
            onWallvel = 0.5 * (U0.m_U[0] + U1.m_U[0]);
            slope[2] = EBM_minmod(2.0*(U1.m_U[0]-onWallvel)/dz, (U2.m_U[0]-U1.m_U[0])/dz);
        }
        else if (bNoBoundary[5] == 3) // reflect
        {
            U0 = cell_center[index_nb[4]].m_state;
            U2 = cell_center[index_nb[5]].m_state;
            onWallvel = 0.5 * (U1.m_U[0] + U2.m_U[0]);
            slope[2] = EBM_minmod(2.0*(onWallvel-U1.m_U[0])/dz, (U1.m_U[0]-U0.m_U[0])/dz);
        }
        else assert(false);
    }

    //v-faces
    else if (xyz==COORD_Y)
    {
        //v_x
        if (!bNoBoundary[2] && !bNoBoundary[3])
        {
            U0 = cell_center[index_nb[0]].m_state;
            U2 = cell_center[index_nb[1]].m_state;
            slope[0] = EBM_minmod((U1.m_U[1]-U0.m_U[1])/dx, (U2.m_U[1]-U1.m_U[1])/dx);
        }
        else if (bNoBoundary[2] == 3)
        {
            U0 = cell_center[index_nb[0]].m_state;
            U2 = cell_center[index_nb[1]].m_state;
            onWallvel = 0.5 * (U0.m_U[1] + U1.m_U[1]);
            slope[0] = EBM_minmod(2.0*(U1.m_U[1]-onWallvel)/dx, (U2.m_U[1]-U1.m_U[1])/dx);
        }
        else if (bNoBoundary[3] == 3)
        {
            U0 = cell_center[index_nb[0]].m_state;
            U2 = cell_center[index_nb[1]].m_state;
            onWallvel = 0.5 * (U1.m_U[1] + U2.m_U[1]);
            slope[0] = EBM_minmod(2.0*(onWallvel-U1.m_U[1])/dx, (U1.m_U[1]-U0.m_U[1])/dx);
        }
        else assert(false);

        //v_y
        if ((!bNoBoundary[2] && !bNoBoundary[3]) || bNoBoundary[2] == 3 || bNoBoundary[3] == 3) // PERIODIC, INTERIOR, REFLECT
        {
            U0 = cell_center[index_nb[2]].m_state;
            U2 = cell_center[index_nb[3]].m_state;
            slope[1] = EBM_minmod((U1.m_U[1]-U0.m_U[1])/dy, (U2.m_U[1]-U1.m_U[1])/dy);
            //TODO && FIXME: remove printout should be zero
            /*
            if (bNoBoundary[2] == 3)
                printf("U0.m_U[1] = %24.24f\n", U0.m_U[1]);
            if (bNoBoundary[3] == 3)
                printf("U1.m_U[1] = %24.24f\n", U1.m_U[1]);
            */
        }
        else assert(false);

        //v_z
        if (!bNoBoundary[4] && !bNoBoundary[5])
        {
            U0 = cell_center[index_nb[4]].m_state;
            U2 = cell_center[index_nb[5]].m_state;
            slope[2] = EBM_minmod((U1.m_U[1]-U0.m_U[1])/dz, (U2.m_U[1]-U1.m_U[1])/dz);
        }
        else if (bNoBoundary[4] == 2) //cells on LOWER boundary
        {
            U2 = cell_center[index_nb[5]].m_state;
            slope[2] = EBM_minmod((2.0*U1.m_U[1])/dz, (U2.m_U[1]-U1.m_U[1])/dz);
        }
        else if (bNoBoundary[5] == 2) //cells on UPPER boundary
        {
            U0 = cell_center[index_nb[4]].m_state;
            slope[2] = EBM_minmod((U1.m_U[1]-U0.m_U[1])/dz, (-2.0*U1.m_U[1])/dz);
        }
        else if (bNoBoundary[4] == 3)
        {
            U0 = cell_center[index_nb[4]].m_state;
            U2 = cell_center[index_nb[5]].m_state;
            onWallvel = 0.5 * (U0.m_U[1] + U1.m_U[1]);
            slope[2] = EBM_minmod(2.0*(U1.m_U[1]-onWallvel)/dz, (U2.m_U[1]-U1.m_U[1])/dz);
        }
        else if (bNoBoundary[5] == 3)
        {
            U0 = cell_center[index_nb[4]].m_state;
            U2 = cell_center[index_nb[5]].m_state;
            onWallvel = 0.5 * (U1.m_U[1] + U2.m_U[1]);
            slope[2] = EBM_minmod(2.0*(onWallvel-U1.m_U[1])/dz, (U1.m_U[1]-U0.m_U[1])/dz);
        }
        else assert(false);
    }

    //w-faces
    else if (xyz==COORD_Z)
    {
        //w_x
        if (!bNoBoundary[0] && !bNoBoundary[1])
        {
            U0 = cell_center[index_nb[0]].m_state;
            U2 = cell_center[index_nb[1]].m_state;
            slope[0] = EBM_minmod((U1.m_U[2]-U0.m_U[2])/dx, (U2.m_U[2]-U1.m_U[2])/dx);
        }
        else if (bNoBoundary[0]==3) // WEST REFLECT
        {
            U0 = cell_center[index_nb[0]].m_state;
            U2 = cell_center[index_nb[1]].m_state;
            onWallvel = 0.5 * (U0.m_U[2]+U1.m_U[2]);
            slope[0] = EBM_minmod(2.0*(U1.m_U[2]-onWallvel)/dz, (U2.m_U[2]-U1.m_U[2])/dz);
        }
        else if (bNoBoundary[1]==3) // EAST REFLECT
        {
            U0 = cell_center[index_nb[0]].m_state;
            U2 = cell_center[index_nb[1]].m_state;
            onWallvel = 0.5 * (U1.m_U[2]+U2.m_U[2]);
            slope[0] = EBM_minmod(2.0*(onWallvel-U1.m_U[2])/dz, (U1.m_U[2]-U0.m_U[2])/dz);
        }
        else assert(false);

        //w_y
        if (!bNoBoundary[2] && !bNoBoundary[3])
        {
            U0 = cell_center[index_nb[2]].m_state;
            U2 = cell_center[index_nb[3]].m_state;
            slope[1] = EBM_minmod((U1.m_U[2]-U0.m_U[2])/dy, (U2.m_U[2]-U1.m_U[2])/dy);
        }
        else if (bNoBoundary[2]==3) // SOUTH REFLECT
        {
            U0 = cell_center[index_nb[2]].m_state;
            U2 = cell_center[index_nb[3]].m_state;
            onWallvel = 0.5 *(U0.m_U[2] + U1.m_U[2]);
            slope[1] = EBM_minmod(2.0*(U1.m_U[2]-onWallvel)/dy, (U2.m_U[2]-U1.m_U[2])/dy);
        }
        else if (bNoBoundary[3]==3) // NORTH REFLECT
        {
            U0 = cell_center[index_nb[2]].m_state;
            U2 = cell_center[index_nb[3]].m_state;
            onWallvel = 0.5 *(U1.m_U[2] + U2.m_U[2]);
            slope[1] = EBM_minmod(2.0*(onWallvel-U1.m_U[2])/dy, (U1.m_U[2]-U0.m_U[2])/dy);
        }
        else assert(false);

        //w_z
        if ((!bNoBoundary[4] && !bNoBoundary[5]) || bNoBoundary[4] == 3 || bNoBoundary[5] == 3)
        {
            U0 = cell_center[index_nb[4]].m_state;
            U2 = cell_center[index_nb[5]].m_state;
            slope[2] = EBM_minmod((U1.m_U[2]-U0.m_U[2])/dz, (U2.m_U[2]-U1.m_U[2])/dz);
            /*
            if (bNoBoundary[4] == 3)
                printf("U0.m_U[2] = %24.24f\n",U0.m_U[2]);
            if (bNoBoundary[5] == 3)
                printf("U1.m_U[2] = %24.24f\n",U1.m_U[2]);
                */
        }
        else if (bNoBoundary[4] == 2) //cells on LOWER boundary
        {
            U2 = cell_center[index_nb[5]].m_state;
            slope[2] = EBM_minmod((U1.m_U[2]-0.0)/dz, (U2.m_U[2]-U1.m_U[2])/dz);
        }
        else if (bNoBoundary[5] == 2) //cells on UPPER boundary
        {
            U0 = cell_center[index_nb[4]].m_state;
            slope[2] = EBM_minmod((0.0-U0.m_U[2])/dz, (-U0.m_U[2]-0.0)/dz);
        }
        else assert(false);
    }
    else assert(false);
} /* end getLimitedSlope_Velocity_MAC_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getCenterDifference_vd(
        int *icoords,
        EBM_COORD xyz,
        double slope[5])
{
    double dh, dh0, dh1;
    L_STATE U0, U1, U2;

    U1 = cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_state;

    bool bNoBoundary[2];

    if(xyz==COORD_X)
    {
        dh = top_h[0];

        bNoBoundary[0] = getNeighborOrBoundaryState_vd(icoords,WEST,U0,m_t_old);
        if(bNoBoundary[0])
            dh0 = top_h[0];
        else
            dh0 = top_h[0]/2;

        bNoBoundary[1] = getNeighborOrBoundaryState_vd(icoords,EAST,U2,m_t_old);
        if(bNoBoundary[1])
            dh1 = top_h[0];
        else
            dh1 = top_h[0]/2;
    }
    else if(xyz==COORD_Y)
    {
        dh = top_h[1];

        bNoBoundary[0] = getNeighborOrBoundaryState_vd(icoords,SOUTH,U0,m_t_old);
        if(bNoBoundary[0])
            dh0 = top_h[1];
        else
            dh0 = top_h[1]/2;

        bNoBoundary[1] = getNeighborOrBoundaryState_vd(icoords,NORTH,U2,m_t_old);
        if(bNoBoundary[1])
            dh1 = top_h[1];
        else
            dh1 = top_h[1]/2;
    }
    else  //xyz == COORD_Z
    {
        dh = top_h[2];

        bNoBoundary[0] = getNeighborOrBoundaryState_vd(icoords,LOWER,U0,m_t_old);
        if(bNoBoundary[0])
            dh0 = top_h[2];
        else
            dh0 = top_h[2]/2;

        bNoBoundary[1] = getNeighborOrBoundaryState_vd(icoords,UPPER,U2,m_t_old);
        if(bNoBoundary[1])
            dh1 = top_h[2];
        else
            dh1 = top_h[2]/2;
    }
    slope[0] = ((U1.m_U[0]-U0.m_U[0])/dh0 - (U2.m_U[0]-U1.m_U[0])/dh1)/dh;
    slope[1] = ((U1.m_U[1]-U0.m_U[1])/dh0 - (U2.m_U[1]-U1.m_U[1])/dh1)/dh;
    slope[2] = ((U1.m_U[2]-U0.m_U[2])/dh0 - (U2.m_U[2]-U1.m_U[2])/dh1)/dh;
    slope[3] = ((U1.m_rho-U0.m_rho)/dh0 - (U2.m_rho-U1.m_rho)/dh1)/dh;
    slope[4] = ((U1.m_c-U0.m_c)/dh0 - (U2.m_c-U1.m_c)/dh1)/dh;
} /* end getCenterDifference_vd */


bool Incompress_Solver_Smooth_3D_Cartesian::getNeighborOrBoundaryState_vd(
        int icoords[3],
        GRID_DIRECTION dir,
        L_STATE &state,
        double t)
{
    double crx_coords[MAXD];
    static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};
    POINTER intfc_state;
    HYPER_SURF *hs;

    int index_nb;
    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    int comp = cell_center[index].comp;

    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir,
            comp,&intfc_state,&hs,crx_coords,t) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
    {
        if (wave_type(hs) == DIRICHLET_BOUNDARY) //for LOWER & UPPER bdry only
        {
            state.m_U[0] = 0.0;
            state.m_U[1] = 0.0;
            state.m_U[2] = 0.0;
            if (dir == LOWER)
            {
                state.m_rho = m_rho[0];
                state.m_rho_old = m_rho[0];
                state.m_c = m_c[0];
            }
            if (dir == UPPER)
            {
                state.m_rho = m_rho[1];
                state.m_rho_old = m_rho[1];
                state.m_c = m_c[1];
            }
            return false;
        }
        else if (wave_type(hs) == NEUMANN_BOUNDARY) // for LOWER & UPPER bdry only
        {
            state.m_U[0] = 0.0;
            state.m_U[1] = 0.0;
            state.m_U[2] = 0.0;
            state.m_rho = cell_center[index].m_state.m_rho;
            state.m_rho_old = cell_center[index].m_state.m_rho_old;
            state.m_c = cell_center[index].m_state.m_c;
            return false;
        }
        else
        {
            switch(dir)
            {
            case WEST:
                index_nb = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
                break;
            case EAST:
                index_nb = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
                break;
            case SOUTH:
                index_nb = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
                break;
            case NORTH:
                index_nb = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
                break;
            case LOWER:
                index_nb = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
                break;
            case UPPER:
                index_nb = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);
                break;
            default:
                assert(false);
            }
            state = cell_center[index_nb].m_state;
            return true;
        }
    }
    else
    {
        switch(dir)
        {
        case WEST:
            index_nb = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
            break;
        case EAST:
            index_nb = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
            break;
        case SOUTH:
            index_nb = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
            break;
        case NORTH:
            index_nb = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
            break;
        case LOWER:
            index_nb = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
            break;
        case UPPER:
            index_nb = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);
            break;
        default:
            assert(false);
        }
        state = cell_center[index_nb].m_state;
        return true;
    }
} /* end getNeighborOrBoundaryState_vd */

//function return type from bool to int
int Incompress_Solver_Smooth_3D_Cartesian::getNeighborOrBoundaryScalar_MAC_vd(
        int icoords[3],
        GRID_DIRECTION dir,
        L_STATE &state,
        double t)
{
    double crx_coords[MAXD];
    POINTER intfc_state;
    HYPER_SURF *hs;
    int index_nb;
    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    COMPONENT comp = top_comp[index];
    INTERFACE *intfc = front->interf;
    int dirr, side;

    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir,
            comp,&intfc_state,&hs,crx_coords,t) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
    {
        if (wave_type(hs) == DIRICHLET_BOUNDARY)
        {
/*
            state.m_rho = cell_center[index].m_state.m_rho;
            state.m_rho_old = cell_center[index].m_state.m_rho_old;
            state.m_c = cell_center[index].m_state.m_c;
            state.m_Dcoef = cell_center[index].m_state.m_Dcoef;
*/
            state = cell_center[index].m_state;
            clean_up(ERROR);
            return 1;
        }
        else if (wave_type(hs) == NEUMANN_BOUNDARY)
        {
            state = cell_center[index].m_state;
            return 2;
        }
        else //we don't consider other BC types here
        {
            assert (false);
            return 0;
        }
    }
    else
    {
        convertGridDirectionToDirSide(dir, &dirr, &side);
        if (rect_boundary_type(intfc,dirr, side) == REFLECTION_BOUNDARY && FT_Reflect(icoords, dirr, side))
            return 3;
        switch(dir)
        {
        case WEST:
            index_nb = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
            break;
        case EAST:
            index_nb = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
            break;
        case SOUTH:
            index_nb = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
            break;
        case NORTH:
            index_nb = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
            break;
        case LOWER:
            index_nb = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
            break;
        case UPPER:
            index_nb = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);
            break;
        default:
            assert(false);
        }
        state = cell_center[index_nb].m_state;
        return 0;
    }
} /* end getNeighborOrBoundaryScalar_MAC_vd */


//split scheme, use constant extrapolations for ghost cell states rho and c
bool Incompress_Solver_Smooth_3D_Cartesian::getNeighborOrBoundaryScalar_MAC_GhostCell_vd(
        int icoords[3],
        GRID_DIRECTION dir,
        L_STATE &state,
        double t)
{
    double crx_coords[MAXD];
    POINTER intfc_state;
    HYPER_SURF *hs;

    int index_nb;
    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    COMPONENT comp = top_comp[index];

    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir,
            comp,&intfc_state,&hs,crx_coords,t) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
    {
        if (wave_type(hs) == DIRICHLET_BOUNDARY ||
            wave_type(hs) == NEUMANN_BOUNDARY)
        {
/*
            state.m_rho = cell_center[index].m_state.m_rho;
            state.m_rho_old = cell_center[index].m_state.m_rho_old;
            state.m_c = cell_center[index].m_state.m_c;
            state.m_Dcoef = cell_center[index].m_state.m_Dcoef;
*/
            state = cell_center[index].m_state;
            return false;
        }
        else //we don't consider other BC types here
        {
            assert (false);
            return true;
        }
    }
    else
    {
        switch(dir)
        {
        case WEST:
            index_nb = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
            if (top_comp[index] == top_comp[index_nb])
                state = cell_center[index_nb].m_state;
            else
                state = cell_center[index].m_state;
            break;
        case EAST:
            index_nb = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
            if (top_comp[index] == top_comp[index_nb])
                state = cell_center[index_nb].m_state;
            else
                state = cell_center[index].m_state;
            break;
        case SOUTH:
            index_nb = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
            if (top_comp[index] == top_comp[index_nb])
                state = cell_center[index_nb].m_state;
            else
                state = cell_center[index].m_state;
            break;
        case NORTH:
            index_nb = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
            if (top_comp[index] == top_comp[index_nb])
                state = cell_center[index_nb].m_state;
            else
                state = cell_center[index].m_state;
            break;
        case LOWER:
            index_nb = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
            if (top_comp[index] == top_comp[index_nb])
                state = cell_center[index_nb].m_state;
            else
                state = cell_center[index].m_state;
            break;
        case UPPER:
            index_nb = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);
            if (top_comp[index] == top_comp[index_nb])
                state = cell_center[index_nb].m_state;
            else
                state = cell_center[index].m_state;
            break;
        default:
            assert(false);
        }
        return true;
    }
} /* end getNeighborOrBoundaryScalar_MAC_GhostCell_vd */


bool Incompress_Solver_Smooth_3D_Cartesian::getNeighborOrBoundaryState_tmp_vd(
        int icoords[3],
        GRID_DIRECTION dir,
        L_STATE &state,
        double t)
{
    return false;
} /* end getNeighborOrBoundaryState_tmp_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getRiemannSolution_Velocity_vd(
        EBM_COORD xyz,
        L_STATE &state_left,
        L_STATE &state_right,
        L_STATE &ans,
        int *icoords,
        int *ICoords)
{
    L_STATE sl, sr;
    int index, index_nb;
    double tmp;

    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    index_nb = d_index3d(ICoords[0],ICoords[1],ICoords[2],top_gmax);

    // rotate state
    if(xyz==COORD_X)
    {
        sl.m_U[0] = state_left.m_U[0];
        sl.m_U[1] = state_left.m_U[1];
        sl.m_U[2] = state_left.m_U[2];

        sr.m_U[0] = state_right.m_U[0];
        sr.m_U[1] = state_right.m_U[1];
        sr.m_U[2] = state_right.m_U[2];

        tmp = cell_center[index].m_state.m_U[0] + cell_center[index_nb].m_state.m_U[0];
    }
    else if(xyz==COORD_Y)
    {
        sl.m_U[0] = state_left.m_U[1];
        sl.m_U[1] = state_left.m_U[0];
        sl.m_U[2] = state_left.m_U[2];

        sr.m_U[0] = state_right.m_U[1];
        sr.m_U[1] = state_right.m_U[0];
        sr.m_U[2] = state_right.m_U[2];

        tmp = cell_center[index].m_state.m_U[1] + cell_center[index_nb].m_state.m_U[1];
    }
    else
    {
        sl.m_U[0] = state_left.m_U[2];
        sl.m_U[1] = state_left.m_U[1];
        sl.m_U[2] = state_left.m_U[0];

        sr.m_U[0] = state_right.m_U[2];
        sr.m_U[1] = state_right.m_U[1];
        sr.m_U[2] = state_right.m_U[0];

        tmp = cell_center[index].m_state.m_U[2] + cell_center[index_nb].m_state.m_U[2];
    }

    // calculate the Riemann solution
    double uL = sl.m_U[0];
    double uR = sr.m_U[0];

    // BCG, JCP 85, 257-283 (1989)
    // u_t + u*u_x = 0
    if (tmp > 0.0)
        ans.m_U[0] = uL;
    else
        ans.m_U[0] = uR;

    // v_t + u*v_x = 0
    // w_t + u*w_x = 0
    if (ans.m_U[0] > 0.0)
    {
        ans.m_U[1] = sl.m_U[1];
        ans.m_U[2] = sl.m_U[2];
    }
    else if (ans.m_U[0] < 0.0)
    {
        ans.m_U[1] = sr.m_U[1];
        ans.m_U[2] = sr.m_U[2];
    }
    else
    {
        ans.m_U[1] = 0.5*(sl.m_U[1] + sr.m_U[1]);
        ans.m_U[2] = 0.5*(sl.m_U[2] + sr.m_U[2]);
    }

    // rotate state
    if(xyz==COORD_X)
        ; // do nothing
    else if(xyz==COORD_Y)
        std::swap(ans.m_U[0],ans.m_U[1]);
    else
        std::swap(ans.m_U[0],ans.m_U[2]);
} /* end getRiemannSolution_Velocity_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getRiemannSolution_MAC_CenterVelocity_vd(
        EBM_COORD xyz,
        L_STATE &state_left,
        L_STATE &state_right,
        L_STATE &ans,
        int *icoords)
{
    int i, j, k;
    bool bNoBoundary[6];
    int index, index_nb[6];
    double tmp, uL, uR;
    COMPONENT comp;
    double crx_coords[MAXD];
    POINTER intfc_state;
    HYPER_SURF *hs;

    i = icoords[0];
    j = icoords[1];
    k = icoords[2];
    index = d_index3d(i,j,k,top_gmax);
    comp = top_comp[index];

    //6 neighbours of the center cell
    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    // 4 directions
    bNoBoundary[0] = YES;
    bNoBoundary[1] = YES;
    bNoBoundary[2] = YES;
    bNoBoundary[3] = YES;

    // LOWER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[4] = NO;
    else
        bNoBoundary[4] = YES;

    // UPPER
    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        bNoBoundary[5] = NO;
    else
        bNoBoundary[5] = YES;

    uL = state_left.m_U[xyz];
    uR = state_right.m_U[xyz];

    if (!bNoBoundary[4] && xyz==2) //for w on LOWER boundary cells
        tmp = 0.5*(cell_center[index].m_state.m_U[xyz] + 0.0);
    else if (!bNoBoundary[5] && xyz==2) //for w on UPPER boundary cells
        tmp = 0.5*(0.0 + cell_center[index_nb[2*xyz]].m_state.m_U[xyz]);
    else
        tmp = 0.5*(cell_center[index].m_state.m_U[xyz] +
                   cell_center[index_nb[2*xyz]].m_state.m_U[xyz]);

    //calculate the Riemann solution
    if (tmp > 0.0)
        ans.m_U[xyz] = uL;
    else
        ans.m_U[xyz] = uR;
} /* end getRiemannSolution_MAC_CenterVelocity_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getRiemannSolution_MAC_EdgeVelocity_vd(
        EBM_COORD xyz,
        EBM_COORD XYZ,
        L_STATE &state_left,
        L_STATE &state_right,
        L_STATE &state_LEFT,
        L_STATE &state_RIGHT,
        L_STATE &ans)
{
    double uL,uR,vL,vR;

    ans.m_U[0] = 0.0;
    ans.m_U[1] = 0.0;
    ans.m_U[2] = 0.0;

    switch (xyz + XYZ)
    {
    case 1: //(u,v)
        uL = state_left.m_U[0];
        uR = state_right.m_U[0];
        vL = state_LEFT.m_U[1];
        vR = state_RIGHT.m_U[1];
        break;
    case 2: //(w,u)
        uL = state_left.m_U[2];
        uR = state_right.m_U[2];
        vL = state_LEFT.m_U[0];
        vR = state_RIGHT.m_U[0];
        break;
    case 3: //(v,w)
        uL = state_left.m_U[1];
        uR = state_right.m_U[1];
        vL = state_LEFT.m_U[2];
        vR = state_RIGHT.m_U[2];
        break;
    default:
        assert(false);
    }

    //calculate the Riemann solution
    if ((uL*uR < 0.0) && (vL*vR < 0.0)) //no pair has the same sign
    {
        ans.m_U[0] = 0.5*(uL + uR);
        ans.m_U[1] = 0.5*(vL + vR);
    }
    else if ((uL*uR >= 0.0) && (vL*vR < 0.0)) //only one pair has the same sign
    {
        // ambiguities occur when uL = uR = 0 ???
        if ((uL == 0.0) && (uR == 0.0)) //all zero's
            ;
        if ((uL >= 0.0) && (uR >= 0.0)) //all non-negative
            ans.m_U[1] = vL;
        if ((uL <= 0.0) && (uR <= 0.0)) //all non-positive
            ans.m_U[1] = vR;

        if (ans.m_U[1] > 0.0)
            ans.m_U[0] = uL;
        else if (ans.m_U[1] < 0.0)
            ans.m_U[0] = uR;
        else //ans.m_U[1] == 0.0
            ans.m_U[0] = 0.5*(uL + uR);
    }
    else //vL*vR >= 0.0, either only one pair or both pairs have the same sign
    {
        // ambiguities occur when vL = vR = 0 ???
        if ((vL == 0.0) && (vR == 0.0)) //all zero's
            ;
        if ((vL >= 0.0) && (vR >= 0.0)) //all non-negative
            ans.m_U[0] = uL;
        if ((vL <= 0.0) && (vR <= 0.0)) //all non-positive
            ans.m_U[0] = uR;

        if (ans.m_U[0] > 0.0)
            ans.m_U[1] = vL;
        else if (ans.m_U[0] < 0.0)
            ans.m_U[1] = vR;
        else //ans.m_U[0] == 0.0
            ans.m_U[1] = 0.5*(vL + vR);
    }

    //rotate coordinates
    switch (xyz + XYZ)
    {
    case 1: //(u,v)
        break;
    case 2: //(w,u)
        ans.m_U[2] = ans.m_U[0];
        ans.m_U[0] = ans.m_U[1];
        break;
    case 3: //(v,w)
        ans.m_U[2] = ans.m_U[1];
        ans.m_U[1] = ans.m_U[0];
        break;
    default:
        assert(false);
    }
} /* end getRiemannSolution_MAC_EdgeVelocity_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getRiemannSolution_Scalar_vd(
        L_STATE &state_left,
        L_STATE &state_right,
        L_STATE &ans,
        int *icoords,
        GRID_DIRECTION dir)
{
} /* end getRiemannSolution_Scalar_vd */


void Incompress_Solver_Smooth_3D_Cartesian::getRiemannSolution_MAC_Scalar_vd(
        L_STATE &state_left,
        L_STATE &state_right,
        L_STATE &ans,
        int *icoords,
        GRID_DIRECTION dir)
{
    L_STATE sl, sr;
    int index, index_nb[6];
    double tmp; //U_face_bar after MAC projection

    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    index_nb[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
    index_nb[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
    index_nb[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
    index_nb[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
    index_nb[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
    index_nb[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);

    sl.m_rho  = state_left.m_rho;
    sl.m_c    = state_left.m_c;
    sr.m_rho  = state_right.m_rho;
    sr.m_c    = state_right.m_c;

    switch(dir)
    {
    case WEST:
        tmp = cell_center[index_nb[0]].m_state.m_U_face_bar[0];
        break;
    case EAST:
        tmp = cell_center[index].m_state.m_U_face_bar[0];
        break;
    case SOUTH:
        tmp = cell_center[index_nb[2]].m_state.m_U_face_bar[1];
        break;
    case NORTH:
        tmp = cell_center[index].m_state.m_U_face_bar[1];
        break;
    case LOWER:
        tmp = cell_center[index_nb[4]].m_state.m_U_face_bar[2];
        break;
    case UPPER:
        tmp = cell_center[index].m_state.m_U_face_bar[2];
        break;
    default:
        assert(false);
    }

    // rho_t + u*rho_x = 0
    // c_t + u*c_x = 0
    if (tmp > 0.0)
    {
        ans.m_rho  = sl.m_rho;
        ans.m_c    = sl.m_c;
    }
    else if (tmp < 0.0)
    {
        ans.m_rho  = sr.m_rho;
        ans.m_c    = sr.m_c;
    }
    else
    {
        ans.m_rho  = 0.5*(sl.m_rho + sr.m_rho);
        ans.m_c    = 0.5*(sl.m_c + sr.m_c);
    }
} /* end getRiemannSolution_MAC_Scalar_vd */


void Incompress_Solver_Smooth_3D_Cartesian::setInitialCondition_vd(LEVEL_FUNC_PACK *level_func_pack)
{
        COMPONENT comp;
        double coords[MAXD];
        int i,j,k,l,index,dim,sign;
        double grad_rho[MAXD];
        int index_nb[6],icoords[MAXD];
        double rho,Dcoef,speed,rho0,rho1;
        double max_speed=0;
        double dist=0, depth_diff_layer=0;
        double **vel = iFparams->field->vel;
        FOURIER_POLY *pert = (FOURIER_POLY*)level_func_pack->func_params;

        if (iFparams->width_idl != 0)
            depth_diff_layer = iFparams->width_idl;
        else
            depth_diff_layer = top_h[MAXD-1];
        printf("\nThe width of initial diffusion layer is: %.16g\n",depth_diff_layer);

        FT_MakeGridIntfc(front);
        setDomain_vd();
        setComponent();

        setGlobalIndex();
        setIndexMap();

        dim = front->rect_grid->dim;
        m_rho[0] = m_rho_old[0] = iFparams->rho1;
        m_rho[1] = m_rho_old[1] = iFparams->rho2;
        m_c[0] = iFparams->c1;
        m_c[1] = iFparams->c2;
        m_mu[0] = iFparams->mu1;
        m_mu[1] = iFparams->mu2;
        m_Dcoef[0] = iFparams->Dcoef1;
        m_Dcoef[1] = iFparams->Dcoef2;
        m_comp[0] = iFparams->m_comp1;
        m_comp[1] = iFparams->m_comp2;
        m_smoothing_radius = iFparams->smoothing_radius;
        m_sigma = iFparams->surf_tension;
        mu_min = rho_min = Dcoef_min = c_min = HUGE;
        c_max = -HUGE;
        bGhostCell = NO;
        z0 = pert->z0;

        fprintf(stdout, "PRINT-HK z0 %e depth_diff_layer %e\n",z0,depth_diff_layer);
        fflush(NULL);

/*
        static boolean first = YES;
        if (first)
        {
            first = NO;
            printf("m_rho[0] = %lf, m_rho[1] = %lf\n", m_rho[0],m_rho[1]);
        }
*/

        for (i = 0; i < 2; ++i)
        {
            if (ifluid_comp(m_comp[i]))
            {
                mu_min = std::min(mu_min,m_mu[i]);
                rho_min = std::min(rho_min,m_rho[i]);
                Dcoef_min = std::min(Dcoef_min,m_Dcoef[i]);
                c_min = std::min(c_min,m_c[i]);
                c_max = std::max(c_max,m_c[i]);
            }
        }

        // Initialize state at cell_center
        for (i = 0; i < cell_center.size(); i++)
        {
            getRectangleCenter(i, coords);
            cell_center[i].m_state.setZero(); //set U, P, q and grad_q to be zero
            comp = top_comp[i];
            if (getInitialState != NULL)
                (*getInitialState)(comp,coords,cell_center[i].m_state,dim,iFparams);
                                              //set U to be zero
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            comp = cell_center[index].comp;
            if (!ifluid_comp(comp)) continue;

            cell_center[index].m_state.m_Dcoef = m_Dcoef[0];

            getRectangleCenter(index, coords);
            dist = level_wave_func(level_func_pack->func_params, coords);

            rho0 = std::min(m_rho[0], m_rho[1]);
            rho1 = std::max(m_rho[0], m_rho[1]);
            cell_center[index].m_state.m_rho = 0.5*(m_rho[0]+m_rho[1]) + 0.5*(m_rho[0]-m_rho[1])*
                                               erf(dist/depth_diff_layer*4.0);
            cell_center[index].m_state.m_rho_old = cell_center[index].m_state.m_rho;

            //conservation of volume, c is the mass fraction of light fluid
            rho = cell_center[index].m_state.m_rho;
            cell_center[index].m_state.m_c = rho0*(rho1-rho)/rho/(rho1-rho0);
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_Dcoef;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_Dcoef = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_rho;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_rho = array[index];
            cell_center[index].m_state.m_rho_old = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_c;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_c = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_P;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_P = array[index];
            cell_center[index].m_state.m_q = array[index];
        }

        for (l = 0; l < dim; l++)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }

	//calc div_U before adding diffusion velo to U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l];
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);
        }
        FT_ParallelExchGridArrayBuffer(source,front);
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        if (debugging("step_size"))
        {
            double value,sum_div;
            value = sum_div = max_value = 0.0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                sum_div += cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U_0 without diffusion velocity added is %.16g\n",
                                        sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U_0 without diffusion velocity added is %.16g\n",
                                        max_value);
            max_value = 0.0;
        }

        // ***************************************************************
        //  Add diffusion velocity according to the divergence constraint
        // ***************************************************************
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_rho;
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);

            Dcoef = cell_center[index].m_state.m_Dcoef;
            computeFieldPointGrad_MAC_vd(icoords,array,grad_rho);

            rho = 2.0/(1.0/cell_center[index].m_state.m_rho +
                       1.0/cell_center[index_nb[1]].m_state.m_rho);
            cell_center[index].m_state.m_U[0] -= Dcoef*grad_rho[0]/rho;

            rho = 2.0/(1.0/cell_center[index].m_state.m_rho +
                       1.0/cell_center[index_nb[3]].m_state.m_rho);
            cell_center[index].m_state.m_U[1] -= Dcoef*grad_rho[1]/rho;

            if (ijk_to_I[i][j][k+1] < 0) //cells on UPPER bdry
                rho = cell_center[index].m_state.m_rho;
            else //other cells
                rho = 2.0/(.0/cell_center[index].m_state.m_rho +
                           1.0/cell_center[index_nb[5]].m_state.m_rho);
            cell_center[index].m_state.m_U[2] -= Dcoef*grad_rho[2]/rho;

            speed = fabs(cell_center[index].m_state.m_U[0]) +
                    fabs(cell_center[index].m_state.m_U[1]) +
                    fabs(cell_center[index].m_state.m_U[2]);
            if (speed > max_speed)
                max_speed = speed;
        }
        pp_global_max(&max_speed,1);

        if (debugging("step_size"))
            (void) printf("\nmax_speed of U^0 WITH diffusion velocity added is %24.18g\n",max_speed);

        for (l = 0; l < dim; l++)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }

        //calc div_U after adding diffusion velo to U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l];
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);
        }
        FT_ParallelExchGridArrayBuffer(source,front);
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        if (debugging("step_size"))
        {
            double value,sum_div;
            value = sum_div = max_value = 0.0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                sum_div += cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U^0 WITH diffusion velocity added is %.16g\n",
                                        sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U^0 WITH diffusion velocity added is %.16g\n",
                                        max_value);
            max_value = 0.0;
        }

        if (debugging("step_size"))
        {
            double value, sum_div, diffusion;
            value = max_value = sum_div = 0.0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;
                getDivU_MAC_vd(icoords,&diffusion,0,bGhostCell); //get the divergence constraint S^0
                source[index] = diffusion;
            }
            FT_ParallelExchGridArrayBuffer(source,front);

            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.div_U_tmp = source[index];
            }

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U_tmp);
                sum_div += cell_center[index].m_state.div_U_tmp;
                if(value > max_value)
                    max_value = value;
            }

            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of S^0 is %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of S^0 is %.16g\n",max_value);
            (void) printf("\n");
            max_value = 0.0;
        }

        computeGradientQ_MAC_vd();
        copyMeshStates_vd();
        setAdvectionDt_vd();
} /* end setInitialCondition_vd */


void Incompress_Solver_Smooth_3D_Cartesian::setInitialCondition_RSRV_vd(LEVEL_FUNC_PACK *level_func_pack)
{
        COMPONENT comp;
        double coords[MAXD];
        int i,j,k,l,index,index_nb[2*MAXD],I,I_nb[2*MAXD],dim,sign,icoords[MAXD];
        double grad_rho[MAXD];
        double rho,mu,speed, max_speed=0;
        double x_face=0, z_face=0;
        double z_intfc_pert=0, initial_diff_layer=0;
        double **vel = iFparams->field->vel;
        FOURIER_POLY *pert = (FOURIER_POLY*)level_func_pack->func_params;
        RECT_GRID *gr = front->rect_grid;

	if (iFparams->width_idl != 0)
	    initial_diff_layer = iFparams->width_idl;
	else
	    initial_diff_layer = top_h[MAXD-1];
        printf("\nThe depth of initial diffusion layer for RSRV case is: %.16g\n", initial_diff_layer);
        fprintf(stdout, "PRINT-HK z0 %e initial_diff_layer %e\n",pert->z0,initial_diff_layer);
        fflush(NULL);


        FT_MakeGridIntfc(front);
        setDomain_vd();
        setComponent();

        setGlobalIndex();
        setIndexMap();

        dim = front->rect_grid->dim;
        m_rho[0] = m_rho_old[0] = iFparams->rho1;
        m_rho[1] = m_rho_old[1] = iFparams->rho2;
        m_c[0] = iFparams->c1;
        m_c[1] = iFparams->c2;
        m_mu[0] = iFparams->mu1;
        m_mu[1] = iFparams->mu2;
        m_Dcoef[0] = iFparams->Dcoef1;
        m_Dcoef[1] = iFparams->Dcoef2;
        m_comp[0] = iFparams->m_comp1;
        m_comp[1] = iFparams->m_comp2;
        m_smoothing_radius = iFparams->smoothing_radius;
        m_sigma = iFparams->surf_tension;
        mu_min = rho_min = Dcoef_min = c_min = HUGE;
        c_max = -HUGE;
        bGhostCell = NO;
        z0 = pert->z0;

        for (i = 0; i < 2; ++i)
        {
            if (ifluid_comp(m_comp[i]))
            {
                mu_min = std::min(mu_min,m_mu[i]);
                rho_min = std::min(rho_min,m_rho[i]);
                Dcoef_min = std::min(Dcoef_min,m_Dcoef[i]);
                c_min = std::min(c_min,m_c[i]);
                c_max = std::max(c_max,m_c[i]);
            }
        }

        //Initialize state at cell_center
        for (i = 0; i < cell_center.size(); ++i)
        {
            getRectangleCenter(i, coords);
            cell_center[i].m_state.setZero(); //set U, P, q, D_t, mu_t to be zero
            comp = top_comp[i];
            if (getInitialState != NULL)
                (*getInitialState)(comp,coords,cell_center[i].m_state,dim,iFparams);
                                              //set U[i] to be zero
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            comp = cell_center[index].comp;
            if (!ifluid_comp(comp)) continue;

            cell_center[index].m_state.m_Dcoef = m_Dcoef[0];

            getRectangleCenter(index, coords);


            for (l=0; l<dim; ++l) {
                if (coords[l] < gr->GL[l]  ||  coords[l] > gr->GU[l]) {
                    printf("ERROR: (%lf, %lf, %lf) out of range in direction %d\n", coords[0],coords[1],coords[2],l);
                    clean_up(ERROR);
                }
            }

            z_intfc_pert = pert_interface_vd(pert,coords,0,dim); //TODO: should call "pert_height_vd"
            //For SY
            //double dist = level_wave_func_cylindrical_init(level_func_pack->func_params, coords);

            cell_center[index].m_state.m_rho = 0.5*(m_rho[0]+m_rho[1]) + 0.5*(m_rho[0]-m_rho[1])*
                                               erf((coords[dim-1]-z0 + z_intfc_pert-z0)/initial_diff_layer*2.0);
//                                               erf((coords[dim-1]-z_intfc_pert)/initial_diff_layer*2.0)
            //For SY
            //cell_center[index].m_state.m_rho = 0.5*(m_rho[0]+m_rho[1]) + 0.5*(m_rho[0]-m_rho[1])*erf(dist/initial_diff_layer*4.0);
            cell_center[index].m_state.m_rho_old = cell_center[index].m_state.m_rho;
            //define dynamic mu = rho*(mu0+mu1)/(rho0+rho1) as eqn. (3.3) of Mueschke thesis
            mu = cell_center[index].m_state.m_rho*(m_mu[0]+m_mu[1])/(m_rho[0]+m_rho[1]);
            cell_center[index].m_state.m_mu = cell_center[index].m_state.m_mu_old = mu;

            //m_c: mass fraction
            //see Eq.(9b) of Transition stages of Rayleigh-Taylor instability between miscible fluids (2001)
            rho = cell_center[index].m_state.m_rho;
            cell_center[index].m_state.m_c = m_rho[0]*(m_rho[1]-rho)/rho/(m_rho[1]-m_rho[0]);
            if (cell_center[index].m_state.m_c < 0.0)	cell_center[index].m_state.m_c = 0.0;
            if (cell_center[index].m_state.m_c > 1.0)	cell_center[index].m_state.m_c = 1.0;


            //for MAC grid
            x_face = coords[0] + 0.5*top_h[0];
            if (pert->pert_bdry_type[0] == PERIODIC)
            {
                while (x_face < gr->GL[0])
                    x_face += (gr->GU[0] - gr->GL[0]);
                while (x_face > gr->GU[0])
                    x_face -= (gr->GU[0] - gr->GL[0]);
            }
            else
	    {
                printf("pert->pert_bdry_type[0] == Other BC\n");
                clean_up(ERROR);
            }
            z_face = coords[dim-1] + 0.5*top_h[dim-1];

            for (l = 0; l < pert->num_modes_x; ++l)
            {
                sign = (coords[dim-1] <= z0)? -1 : 1;
                cell_center[index].m_state.m_U[0] += sign*pert->vel_A[l]*cos(pert->nu[l][0]*x_face-pert->vel_phase[l])*
                                                     exp(-1.0*pert->nu[l][0]*fabs(coords[dim-1]-z0));
                cell_center[index].m_state.m_U[dim-1] -= pert->vel_A[l]*sin(pert->nu[l][0]*coords[0]-pert->vel_phase[l])*
                                                         exp(-1.0*pert->nu[l][0]*fabs(z_face-z0));
            }
        }

        //scatter states
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_Dcoef;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_Dcoef = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_mu;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_mu = array[index];
            cell_center[index].m_state.m_mu_old = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_rho;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_rho = array[index];
            cell_center[index].m_state.m_rho_old = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_c;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_c = array[index];
        }

        for (l = 0; l < dim; ++l)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }

        for (l = 0; l < dim; ++l)
        {
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U_velo_var[l] = 0;
            }
        }

        for (l = 0; l < dim+1; ++l)
        {
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_mu_turbulent[l] = 0;
                cell_center[index].m_state.m_Dcoef_turbulent[l] = 0;
            }
        }

        //calc div_U before adding diffusion velo to U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l];
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);
        }
        FT_ParallelExchGridArrayBuffer(source,front);
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        if (debugging("step_size"))
        {
            double value,sum_div;
            value = sum_div = max_value = max_speed = 0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;
                index = d_index3d(i,j,k,top_gmax);
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
                        fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                    max_speed = speed;
            }
            pp_global_max(&max_speed,1);
            (void) printf("\nmax_speed of U_0 without diffusion velocity added is %.16g\n",max_speed);

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                sum_div += cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U_0 without diffusion velocity added is %.16g\n",
                                        sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U_0 without diffusion velocity added is %.16g\n",
                                        max_value);
            max_value = 0;
        }

        // ***************************************************************
        //     Add diffusion velocity due to the divergence constraint
        // ***************************************************************

        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_rho;
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index = d_index3d(i,j,k,top_gmax);
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
            computeFieldPointGradRho_MAC_vd(icoords,array,grad_rho);
            for (l = 0; l < dim; ++l) {
                if (I_nb[5] >= 0)
                    rho = 2.0/(1.0/cell_center[index].m_state.m_rho + 1.0/cell_center[index_nb[2*l+1]].m_state.m_rho);
                else //UPPER bdry
                    rho = cell_center[index].m_state.m_rho;

                cell_center[index].m_state.m_U[l] -= cell_center[index].m_state.m_Dcoef*grad_rho[l]/rho;
            }

            speed = fabs(cell_center[index].m_state.m_U[0]) +
                    fabs(cell_center[index].m_state.m_U[1]) +
                    fabs(cell_center[index].m_state.m_U[2]);
            if (speed > max_speed)
                max_speed = speed;
        }
        pp_global_max(&max_speed,1);

        if (debugging("step_size"))
            (void) printf("\nmax_speed of U_0 with diffusion velocity added is %.16g\n",max_speed);

        for (l = 0; l < dim; ++l)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }

        //calc div_U after adding diffusion velo to U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l];
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);
        }
        FT_ParallelExchGridArrayBuffer(source,front);
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        if (debugging("step_size"))
        {
            double value,sum_div;
            value = sum_div = max_value = 0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                sum_div += cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U_0 with diffusion velocity added is %.16g\n",
                                        sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U_0 with diffusion velocity added is %.16g\n",
                                        max_value);
            max_value = 0;
        }

        if (debugging("step_size"))
        {
            double value, sum_div, diffusion;
            value = max_value = sum_div = 0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;
                getDivU_MAC_vd(icoords,&diffusion,0,NO); //get the divergence constraint S^0
                source[index] = diffusion;
            }
            FT_ParallelExchGridArrayBuffer(source,front);
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.div_U_tmp = source[index];
            }

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U_tmp);
                sum_div += cell_center[index].m_state.div_U_tmp;
                if(value > max_value)
                    max_value = value;
            }

            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of S^0 is %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of S^0 is %.16g\n",max_value);
            (void) printf("\n");
            max_value = 0;
        }

        fflush(stdout);


        // ***************************************************************
        //    Initialize pressure field by laplace(p) = div(rho*g)
        //                                 dp/dn = rho*g on bdry
        // ***************************************************************
        double P_max,P_min;
        double gz = iFparams->gravity[dim-1];

        for (l = 0; l < dim; l++)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_rho*iFparams->gravity[l];
            //nonhomogeneous Neumann BC for top/bottom bdry
            diff_coeff[index] = cell_center[index].m_state.m_rho*gz;
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            if (I == -1) continue;

            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            //rhs = div(rho*g), vel stands for g.
            source[index] = computeFieldPointDiv_Neumann_vd(icoords,vel);
        }
        FT_ParallelExchGridArrayBuffer(source,front);
        poisson_solver3d_P0_vd(front,ilower,iupper,ijk_to_I,source,diff_coeff,
                        array,&P_max,&P_min);

        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_P = array[index];
            cell_center[index].m_state.m_q = array[index];
        }


        computeGradientQ_MAC_vd();
        copyMeshStates_vd();
        setAdvectionDt_vd();
} /* end setInitialCondition_RSRV_vd */


//setInitialCondition_RSSY_MAC_vd
void Incompress_Solver_Smooth_3D_Cartesian::setInitialCondition_RSSY_vd(LEVEL_FUNC_PACK *level_func_pack)
{
        COMPONENT comp;
        double coords[MAXD];
        int i,j,k,l,index,dim,sign;
        double grad_rho[MAXD];
        int index_nb[6],icoords[MAXD];
        double rho,Dcoef,speed,rho0,rho1;
        double max_speed=0;
        double dist=0, depth_diff_layer=0;
        double **vel = iFparams->field->vel;
        FOURIER_POLY *pert = (FOURIER_POLY*)level_func_pack->func_params;

        //depth_diff_layer = top_h[MAXD-1];
        depth_diff_layer = iFparams->width_idl;
        printf("\nThe width of initial diffusion layer for RSSY case is: %.16g\n",depth_diff_layer);

        FT_MakeGridIntfc(front);
        setDomain_vd();
        setComponent();

        setGlobalIndex();
        setIndexMap();

        dim = front->rect_grid->dim;
        m_rho[0] = m_rho_old[0] = iFparams->rho1;
        m_rho[1] = m_rho_old[1] = iFparams->rho2;
        m_c[0] = iFparams->c1;
        m_c[1] = iFparams->c2;
        m_mu[0] = iFparams->mu1;
        m_mu[1] = iFparams->mu2;
        m_Dcoef[0] = iFparams->Dcoef1;
        m_Dcoef[1] = iFparams->Dcoef2;
        m_comp[0] = iFparams->m_comp1;
        m_comp[1] = iFparams->m_comp2;
        m_smoothing_radius = iFparams->smoothing_radius;
        m_sigma = iFparams->surf_tension;
        mu_min = rho_min = Dcoef_min = c_min = HUGE;
        c_max = -HUGE;
        bGhostCell = NO;
        z0 = pert->z0;

        for (i = 0; i < 2; ++i)
        {
            if (ifluid_comp(m_comp[i]))
            {
                mu_min = std::min(mu_min,m_mu[i]);
                rho_min = std::min(rho_min,m_rho[i]);
                Dcoef_min = std::min(Dcoef_min,m_Dcoef[i]);
                c_min = std::min(c_min,m_c[i]);
                c_max = std::max(c_max,m_c[i]);
            }
        }

        // Initialize state at cell_center
        for (i = 0; i < cell_center.size(); i++)
        {
            getRectangleCenter(i, coords);
            cell_center[i].m_state.setZero(); //set U, P, q and grad_q to be zero
            comp = top_comp[i];
            if (getInitialState != NULL)
                (*getInitialState)(comp,coords,cell_center[i].m_state,dim,iFparams);
                                              //set U[i] to be zero
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            comp = cell_center[index].comp;
            if (!ifluid_comp(comp)) continue;

            cell_center[index].m_state.m_Dcoef = m_Dcoef[0];

            getRectangleCenter(index, coords);
            dist = level_wave_func_Meniscus(level_func_pack->func_params, coords);

            rho0 = std::min(m_rho[0], m_rho[1]);
            rho1 = std::max(m_rho[0], m_rho[1]);
            cell_center[index].m_state.m_rho = 0.5*(m_rho[1]+m_rho[0]) + 0.5*(m_rho[0]-m_rho[1])*
                                               erf(dist/depth_diff_layer*2.0);
            cell_center[index].m_state.m_rho_old = cell_center[index].m_state.m_rho;

            //conservation of volume
            rho = cell_center[index].m_state.m_rho;
            cell_center[index].m_state.m_c = rho0*(rho1-rho)/rho/(rho1-rho0);

/* ??????????????????????????????????????????
            for(l = 0; l < 25; ++l)
            {
                cell_center[index].m_state.m_U[0] += pert->vel_A[l]*exp(-1.0*(pert->nu[l][0])*fabs(dist))*
                                                     cos(pert->nu[l][0]*coords[0]-pert->vel_phase[l]);
                if(dist <= 0)
                    sign = -1;
                else
                    sign = 1;
                cell_center[index].m_state.m_U[dim-1] += -1.0*sign*pert->vel_A[l]*exp(-1.0*(pert->nu[l][0])*fabs(dist))*
                                                         sin(pert->nu[l][0]*coords[0]-pert->vel_phase[l]);
            }
*/
        }

        /*
         * removal tag: HAOZ
         * REFLECTION BOUNDARY CONDITION COMMENTS. WILL BE REMOVED LATER.
         * Search for HAOZ as tag for all possible comment spots
         * scatMeshArray() mirror
         * scalar field:
         * m_Dcoef
         * m_rho
         * m_rho_old
         * m_c
         * m_P
         * vector field:
         * m_U
         *
        */
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_Dcoef;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_Dcoef = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_rho;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_rho = array[index];
            cell_center[index].m_state.m_rho_old = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_c;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_c = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_P;
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_P = array[index];
            cell_center[index].m_state.m_q = array[index];
        }

        for (l = 0; l < dim; l++)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();//TODO && FIXME: REFLECTION BOUNDARY CONDITION REQUIRES SPECIAL VELOCITY TREATMENT. removal tag HAOZ
            //Solute_Reflect(l, array); // removal tag HAOZ
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }

        //calc div_U before adding diffusion velo to U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l];
        }
        /*
         * removal tag: HAOZ
         *
         * before computation starts, the vector field need further REFLECTION BOUNDARY CONDITION TREATMENT.
         * concern: match field->vel and cell_center().m_state.m_U
         *
        */
        enforceReflectionState(vel);// on m_U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U[l] = vel[l][index];
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            // removal tag: HAOZ, REFLECTION B.C.
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);//removal tag: adjust Reflection B.C. afterwards
        }
        // removal tag: HAOZ adjust
        FT_ParallelExchGridArrayBuffer(source,front);
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        if (debugging("step_size"))
        {
            double value,sum_div;
            value = sum_div = max_value = 0.0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                sum_div += cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U^0 WITHOUT diffusion velocity added is %.16g\n",
                                        sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U^0 WITHOUT diffusion velocity added is %.16g\n",
                                        max_value);
            max_value = 0.0;
        }

        // ***************************************************************
        //  Add diffusion velocity according to the divergence constraint
        // ***************************************************************
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_rho;
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);

            Dcoef = cell_center[index].m_state.m_Dcoef;
            computeFieldPointGrad_MAC_vd(icoords,array,grad_rho);

            rho = 2.0/(1.0/cell_center[index].m_state.m_rho +
                       1.0/cell_center[index_nb[1]].m_state.m_rho);
            cell_center[index].m_state.m_U[0] -= Dcoef*grad_rho[0]/rho;

            rho = 2.0/(1.0/cell_center[index].m_state.m_rho +
                       1.0/cell_center[index_nb[3]].m_state.m_rho);
            cell_center[index].m_state.m_U[1] -= Dcoef*grad_rho[1]/rho;

            if (ijk_to_I[i][j][k+1] < 0) //cells on UPPER bdry
                rho = cell_center[index].m_state.m_rho;
            else //other cells
                rho = 2.0/(1.0/cell_center[index].m_state.m_rho +
                           1.0/cell_center[index_nb[5]].m_state.m_rho);
            cell_center[index].m_state.m_U[2] -= Dcoef*grad_rho[2]/rho;

            speed = fabs(cell_center[index].m_state.m_U[0]) +
                    fabs(cell_center[index].m_state.m_U[1]) +
                    fabs(cell_center[index].m_state.m_U[2]);
            if (speed > max_speed)
                max_speed = speed;
        }
        pp_global_max(&max_speed,1);

        if (debugging("step_size"))
            (void) printf("\nmax_speed of U^0 WITH diffusion velocity added is %24.18g\n",max_speed);

        for (l = 0; l < dim; ++l)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();//TODO && FIXME: REFLECTION BOUNDARY CONDITION
            //Solute_Reflect(l, array);
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }

        //calc div_U after adding diffusion velo to U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l];
        }
        // removal tag: HAOZ
        // REFLECTION installed
        enforceReflectionState(vel);// on m_U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U[l] = vel[l][index];
        }
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index3d(i,j,k,top_gmax);
            source[index] = computeFieldPointDiv_MAC_vd(icoords,vel);//removal tag: adjust Reflection B.C. afterwards
        }
        // removal tag: HAOZ adjust
        FT_ParallelExchGridArrayBuffer(source,front);
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.div_U = source[index];
        }

        if (debugging("step_size"))
        {
            double value,sum_div;
            value = sum_div = max_value = 0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                //TODO && FIXME: a better way to understand divergence of velocity.
                sum_div += cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U^0 WITH diffusion velocity added is %.16g\n",
                                        sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U^0 WITH diffusion velocity added is %.16g\n",
                                        max_value);
            // TODO && FIXME: If it's not divergence free, call clean_up(ERROR)
            /*
            if (fabs(sum_div) > 1e-12)
            {
                printf("It's not divergence free!\n");
                clean_up(ERROR);
            }
            */
            max_value = 0;
        }

        if (debugging("step_size"))
        {
            double value, sum_div, diffusion;
            value = max_value = sum_div = 0;

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;
                // removal tag: HAOZ
                // scalar function
                // wrt REFLECTION
                getDivU_MAC_vd(icoords,&diffusion,0,bGhostCell); //get the divergence constraint S^0
                source[index] = diffusion;
            }
            FT_ParallelExchGridArrayBuffer(source,front);

            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.div_U_tmp = source[index];
            }

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                value = fabs(cell_center[index].m_state.div_U_tmp);
                sum_div += cell_center[index].m_state.div_U_tmp;
                if(value > max_value)
                    max_value = value;
            }

            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of S^0 is %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of S^0 is %.16g\n",max_value);
            (void) printf("\n");
            max_value = 0;
        }
//TEST INITIALIZE PRESSURE START
//TEST INITIALIZE PRESSURE END
        computeGradientQ_MAC_vd();
        copyMeshStates_vd();
        setAdvectionDt_vd();
} /* end setInitialCondition_RSSY_vd */


void Incompress_Solver_Smooth_3D_Cartesian::printInteriorVelocity(char *out_name)
{
        int   i,j,k,l,index,totalpoints;
        double coord_x,coord_y,coord_z;
        double vel_x, vel_y, vel_z;
        int first;
        int step = front->step;
        char name[200];
        char dirname[256];
        static char   *filename = NULL;
        static size_t filename_len = 0;
        FILE *outfile;

        if (pp_numnodes() > 1)
            sprintf(dirname,"%s/P-%s",out_name,right_flush(pp_mynode(),4));
        else
            sprintf(dirname,"%s",out_name);
        sprintf(dirname,"%s/state",dirname,right_flush(pp_mynode(),4));
        if (!create_directory(dirname,YES))
        {
            screen("Cannot create directory %s\n",dirname);
            clean_up(ERROR);
        }

        sprintf(name,"state.t%s",right_flush(step,7));
        if (pp_numnodes() > 1)
            sprintf(name,"%s-p%s",name,right_flush(pp_mynode(),4));

        filename = get_vtk_file_name(filename,dirname,name,&filename_len);

        if (iFparams->movie_option->plot_velo ||
            iFparams->movie_option->plot_pres ||
            iFparams->movie_option->plot_dens ||
            iFparams->movie_option->plot_conc ||
            iFparams->movie_option->plot_comp ||
            iFparams->movie_option->plot_vort)
        {
	    outfile = fopen(filename,"w");
            fprintf(outfile,"# vtk DataFile Version 3.0\n");
            fprintf(outfile,"States of the whole computational domain\n");
            fprintf(outfile,"ASCII\n");
            fprintf(outfile,"DATASET STRUCTURED_GRID\n");
            int pointsx, pointsy, pointsz;

	    pointsz = top_gmax[2] + 1;
            pointsy = top_gmax[1] + 1;
            pointsx = top_gmax[0] + 1;
            totalpoints = pointsx*pointsy*pointsz;
            fprintf(outfile, "DIMENSIONS %d %d %d\n", pointsz, pointsy, pointsx);
            fprintf(outfile, "POINTS %d double\n", totalpoints);
            for(i = 0; i <= top_gmax[0]; ++i)
            for(j = 0; j <= top_gmax[1]; ++j)
            for(k = 0; k <= top_gmax[2]; ++k)
            {
                index = d_index3d(i,j,k,top_gmax);
                coord_x = cell_center[index].m_coords[0];
                coord_y = cell_center[index].m_coords[1];
                coord_z = cell_center[index].m_coords[2];
                fprintf(outfile, "%.16g %.16g %.16g\n", coord_x, coord_y, coord_z);
            }
            fprintf(outfile, "POINT_DATA %i\n", totalpoints);
            if(iFparams->movie_option->plot_velo)
            {
                fprintf(outfile, "VECTORS velocity double\n");
                for(i = 0; i <= top_gmax[0]; ++i)
                for(j = 0; j <= top_gmax[1]; ++j)
                for(k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    vel_x = cell_center[index].m_state.m_U[0];
                    vel_y = cell_center[index].m_state.m_U[1];
                    vel_z = cell_center[index].m_state.m_U[2];
                    fprintf(outfile, "%.16g %.16g %.16g\n",vel_x, vel_y, vel_z);
                }
            }
            if(iFparams->movie_option->plot_velo)
            {
                fprintf(outfile, "SCALARS u double\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
                for(i = 0; i <= top_gmax[0]; ++i)
                for(j = 0; j <= top_gmax[1]; ++j)
                for(k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    vel_x = cell_center[index].m_state.m_U[0];
                    fprintf(outfile, "%.16g\n",vel_x);
                }
            }
            if(iFparams->movie_option->plot_velo)
            {
                fprintf(outfile, "SCALARS v double\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
                for(i = 0; i <= top_gmax[0]; ++i)
                for(j = 0; j <= top_gmax[1]; ++j)
                for(k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    vel_y = cell_center[index].m_state.m_U[1];
                    fprintf(outfile, "%.16g\n",vel_y);
                }
            }
            if(iFparams->movie_option->plot_velo)
            {
                fprintf(outfile, "SCALARS w double\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
                for(i = 0; i <= top_gmax[0]; ++i)
                for(j = 0; j <= top_gmax[1]; ++j)
                for(k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    vel_z = cell_center[index].m_state.m_U[2];
                    fprintf(outfile, "%.16g\n",vel_z);
                }
            }

            if(iFparams->movie_option->plot_velo)
            {
                fprintf(outfile, "SCALARS DivU double\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
                for(i = 0; i <= top_gmax[0]; ++i)
                for(j = 0; j <= top_gmax[1]; ++j)
                for(k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    fprintf(outfile, "%.16g\n",cell_center[index].m_state.div_U);
                }
            }

            if(iFparams->movie_option->plot_pres)
            {
                fprintf(outfile, "SCALARS pressure double\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
                for(i = 0; i <= top_gmax[0]; ++i)
                for(j = 0; j <= top_gmax[1]; ++j)
                for(k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    fprintf(outfile,"%.16g\n",cell_center[index].m_state.m_P);
                }
            }
            //DEBUG REFLECTION BOUNDARY CONDITION. TODO && FIXME:
            {
                fprintf(outfile, "SCALARS Dcoef double\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
                for(i = 0; i <= top_gmax[0]; ++i)
                for(j = 0; j <= top_gmax[1]; ++j)
                for(k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    fprintf(outfile,"%.16g\n",cell_center[index].m_state.m_Dcoef);
                }
            }
            {
                fprintf(outfile, "SCALARS Phi double\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
                for(i = 0; i <= top_gmax[0]; ++i)
                for(j = 0; j <= top_gmax[1]; ++j)
                for(k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    fprintf(outfile,"%.16g\n",cell_center[index].m_state.m_phi);
                }
            }
            if(iFparams->movie_option->plot_dens)
            {
                fprintf(outfile, "SCALARS density double\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
                for(i = 0; i <= top_gmax[0]; ++i)
                for(j = 0; j <= top_gmax[1]; ++j)
                for(k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    fprintf(outfile,"%.16g\n",cell_center[index].m_state.m_rho);
                }
            }
            if(iFparams->movie_option->plot_conc)
            {
                fprintf(outfile, "SCALARS concentration double\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
                for(i = 0; i <= top_gmax[0]; ++i)
                for(j = 0; j <= top_gmax[1]; ++j)
                for(k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    fprintf(outfile,"%.16g\n",cell_center[index].m_state.m_c);
                }
            }
            if(iFparams->movie_option->plot_comp)
            {
                fprintf(outfile, "SCALARS component int\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
                for(i = 0; i <= top_gmax[0]; ++i)
                for(j = 0; j <= top_gmax[1]; ++j)
                for(k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    fprintf(outfile,"%d\n",top_comp[index]);
                }
            }
            if(iFparams->movie_option->plot_vort)
            {
            }
            fclose(outfile);
        }
}  //end printInteriorVelocity_vd


/*
void Incompress_Solver_Smooth_3D_Cartesian::printInteriorVelocity(char *out_name)
{
        int   i,j,k,l,index,totalpoints;
        double coord_x,coord_y,coord_z;
        double vel_x, vel_y, vel_z;
        int first;
        int step = front->step;
        char name[200];
        char dirname[256];
        static char   *filename = NULL;
        static size_t filename_len = 0;
        FILE *outfile;

        if (pp_numnodes() > 1)
            sprintf(dirname,"%s/P-%s",out_name,right_flush(pp_mynode(),4));
        else
            sprintf(dirname,"%s",out_name);
        sprintf(dirname,"%s/state",dirname,right_flush(pp_mynode(),4));
        if (!create_directory(dirname,YES))
        {
            screen("Cannot create directory %s\n",dirname);
            clean_up(ERROR);
        }

        sprintf(name,"state.t%s",right_flush(step,7));
        if (pp_numnodes() > 1)
            sprintf(name,"%s-p%s",name,right_flush(pp_mynode(),4));

        filename = get_vtk_file_name(filename,dirname,name,&filename_len);

        if (iFparams->movie_option->plot_velo ||
            iFparams->movie_option->plot_pres ||
            iFparams->movie_option->plot_dens ||
            iFparams->movie_option->plot_comp ||
            iFparams->movie_option->plot_vort)
        {
            outfile = fopen(filename,"w");
            fprintf(outfile,"# vtk DataFile Version 3.0\n");
            fprintf(outfile,"States of the whole computational domain\n");
            fprintf(outfile,"ASCII\n");
            fprintf(outfile,"DATASET STRUCTURED_GRID\n");
            int pointsx, pointsy, pointsz;

            pointsz = kmax - kmin + 2;
            pointsy = jmax - jmin + 2;
            pointsx = imax - imin + 2;
            totalpoints = pointsx*pointsy*pointsz;
            fprintf(outfile, "DIMENSIONS %d %d %d\n", pointsz, pointsy, pointsx);
            fprintf(outfile, "POINTS %d double\n", totalpoints);
            for(i = 0; i < pointsx; ++i)
            for(j = 0; j < pointsy; ++j)
            for(k = 0; k < pointsz; ++k)
            {
                index = d_index3d(i,j,k,top_gmax);
                coord_x = cell_center[index].m_coords[0];
                coord_y = cell_center[index].m_coords[1];
                coord_z = cell_center[index].m_coords[2];
                fprintf(outfile, "%.16g %.16g %.16g\n", coord_x, coord_y, coord_z);
            }
            fprintf(outfile, "POINT_DATA %i\n", totalpoints);
            if(iFparams->movie_option->plot_velo)
            {
                fprintf(outfile, "VECTORS velocity double\n");
                for(i = 0; i < pointsx; ++i)
                for(j = 0; j < pointsy; ++j)
                for(k = 0; k < pointsz; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    vel_x = cell_center[index].m_state.m_U[0];
                    vel_y = cell_center[index].m_state.m_U[1];
                    vel_z = cell_center[index].m_state.m_U[2];
                    fprintf(outfile, "%.16g %.16g %.16g\n",vel_x, vel_y, vel_z);
                }
            }
            if(iFparams->movie_option->plot_pres)
            {
                fprintf(outfile, "SCALARS pressure double\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
                for(i = 0; i < pointsx; ++i)
                for(j = 0; j < pointsy; ++j)
                for(k = 0; k < pointsz; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    fprintf(outfile,"%.16g\n",cell_center[index].m_state.m_P);
                }
            }
            if(iFparams->movie_option->plot_dens)
            {
                fprintf(outfile, "SCALARS density double\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
                for(i = 0; i < pointsx; ++i)
                for(j = 0; j < pointsy; ++j)
                for(k = 0; k < pointsz; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    fprintf(outfile,"%.16g\n",cell_center[index].m_state.m_rho);
                }
            }
            if(iFparams->movie_option->plot_comp)
            {
                fprintf(outfile, "SCALARS component int\n");
                fprintf(outfile, "LOOKUP_TABLE default\n");
                for(i = 0; i < pointsx; ++i)
                for(j = 0; j < pointsy; ++j)
                for(k = 0; k < pointsz; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    fprintf(outfile,"%d\n",top_comp[index]);
                }
            }
            if(iFparams->movie_option->plot_vort)
            {
            }
            fclose(outfile);
        }
} // end printInteriorVelocity
*/


void Incompress_Solver_Smooth_3D_Cartesian::printExpandedMesh(char *out_name, bool binary)
{
    if(binary == YES)
    {
        if(hardware_is_little_endian())
            printExpandedMesh_little_endian(out_name);
        else
            printExpandedMesh_big_endian(out_name);
    }
    else
        printExpandedMesh_ascii(out_name);
}


void Incompress_Solver_Smooth_3D_Cartesian::printExpandedMesh_big_endian(char *out_name)
{
}


void Incompress_Solver_Smooth_3D_Cartesian::printExpandedMesh_little_endian(char *out_name)
{
}


void Incompress_Solver_Smooth_3D_Cartesian::printExpandedMesh_ascii(char *out_name)
{
        int   i,j,k,l,index,totalpoints;
        double coord_x,coord_y,coord_z;
        int first;
        //char filename[200];
        //FILE *outfile;
        char dirname[256];
        char filename1[200],filename2[200],filename3[200],filename4[200],filename5[200],filename6[200];
        FILE *outfile1,*outfile2,*outfile3,*outfile4,*outfile5,*outfile6;
        INTERFACE *intfc = front->interf;
        RECT_GRID *gr = &topological_grid(intfc);
        int *ppgmax = front->pp_grid->gmax;
        int pointsx, pointsy, pointsz;
        int gmax[3];
        gmax[0] = front->gmax[0];
        gmax[1] = front->gmax[1];
        gmax[2] = front->gmax[2];

        if(pp_mynode() == 0)
        {
            sprintf(dirname,"%s/domain-bdry",out_name);
           if (!create_directory(dirname,YES))
            {
                screen("Cannot create directory %s\n",dirname);
                clean_up(ERROR);
            }

            sprintf(filename1,"%s/surface1.vtk",dirname);
            outfile1 = fopen(filename1,"w");
            fprintf(outfile1,"# vtk DataFile Version 3.0\n");
            fprintf(outfile1,"Mesh of the topological grid\n");
            fprintf(outfile1,"ASCII\n");
            fprintf(outfile1,"DATASET STRUCTURED_GRID\n");
            pointsz = gmax[2] + 1;
            pointsy = gmax[1] + 1;
            pointsx = 1;
            totalpoints = pointsx*pointsy*pointsz;
            fprintf(outfile1, "DIMENSIONS %d %d %d\n", pointsz, pointsy, pointsx);
            fprintf(outfile1, "POINTS %d double\n", totalpoints);
            i = 0;
            for(j = 0; j <= pointsy-1; ++j)
            for(k = 0; k <= pointsz-1; ++k)
            {
                coord_x = gr->GL[0] + i*top_h[0];
                coord_y = gr->GL[1] + j*top_h[1];
                coord_z = gr->GL[2] + k*top_h[2];
                fprintf(outfile1, "%.16g %.16g %.16g\n", coord_x, coord_y, coord_z);
            }
            fclose(outfile1);

            sprintf(filename2,"%s/surface2.vtk",dirname);
            outfile2 = fopen(filename2,"w");
            fprintf(outfile2,"# vtk DataFile Version 3.0\n");
            fprintf(outfile2,"Mesh of the topological grid\n");
            fprintf(outfile2,"ASCII\n");
            fprintf(outfile2,"DATASET STRUCTURED_GRID\n");
            pointsz = gmax[2] + 1;
            pointsy = 1;
            pointsx = gmax[0] + 1;
            totalpoints = pointsx*pointsy*pointsz;
            fprintf(outfile2, "DIMENSIONS %d %d %d\n", pointsz, pointsy, pointsx);
            fprintf(outfile2, "POINTS %d double\n", totalpoints);
            j = 0;
            for(i = 0; i <= pointsx-1; ++i)
            for(k = 0; k <= pointsz-1; ++k)
            {
                coord_x = gr->GL[0] + i*top_h[0];
                coord_y = gr->GL[1] + j*top_h[1];
                coord_z = gr->GL[2] + k*top_h[2];
                fprintf(outfile2, "%.16g %.16g %.16g\n", coord_x, coord_y, coord_z);
            }
            fclose(outfile2);

            sprintf(filename3,"%s/surface3.vtk",dirname);
            outfile3 = fopen(filename3,"w");
            fprintf(outfile3,"# vtk DataFile Version 3.0\n");
            fprintf(outfile3,"Mesh of the topological grid\n");
            fprintf(outfile3,"ASCII\n");
            fprintf(outfile3,"DATASET STRUCTURED_GRID\n");
            pointsz = 1;
            pointsy = gmax[1] + 1;
            pointsx = gmax[0] + 1;
            totalpoints = pointsx*pointsy*pointsz;
            fprintf(outfile3, "DIMENSIONS %d %d %d\n", pointsz, pointsy, pointsx);
            fprintf(outfile3, "POINTS %d double\n", totalpoints);
            k = 0;
            for(i = 0; i <= pointsx-1; ++i)
            for(j = 0; j <= pointsy-1; ++j)
            {
                coord_x = gr->GL[0] + i*top_h[0];
                coord_y = gr->GL[1] + j*top_h[1];
                coord_z = gr->GL[2] + k*top_h[2];
                fprintf(outfile3, "%.16g %.16g %.16g\n", coord_x, coord_y, coord_z);
            }
            fclose(outfile3);

            sprintf(filename4,"%s/surface4.vtk",dirname);
            outfile4 = fopen(filename4,"w");
            fprintf(outfile4,"# vtk DataFile Version 3.0\n");
            fprintf(outfile4,"Mesh of the topological grid\n");
            fprintf(outfile4,"ASCII\n");
            fprintf(outfile4,"DATASET STRUCTURED_GRID\n");
            pointsz = gmax[2] + 1;
            pointsy = gmax[1] + 1;
            pointsx = 1;
            totalpoints = pointsx*pointsy*pointsz;
            fprintf(outfile4, "DIMENSIONS %d %d %d\n", pointsz, pointsy, pointsx);
            fprintf(outfile4, "POINTS %d double\n", totalpoints);
            i = gmax[0];
            for(j = 0; j <= pointsy-1; ++j)
            for(k = 0; k <= pointsz-1; ++k)
            {
                coord_x = gr->GL[0] + i*top_h[0];
                coord_y = gr->GL[1] + j*top_h[1];
                coord_z = gr->GL[2] + k*top_h[2];
                fprintf(outfile4, "%.16g %.16g %.16g\n", coord_x, coord_y, coord_z);
            }
            fclose(outfile4);

            sprintf(filename5,"%s/surface5.vtk",dirname);
            outfile5 = fopen(filename5,"w");
            fprintf(outfile5,"# vtk DataFile Version 3.0\n");
            fprintf(outfile5,"Mesh of the topological grid\n");
            fprintf(outfile5,"ASCII\n");
            fprintf(outfile5,"DATASET STRUCTURED_GRID\n");
            pointsz = gmax[2] + 1;
            pointsy = 1;
            pointsx = gmax[0] + 1;
            totalpoints = pointsx*pointsy*pointsz;
            fprintf(outfile5, "DIMENSIONS %d %d %d\n", pointsz, pointsy, pointsx);
            fprintf(outfile5, "POINTS %d double\n", totalpoints);
            j = gmax[1];
            for(i = 0; i <= pointsx-1; ++i)
            for(k = 0; k <= pointsz-1; ++k)
            {
                coord_x = gr->GL[0] + i*top_h[0];
                coord_y = gr->GL[1] + j*top_h[1];
                coord_z = gr->GL[2] + k*top_h[2];
                fprintf(outfile5, "%.16g %.16g %.16g\n", coord_x, coord_y, coord_z);
            }
            fclose(outfile5);

            sprintf(filename6,"%s/surface6.vtk",dirname);
            outfile6 = fopen(filename6,"w");
            fprintf(outfile6,"# vtk DataFile Version 3.0\n");
            fprintf(outfile6,"Mesh of the topological grid\n");
            fprintf(outfile6,"ASCII\n");
            fprintf(outfile6,"DATASET STRUCTURED_GRID\n");
            pointsz = 1;
            pointsy = gmax[1] + 1;
            pointsx = gmax[0] + 1;
            totalpoints = pointsx*pointsy*pointsz;
            fprintf(outfile6, "DIMENSIONS %d %d %d\n", pointsz, pointsy, pointsx);
            fprintf(outfile6, "POINTS %d double\n", totalpoints);
            k = gmax[2];
            for(i = 0; i <= pointsx-1; ++i)
            for(j = 0; j <= pointsy-1; ++j)
            {
                coord_x = gr->GL[0] + i*top_h[0];
                coord_y = gr->GL[1] + j*top_h[1];
                coord_z = gr->GL[2] + k*top_h[2];
                fprintf(outfile6, "%.16g %.16g %.16g\n", coord_x, coord_y, coord_z);
            }
            fclose(outfile6);
        }
} /* end printExpandedMesh_ascii */


void Incompress_Solver_Smooth_3D_Cartesian::computeAdvection(void)
{
	int i,j,k,l;
	int index,index00,index01,index10,index11,index20,index21,size;
	L_STATE state;
	COMPONENT comp;
	double speed;
	double *u, *v, *w;
	double u0,u00,u01,u10,u11,u20,u21;
	double v0,v00,v01,v10,v11,v20,v21;
	double w0,w00,w01,w10,w11,w20,w21;
	double crx_coords[MAXD];
	int icoords[MAXD];
	POINTER intfc_state;
	HYPER_SURF *hs;

	max_speed = 0.0;

	size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
	FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&w,size,sizeof(double));

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    u[index] = cell_center[index].m_state.m_U[0];
	    v[index] = cell_center[index].m_state.m_U[1];
	    w[index] = cell_center[index].m_state.m_U[2];
	}

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    index = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    if (comp == SOLID_COMP)
	    {
	    	cell_center[index].m_state.m_U[0] = 0.0;
	    	cell_center[index].m_state.m_U[1] = 0.0;
	    	cell_center[index].m_state.m_U[2] = 0.0;
		continue;
	    }
	    u0 = u[index];
	    v0 = v[index];
	    w0 = w[index];
	    index00 = d_index3d(i-1,j,k,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,WEST,comp,
			        &intfc_state,&hs,crx_coords,m_t_old) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u00 = getStateXvel(intfc_state);
		v00 = getStateYvel(intfc_state);
		w00 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u00 = u[index00];
		v00 = v[index00];
		w00 = w[index00];
	    }
	    index01 = d_index3d(i+1,j,k,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,EAST,comp,
			        &intfc_state,&hs,crx_coords,m_t_old) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u01 = getStateXvel(intfc_state);
		v01 = getStateYvel(intfc_state);
		w01 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u01 = u[index01];
		v01 = v[index01];
		w01 = w[index01];
	    }
	    index10 = d_index3d(i,j-1,k,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,SOUTH,comp,
			        &intfc_state,&hs,crx_coords,m_t_old) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u10 = getStateXvel(intfc_state);
		v10 = getStateYvel(intfc_state);
		w10 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u10 = u[index10];
		v10 = v[index10];
		w10 = w[index10];
	    }
	    index11 = d_index3d(i,j+1,k,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,NORTH,comp,
			        &intfc_state,&hs,crx_coords,m_t_old) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u11 = getStateXvel(intfc_state);
		v11 = getStateYvel(intfc_state);
		w11 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u11 = u[index11];
		v11 = v[index11];
		w11 = w[index11];
	    }
	    index20 = d_index3d(i,j,k-1,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,comp,
			        &intfc_state,&hs,crx_coords,m_t_old) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u20 = getStateXvel(intfc_state);
		v20 = getStateYvel(intfc_state);
		w20 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u20 = u[index20];
		v20 = v[index20];
		w20 = w[index20];
	    }
	    index21 = d_index3d(i,j,k+1,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,comp,
			        &intfc_state,&hs,crx_coords,m_t_old) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u21 = getStateXvel(intfc_state);
		v21 = getStateYvel(intfc_state);
		w21 = getStateZvel(intfc_state);
	    }
	    else
	    {
		u21 = u[index21];
		v21 = v[index21];
		w21 = w[index21];
	    }

	    cell_center[index].m_state.m_U[0] += -m_dt*(
				burger_flux(u00,u0,u01)/top_h[0] +
				linear_flux(v0,u10,u0,u11)/top_h[1] +
				linear_flux(w0,u20,u0,u21)/top_h[2]);

	    cell_center[index].m_state.m_U[1] += - m_dt*(
				linear_flux(u0,v00,v0,v01)/top_h[0] +
			 	burger_flux(v10,v0,v11)/top_h[1] +
				linear_flux(w0,v20,v0,v21)/top_h[2]);

	    cell_center[index].m_state.m_U[2] += - m_dt*(
				linear_flux(u0,w00,w0,w01)/top_h[0] +
				linear_flux(v0,w10,w0,w11)/top_h[1] +
			 	burger_flux(w20,w0,w21)/top_h[2]);

	    speed = fabs(cell_center[index].m_state.m_U[0]) +
		    fabs(cell_center[index].m_state.m_U[1]) +
		    fabs(cell_center[index].m_state.m_U[2]);
	    if (speed > max_speed)
		max_speed = speed;
	}
	for (l = 0; l < 3; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
	    	array[index] = cell_center[index].m_state.m_U[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].m_state.m_U[l] = array[index];
	    }
	}
	pp_global_max(&max_speed,1);
	FT_FreeThese(3,u,v,w);
}	/* end computeAdvection3d */



void Incompress_Solver_Smooth_3D_Cartesian::
	compDiffWithSmoothProperty_1st_coupled(void)
{
        COMPONENT comp;
        int index,index_nb[18],size;
        int I,I_nb[18];
        double coords[MAXD],crx_coords[MAXD];
	double coeff[18],mu[6],mu_edge[6],mu0,rho,rhs,U0_nb[18],U1_nb[18],U2_nb[18];
	int flag[6]; //denote whether this is dirichlet or neumann boundary
        L_STATE state;
        int i,j,k,l,nb,icoords[MAXD];
        INTERFACE *intfc = front->interf;
        double speed;
	double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 15, 15);
	// 7u + 4v + 4w for the first equation
	// 7v + 4u + 4w for the second equation
	// 7w + 4u + 4v for the third equation

	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index = d_index3d(i,j,k,top_gmax);
	//6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	//xy cut neighbours
	    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
	    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
	    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
	    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);
	//yz cut neighbours
	    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
	    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
	    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
	    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);
	//xz cut neighbours
	    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
	    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
	    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
	    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);


	//6 neighbours of the center cell
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];

	//xy cut neighbours
            I_nb[6] = ijk_to_I[i-1][j-1][k];
	    I_nb[7] = ijk_to_I[i+1][j-1][k];
	    I_nb[8] = ijk_to_I[i+1][j+1][k];
	    I_nb[9] = ijk_to_I[i-1][j+1][k];
	//yz cut neighbours
	    I_nb[10] = ijk_to_I[i][j-1][k-1];
	    I_nb[11] = ijk_to_I[i][j+1][k-1];
	    I_nb[12] = ijk_to_I[i][j+1][k+1];
	    I_nb[13] = ijk_to_I[i][j-1][k+1];
	//xz cut neighbours
	    I_nb[14] = ijk_to_I[i-1][j][k-1];
	    I_nb[15] = ijk_to_I[i+1][j][k-1];
	    I_nb[16] = ijk_to_I[i+1][j][k+1];
	    I_nb[17] = ijk_to_I[i-1][j][k+1];

	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    comp = top_comp[index];

	    mu0 = cell_center[index].m_state.m_mu;
	    rho = cell_center[index].m_state.m_rho;
            for (nb = 0; nb < 6; nb++)
            {
                if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords,m_t_new) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	       	{
		    flag[nb] = 1;
		    U0_nb[nb] = getStateVel[0](intfc_state);
		    U1_nb[nb] = getStateVel[1](intfc_state);
		    U2_nb[nb] = getStateVel[2](intfc_state);

		    if (wave_type(hs) == DIRICHLET_BOUNDARY ||
			wave_type(hs) == NEUMANN_BOUNDARY)
		    {
			mu[nb] = mu0;
			mu_edge[nb] = mu0;
		    }

		    else
		    {
			    mu[nb] = 0.5*(mu0 + cell_center[index_nb[nb]].m_state.m_mu);
			    mu_edge[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		    }
		}
                    else
		    {
			flag[nb] = 0;
                    	U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
			U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
			U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
			mu_edge[nb] = 0.5*(mu0 + cell_center[index_nb[nb]].m_state.m_mu);
			mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		    }
	    }

	    //xy cut neighbour cells

	    if (flag[0] == 1 && flag[2] == 1) //neighbourcell 6
	    {
		U0_nb[6] = 0.5*(U0_nb[0] + U0_nb[2]);
		U1_nb[6] = 0.5*(U1_nb[0] + U1_nb[2]);
		U2_nb[6] = 0.5*(U2_nb[0] + U2_nb[2]);
	    }
	    else
	    {
		U0_nb[6] = cell_center[index_nb[6]].m_state.m_U[0];
		U1_nb[6] = cell_center[index_nb[6]].m_state.m_U[1];
		U2_nb[6] = cell_center[index_nb[6]].m_state.m_U[2];
	    }
	    if (flag[1] == 1 && flag[2] == 1) //neighbourcell 7
	    {
		U0_nb[7] = 0.5*(U0_nb[1] + U0_nb[2]);
		U1_nb[7] = 0.5*(U1_nb[1] + U1_nb[2]);
		U2_nb[7] = 0.5*(U2_nb[1] + U2_nb[2]);
	    }
	    else
	    {
		U0_nb[7] = cell_center[index_nb[7]].m_state.m_U[0];
		U1_nb[7] = cell_center[index_nb[7]].m_state.m_U[1];
		U2_nb[7] = cell_center[index_nb[7]].m_state.m_U[2];
	    }
	    if (flag[1] == 1 && flag[3] == 1) //neighbourcell 8
	    {
		U0_nb[8] = 0.5*(U0_nb[1] + U0_nb[3]);
		U1_nb[8] = 0.5*(U1_nb[1] + U1_nb[3]);
		U2_nb[8] = 0.5*(U2_nb[1] + U2_nb[3]);
	    }
	    else
	    {
		U0_nb[8] = cell_center[index_nb[8]].m_state.m_U[0];
		U1_nb[8] = cell_center[index_nb[8]].m_state.m_U[1];
		U2_nb[8] = cell_center[index_nb[8]].m_state.m_U[2];
	    }
	    if (flag[0] == 1 && flag[3] == 1) //neighbourcell 9
	    {
		U0_nb[9] = 0.5*(U0_nb[0] + U0_nb[3]);
		U1_nb[9] = 0.5*(U1_nb[0] + U1_nb[3]);
		U2_nb[9] = 0.5*(U2_nb[0] + U2_nb[3]);
	    }
	    else
	    {
		U0_nb[9] = cell_center[index_nb[9]].m_state.m_U[0];
		U1_nb[9] = cell_center[index_nb[9]].m_state.m_U[1];
		U2_nb[9] = cell_center[index_nb[9]].m_state.m_U[2];
	    }
	    //yz cut neighbours
	    if (flag[2] == 1 && flag[4] == 1) //neighbourcell 10
	    {
		U0_nb[10] = 0.5*(U0_nb[2] + U0_nb[4]);
		U1_nb[10] = 0.5*(U1_nb[2] + U1_nb[4]);
		U2_nb[10] = 0.5*(U2_nb[2] + U2_nb[4]);
	    }
	    else
	    {
		U0_nb[10] = cell_center[index_nb[10]].m_state.m_U[0];
		U1_nb[10] = cell_center[index_nb[10]].m_state.m_U[1];
		U2_nb[10] = cell_center[index_nb[10]].m_state.m_U[2];
	    }
	    if (flag[3] == 1 && flag[4] == 1) //neighbourcell 11
	    {
		U0_nb[11] = 0.5*(U0_nb[3] + U0_nb[4]);
		U1_nb[11] = 0.5*(U1_nb[3] + U1_nb[4]);
		U2_nb[11] = 0.5*(U2_nb[3] + U2_nb[4]);
	    }
	    else
	    {
		U0_nb[11] = cell_center[index_nb[11]].m_state.m_U[0];
		U1_nb[11] = cell_center[index_nb[11]].m_state.m_U[1];
		U2_nb[11] = cell_center[index_nb[11]].m_state.m_U[2];
	    }
	    if (flag[3] == 1 && flag[5] == 1) //neighbourcell 12
	    {
		U0_nb[12] = 0.5*(U0_nb[3] + U0_nb[5]);
		U1_nb[12] = 0.5*(U1_nb[3] + U1_nb[5]);
		U2_nb[12] = 0.5*(U2_nb[3] + U2_nb[5]);
	    }
	    else
	    {
		U0_nb[12] = cell_center[index_nb[12]].m_state.m_U[0];
		U1_nb[12] = cell_center[index_nb[12]].m_state.m_U[1];
		U2_nb[12] = cell_center[index_nb[12]].m_state.m_U[2];
	    }
	    if (flag[2] == 1 && flag[5] == 1) //neighbourcell 13
	    {
		U0_nb[13] = 0.5*(U0_nb[2] + U0_nb[5]);
		U1_nb[13] = 0.5*(U1_nb[2] + U1_nb[5]);
		U2_nb[13] = 0.5*(U2_nb[2] + U2_nb[5]);
	    }
	    else
	    {
		U0_nb[13] = cell_center[index_nb[13]].m_state.m_U[0];
		U1_nb[13] = cell_center[index_nb[13]].m_state.m_U[1];
		U2_nb[13] = cell_center[index_nb[13]].m_state.m_U[2];
	    }
	    //xz cut neighbours
	    if (flag[0] == 1 && flag[4] == 1) //neighbourcell 14
	    {
		U0_nb[14] = 0.5*(U0_nb[0] + U0_nb[4]);
		U1_nb[14] = 0.5*(U1_nb[0] + U1_nb[4]);
		U2_nb[14] = 0.5*(U2_nb[0] + U2_nb[4]);
	    }
	    else
	    {
		U0_nb[14] = cell_center[index_nb[14]].m_state.m_U[0];
		U1_nb[14] = cell_center[index_nb[14]].m_state.m_U[1];
		U2_nb[14] = cell_center[index_nb[14]].m_state.m_U[2];
	    }
	    if (flag[1] == 1 && flag[4] == 1) //neighbourcell 15
	    {
		U0_nb[15] = 0.5*(U0_nb[1] + U0_nb[4]);
		U1_nb[15] = 0.5*(U1_nb[1] + U1_nb[4]);
		U2_nb[15] = 0.5*(U2_nb[1] + U2_nb[4]);
	    }
	    else
	    {
		U0_nb[15] = cell_center[index_nb[15]].m_state.m_U[0];
		U1_nb[15] = cell_center[index_nb[15]].m_state.m_U[1];
		U2_nb[15] = cell_center[index_nb[15]].m_state.m_U[2];
	    }
	    if (flag[1] == 1 && flag[5] == 1) //neighbourcell 16
	    {
		U0_nb[16] = 0.5*(U0_nb[1] + U0_nb[5]);
		U1_nb[16] = 0.5*(U1_nb[1] + U1_nb[5]);
		U2_nb[16] = 0.5*(U2_nb[1] + U2_nb[5]);
	    }
	    else
	    {
		U0_nb[16] = cell_center[index_nb[16]].m_state.m_U[0];
		U1_nb[16] = cell_center[index_nb[16]].m_state.m_U[1];
		U2_nb[16] = cell_center[index_nb[16]].m_state.m_U[2];
	    }
	    if (flag[0] == 1 && flag[5] == 1) //neighbourcell 17
	    {
		U0_nb[17] = 0.5*(U0_nb[0] + U0_nb[5]);
		U1_nb[17] = 0.5*(U1_nb[0] + U1_nb[5]);
		U2_nb[17] = 0.5*(U2_nb[0] + U2_nb[5]);
	    }
	    else
	    {
		U0_nb[17] = cell_center[index_nb[17]].m_state.m_U[0];
		U1_nb[17] = cell_center[index_nb[17]].m_state.m_U[1];
		U2_nb[17] = cell_center[index_nb[17]].m_state.m_U[2];
	    }



            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);


    //Setting the coeffecients for the first equation


	    coeff[0] = m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff[1] = m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff[4] = 0.5*m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    coeff[6] =  0.5*m_dt/rho * mu[2]/(4*top_h[0]*top_h[1]);
	    coeff[7] = -0.5*m_dt/rho * mu[2]/(4*top_h[0]*top_h[1]);
	    coeff[8] =  0.5*m_dt/rho * mu[3]/(4*top_h[0]*top_h[1]);
	    coeff[9] = -0.5*m_dt/rho * mu[3]/(4*top_h[0]*top_h[1]);

	    coeff[14] =  0.5*m_dt/rho * mu[4]/(4*top_h[0]*top_h[2]);
	    coeff[15] = -0.5*m_dt/rho * mu[4]/(4*top_h[0]*top_h[2]);
	    coeff[16] =  0.5*m_dt/rho * mu[5]/(4*top_h[0]*top_h[2]);
	    coeff[17] = -0.5*m_dt/rho * mu[5]/(4*top_h[0]*top_h[2]);

	    solver.Set_A(I*3,I*3,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*cell_center[index].m_state.m_U[0];

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3,I_nb[nb]*3,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }
	    for (nb = 6; nb < 10; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3,I_nb[nb]*3+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }
	    for (nb = 14; nb < 18; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3,I_nb[nb]*3+2,-coeff[nb]);
		    rhs += coeff[nb]*U2_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U2_nb[nb];
	    }

	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;

	    solver.Set_b(I*3, rhs);

	    /************************************************************************/


    //Setting the coeffecients for the second equation


	    coeff[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff[2] = m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff[3] = m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff[4] = 0.5*m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff[5] = 0.5*m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    coeff[6] =  0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[1]);
	    coeff[7] = -0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[1]);
	    coeff[8] =  0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[1]);
	    coeff[9] = -0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[1]);

	    coeff[10] =  0.5*m_dt/rho * mu[4]/(4*top_h[1]*top_h[2]);
	    coeff[11] = -0.5*m_dt/rho * mu[4]/(4*top_h[1]*top_h[2]);
	    coeff[12] =  0.5*m_dt/rho * mu[5]/(4*top_h[1]*top_h[2]);
	    coeff[13] = -0.5*m_dt/rho * mu[5]/(4*top_h[1]*top_h[2]);

	    solver.Set_A(I*3+1,I*3+1,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*cell_center[index].m_state.m_U[1];

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+1,I_nb[nb]*3+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }
	    for (nb = 6; nb < 10; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+1,I_nb[nb]*3,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }
	    for (nb = 10; nb < 14; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+1,I_nb[nb]*3+2,-coeff[nb]);
		    rhs += coeff[nb]*U2_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U2_nb[nb];
	    }

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;

	    solver.Set_b(I*3+1, rhs);

	    /************************************************************************/

    //Setting the coeffecients for the third equation


	    coeff[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff[4] = m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff[5] = m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    coeff[14] =  0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[2]);
	    coeff[15] = -0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[2]);
	    coeff[16] =  0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[2]);
	    coeff[17] = -0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[2]);

	    coeff[10] =  0.5*m_dt/rho * mu[2]/(4*top_h[1]*top_h[2]);
	    coeff[11] = -0.5*m_dt/rho * mu[3]/(4*top_h[1]*top_h[2]);
	    coeff[12] =  0.5*m_dt/rho * mu[3]/(4*top_h[1]*top_h[2]);
	    coeff[13] = -0.5*m_dt/rho * mu[2]/(4*top_h[1]*top_h[2]);

	    solver.Set_A(I*3+2,I*3+2,1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5]);
	    rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*cell_center[index].m_state.m_U[2];

	    for (nb = 0; nb < 6; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+2,I_nb[nb]*3+2,-coeff[nb]);
		    rhs += coeff[nb]*U2_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U2_nb[nb];
	    }
	    for (nb = 14; nb < 18; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+2,I_nb[nb]*3,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }
	    for (nb = 10; nb < 14; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*3+2,I_nb[nb]*3+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }

	    rhs += m_dt*state.m_U[2];
	    rhs += m_dt*cell_center[index].m_state.f_surf[2];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[2]/rho;

	    solver.Set_b(I*3+2, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-10);

	start_clock("Before Petsc Solve");
        solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);

	//if (rel_residual > 1)
	//{
	//    printf("\nThe solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
	//    solver.Reset_x();
	//    solver.Solve_GMRES();
	//    solver.GetNumIterations(&num_iter);
	//    solver.GetFinalRelativeResidualNorm(&rel_residual);
	//}

	stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_1st_coupled: "
                        "num_iter = %d, rel_residual = %le. \n",
                        num_iter,rel_residual);


	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
                cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
		cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
			fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                    max_speed = speed;
            }
            else
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
		cell_center[index].m_state.m_U[2] = 0.0;
            }
        }
        for (l = 0; l < 3; ++l)
        {
	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }
        pp_global_max(&max_speed,1);

        FT_FreeThese(1,x);
}       /* end compDiffWithSmoothProperty3d_1st_coupled */


void Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_2nd_coupled(void)
{
        COMPONENT comp;
        int index,index_nb[18],size;
        int I,I_nb[18];
        double coords[MAXD],corner_coords[MAXD],crx_coords[MAXD];
	double coeff0[6],coeff1[6],coeff2[6],coeff_temp;
	double mu[6],mu_edge[6],mu0,rho,rhs;
	double U0_nb[18],U1_nb[18],U2_nb[18],U0_center,U1_center,U2_center;
        double U0_nb_new[6],U1_nb_new[6],U2_nb_new[6];
	int flag[6]; //denote whether this is dirichlet or neumann boundary
        L_STATE state,corner_state,corner_state_new;
        int i,j,k,l,nb,icoords[MAXD];
        INTERFACE *intfc = front->interf;
        double speed;
	double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,3*size,sizeof(double));

        PETSc solver;
        solver.Create(3*ilower, 3*iupper-1, 30, 30);
	// 7u + 9v + 9w for the first equation
	// 7v + 9u + 9w for the second equation
	// 7w + 9u + 9v for the third equation

	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index = d_index3d(i,j,k,top_gmax);
	//6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	//xy cut neighbours
	    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
	    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
	    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
	    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);
	//yz cut neighbours
	    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
	    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
	    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
	    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);
	//xz cut neighbours
	    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
	    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
	    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
	    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);


	//6 neighbours of the center cell
            I_nb[0] = ijk_to_I[i-1][j][k];
            I_nb[1] = ijk_to_I[i+1][j][k];
            I_nb[2] = ijk_to_I[i][j-1][k];
            I_nb[3] = ijk_to_I[i][j+1][k];
            I_nb[4] = ijk_to_I[i][j][k-1];
            I_nb[5] = ijk_to_I[i][j][k+1];

	//xy cut neighbours
            I_nb[6] = ijk_to_I[i-1][j-1][k];
	    I_nb[7] = ijk_to_I[i+1][j-1][k];
	    I_nb[8] = ijk_to_I[i+1][j+1][k];
	    I_nb[9] = ijk_to_I[i-1][j+1][k];
	//yz cut neighbours
	    I_nb[10] = ijk_to_I[i][j-1][k-1];
	    I_nb[11] = ijk_to_I[i][j+1][k-1];
	    I_nb[12] = ijk_to_I[i][j+1][k+1];
	    I_nb[13] = ijk_to_I[i][j-1][k+1];
	//xz cut neighbours
	    I_nb[14] = ijk_to_I[i-1][j][k-1];
	    I_nb[15] = ijk_to_I[i+1][j][k-1];
	    I_nb[16] = ijk_to_I[i+1][j][k+1];
	    I_nb[17] = ijk_to_I[i-1][j][k+1];

	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    comp = top_comp[index];

	    mu0 = cell_center[index].m_state.m_mu;
	    rho = cell_center[index].m_state.m_rho;
	    U0_center = cell_center[index].m_state.m_U[0];
	    U1_center = cell_center[index].m_state.m_U[1];
	    U2_center = cell_center[index].m_state.m_U[2];

            for (nb = 0; nb < 6; nb++)
            {
                if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	       	{
		    flag[nb] = 1;
                    // old boundary condition
		    U0_nb[nb] = getStateVel[0](intfc_state);
		    U1_nb[nb] = getStateVel[1](intfc_state);
		    U2_nb[nb] = getStateVel[2](intfc_state);
                    // new boundary condition
                    FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords,m_t_new);
                    U0_nb_new[nb] = getStateVel[0](intfc_state);
                    U1_nb_new[nb] = getStateVel[1](intfc_state);
                    U2_nb_new[nb] = getStateVel[2](intfc_state);

		    if (wave_type(hs) == DIRICHLET_BOUNDARY ||
			wave_type(hs) == NEUMANN_BOUNDARY)
		    {
			mu[nb] = mu0;
			mu_edge[nb] = mu0;
		    }
		    else
		    {
			mu[nb] = 0.5*(mu0 + cell_center[index_nb[nb]].m_state.m_mu);
			mu_edge[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		    }
		}
		else
		{
		    flag[nb] = 0;
                    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
		    mu_edge[nb] = 0.5*(mu0 + cell_center[index_nb[nb]].m_state.m_mu);
		    mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		}
	    }

	    for (nb = 6; nb < 18; nb++) //cut corner values, interior
	    {
		if (I_nb[nb] != -1)
		{
		    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    U2_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[2];
		}
	    }


	    // source term
            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);


            /************************************************************************/
	    //Setting the coeffecients for U0 in the first equation

	    coeff0[0] = m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff0[1] = m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff0[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff0[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff0[4] = 0.5*m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff0[5] = 0.5*m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    solver.Add_A(I*3,I*3,1.0);
	    rhs = U0_center;

	    for (nb = 0; nb < 6; nb++)
	    {
		//if (I_nb[nb] != -1)
                if (flag[nb] == 0)
		{
		    solver.Add_A(I*3,I_nb[nb]*3,-coeff0[nb]);
		    rhs += coeff0[nb]*U0_nb[nb];
		}
                else if (flag[nb] == 1)
		{
		    coeff0[nb] = 2.0*coeff0[nb];
		    //rhs += 2.0*coeff0[nb]*U0_nb[nb];
                    rhs += coeff0[nb]*(U0_nb_new[nb] + U0_nb[nb]);
		}
		solver.Add_A(I*3,I*3,coeff0[nb]);
		rhs -= coeff0[nb]*U0_center;
	    }

	    //set the coefficients for U1 in the first equation (mu*v_x)_y
	    //traverse the four corners
	    //(i-1/2,j-1/2,k),(i+1/2,j-1/2,k),(i+1/2,j+1/2,k),(i-1/2,j+1/2,k)

	    //corner (i-1/2,j-1/2,k)

	    coeff_temp = mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
	    if (flag[0] == 1 && flag[2] == 0)
		rhs += coeff_temp*U1_nb[0];
	    else if(flag[0] == 0 && flag[2] == 1)
		rhs += coeff_temp*U1_nb[2];
	    else if(flag[0] == 1 && flag[2] == 1)
		rhs += coeff_temp*U1_nb[0];
	    else {
		rhs += coeff_temp*(U1_nb[0]+U1_nb[2]+U1_nb[6]+U1_center)/8.0;

		solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[0]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[2]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[6]*3+1, -coeff_temp/8.0);
	    }  */
            if (flag[0] == 0 && flag[2] == 0)
            {
                rhs += coeff_temp*(U1_nb[0]+U1_nb[2]+U1_nb[6]+U1_center)/8.0;
                solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[0]*3+1, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[2]*3+1, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[6]*3+1, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                corner_coords[2] = coords[2] - 0.0*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }

	    //corner (i+1/2,j-1/2,k)

	    coeff_temp = -mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
	    if (flag[2] == 1 && flag[1] == 0)
		rhs += coeff_temp*U1_nb[2];
	    else if(flag[2] == 0 && flag[1] == 1)
		rhs += coeff_temp*U1_nb[1];
	    else if(flag[2] == 1 && flag[1] == 1)
		rhs += coeff_temp*U1_nb[2];
	    else {
		rhs += coeff_temp*(U1_nb[2]+U1_nb[1]+U1_nb[7]+U1_center)/8.0;

		solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[2]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[1]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[7]*3+1, -coeff_temp/8.0);
	    }  */
            if (flag[2] == 0 && flag[1] == 0)
            {
                rhs += coeff_temp*(U1_nb[2]+U1_nb[1]+U1_nb[7]+U1_center)/8.0;
                solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[2]*3+1, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[1]*3+1, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[7]*3+1, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                corner_coords[2] = coords[2] - 0.0*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }

	    //corner (i+1/2,j+1/2,k)

	    coeff_temp = mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
	    if (flag[1] == 1 && flag[3] == 0)
		rhs += coeff_temp*U1_nb[1];
	    else if(flag[1] == 0 && flag[3] == 1)
		rhs += coeff_temp*U1_nb[3];
	    else if(flag[1] == 1 && flag[3] == 1)
		rhs += coeff_temp*U1_nb[1];
	    else {
		rhs += coeff_temp*(U1_nb[1]+U1_nb[3]+U1_nb[8]+U1_center)/8.0;

		solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[1]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[3]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[8]*3+1, -coeff_temp/8.0);
	    } */

            if (flag[1] == 0 && flag[3] == 0)
            {
                rhs += coeff_temp*(U1_nb[1]+U1_nb[3]+U1_nb[8]+U1_center)/8.0;
                solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[1]*3+1, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[3]*3+1, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[8]*3+1, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                corner_coords[2] = coords[2] - 0.0*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }

	    //corner (i-1/2,j+1/2,k)

	    coeff_temp = -mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
	    if (flag[3] == 1 && flag[0] == 0)
		rhs += coeff_temp*U1_nb[3];
	    else if(flag[3] == 0 && flag[0] == 1)
		rhs += coeff_temp*U1_nb[0];
	    else if(flag[3] == 1 && flag[0] == 1)
		rhs += coeff_temp*U1_nb[3];
	    else {
		rhs += coeff_temp*(U1_nb[3]+U1_nb[0]+U1_nb[9]+U1_center)/8.0;

		solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[3]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[0]*3+1, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[9]*3+1, -coeff_temp/8.0);
	    } */

            if (flag[3] == 0 && flag[0] == 0)
            {
                rhs += coeff_temp*(U1_nb[3]+U1_nb[0]+U1_nb[9]+U1_center)/8.0;
                solver.Add_A(I*3,I*3+1,       -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[3]*3+1, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[0]*3+1, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[9]*3+1, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                corner_coords[2] = coords[2] - 0.0*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }


	    //set the coefficients for U2 in the first equation (mu*w_x)_z
	    //traverse the four corners
	    //(i-1/2,j,k-1/2),(i+1/2,j,k-1/2),(i+1/2,j,k+1/2),(i-1/2,j,k+1/2)

	    //corner (i-1/2,j,k-1/2)

	    coeff_temp = mu_edge[4]*m_dt/(top_h[0]*top_h[2])/rho;
            /*
	    if (flag[0] == 1 && flag[4] == 0)
		rhs += coeff_temp*U2_nb[0];
	    else if(flag[0] == 0 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[4];
	    else if(flag[0] == 1 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[0];
	    else {
		rhs += coeff_temp*(U2_nb[0]+U2_nb[4]+U2_nb[14]+U2_center)/8.0;

		solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[0]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[4]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[14]*3+2, -coeff_temp/8.0);
	    } */

            if (flag[0] == 0 && flag[4] == 0)
            {
                rhs += coeff_temp*(U2_nb[0]+U2_nb[4]+U2_nb[14]+U2_center)/8.0;

                solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[0]*3+2, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[4]*3+2, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[14]*3+2, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.0*top_h[1];
                corner_coords[2] = coords[2] - 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[2]+
                    corner_state_new.m_U[2]);
            }

	    //corner (i+1/2,j,k-1/2)

	    coeff_temp = -mu_edge[4]*m_dt/(top_h[0]*top_h[2])/rho;
            /*
	    if (flag[1] == 1 && flag[4] == 0)
		rhs += coeff_temp*U2_nb[1];
	    else if(flag[1] == 0 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[4];
	    else if(flag[1] == 1 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[1];
	    else {
		rhs += coeff_temp*(U2_nb[1]+U2_nb[4]+U2_nb[15]+U2_center)/8.0;

		solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[1]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[4]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[15]*3+2, -coeff_temp/8.0);
	    } */

            if (flag[1] == 0 && flag[4] == 0)
            {
                rhs += coeff_temp*(U2_nb[1]+U2_nb[4]+U2_nb[15]+U2_center)/8.0;

                solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[1]*3+2, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[4]*3+2, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[15]*3+2, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.0*top_h[1];
                corner_coords[2] = coords[2] - 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[2]+
                    corner_state_new.m_U[2]);
            }

	    //corner (i+1/2,j,k+1/2)

	    coeff_temp = mu_edge[5]*m_dt/(top_h[0]*top_h[2])/rho;
            /*
	    if (flag[1] == 1 && flag[5] == 0)
		rhs += coeff_temp*U2_nb[1];
	    else if(flag[1] == 0 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[5];
	    else if(flag[1] == 1 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[1];
	    else {
		rhs += coeff_temp*(U2_nb[1]+U2_nb[5]+U2_nb[16]+U2_center)/8.0;

		solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[1]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[5]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[16]*3+2, -coeff_temp/8.0);
	    } */

            if (flag[1] == 0 && flag[5] == 0)
            {
                rhs += coeff_temp*(U2_nb[1]+U2_nb[5]+U2_nb[16]+U2_center)/8.0;

                solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[1]*3+2, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[5]*3+2, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[16]*3+2, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.0*top_h[1];
                corner_coords[2] = coords[2] + 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[2]+
                    corner_state_new.m_U[2]);
            }

	    //corner (i-1/2,j,k+1/2)

	    coeff_temp = -mu_edge[5]*m_dt/(top_h[0]*top_h[2])/rho;
            /*
	    if (flag[5] == 1 && flag[0] == 0)
		rhs += coeff_temp*U2_nb[5];
	    else if(flag[5] == 0 && flag[0] == 1)
		rhs += coeff_temp*U2_nb[0];
	    else if(flag[5] == 1 && flag[0] == 1)
		rhs += coeff_temp*U2_nb[5];
	    else {
		rhs += coeff_temp*(U2_nb[5]+U2_nb[0]+U2_nb[17]+U2_center)/8.0;

		solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[5]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[0]*3+2, -coeff_temp/8.0);
		solver.Add_A(I*3,I_nb[17]*3+2, -coeff_temp/8.0);
	    } */

            if (flag[5] == 0 && flag[0] == 0)
            {
                rhs += coeff_temp*(U2_nb[5]+U2_nb[0]+U2_nb[17]+U2_center)/8.0;

                solver.Add_A(I*3,I*3+2,       -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[5]*3+2, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[0]*3+2, -coeff_temp/8.0);
                solver.Add_A(I*3,I_nb[17]*3+2, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.0*top_h[1];
                corner_coords[2] = coords[2] + 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[2]+
                    corner_state_new.m_U[2]);
            }

	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[0];

	    solver.Add_b(I*3, rhs);

	    /************************************************************************/


	    //Setting the coeffecients for U1 in the second equation


	    coeff1[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff1[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff1[2] = m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff1[3] = m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff1[4] = 0.5*m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff1[5] = 0.5*m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    solver.Add_A(I*3+1,I*3+1,1.0);
	    rhs = U1_center;

	    for (nb = 0; nb < 6; nb++)
	    {
		//if (I_nb[nb] != -1)
                if (flag[nb] == 0)
		{
		    solver.Add_A(I*3+1,I_nb[nb]*3+1,-coeff1[nb]);
		    rhs += coeff1[nb]*U1_nb[nb];
		}
                else if (flag[nb] == 1)
		{
		    coeff1[nb] = 2.0*coeff1[nb];
		    rhs += coeff1[nb]*(U1_nb[nb] + U1_nb_new[nb]);
		}
		solver.Add_A(I*3+1,I*3+1,coeff1[nb]);
		rhs -= coeff1[nb]*U1_center;
	    }

	    //set the coefficients for U0 in the second equation (mu*u_y)_x
	    //traverse the four corners
	    //(i-1/2,j-1/2,k),(i+1/2,j-1/2,k),(i+1/2,j+1/2,k),(i-1/2,j+1/2,k)

	    //corner (i-1/2,j-1/2,k)

	    coeff_temp = mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
	    if (flag[0] == 1 && flag[2] == 0)
		rhs += coeff_temp*U0_nb[0];
	    else if (flag[0] == 0 && flag[2] == 1)
		rhs += coeff_temp*U0_nb[2];
	    else if (flag[0] == 1 && flag[2] == 1)
		rhs += coeff_temp*U0_nb[0];
	    else {
		rhs += coeff_temp*(U0_nb[0]+U0_nb[2]+U0_nb[6]+U0_center)/8.0;

		solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[0]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[2]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[6]*3,  -coeff_temp/8.0);
	    } */

            if (flag[0] == 0 && flag[2] == 0)
            {
                rhs += coeff_temp*(U0_nb[0]+U0_nb[2]+U0_nb[6]+U0_center)/8.0;

                solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[0]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[2]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[6]*3,  -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                corner_coords[2] = coords[2] - 0.0*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }

	    //corner (i+1/2,j-1/2,k)

	    coeff_temp = -mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
	    if (flag[1] == 1 && flag[2] == 0)
		rhs += coeff_temp*U0_nb[1];
	    else if (flag[1] == 0 && flag[2] == 1)
		rhs += coeff_temp*U0_nb[2];
	    else if (flag[1] == 1 && flag[2] == 1)
		rhs += coeff_temp*U0_nb[1];
	    else {
		rhs += coeff_temp*(U0_nb[1]+U0_nb[2]+U0_nb[7]+U0_center)/8.0;

		solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[1]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[2]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[7]*3,  -coeff_temp/8.0);
	    } */

            if (flag[1] == 0 && flag[2] == 0)
            {
                rhs += coeff_temp*(U0_nb[1]+U0_nb[2]+U0_nb[7]+U0_center)/8.0;

                solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[1]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[2]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[7]*3,  -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                corner_coords[2] = coords[2] - 0.0*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }

	    //corner (i+1/2,j+1/2,k)

	    coeff_temp = mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
	    if (flag[1] == 1 && flag[3] == 0)
		rhs += coeff_temp*U0_nb[1];
	    else if (flag[1] == 0 && flag[3] == 1)
		rhs += coeff_temp*U0_nb[3];
	    else if (flag[1] == 1 && flag[3] == 1)
		rhs += coeff_temp*U0_nb[1];
	    else {
		rhs += coeff_temp*(U0_nb[1]+U0_nb[3]+U0_nb[8]+U0_center)/8.0;

		solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[1]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[3]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[8]*3,  -coeff_temp/8.0);
	    } */

            if (flag[1] == 0 && flag[3] == 0)
            {
                rhs += coeff_temp*(U0_nb[1]+U0_nb[3]+U0_nb[8]+U0_center)/8.0;

                solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[1]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[3]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[8]*3,  -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                corner_coords[2] = coords[2] - 0.0*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }

	    //corner (i-1/2,j+1/2,k)

	    coeff_temp = -mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
	    if (flag[0] == 1 && flag[3] == 0)
		rhs += coeff_temp*U0_nb[0];
	    else if (flag[0] == 0 && flag[3] == 1)
		rhs += coeff_temp*U0_nb[3];
	    else if (flag[0] == 1 && flag[3] == 1)
		rhs += coeff_temp*U0_nb[0];
	    else {
		rhs += coeff_temp*(U0_nb[0]+U0_nb[3]+U0_nb[9]+U0_center)/8.0;

		solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[0]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[3]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[9]*3,  -coeff_temp/8.0);
	    } */

            if (flag[0] == 0 && flag[3] == 0)
            {
                rhs += coeff_temp*(U0_nb[0]+U0_nb[3]+U0_nb[9]+U0_center)/8.0;

                solver.Add_A(I*3+1,I*3,        -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[0]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[3]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[9]*3,  -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                corner_coords[2] = coords[2] - 0.0*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }


	    //set the coefficients for U2 in the second equation (mu*w_y)_z
	    //traverse the four corners
	    //(i,j-1/2,k-1/2),(i,j+1/2,k-1/2),(i,j+1/2,k+1/2),(i,j-1/2,k+1/2)

	    //corner (i,j-1/2,k-1/2)

	    coeff_temp = mu_edge[4]*m_dt/(top_h[1]*top_h[2])/rho;
            /*
	    if (flag[2] == 1 && flag[4] == 0)
		rhs += coeff_temp*U2_nb[2];
	    else if (flag[2] == 0 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[4];
	    else if (flag[2] == 1 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[2];
	    else {
		rhs += coeff_temp*(U2_nb[2]+U2_nb[4]+U2_nb[10]+U2_center)/8.0;

		solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[2]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[4]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[10]*3+2,  -coeff_temp/8.0);
	    } */

            if (flag[2] == 0 && flag[4] == 0)
            {
                rhs += coeff_temp*(U2_nb[2]+U2_nb[4]+U2_nb[10]+U2_center)/8.0;

                solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[2]*3+2,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[4]*3+2,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[10]*3+2,  -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.0*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                corner_coords[2] = coords[2] - 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[2]+
                    corner_state_new.m_U[2]);
            }

	    //corner (i,j+1/2,k-1/2)

	    coeff_temp = -mu_edge[4]*m_dt/(top_h[1]*top_h[2])/rho;
            /*
	    if (flag[3] == 1 && flag[4] == 0)
		rhs += coeff_temp*U2_nb[3];
	    else if (flag[3] == 0 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[4];
	    else if (flag[3] == 1 && flag[4] == 1)
		rhs += coeff_temp*U2_nb[3];
	    else {
		rhs += coeff_temp*(U2_nb[3]+U2_nb[4]+U2_nb[11]+U2_center)/8.0;

		solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[3]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[4]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[11]*3+2,  -coeff_temp/8.0);
	    } */

            if (flag[3] == 0 && flag[4] == 0)
            {
                rhs += coeff_temp*(U2_nb[3]+U2_nb[4]+U2_nb[11]+U2_center)/8.0;

                solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[3]*3+2,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[4]*3+2,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[11]*3+2,  -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.0*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                corner_coords[2] = coords[2] - 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[2]+
                    corner_state_new.m_U[2]);
            }

	    //corner (i,j+1/2,k+1/2)

	    coeff_temp = mu_edge[5]*m_dt/(top_h[1]*top_h[2])/rho;
            /*
	    if (flag[3] == 1 && flag[5] == 0)
		rhs += coeff_temp*U2_nb[3];
	    else if (flag[3] == 0 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[5];
	    else if (flag[3] == 1 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[3];
	    else {
		rhs += coeff_temp*(U2_nb[3]+U2_nb[5]+U2_nb[12]+U2_center)/8.0;

		solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[3]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[5]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[12]*3+2,  -coeff_temp/8.0);
	    } */

            if (flag[3] == 0 && flag[5] == 0)
            {
                rhs += coeff_temp*(U2_nb[3]+U2_nb[5]+U2_nb[12]+U2_center)/8.0;

                solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[3]*3+2,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[5]*3+2,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[12]*3+2,  -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.0*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                corner_coords[2] = coords[2] + 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[2]+
                    corner_state_new.m_U[2]);
            }

	    //corner (i,j-1/2,k+1/2)

	    coeff_temp = -mu_edge[5]*m_dt/(top_h[1]*top_h[2])/rho;
            /*
	    if (flag[2] == 1 && flag[5] == 0)
		rhs += coeff_temp*U2_nb[2];
	    else if (flag[2] == 0 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[5];
	    else if (flag[2] == 1 && flag[5] == 1)
		rhs += coeff_temp*U2_nb[2];
	    else {
		rhs += coeff_temp*(U2_nb[2]+U2_nb[5]+U2_nb[13]+U2_center)/8.0;

		solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[2]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[5]*3+2,  -coeff_temp/8.0);
		solver.Add_A(I*3+1,I_nb[13]*3+2,  -coeff_temp/8.0);
	    } */

            if (flag[2] == 0 && flag[5] == 0)
            {
                rhs += coeff_temp*(U2_nb[2]+U2_nb[5]+U2_nb[13]+U2_center)/8.0;

                solver.Add_A(I*3+1,I*3+2,        -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[2]*3+2,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[5]*3+2,  -coeff_temp/8.0);
                solver.Add_A(I*3+1,I_nb[13]*3+2,  -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.0*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                corner_coords[2] = coords[2] + 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[2]+
                    corner_state_new.m_U[2]);
            }

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[1];

	    solver.Add_b(I*3+1, rhs);

	    /************************************************************************/

	    //Setting the coeffecients of U2 for the third equation


	    coeff2[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);
	    coeff2[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);
	    coeff2[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);
	    coeff2[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);
	    coeff2[4] = m_dt/rho * mu_edge[4]/(top_h[2]*top_h[2]);
	    coeff2[5] = m_dt/rho * mu_edge[5]/(top_h[2]*top_h[2]);

	    solver.Add_A(I*3+2,I*3+2,1.0);
	    rhs = U2_center;

	    for (nb = 0; nb < 6; nb++)
	    {
		//if (I_nb[nb] != -1)
                if (flag[nb] == 0)
		{
		    solver.Add_A(I*3+2,I_nb[nb]*3+2,-coeff2[nb]);
		    rhs += coeff2[nb]*U2_nb[nb];
		}
                else if (flag[nb] == 1)
		{
		    coeff2[nb] = 2.0*coeff2[nb];
		    rhs += coeff2[nb]*(U2_nb[nb] + U2_nb_new[nb]);
		}
		solver.Add_A(I*3+2,I*3+2,coeff2[nb]);
		rhs -= coeff2[nb]*U2_center;
	    }

	    //set the coefficients for U0 in the thrid equation (mu*u_z)_x
	    //traverse the four corners
	    //(i-1/2,j,k-1/2),(i+1/2,j,k-1/2),(i+1/2,j,k+1/2),(i-1/2,j,k+1/2)

	    //corner (i-1/2,j,k-1/2)

	    coeff_temp = mu_edge[0]*m_dt/(top_h[0]*top_h[2])/rho;
            /*
	    if (flag[0] == 1 && flag[4] == 0)
		rhs += coeff_temp*U0_nb[0];
	    else if (flag[0] == 0 && flag[4] == 1)
		rhs += coeff_temp*U0_nb[4];
	    else if (flag[0] == 1 && flag[4] == 1)
		rhs += coeff_temp*U0_nb[0];
	    else {
		rhs += coeff_temp*(U0_nb[0]+U0_nb[4]+U0_nb[14]+U0_center)/8.0;

		solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[0]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[4]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[14]*3, -coeff_temp/8.0);
	    } */

            if (flag[0] == 0 && flag[4] == 0)
            {
                rhs += coeff_temp*(U0_nb[0]+U0_nb[4]+U0_nb[14]+U0_center)/8.0;

                solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[0]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[4]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[14]*3, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.0*top_h[1];
                corner_coords[2] = coords[2] - 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }

	    //corner (i+1/2,j,k-1/2)

	    coeff_temp = -mu_edge[1]*m_dt/(top_h[0]*top_h[2])/rho;
            /*
	    if (flag[1] == 1 && flag[4] == 0)
		rhs += coeff_temp*U0_nb[1];
	    else if (flag[1] == 0 && flag[4] == 1)
		rhs += coeff_temp*U0_nb[4];
	    else if (flag[1] == 1 && flag[4] == 1)
		rhs += coeff_temp*U0_nb[1];
	    else {
		rhs += coeff_temp*(U0_nb[1]+U0_nb[4]+U0_nb[15]+U0_center)/8.0;

		solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[1]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[4]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[15]*3, -coeff_temp/8.0);
	    } */

            if (flag[1] == 0 && flag[4] == 0)
            {
                rhs += coeff_temp*(U0_nb[1]+U0_nb[4]+U0_nb[15]+U0_center)/8.0;

                solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[1]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[4]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[15]*3, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.0*top_h[1];
                corner_coords[2] = coords[2] - 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }

	    //corner (i+1/2,j,k+1/2)

	    coeff_temp = mu_edge[1]*m_dt/(top_h[0]*top_h[2])/rho;
            /*
	    if (flag[1] == 1 && flag[5] == 0)
		rhs += coeff_temp*U0_nb[1];
	    else if (flag[1] == 0 && flag[5] == 1)
		rhs += coeff_temp*U0_nb[5];
	    else if (flag[1] == 1 && flag[5] == 1)
		rhs += coeff_temp*U0_nb[1];
	    else {
		rhs += coeff_temp*(U0_nb[1]+U0_nb[5]+U0_nb[16]+U0_center)/8.0;

		solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[1]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[5]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[16]*3, -coeff_temp/8.0);
	    } */

            if (flag[1] == 0 && flag[5] == 0)
            {
                rhs += coeff_temp*(U0_nb[1]+U0_nb[5]+U0_nb[16]+U0_center)/8.0;

                solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[1]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[5]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[16]*3, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.0*top_h[1];
                corner_coords[2] = coords[2] + 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }

	    //corner (i-1/2,j,k+1/2)

	    coeff_temp = -mu_edge[0]*m_dt/(top_h[0]*top_h[2])/rho;
            /*
	    if (flag[0] == 1 && flag[5] == 0)
		rhs += coeff_temp*U0_nb[0];
	    else if (flag[0] == 0 && flag[5] == 1)
		rhs += coeff_temp*U0_nb[5];
	    else if (flag[0] == 1 && flag[5] == 1)
		rhs += coeff_temp*U0_nb[0];
	    else {
		rhs += coeff_temp*(U0_nb[0]+U0_nb[5]+U0_nb[17]+U0_center)/8.0;

		solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[0]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[5]*3,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[17]*3, -coeff_temp/8.0);
	    } */

            if (flag[0] == 0 && flag[5] == 0)
            {
                rhs += coeff_temp*(U0_nb[0]+U0_nb[5]+U0_nb[17]+U0_center)/8.0;

                solver.Add_A(I*3+2,I*3,        -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[0]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[5]*3,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[17]*3, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.0*top_h[1];
                corner_coords[2] = coords[2] + 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }


	    //set the coefficients for U1 in the thrid equation (mu*v_z)_y
	    //traverse the four corners
	    //(i,j-1/2,k-1/2),(i,j+1/2,k-1/2),(i,j+1/2,k+1/2),(i,j-1/2,k+1/2)

	    //corner (i,j-1/2,k-1/2)

	    coeff_temp = mu_edge[2]*m_dt/(top_h[1]*top_h[2])/rho;
            /*
	    if (flag[2] == 1 && flag[4] == 0)
		rhs += coeff_temp*U1_nb[2];
	    else if (flag[2] == 0 && flag[4] == 1)
		rhs += coeff_temp*U1_nb[4];
	    else if (flag[2] == 1 && flag[4] == 1)
		rhs += coeff_temp*U1_nb[2];
	    else {
		rhs += coeff_temp*(U1_nb[2]+U1_nb[4]+U1_nb[10]+U1_center)/8.0;

		solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[2]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[4]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[10]*3+1, -coeff_temp/8.0);
	    } */

            if (flag[2] == 0 && flag[4] == 0)
            {
                rhs += coeff_temp*(U1_nb[2]+U1_nb[4]+U1_nb[10]+U1_center)/8.0;

                solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[2]*3+1,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[4]*3+1,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[10]*3+1, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.0*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                corner_coords[2] = coords[2] - 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }

	    //corner (i,j+1/2,k-1/2)

	    coeff_temp = -mu_edge[3]*m_dt/(top_h[1]*top_h[2])/rho;
            /*
	    if (flag[3] == 1 && flag[4] == 0)
		rhs += coeff_temp*U1_nb[3];
	    else if (flag[3] == 0 && flag[4] == 1)
		rhs += coeff_temp*U1_nb[4];
	    else if (flag[3] == 1 && flag[4] == 1)
		rhs += coeff_temp*U1_nb[3];
	    else {
		rhs += coeff_temp*(U1_nb[3]+U1_nb[4]+U1_nb[11]+U1_center)/8.0;

		solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[3]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[4]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[11]*3+1, -coeff_temp/8.0);
	    } */

            if (flag[3] == 0 && flag[4] == 0)
            {
                rhs += coeff_temp*(U1_nb[3]+U1_nb[4]+U1_nb[11]+U1_center)/8.0;

                solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[3]*3+1,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[4]*3+1,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[11]*3+1, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.0*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                corner_coords[2] = coords[2] - 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }

	    //corner (i,j+1/2,k+1/2)

	    coeff_temp = mu_edge[3]*m_dt/(top_h[1]*top_h[2])/rho;
            /*
	    if (flag[3] == 1 && flag[5] == 0)
		rhs += coeff_temp*U1_nb[3];
	    else if (flag[3] == 0 && flag[5] == 1)
		rhs += coeff_temp*U1_nb[5];
	    else if (flag[3] == 1 && flag[5] == 1)
		rhs += coeff_temp*U1_nb[3];
	    else {
		rhs += coeff_temp*(U1_nb[3]+U1_nb[5]+U1_nb[12]+U1_center)/8.0;

		solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[3]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[5]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[12]*3+1, -coeff_temp/8.0);
	    } */

            if (flag[3] == 0 && flag[5] == 0)
            {
                rhs += coeff_temp*(U1_nb[3]+U1_nb[5]+U1_nb[12]+U1_center)/8.0;

                solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[3]*3+1,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[5]*3+1,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[12]*3+1, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.0*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                corner_coords[2] = coords[2] + 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }

	    //corner (i,j-1/2,k+1/2)

	    coeff_temp = -mu_edge[2]*m_dt/(top_h[1]*top_h[2])/rho;
            /*
	    if (flag[2] == 1 && flag[5] == 0)
		rhs += coeff_temp*U1_nb[2];
	    else if (flag[2] == 0 && flag[5] == 1)
		rhs += coeff_temp*U1_nb[5];
	    else if (flag[2] == 1 && flag[5] == 1)
		rhs += coeff_temp*U1_nb[2];
	    else {
		rhs += coeff_temp*(U1_nb[2]+U1_nb[5]+U1_nb[13]+U1_center)/8.0;

		solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[2]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[5]*3+1,  -coeff_temp/8.0);
		solver.Add_A(I*3+2,I_nb[13]*3+1, -coeff_temp/8.0);
	    } */

            if (flag[2] == 0 && flag[5] == 0)
            {
                rhs += coeff_temp*(U1_nb[2]+U1_nb[5]+U1_nb[13]+U1_center)/8.0;

                solver.Add_A(I*3+2,I*3+1,        -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[2]*3+1,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[5]*3+1,  -coeff_temp/8.0);
                solver.Add_A(I*3+2,I_nb[13]*3+1, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.0*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                corner_coords[2] = coords[2] + 0.5*top_h[2];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }

	    rhs += m_dt*state.m_U[2];
	    rhs += m_dt*cell_center[index].m_state.f_surf[2];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[2]/rho;
	    rhs -= m_dt*cell_center[index].m_state.m_adv[2];

	    solver.Add_b(I*3+2, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-10);

	start_clock("Before Petsc Solve");
        solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);

	//if (rel_residual > 1)
	//{
	//    printf("\nThe solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
	//    solver.Reset_x();
	//    solver.Solve_GMRES();
	//    solver.GetNumIterations(&num_iter);
	//    solver.GetFinalRelativeResidualNorm(&rel_residual);
	//}

	stop_clock("After Petsc Solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_2nd_coupled: "
                        "num_iter = %d, rel_residual = %le. \n",
                        num_iter,rel_residual);


	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*3-ilower*3];
                cell_center[index].m_state.m_U[1] = x[I*3-ilower*3+1];
		cell_center[index].m_state.m_U[2] = x[I*3-ilower*3+2];
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]) +
			fabs(cell_center[index].m_state.m_U[2]);
                if (speed > max_speed)
                    max_speed = speed;
            }
            else
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
		cell_center[index].m_state.m_U[2] = 0.0;
            }
        }
        for (l = 0; l < 3; ++l)
        {
	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)

            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)

            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }
        pp_global_max(&max_speed,1);

        FT_FreeThese(1,x);
}       /* end compDiffWithSmoothProperty3d_2nd_coupled */


void Incompress_Solver_Smooth_3D_Cartesian::computeProjection(void)
{
	int index, index_nb[6], size;
	double rhs, coeff[6], rho[6], rho0;
	int I,I_nb[6];
	int i,j,k,l,icoords[MAXD];
	INTERFACE *intfc = front->interf;
	double P_max,P_min;
	int icrds_Pmax[MAXD],icrds_Pmin[MAXD];
	COMPONENT comp;
	double aII,rho_nb[6];
	double coords[MAXD],crx_coords[MAXD];
	double **vel = iFparams->field->vel;
	int num_nb;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	boolean use_neumann_solver = YES;
	max_value = 0.0;
	double value;
	double sum_div;
	double sum_phi;
	sum_div = 0.0;
	sum_phi = 0.0;
	PetscInt num_iter = 0;
	double rel_residual = 0.0;

	PETSc solver;
	solver.Create(ilower, iupper-1, 7, 7);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();
	size = iupper - ilower;

	setIndexMap();

	for (l = 0; l < dim; ++l)
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    vel[l][index] = cell_center[index].m_state.m_U[l];
	}

	/* Compute velocity divergence */
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = computeFieldPointDiv(icoords,vel);
	}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.div_U = array[index];
	}


	if(debugging("step_size"))
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
	        value = fabs(cell_center[index].m_state.div_U);
		sum_div = sum_div + cell_center[index].m_state.div_U;
	        if(value > max_value)
		    max_value = value;
	    }
	    pp_global_sum(&sum_div,1);
	    printf("\nThe summation of divergence of U is %.16g\n",sum_div);
	    pp_global_max(&max_value,1);
	    printf("\nThe max value of divergence of U is %.16g\n",max_value);
	    max_value = 0.0;
	}


	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
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

	    rho0   = cell_center[index].m_state.m_rho;
	    num_nb = 0;
	    for (l = 0; l < 6; ++l)
	    {
		if (I_nb[l] == -1)
		    index_nb[l] = index;
		else num_nb++;
		rho[l] = 1.0/2*(rho0 + cell_center[index_nb[l]].m_state.m_rho);
		coeff[l] = 1.0/rho[l]/(top_h[l/2]*top_h[l/2]);
	    }

	    rhs = cell_center[index].m_state.div_U/accum_dt;

	    aII = 0.0;
	    for (l = 0; l < 6; ++l)
	    {
		if (num_nb <= 1) break;
	    	if (I_nb[l] != -1)
		{
		    solver.Set_A(I,I_nb[l],coeff[l]);
                    aII += -coeff[l];
		}
		else
                {
                    if (isDirichletPresetBdry(front,icoords,dir[l],comp))
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
            if (num_nb > 1)
	    {
                solver.Set_A(I,I,aII);
	    }
            else
            {
		printf("\nUnsolvable psudo-pressure, using original pressure!\n");
                solver.Set_A(I,I,1.0);
                rhs = cell_center[index].m_state.m_P;
            }
            solver.Set_b(I,rhs);
	}
	use_neumann_solver = pp_min_status(use_neumann_solver);

	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);

	start_clock("Before Petsc Solver in Projection step");
	if (use_neumann_solver)
	{
	    printf("\nUsing Neumann Solver!\n");
	    solver.Solve_withPureNeumann();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	    if(rel_residual > 1)
	    {
		printf("\n The solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_withPureNeumann_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	else
	{
	    printf("\nUsing non-Neumann Solver!\n");
	    solver.Solve();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);

	    if(rel_residual > 1)
	    {
		printf("\n The solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }
	}
	stop_clock("After Petsc Solver in Projection step");

	double *x;
	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

	if (debugging("PETSc"))
	    (void) printf("Incompress_Solver_Smooth_3D_Cartesian::"
			"computeProjection: "
	       		"num_iter = %d, rel_residual = %le \n",
			num_iter, rel_residual);

	P_max = -HUGE;
        P_min = HUGE;
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    I = ijk_to_I[i][j][k];
	    array[index] = x[I-ilower];
	}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_phi = array[index];
	}

	if(debugging("step_size"))
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		value = fabs(cell_center[index].m_state.m_phi);
		sum_phi = sum_phi + cell_center[index].m_state.m_phi*top_h[0]*top_h[1]*top_h[2];
		if (value > max_value)
		    max_value = value;
	    }
	    pp_global_sum(&sum_phi,1);
	    printf("\nThe integration of phi is %.16g\n", sum_phi);
	    pp_global_max(&max_value,1);
	    printf("\nThe max value of phi is %.16g\n",max_value);
	}

	FT_FreeThese(1,x);
} /* end computeProjection3d */


void Incompress_Solver_Smooth_3D_Cartesian::computeNewVelocity(void)
{
	int i, j, k, l, index;
	double grad_phi[3], rho;
	COMPONENT comp;
	double speed;
	int icoords[MAXD];

	max_speed = 0.0;

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = cell_center[index].m_state.m_phi;
	}
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    if (!ifluid_comp(comp))
	    {
		for (l = 0; l < 3; ++l)
		    cell_center[index].m_state.m_U[l] = 0.0;
		continue;
	    }
	    rho = cell_center[index].m_state.m_rho;
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    computeFieldPointGradPhi(icoords,array,grad_phi);
	    speed = 0.0;
	    for (l = 0; l < 3; ++l)
	    {
	    	cell_center[index].m_state.m_U[l] -= accum_dt/rho*grad_phi[l];
		speed += fabs(cell_center[index].m_state.m_U[l]);
	    }
	    if (speed > max_speed)
		max_speed = speed;
	}
	for (l = 0; l < 3; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
	    	array[index] = cell_center[index].m_state.m_U[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].m_state.m_U[l] = array[index];
	    }
	}
	pp_global_max(&max_speed,1);
}	/* end computeNewVelocity3d */

void Incompress_Solver_Smooth_3D_Cartesian::computeSourceTerm_Adv(double *coords, L_STATE &state)
{
	int i;
	for (i = 0; i < dim; ++i)
	    state.m_U[i] = iFparams->gravity[i];

	state.m_P = HUGE_VAL;
}


void Incompress_Solver_Smooth_3D_Cartesian::computeSourceTerm(double *coords, L_STATE &state)
{
	int i;
	for (i = 0; i < dim; ++i)
	    state.m_U[i] = iFparams->gravity[i];

	state.m_P = HUGE_VAL;
}
void Incompress_Solver_Smooth_3D_Cartesian::computeSourceTerm(double *coords, double t, L_STATE &state)
{
	computeSourceTerm(coords, state);
}

// for initial condition:
// 		setInitialCondition();
// this function should be called before solve()
// for the source term of the momentum equation:
// 		computeSourceTerm();
void Incompress_Solver_Smooth_3D_Cartesian::solve(double dt)
{

        printf("\nEntering Incompress Solve 3D!! The dt for this step is %.16g\n",dt);

	m_t_old = front->time;
	m_t_int = front->time + front->dt/2.0;
	m_t_new = front->time + front->dt;

	static boolean first = YES;
	if (first)
	{
	    accum_dt = 0.0;
	    first = NO;
	}
	m_dt = dt;
	max_speed = 0.0;

	start_clock("solve");
	setDomain();

	setComponent();
	if (debugging("trace"))
	    printf("Passed setComponent()\n");
	setGlobalIndex();
	if (debugging("trace"))
	    printf("Passed setGlobalIndex()\n");
	start_clock("setSmoothedProperties");
	setSmoothedProperties();
	stop_clock("setSmoothedProperties");
	if (debugging("trace"))
	    printf("Passed setSmoothedProperties()\n");

	// 1) solve for intermediate velocity
	start_clock("computeAdvection");
	computeAdvection();
	stop_clock("computeAdvection");
	if (debugging("trace"))
	    printf("max_speed after computeAdvection(): %20.14f\n",
				max_speed);
	if (debugging("sample_velocity"))
	    sampleVelocity();

	start_clock("compDiffWithSmoothProperty");
	compDiffWithSmoothProperty_1st_decoupled();
	//compDiffWithSmoothProperty();
	stop_clock("compDiffWithSmoothProperty");
	if (debugging("sample_velocity"))
	    sampleVelocity();

        start_clock("compSGS");
        //compSGS();	//Subgrid model by Hyunkyun Lim
        stop_clock("compSGS");

	if (debugging("trace"))
	    printf("max_speed after compDiffWithSmoothProperty(): %20.14f\n",
				max_speed);

	// 2) projection step
	accum_dt += m_dt;
	if (accum_dt >= min_dt)
	{
	    start_clock("computeProjection");
	    computeProjection();
	    stop_clock("computeProjection");

	    start_clock("computePressure");
	    computePressure();
	    stop_clock("computePressure");

	    start_clock("computeNewVelocity");
	    computeNewVelocity();
	    stop_clock("computeNewVelocity");
	    accum_dt = 0.0;
	}

	if (debugging("sample_velocity"))
	    sampleVelocity();

	if (debugging("trace"))
	    printf("max_speed after computeNewVelocity(): %20.14f\n",
				max_speed);

	start_clock("copyMeshStates");
	copyMeshStates();
	stop_clock("copyMeshStates");

	setAdvectionDt();
	stop_clock("solve");
}	/* end solve */


//solve_MAC_vd
void Incompress_Solver_Smooth_3D_Cartesian::solve_vd(double dt)
{
        (void) printf("\nEntering Incompress Solve_vd 3D! The dt for this time step is %.16g\n\n",dt);

        m_t_old = front->time;
        m_t_int = front->time + dt/2.0;
        m_t_new = front->time + dt;

        static boolean first = YES;
        bGhostCell = YES;
        (void) printf("useGhostCell = %d for vd problem!\n", bGhostCell);

        if (first) {
            accum_dt = 0;
            first = NO;
        }
        m_dt = dt;
        max_speed = 0;
        max_density = 0;

        //solve_MAC_vd
        start_clock("solve_vd");

        setDomain_vd();
        setComponent();

        setGlobalIndex();
        setIndexMap();

        //iterate for 1st step to get a good approx. for pressure
        int i,j,k,l,index,ite_num,tol_num;

        if (front->step==0 && dt==0) //zeroth time-step
            tol_num = 1;
        else if (front->step==0 && dt!=0) //1st time-step
            tol_num = 4;
        else //front->step>0, other time-steps
            tol_num = 1;

        ite_num = 0;
        while (ite_num < tol_num)
        {
            // 1) solve for intermediate U, density and concentration
            //solve for approximated adv terms using lagged source term
            start_clock("compAdvectionTerm_MAC_decoupled_vd(0)");
            compAdvectionTerm_MAC_decoupled_vd(0, bGhostCell);
            stop_clock("compAdvectionTerm_MAC_decoupled_vd(0)");

            //computeNewDensity_MAC_constRho_vd();

            //solve for estimated density explicitly
            start_clock("compNewDensity_vd(0)");
            computeNewDensity_vd(0);
            stop_clock("compNewDensity_vd(0)");
            if (debugging("step_size"))
            {
                (void) printf("max_density after computeNewDensity_vd(0): %.16g\n",
                              max_density);
                (void) printf("min_density after computeNewDensity_vd(0): %.16g\n",
                              min_density);
            }

            //solve for accurate adv terms
            start_clock("compAdvectionTerm_MAC_decoupled_vd(1)");
            compAdvectionTerm_MAC_decoupled_vd(1, bGhostCell);
            stop_clock("compAdvectionTerm_MAC_decoupled_vd(1)");

            //solve for accurate density explicitly, and calc dynamic viscosity
            start_clock("compNewDensity_vd(1)");
            computeNewDensity_vd(1);
            stop_clock("compNewDensity_vd(1)");
            if (debugging("step_size"))
            {
                (void) printf("max_density after computeNewDensity_vd(1): %.16g\n",
                              max_density);
                (void) printf("min_density after computeNewDensity_vd(1): %.16g\n",
                              min_density);
            }

/*
            if (ite_num == (tol_num-1)) //in the last iteration step, update concentrations
            {
                //solve for concentration implicitly
                start_clock("compNewConcentration_vd");
                computeNewConcentration_vd();
                stop_clock("compNewConcentration_vd");
                if (debugging("step_size"))
                {
                    (void) printf("max_concentration after computeNewConcentration_vd: %.16g\n",
                                   max_concentration);
                    (void) printf("min_concentration after computeNewConcentration_vd: %.16g\n",
                                   min_concentration);
                }
            }
*/

            // 2) projection step
            //solve for U^{star} implicitly
            start_clock("compDiffWithSmoothProperty_velocity_MAC_vd");
            //for D=0 and constant mu
            //compDiffWithSmoothProperty_velocity_MAC_decoupled_vd();

            //valid for variable dynamic viscosity
            compDiffWithSmoothProperty_velocity_MAC_coupled_vd();
            stop_clock("compDiffWithSmoothProperty_velocity_MAC_vd");
            if (debugging("step_size"))
            {
                if (dt == 0)
                {
                    (void) printf("max_speed after compDiffWithSmoothProperty_velocity_MAC_vd() "
                                  "after #ite %d in time-step 0 is: %.16g\n",ite_num+1,max_speed);
                }
                else
                {
                    (void) printf("max_speed after compDiffWithSmoothProperty_velocity_MAC_vd() "
                                  "after #ite %d in time-step %d is: %.16g\n",
                                  ite_num+1,front->step+1,max_speed);
                }
            }

            start_clock("compSGS");
            //Subgrid-scale model for vd
            //TODO && FIXME: Smeeton Young 105s don't need SGS term.
            computeSubgridModel_vd();
            stop_clock("compSGS");

            //projection and updating velocities & pressures
            accum_dt += m_dt;
            if (accum_dt >= min_dt)
            {
                start_clock("computeProjection_MAC_vd");
                computeProjection_MAC_vd();
                stop_clock("computeProjection_MAC_vd");

                start_clock("computePressure_MAC_vd");
                //Could PmII and PmIII be modified to be applied to vd?????????
                computePressure_MAC_vd();
                stop_clock("computePressure_MAC_vd");

                if (ite_num < (tol_num-1)) //set U to be U^n, rho to be rho^n
                {
                    for (k = 0; k <= top_gmax[2]; k++)
                    for (j = 0; j <= top_gmax[1]; j++)
                    for (i = 0; i <= top_gmax[0]; i++)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        for (l = 0; l < 3; l++)
                            cell_center[index].m_state.m_U[l] = cell_center[index].m_state.m_U_tmp[l];
                    }

                    for (k = 0; k <= top_gmax[2]; k++)
                    for (j = 0; j <= top_gmax[1]; j++)
                    for (i = 0; i <= top_gmax[0]; i++)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        cell_center[index].m_state.m_rho = cell_center[index].m_state.m_rho_old;
                        cell_center[index].m_state.m_rho_old = cell_center[index].m_state.m_rho;
                        double nu = (m_mu[0]+m_mu[1])/(m_rho[0]+m_rho[1]);
                        cell_center[index].m_state.m_mu = cell_center[index].m_state.m_rho*nu;
                        cell_center[index].m_state.m_mu_old = cell_center[index].m_state.m_rho_old*nu;
                    }
                }
                else //the last iteration step, update velocities
                {
                    start_clock("computeNewVelocity_fullMAC_vd");
                    computeNewVelocity_fullMAC_vd();
                    stop_clock("computeNewVelocity_fullMAC_vd");

		    for (l = 0; l < dim; ++l)
		    for (k = 0; k <= top_gmax[2]; k++)
		    for (j = 0; j <= top_gmax[1]; j++)
		    for (i = 0; i <= top_gmax[0]; i++)
		    {
			index = d_index3d(i,j,k,top_gmax);
			cell_center[index].m_state.m_U_velo_var[l] += cell_center[index].m_state.m_U[l];
		    }


                    if (debugging("step_size"))
                    {
                        (void) printf("max_speed after computeNewVelocity_fullMAC_vd() "
                                      "in time-step %d is: %.16g\n", front->step+1,max_speed);
                    }
                }
                accum_dt = 0;
            }
            ++ite_num;
        } //while-loop ends

        if (debugging("sample_velocity"))
            sampleVelocity();

        start_clock("copyMeshStates_vd");
        copyMeshStates_vd();
        stop_clock("copyMeshStates_vd");
        setAdvectionDt_vd();

        stop_clock("solve_vd");
} /* end solve_vd */

//calc parameters for RT simulation
void Incompress_Solver_Smooth_3D_Cartesian::computeRTParameters(double dt, char *out_name, int nStep_velo_var)
{
        INTERFACE *intfc = front->interf;
	const RECT_GRID *rgr = front->rect_grid;
	const double *GL = rgr->GL;
        const double *GU = rgr->GU;
        POINT   *p;
        SURFACE **s;
        TRI     *tri;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        boolean  bdry = NO;
        int i,j,k,l,index,count,ucount;
        double *crds, coords[MAXD];
        double zmin=HUGE, zmax=-HUGE;
        double A,g,t2,Agt2,tau;
	double h_bubble_intfc,h_spike_intfc; //bubble/spike height of contact interface
        double h_bubble_vf,h_spike_vf; //bubble/spike height of n% volume fraction contour
        //H = 32cm is the vertical height of the water channel
        //refer to Eq. (18) of Mueschke's paper in 2009
        double H = 32.0;
        FILE *outfile;
        double z_coord,max_dens,min_dens,Dens;
        double sum_vf1,sum_vf2,sum_vf1vf2,mean_vf1,mean_vf2,mean_vf1vf2;
        double usum_vf1,usum_vf2,usum_vf1vf2,umean_vf1,umean_vf2,umean_vf1vf2;
        double vf1,vf2,vf1vf2,uvf1,uvf2,uvf1vf2;

	z0 = (GL[dim-1] + GU[dim-1])/2.0;
        A = fabs(m_rho[1]-m_rho[0])/(m_rho[1]+m_rho[0]);
        g = fabs(iFparams->gravity[dim-1]);
        t2 = front->time*front->time;
        Agt2 = A*g*t2;
        tau = sqrt(A*g/H)*front->time;
        max_dens = (m_rho[0]>m_rho[1])? m_rho[0] : m_rho[1];
        min_dens = (m_rho[0] + m_rho[1]) - max_dens;

	//0. print intfc for selected tau's
        if ( fabs(tau-0.21) <= 0.5*sqrt(A*g/H)*front->dt ||
             fabs(tau-0.50) <= 0.5*sqrt(A*g/H)*front->dt ||
             fabs(tau-1.01) <= 0.5*sqrt(A*g/H)*front->dt ||
             fabs(tau-1.52) <= 0.5*sqrt(A*g/H)*front->dt )
	{
            //print P-node/intfc files for making plot
            FT_AddMovieFrame(front,out_name,NO);
	}


	//1. get bubble/spike height of interface

	bdry = NO;
        for (s = intfc->surfaces; s && *s; ++s)
        {
            for (tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
            {
                for (k = 0; k < dim; ++k)
                    Index_of_point(Point_of_tri(tri)[k]) = -1;
            }
        }

        for (s = intfc->surfaces; s && *s; ++s)
        {
            if (bdry == YES  &&  !Boundary(*s))
                continue;
            if (bdry == NO  &&  Boundary(*s))
                continue;
            else
            {
                for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s);
                     tri = tri->next)
                {
                    for (k = 0; k < dim; ++k)
                    {
                        p = Point_of_tri(tri)[k];
                        if (Index_of_point(p) == -1)
                        {
                            crds = Coords(p);
                            if (crds[dim-1] > zmax)
                                zmax = crds[dim-1];
                            if (crds[dim-1] < zmin)
                                zmin = crds[dim-1];
                        }
                    }
                }
            }
        }
        pp_global_max(&zmax,1);
        pp_global_min(&zmin,1);
	// store zmax and zmin
	zmax_intfc = zmax;
	zmin_intfc = zmin;
        if (debugging("glayer"))
        {
            printf("\n\tzmin_intfc = %lf\n", zmin);
            printf("\tzmax_intfc = %lf\n", zmax);
        }
        h_bubble_intfc = (m_rho[0]<m_rho[1]) ? fabs(zmin-z0):fabs(zmax-z0);
        h_spike_intfc = (m_rho[0]<m_rho[1]) ? fabs(zmax-z0):fabs(zmin-z0);


	//2. get bubble/spike height of volume fraction contour
	//1%~99% or 5%~95% vol frac contour
	//Note that vf = 100% at top bdry and vf = 0 at bottom bdry

        double tol = iFparams->vol_frac_threshold;
	if (tol > 0.5)	tol = 1.0-tol;
        //if (debugging("glayer"))
	{
	    printf("vol_frac_threshold = %lf\n", tol);
	}
	//search for zmin from lowest cell centers in z-direction
        zmin = new_height_at_fraction_vd(GL[dim-1]+0.5*top_h[dim-1],+1,tol);
        //search for zmax from highest cell centers in z-direction
        zmax = new_height_at_fraction_vd(GU[dim-1]-0.5*top_h[dim-1],-1,1.0-tol);
        // store zmax and zmin
        zmax_vf = zmax;
        zmin_vf = zmin;
        if (debugging("glayer"))
        {
            printf("\nAfter calls to height_at_fraction_vd():\n");
            printf("\n\tzmin_vf = %lf\n", zmin);
            printf("\tzmax_vf = %lf\n", zmax);
        }
        h_bubble_vf = (m_rho[0]<m_rho[1]) ? fabs(zmin-z0):fabs(zmax-z0);
        h_spike_vf = (m_rho[0]<m_rho[1]) ? fabs(zmax-z0):fabs(zmin-z0);


        //3. get molecular mixing parameter theta
        //theta = <vf1vf2>/(<vf1><vf2>)

        count = 0;
        sum_vf1 = sum_vf2 = sum_vf1vf2 = 0;

        ucount = 0;
        usum_vf1 = usum_vf2 = usum_vf1vf2 = 0;

        double vf1_midPlane,vf2_midPlane,vf1vf2_midPlane,mean_vf1_midPlane,mean_vf2_midPlane,mean_vf1vf2_midPlane;
        double Dens_nb, Dens_midPlane;
        double sum_vf1_midPlane = 0, sum_vf2_midPlane = 0, sum_vf1vf2_midPlane = 0;

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            Dens = cell_center[index].m_state.m_rho;
            z_coord = cell_center[index].m_coords[dim-1];

            if (z_coord > z0-top_h[dim-1] && z_coord < z0)
            {
                count++;

                //volume/mole fraction
                vf2 = (Dens-max_dens)/(min_dens-max_dens);
                vf1 = (Dens-min_dens)/(max_dens-min_dens);
                vf1vf2 = vf1*vf2;
                sum_vf1 += vf1;
                sum_vf2 += vf2;
                sum_vf1vf2 += vf1vf2;

                Dens_nb = cell_center[d_index3d(i,j,k+1,top_gmax)].m_state.m_rho;
                Dens_midPlane = 0.5*(Dens + Dens_nb);
                vf2_midPlane = (Dens_midPlane-max_dens)/(min_dens-max_dens);
                vf1_midPlane = (Dens_midPlane-min_dens)/(max_dens-min_dens);
                vf1vf2_midPlane = vf1_midPlane*vf2_midPlane;
                sum_vf1_midPlane += vf1_midPlane;
                sum_vf2_midPlane += vf2_midPlane;
                sum_vf1vf2_midPlane += vf1vf2_midPlane;
            }

            if (z_coord > z0 && z_coord < z0+top_h[dim-1])
            {
                ucount++;

                //volume/mole fraction
                uvf2 = (Dens-max_dens)/(min_dens-max_dens);
                uvf1 = (Dens-min_dens)/(max_dens-min_dens);
                uvf1vf2 = uvf1*uvf2;
                usum_vf1 += uvf1;
                usum_vf2 += uvf2;
                usum_vf1vf2 += uvf1vf2;
            }
        }

        if (pp_numnodes() > 1)
        {
            pp_global_isum(&count,1);
            pp_global_sum(&sum_vf1,1);
            pp_global_sum(&sum_vf2,1);
            pp_global_sum(&sum_vf1vf2,1);

            pp_global_isum(&ucount,1);
            pp_global_sum(&usum_vf1,1);
            pp_global_sum(&usum_vf2,1);
            pp_global_sum(&usum_vf1vf2,1);

            pp_global_sum(&sum_vf1_midPlane,1);
            pp_global_sum(&sum_vf2_midPlane,1);
            pp_global_sum(&sum_vf1vf2_midPlane,1);
        }
        mean_vf1 = sum_vf1/count;
        mean_vf2 = sum_vf2/count;
        mean_vf1vf2 = sum_vf1vf2/count;

        umean_vf1 = usum_vf1/ucount;
        umean_vf2 = usum_vf2/ucount;
        umean_vf1vf2 = usum_vf1vf2/ucount;

        mean_vf1_midPlane = sum_vf1_midPlane/count;
        mean_vf2_midPlane = sum_vf2_midPlane/count;
        mean_vf1vf2_midPlane = sum_vf1vf2_midPlane/count;


        //4. get velocity variance
        //Var(X) = E[X*X] - E[X]*E[X]

        double u,v,w,mean_u,mean_v,mean_w,mean_uu,mean_vv,mean_ww;
        double sum_u = 0, sum_v = 0, sum_w = 0, sum_uu = 0, sum_vv = 0, sum_ww = 0;
        double sum_u_raw = 0, sum_v_raw = 0, sum_w_raw = 0;
	double sum_uu_raw = 0, sum_vv_raw = 0, sum_ww_raw = 0;
        double sum_uw_raw = 0, sum_vw_raw = 0;

	double mean_u_raw, mean_v_raw, mean_w_raw, mean_uu_raw, mean_vv_raw, mean_ww_raw, mean_uw_raw, mean_vw_raw;

        // escape velocity for laser sheet, unit: cm/s
        // v_esc = 0.5*(w+d)/(1/f)
        // where w = 532 (nm) is width of laser sheet, d = 13 (micron) is mean particle diamater of Conduct-O-Fil silver-coated hollow glass spheres, f = 30 (Hz) is the sampling rate of 2 lasers.
        double v_esc = 0.040596; //unit: cm/s
        int count_nonEsc = 0;
        count = 0;

        int ii, jj, kk;
        //average to interrogation window, whose size is 0.135*0.135 cm^2 (in uw-plane)
        int NB[3] = {1, 1, 1}; // coarse and medium grid
        //if (fabs(top_h[0] - 0.05) < 1e-12) // fine grid
        //    NB[0] = NB[2] = 2;

        for (k = 0; k <= ((kmax-kmin+1)/NB[2])-1; k++)
        for (j = 0; j <= ((jmax-jmin+1)/NB[1])-1; j++)
        for (i = 0; i <= ((imax-imin+1)/NB[0])-1; i++)
        {
            kk = (NB[2]*k)+kmin;
            jj = (NB[1]*j)+jmin;
            ii = (NB[0]*i)+imin;

            if (NB[0] == 1 && NB[1] == 1 && NB[2] == 1) // coarse and medium grid
            {
                z_coord = cell_center[d_index3d(ii,jj,kk,top_gmax)].m_coords[2];
            }
            else if (NB[0] == 2 && NB[1] == 1 && NB[2] == 2) // fine grid
            {
                z_coord = (cell_center[d_index3d(ii,jj,kk,top_gmax)].m_coords[2]+
                           cell_center[d_index3d(ii,jj,kk+1,top_gmax)].m_coords[2])/2;
            }

            if (z_coord > z0-top_h[dim-1]*NB[2] && z_coord < z0)
            {
                count++;

                // cell center of (i, j, k)
                if (NB[0] == 1 && NB[1] == 1 && NB[2] == 1) // coarse and medium grid
                {
                    u = (cell_center[d_index3d(ii,jj,kk,top_gmax)].m_state.m_U_velo_var[0]+
                         cell_center[d_index3d(ii-1,jj,kk,top_gmax)].m_state.m_U_velo_var[0])/2;
                    v = (cell_center[d_index3d(ii,jj,kk,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii,jj-1,kk,top_gmax)].m_state.m_U_velo_var[1])/2;
                    w = (cell_center[d_index3d(ii,jj,kk,top_gmax)].m_state.m_U_velo_var[2]+
                         cell_center[d_index3d(ii,jj,kk-1,top_gmax)].m_state.m_U_velo_var[2])/2;
                }
                else if (NB[0] == 2 && NB[1] == 1 && NB[2] == 2) // fine grid
                {
                    u = (cell_center[d_index3d(ii,jj,kk,top_gmax)].m_state.m_U_velo_var[0]+
                         cell_center[d_index3d(ii,jj,kk+1,top_gmax)].m_state.m_U_velo_var[0])/2;
                    v = (cell_center[d_index3d(ii,jj,kk,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii,jj,kk+1,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii+1,jj,kk+1,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii+1,jj,kk,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii,jj-1,kk,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii,jj-1,kk+1,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii+1,jj-1,kk+1,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii+1,jj-1,kk,top_gmax)].m_state.m_U_velo_var[1])/8;
                    w = (cell_center[d_index3d(ii,jj,kk,top_gmax)].m_state.m_U_velo_var[2]+
                         cell_center[d_index3d(ii+1,jj,kk,top_gmax)].m_state.m_U_velo_var[2])/2;
                }
                else
                {
                    printf("ERROR: unrecognized grids\n");
                    clean_up(ERROR);
                }

                // temporal-average
                u /= nStep_velo_var;
                v /= nStep_velo_var;
                w /= nStep_velo_var;

                sum_u_raw += u;
                sum_v_raw += v;
                sum_w_raw += w;
                sum_uu_raw += u*u;
                sum_vv_raw += v*v;
                sum_ww_raw += w*w;
                sum_uw_raw += u*w;
                sum_vw_raw += v*w;

                if (fabs(v) <= v_esc)
                {
                    count_nonEsc++;
                    sum_u += u;
                    sum_v += v;
                    sum_w += w;
                    sum_uu += u*u;
                    sum_vv += v*v;
                    sum_ww += w*w;
                }
            }
        }

        // re-initialize m_U_velo_var to be 0
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_velo_var[l] = 0;
        }

        if (pp_numnodes() > 1)
        {
            pp_global_isum(&count,1);
            pp_global_isum(&count_nonEsc,1);
            pp_global_sum(&sum_u,1);
            pp_global_sum(&sum_v,1);
            pp_global_sum(&sum_w,1);
            pp_global_sum(&sum_uu,1);
            pp_global_sum(&sum_vv,1);
            pp_global_sum(&sum_ww,1);

            pp_global_sum(&sum_u_raw,1);
            pp_global_sum(&sum_v_raw,1);
            pp_global_sum(&sum_w_raw,1);
            pp_global_sum(&sum_uu_raw,1);
            pp_global_sum(&sum_vv_raw,1);
            pp_global_sum(&sum_ww_raw,1);
            pp_global_sum(&sum_uw_raw,1);
            pp_global_sum(&sum_vw_raw,1);
        }
        mean_u = sum_u/count_nonEsc;
        mean_v = sum_v/count_nonEsc;
        mean_w = sum_w/count_nonEsc;
        mean_uu = sum_uu/count_nonEsc;
        mean_vv = sum_vv/count_nonEsc;
        mean_ww = sum_ww/count_nonEsc;

        mean_u_raw = sum_u_raw/count;
        mean_v_raw = sum_v_raw/count;
        mean_w_raw = sum_w_raw/count;
        mean_uu_raw = sum_uu_raw/count;
        mean_vv_raw = sum_vv_raw/count;
        mean_ww_raw = sum_ww_raw/count;
        mean_uw_raw = sum_uw_raw/count;
        mean_vw_raw = sum_vw_raw/count;

	//printf("(count, count_nonEsc) = (%d, %d)\n", count, count_nonEsc);
        //printf("raw data, (mean_u, mean_v, mean_w) = (%lf, %lf, %lf)\n", mean_u_raw, mean_v_raw, mean_w_raw);
        //printf("filter, (mean_u, mean_v, mean_w) = (%lf, %lf, %lf)\n", mean_u, mean_v, mean_w);

        //print output file
        char filename[256];
        sprintf(filename,"%s/sim_hc_params.txt",out_name);
        if (pp_mynode() == 0)
        {
            if (front->step == 0)
            {
                outfile = fopen(filename,"w");
                fprintf(outfile,"  ts     time         tau          Agt2       hb_intfc     hs_intfc      hb_vf        hs_vf        theta1       theta2    var(u)/AgH                 var(v)/AgH   var(w)/AgH   raw_variances\n");
            }
            else {
                outfile = fopen(filename,"a");
            }
            //fprintf(outfile,"%4d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", front->step,front->time,tau,Agt2,h_bubble_intfc,h_spike_intfc,h_bubble_vf,h_spike_vf,0.5*(mean_vf1vf2/(mean_vf1*mean_vf2) + umean_vf1vf2/(umean_vf1*umean_vf2)),mean_vf1vf2_midPlane/(mean_vf1_midPlane*mean_vf2_midPlane),0.0,(mean_uu-mean_u*mean_u)/fabs(A*g*H),(mean_vv-mean_v*mean_v)/fabs(A*g*H),(mean_ww-mean_w*mean_w)/fabs(A*g*H),(mean_uu_raw-mean_u_raw*mean_u_raw)/fabs(A*g*H),(mean_vv_raw-mean_v_raw*mean_v_raw)/fabs(A*g*H),(mean_ww_raw-mean_w_raw*mean_w_raw)/fabs(A*g*H));
            fprintf(outfile,"%4d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", front->step,front->time,tau,Agt2,h_bubble_intfc,h_spike_intfc,h_bubble_vf,h_spike_vf,0.5*(mean_vf1vf2/(mean_vf1*mean_vf2) + umean_vf1vf2/(umean_vf1*umean_vf2)),mean_vf1vf2_midPlane/(mean_vf1_midPlane*mean_vf2_midPlane),0.0,(mean_uu_raw-mean_u_raw*mean_u_raw),(mean_vv_raw-mean_v_raw*mean_v_raw),(mean_ww_raw-mean_w_raw*mean_w_raw),(mean_uw_raw-mean_u_raw*mean_w_raw),(mean_vw_raw-mean_v_raw*mean_w_raw));

            fclose(outfile);
        }


        if (debugging("trace"))
            printf("Leave computeRTParameters()\n");
}  /* end computeRTParameters */
/*
//calc parameters for RT simulation
void Incompress_Solver_Smooth_3D_Cartesian::computeRTParameters(double dt, char *out_name, int nStep_velo_var)
{
        INTERFACE *intfc = front->interf;
	const RECT_GRID *rgr = front->rect_grid;
	const double *GL = rgr->GL;
        const double *GU = rgr->GU;
        POINT   *p;
        SURFACE **s;
        TRI     *tri;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        boolean  bdry = NO;
        int i,j,k,l,index,count,ucount;
        double *crds, coords[MAXD];
        double zmin=HUGE, zmax=-HUGE;
        double A,g,t2,Agt2,tau;
	double h_bubble_intfc,h_spike_intfc; //bubble/spike height of contact interface
        double h_bubble_vf,h_spike_vf; //bubble/spike height of n% volume fraction contour
        //H = 32cm is the vertical height of the water channel
        //refer to Eq. (18) of Mueschke's paper in 2009
        double H = 32.0;//TODO: This is hard-wired. NEED TO BE FIXME: Be aware of the unit
        FILE *outfile;
        double z_coord,max_dens,min_dens,Dens;
        double sum_vf1,sum_vf2,sum_vf1vf2,mean_vf1,mean_vf2,mean_vf1vf2;
        double usum_vf1,usum_vf2,usum_vf1vf2,umean_vf1,umean_vf2,umean_vf1vf2;
        double vf1,vf2,vf1vf2,uvf1,uvf2,uvf1vf2;
        double h_bubble_corner, h_spike_corner, h_bubble_edge, h_spike_edge;// This is new measurement based on Glimm's new discovery.
        double zmax_bubble_corner = -HUGE, zmin_spike_corner = HUGE, zmax_bubble_edge = -HUGE, zmin_spike_edge = HUGE;

	z0 = (GL[dim-1] + GU[dim-1])/2.0;
        A = fabs(m_rho[1]-m_rho[0])/(m_rho[1]+m_rho[0]);
        g = fabs(iFparams->gravity[dim-1]);
        t2 = front->time*front->time;
        Agt2 = A*g*t2;
        tau = sqrt(A*g/H)*front->time;
        max_dens = (m_rho[0]>m_rho[1])? m_rho[0] : m_rho[1];
        min_dens = (m_rho[0] + m_rho[1]) - max_dens;

	//0. print intfc for selected tau's
        if ( fabs(tau-0.21) <= 0.5*sqrt(A*g/H)*front->dt ||
             fabs(tau-0.50) <= 0.5*sqrt(A*g/H)*front->dt ||
             fabs(tau-1.01) <= 0.5*sqrt(A*g/H)*front->dt ||
             fabs(tau-1.52) <= 0.5*sqrt(A*g/H)*front->dt )
	{
            //print P-node/intfc files for making plot
            FT_AddMovieFrame(front,out_name,NO);
	}


	//1. get bubble/spike height of interface

	bdry = NO;
        for (s = intfc->surfaces; s && *s; ++s)
        {
            for (tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
            {
                for (k = 0; k < dim; ++k)
                    Index_of_point(Point_of_tri(tri)[k]) = -1;
            }
        }

        for (s = intfc->surfaces; s && *s; ++s)
        {
            if (bdry == YES  &&  !Boundary(*s))
                continue;
            if (bdry == NO  &&  Boundary(*s))
                continue;
            else
            {
                for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s);
                     tri = tri->next)
                {
                    for (k = 0; k < dim; ++k)
                    {
                        p = Point_of_tri(tri)[k];
                        if (Index_of_point(p) == -1)
                        {
                            crds = Coords(p);
                            if (crds[dim-1] > zmax)
                                zmax = crds[dim-1];
                            if (crds[dim-1] < zmin)
                                zmin = crds[dim-1];
                            //TODO: corner value 0 < x < dh[0]
                            if (crds[0] > 0.0  && crds[0] < rgr->h[0])
                            {
                                if (crds[dim-1] > zmax_bubble_corner)
                                    zmax_bubble_corner = crds[dim-1];
                                if (crds[dim-1] < zmin_spike_corner)
                                    zmin_spike_corner = crds[dim-1];
                            }
                            //TODO: edge value 0 < y < dh[1]
                            if (crds[1] > 0.0 && crds[1] < rgr->h[1])
                            {
                                 if (crds[dim-1] > zmax_bubble_edge)
                                     zmax_bubble_edge = crds[dim-1];
                                 if (crds[dim-1] < zmin_spike_edge)
                                     zmin_spike_edge = crds[dim-1];
                            }
                        }
                    }
                }
            }
        }
        pp_global_max(&zmax,1);
        pp_global_min(&zmin,1);
        pp_global_max(&zmax_bubble_corner,1);
        pp_global_max(&zmax_bubble_edge,1);
        pp_global_min(&zmin_spike_corner,1);
        pp_global_min(&zmin_spike_edge,1);
	// store zmax and zmin
	zmax_intfc = zmax;
	zmin_intfc = zmin;
    // store corner value and edge value
        if (debugging("glayer"))
        {
            printf("\n\tzmin_intfc = %lf\n", zmin);
            printf("\tzmax_intfc = %lf\n", zmax);
        }
        h_bubble_intfc = (m_rho[0]<m_rho[1]) ? fabs(zmin-z0):fabs(zmax-z0);
        h_spike_intfc = (m_rho[0]<m_rho[1]) ? fabs(zmax-z0):fabs(zmin-z0);
        //corner value
        h_bubble_corner = (m_rho[0]<m_rho[1]) ? fabs(zmin_spike_corner-z0):fabs(zmax_bubble_corner-z0);
        h_spike_corner = (m_rho[0]<m_rho[1]) ? fabs(zmax_bubble_corner-z0):fabs(zmin_spike_corner-z0);
        //edge value
        h_bubble_edge = (m_rho[0]<m_rho[1]) ? fabs(zmin_spike_edge-z0):fabs(zmax_bubble_edge-z0);
        h_spike_edge = (m_rho[0]<m_rho[1]) ? fabs(zmax_bubble_edge-z0):fabs(zmin_spike_edge-z0);


	//2. get bubble/spike height of volume fraction contour
	//1%~99% or 5%~95% vol frac contour
	//Note that vf = 100% at top bdry and vf = 0 at bottom bdry

        double tol = iFparams->vol_frac_threshold;
	if (tol > 0.5)	tol = 1.0-tol;
        //if (debugging("glayer"))
	{
	    printf("vol_frac_threshold = %lf\n", tol);
	}
	//search for zmin from lowest cell centers in z-direction
        zmin = new_height_at_fraction_vd(GL[dim-1]+0.5*top_h[dim-1],+1,tol);
        //search for zmax from highest cell centers in z-direction
        zmax = new_height_at_fraction_vd(GU[dim-1]-0.5*top_h[dim-1],-1,1.0-tol);
        // store zmax and zmin
        zmax_vf = zmax;
        zmin_vf = zmin;
        if (debugging("glayer"))
        {
            printf("\nAfter calls to height_at_fraction_vd():\n");
            printf("\n\tzmin_vf = %lf\n", zmin);
            printf("\tzmax_vf = %lf\n", zmax);
        }
        h_bubble_vf = (m_rho[0]<m_rho[1]) ? fabs(zmin-z0):fabs(zmax-z0);
        h_spike_vf = (m_rho[0]<m_rho[1]) ? fabs(zmax-z0):fabs(zmin-z0);


        //3. get molecular mixing parameter theta
        //theta = <vf1vf2>/(<vf1><vf2>)

        count = 0;
        sum_vf1 = sum_vf2 = sum_vf1vf2 = 0;

        ucount = 0;
        usum_vf1 = usum_vf2 = usum_vf1vf2 = 0;

        double vf1_midPlane,vf2_midPlane,vf1vf2_midPlane,mean_vf1_midPlane,mean_vf2_midPlane,mean_vf1vf2_midPlane;
        double Dens_nb, Dens_midPlane;
        double sum_vf1_midPlane = 0, sum_vf2_midPlane = 0, sum_vf1vf2_midPlane = 0;

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            Dens = cell_center[index].m_state.m_rho;
            z_coord = cell_center[index].m_coords[dim-1];

            if (z_coord > z0-top_h[dim-1] && z_coord < z0)
            {
                count++;

                //volume/mole fraction
                vf2 = (Dens-max_dens)/(min_dens-max_dens);
                vf1 = (Dens-min_dens)/(max_dens-min_dens);
                vf1vf2 = vf1*vf2;
                sum_vf1 += vf1;
                sum_vf2 += vf2;
                sum_vf1vf2 += vf1vf2;

                Dens_nb = cell_center[d_index3d(i,j,k+1,top_gmax)].m_state.m_rho;
                Dens_midPlane = 0.5*(Dens + Dens_nb);
                vf2_midPlane = (Dens_midPlane-max_dens)/(min_dens-max_dens);
                vf1_midPlane = (Dens_midPlane-min_dens)/(max_dens-min_dens);
                vf1vf2_midPlane = vf1_midPlane*vf2_midPlane;
                sum_vf1_midPlane += vf1_midPlane;
                sum_vf2_midPlane += vf2_midPlane;
                sum_vf1vf2_midPlane += vf1vf2_midPlane;
            }

            if (z_coord > z0 && z_coord < z0+top_h[dim-1])
            {
                ucount++;

                //volume/mole fraction
                uvf2 = (Dens-max_dens)/(min_dens-max_dens);
                uvf1 = (Dens-min_dens)/(max_dens-min_dens);
                uvf1vf2 = uvf1*uvf2;
                usum_vf1 += uvf1;
                usum_vf2 += uvf2;
                usum_vf1vf2 += uvf1vf2;
            }
        }

        if (pp_numnodes() > 1)
        {
            pp_global_isum(&count,1);
            pp_global_sum(&sum_vf1,1);
            pp_global_sum(&sum_vf2,1);
            pp_global_sum(&sum_vf1vf2,1);

            pp_global_isum(&ucount,1);
            pp_global_sum(&usum_vf1,1);
            pp_global_sum(&usum_vf2,1);
            pp_global_sum(&usum_vf1vf2,1);

            pp_global_sum(&sum_vf1_midPlane,1);
            pp_global_sum(&sum_vf2_midPlane,1);
            pp_global_sum(&sum_vf1vf2_midPlane,1);
        }
        mean_vf1 = sum_vf1/count;
        mean_vf2 = sum_vf2/count;
        mean_vf1vf2 = sum_vf1vf2/count;

        umean_vf1 = usum_vf1/ucount;
        umean_vf2 = usum_vf2/ucount;
        umean_vf1vf2 = usum_vf1vf2/ucount;

        mean_vf1_midPlane = sum_vf1_midPlane/count;
        mean_vf2_midPlane = sum_vf2_midPlane/count;
        mean_vf1vf2_midPlane = sum_vf1vf2_midPlane/count;


        //4. get velocity variance
        //Var(X) = E[X*X] - E[X]*E[X]

        double u,v,w,mean_u,mean_v,mean_w,mean_uu,mean_vv,mean_ww;
        double sum_u = 0, sum_v = 0, sum_w = 0, sum_uu = 0, sum_vv = 0, sum_ww = 0;
        double sum_u_raw = 0, sum_v_raw = 0, sum_w_raw = 0;
	double sum_uu_raw = 0, sum_vv_raw = 0, sum_ww_raw = 0;
        double sum_uw_raw = 0, sum_vw_raw = 0;

	double mean_u_raw, mean_v_raw, mean_w_raw, mean_uu_raw, mean_vv_raw, mean_ww_raw, mean_uw_raw, mean_vw_raw;

        // escape velocity for laser sheet, unit: cm/s
        // v_esc = 0.5*(w+d)/(1/f)
        // where w = 532 (nm) is width of laser sheet, d = 13 (micron) is mean particle diamater of Conduct-O-Fil silver-coated hollow glass spheres, f = 30 (Hz) is the sampling rate of 2 lasers.
        double v_esc = 0.040596; //unit: cm/s
        int count_nonEsc = 0;
        count = 0;

        int ii, jj, kk;
        //average to interrogation window, whose size is 0.135*0.135 cm^2 (in uw-plane)
        int NB[3] = {1, 1, 1}; // coarse and medium grid
        //if (fabs(top_h[0] - 0.05) < 1e-12) // fine grid
        //    NB[0] = NB[2] = 2;

        for (k = 0; k <= ((kmax-kmin+1)/NB[2])-1; k++)
        for (j = 0; j <= ((jmax-jmin+1)/NB[1])-1; j++)
        for (i = 0; i <= ((imax-imin+1)/NB[0])-1; i++)
        {
            kk = (NB[2]*k)+kmin;
            jj = (NB[1]*j)+jmin;
            ii = (NB[0]*i)+imin;

            if (NB[0] == 1 && NB[1] == 1 && NB[2] == 1) // coarse and medium grid
            {
                z_coord = cell_center[d_index3d(ii,jj,kk,top_gmax)].m_coords[2];
            }
            else if (NB[0] == 2 && NB[1] == 1 && NB[2] == 2) // fine grid
            {
                z_coord = (cell_center[d_index3d(ii,jj,kk,top_gmax)].m_coords[2]+
                           cell_center[d_index3d(ii,jj,kk+1,top_gmax)].m_coords[2])/2;
            }

            if (z_coord > z0-top_h[dim-1]*NB[2] && z_coord < z0)
            {
                count++;

                // cell center of (i, j, k)
                if (NB[0] == 1 && NB[1] == 1 && NB[2] == 1) // coarse and medium grid
                {
                    u = (cell_center[d_index3d(ii,jj,kk,top_gmax)].m_state.m_U_velo_var[0]+
                         cell_center[d_index3d(ii-1,jj,kk,top_gmax)].m_state.m_U_velo_var[0])/2;
                    v = (cell_center[d_index3d(ii,jj,kk,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii,jj-1,kk,top_gmax)].m_state.m_U_velo_var[1])/2;
                    w = (cell_center[d_index3d(ii,jj,kk,top_gmax)].m_state.m_U_velo_var[2]+
                         cell_center[d_index3d(ii,jj,kk-1,top_gmax)].m_state.m_U_velo_var[2])/2;
                }
                else if (NB[0] == 2 && NB[1] == 1 && NB[2] == 2) // fine grid
                {
                    u = (cell_center[d_index3d(ii,jj,kk,top_gmax)].m_state.m_U_velo_var[0]+
                         cell_center[d_index3d(ii,jj,kk+1,top_gmax)].m_state.m_U_velo_var[0])/2;
                    v = (cell_center[d_index3d(ii,jj,kk,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii,jj,kk+1,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii+1,jj,kk+1,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii+1,jj,kk,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii,jj-1,kk,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii,jj-1,kk+1,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii+1,jj-1,kk+1,top_gmax)].m_state.m_U_velo_var[1]+
                         cell_center[d_index3d(ii+1,jj-1,kk,top_gmax)].m_state.m_U_velo_var[1])/8;
                    w = (cell_center[d_index3d(ii,jj,kk,top_gmax)].m_state.m_U_velo_var[2]+
                         cell_center[d_index3d(ii+1,jj,kk,top_gmax)].m_state.m_U_velo_var[2])/2;
                }
                else
                {
                    printf("ERROR: unrecognized grids\n");
                    clean_up(ERROR);
                }

                // temporal-average
                u /= nStep_velo_var;
                v /= nStep_velo_var;
                w /= nStep_velo_var;

                sum_u_raw += u;
                sum_v_raw += v;
                sum_w_raw += w;
                sum_uu_raw += u*u;
                sum_vv_raw += v*v;
                sum_ww_raw += w*w;
                sum_uw_raw += u*w;
                sum_vw_raw += v*w;

                if (fabs(v) <= v_esc)
                {
                    count_nonEsc++;
                    sum_u += u;
                    sum_v += v;
                    sum_w += w;
                    sum_uu += u*u;
                    sum_vv += v*v;
                    sum_ww += w*w;
                }
            }
        }

        // re-initialize m_U_velo_var to be 0
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U_velo_var[l] = 0;
        }

        if (pp_numnodes() > 1)
        {
            pp_global_isum(&count,1);
            pp_global_isum(&count_nonEsc,1);
            pp_global_sum(&sum_u,1);
            pp_global_sum(&sum_v,1);
            pp_global_sum(&sum_w,1);
            pp_global_sum(&sum_uu,1);
            pp_global_sum(&sum_vv,1);
            pp_global_sum(&sum_ww,1);

            pp_global_sum(&sum_u_raw,1);
            pp_global_sum(&sum_v_raw,1);
            pp_global_sum(&sum_w_raw,1);
            pp_global_sum(&sum_uu_raw,1);
            pp_global_sum(&sum_vv_raw,1);
            pp_global_sum(&sum_ww_raw,1);
            pp_global_sum(&sum_uw_raw,1);
            pp_global_sum(&sum_vw_raw,1);
        }
        mean_u = sum_u/count_nonEsc;
        mean_v = sum_v/count_nonEsc;
        mean_w = sum_w/count_nonEsc;
        mean_uu = sum_uu/count_nonEsc;
        mean_vv = sum_vv/count_nonEsc;
        mean_ww = sum_ww/count_nonEsc;

        mean_u_raw = sum_u_raw/count;
        mean_v_raw = sum_v_raw/count;
        mean_w_raw = sum_w_raw/count;
        mean_uu_raw = sum_uu_raw/count;
        mean_vv_raw = sum_vv_raw/count;
        mean_ww_raw = sum_ww_raw/count;
        mean_uw_raw = sum_uw_raw/count;
        mean_vw_raw = sum_vw_raw/count;

	//printf("(count, count_nonEsc) = (%d, %d)\n", count, count_nonEsc);
        //printf("raw data, (mean_u, mean_v, mean_w) = (%lf, %lf, %lf)\n", mean_u_raw, mean_v_raw, mean_w_raw);
        //printf("filter, (mean_u, mean_v, mean_w) = (%lf, %lf, %lf)\n", mean_u, mean_v, mean_w);

        //print output file
        char filename[256];
        sprintf(filename,"%s/sim_hc_params.txt",out_name);
        if (pp_mynode() == 0)
        {
            if (front->step == 0)
            {
                outfile = fopen(filename,"w");
                //fprintf(outfile,"  ts     time         tau          Agt2       hb_intfc     hs_intfc      hb_vf        hs_vf        theta1       theta2    var(u)/AgH                 var(v)/AgH   var(w)/AgH   raw_variances\n");
                fprintf(outfile, "time        t2      bubblecorner        spikecorner         bubbleedge          spikeedge\n");
            }
            else {
                outfile = fopen(filename,"a");
            }
            //fprintf(outfile,"%4d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", front->step,front->time,tau,Agt2,h_bubble_intfc,h_spike_intfc,h_bubble_vf,h_spike_vf,0.5*(mean_vf1vf2/(mean_vf1*mean_vf2) + umean_vf1vf2/(umean_vf1*umean_vf2)),mean_vf1vf2_midPlane/(mean_vf1_midPlane*mean_vf2_midPlane),0.0,(mean_uu-mean_u*mean_u)/fabs(A*g*H),(mean_vv-mean_v*mean_v)/fabs(A*g*H),(mean_ww-mean_w*mean_w)/fabs(A*g*H),(mean_uu_raw-mean_u_raw*mean_u_raw)/fabs(A*g*H),(mean_vv_raw-mean_v_raw*mean_v_raw)/fabs(A*g*H),(mean_ww_raw-mean_w_raw*mean_w_raw)/fabs(A*g*H));
            //fprintf(outfile,"%4d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", front->step,front->time,tau,Agt2,h_bubble_intfc,h_spike_intfc,h_bubble_vf,h_spike_vf,0.5*(mean_vf1vf2/(mean_vf1*mean_vf2) + umean_vf1vf2/(umean_vf1*umean_vf2)),mean_vf1vf2_midPlane/(mean_vf1_midPlane*mean_vf2_midPlane),0.0,(mean_uu_raw-mean_u_raw*mean_u_raw),(mean_vv_raw-mean_v_raw*mean_v_raw),(mean_ww_raw-mean_w_raw*mean_w_raw),(mean_uw_raw-mean_u_raw*mean_w_raw),(mean_vw_raw-mean_v_raw*mean_w_raw));
            fprintf(outfile, "%e\t%e\t%e\t%e\t%e\t%e\n", front->time, front->time*front->time, h_bubble_corner, h_spike_corner, h_bubble_edge, h_spike_edge);

            fclose(outfile);
        }


        if (debugging("trace"))
            printf("Leave computeRTParameters()\n");
}*/  /* end computeRTParameters */


double Incompress_Solver_Smooth_3D_Cartesian::new_height_at_fraction_vd(
        double                  h0,
        int                     dir,
        double                  fraction)
{
        double            height = h0; //initial guess for height
	const double	  dh = top_h[dim-1];
        static double     old_frac, new_frac;

        if (debugging("glayer"))
	    printf("\nEntering new_height_at_fraction_vd() for %s\n", (dir==-1)?"zmax":"zmin");

        old_frac = 0;
        new_accumulate_fractions_in_layer_vd(height,&old_frac);
	//zmin <-> dir=+1
        //zmax <-> dir=-1
	//mixing layer too close to the bottom, return bottom.
        if (dir == +1  &&  old_frac > fraction)
        {
            if (debugging("glayer"))
                (void) printf("WARNING in new_height_at_fraction(), "
                          "frac > %lf occurred at initial height of %lf.\n",
                          fraction,height);
            return front->rect_grid->GL[dim-1];
        }
        //mixing layer too close to the top, return top.
        else if (dir == -1  &&  old_frac < fraction)
        {
            if (debugging("gfrac"))
                (void) printf("WARNING in new_height_at_fraction(), "
                          "frac < %lf occurred at initial height of %lf.\n",
                          fraction,height);
            return front->rect_grid->GU[dim-1];
        }

        while (height >= front->rect_grid->GL[dim-1] &&
               height <= front->rect_grid->GU[dim-1] )
        {
            new_accumulate_fractions_in_layer_vd(height+dir*dh,&new_frac);
            if ( ((old_frac <= fraction) && (fraction <= new_frac)) ||
                 ((new_frac <= fraction) && (fraction <= old_frac)) )
            {   // we've bracketed the height, so interpolate
                height += dir*dh*(fraction-old_frac)/(new_frac-old_frac);
                break;
            }
            else
            {
                old_frac = new_frac;
                height += dir*dh;
            }
        }

        if (debugging("glayer"))
            printf("Leaving new_height_at_fraction_vd()\n");

        return height;
}               /*end new_height_at_fraction_vd*/


void Incompress_Solver_Smooth_3D_Cartesian::new_accumulate_fractions_in_layer_vd(
        double                  height,
        double                  *frac)
{
	const double		   *origin = NULL;
        const RECT_GRID            *rgrid = front->rect_grid;
        const enum intfc_geometry  geom = (origin == NULL ? planar : spherical);
        int                        i, j, k, index, count=0;
        const double               *L = rgrid->L;
        const double               *U = rgrid->U;
        double                     z_coord=0, Dens=0, upper_fluid_volume=0;

//        debug_print("glayer","Entered new_accumulate_fractions_in_layer_vd(), h = %lf\n",height);

        if (debugging("glayer"))
        {
            (void) printf("proc #%d:\n", pp_mynode());
            (void) printf("h = %lf\n",height);
        }

        switch(geom)
        {
        case planar:
            count = 0;
            upper_fluid_volume = 0;
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
                index = d_index3d(i,j,k,top_gmax);
              	Dens = cell_center[index].m_state.m_rho;
               	z_coord = cell_center[index].m_coords[dim-1];
		//grids centers lie at the same level as height
		if (fabs(z_coord-height) <= 1e-12)
		{
	    	    count++;
                    upper_fluid_volume += (Dens-m_rho[1])/(m_rho[0]-m_rho[1]);
		}
	    }
            if (debugging("glayer"))
            {
                (void) printf("Accumulation of layer fractions "
                              "at h = %lf completed.\n",height);
                (void) printf("count = %d\n", count);
                (void) printf("total upper_fluid_volume = %lf\n",
                              upper_fluid_volume);
            }
            if (pp_numnodes() > 1)
            {
                pp_global_isum(&count,1);
                pp_global_sum(&upper_fluid_volume, 1);
            }
            if (debugging("glayer"))
            {
                (void) printf("After calling pp_global_sum:\n");
		(void) printf("count = %d\n", count);
                (void) printf("total upper_fluid_volume = %lf\n",
                              upper_fluid_volume);
            }
            *frac=(count==0)? 0 : upper_fluid_volume/count;
            break;

        case spherical:
            screen("ERROR in new_accumulate_fractions_in_layer_vd(),\n");
            screen("spherical 1D/2D/3D not supported. In gas/gintext.c, change\n");
            screen("new_accumulate_fraction_in_layer_vd() function calls\n");
            screen("to accumulate_fraction_in_layer_vd() function calls\n");
            clean_up(ERROR);
            break;
        }

        if (debugging("glayer"))
            (void) printf("average upper_fluid_volume = %lf\n",*frac);
//        debug_print("glayer","Leaving new_accumulate_fractions_in_layer_vd()\n");

        return;
}               /*end new_accumulate_fractions_in_layer_vd*/


double Incompress_Solver_Smooth_3D_Cartesian::getVorticityX(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dy,dz;
	double vorticity;

	dy = top_h[1];
	dz = top_h[2];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i,j-1,k,top_gmax);
	index01 = d_index3d(i,j+1,k,top_gmax);
	index10 = d_index3d(i,j,k-1,top_gmax);
	index11 = d_index3d(i,j,k+1,top_gmax);
	v00 = -cell_center[index00].m_state.m_U[2];
	v01 =  cell_center[index01].m_state.m_U[2];
	v10 =  cell_center[index10].m_state.m_U[1];
	v11 = -cell_center[index11].m_state.m_U[1];

	vorticity = (v00 + v01)/2.0/dz + (v10 + v11)/2.0/dy;
	return vorticity;
}	/* end getVorticityX */

double Incompress_Solver_Smooth_3D_Cartesian::getVorticityY(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dz;
	double vorticity;

	dx = top_h[0];
	dz = top_h[2];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i,j,k-1,top_gmax);
	index01 = d_index3d(i,j,k+1,top_gmax);
	index10 = d_index3d(i-1,j,k,top_gmax);
	index11 = d_index3d(i+1,j,k,top_gmax);
	v00 = -cell_center[index00].m_state.m_U[0];
	v01 =  cell_center[index01].m_state.m_U[0];
	v10 =  cell_center[index10].m_state.m_U[2];
	v11 = -cell_center[index11].m_state.m_U[2];

	vorticity = (v00 + v01)/2.0/dx + (v10 + v11)/2.0/dz;
	return vorticity;
}	/* end getVorticityY */

double Incompress_Solver_Smooth_3D_Cartesian::getVorticityZ(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dy;
	double vorticity;

	dx = top_h[0];
	dy = top_h[1];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i-1,j,k,top_gmax);
	index01 = d_index3d(i+1,j,k,top_gmax);
	index10 = d_index3d(i,j-1,k,top_gmax);
	index11 = d_index3d(i,j+1,k,top_gmax);
	v00 = -cell_center[index00].m_state.m_U[1];
	v01 =  cell_center[index01].m_state.m_U[1];
	v10 =  cell_center[index10].m_state.m_U[0];
	v11 = -cell_center[index11].m_state.m_U[0];

	vorticity = (v00 + v01)/2.0/dy + (v10 + v11)/2.0/dx;
	return vorticity;
}	/* end getVorticityZ */


void Incompress_Solver_Smooth_3D_Cartesian::copyMeshStates()
{
	int i,j,k,d,index;
	double **vel = field->vel;
	double *pres = field->pres;
	double *vort = field->vort;
	double **vort3d = field->vort3d;

	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	for (k = kmin; k <= kmax; ++k)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    if (ifluid_comp(top_comp[index]))
	    {
		pres[index] = cell_center[index].m_state.m_P;
	    	vel[0][index] = cell_center[index].m_state.m_U[0];
	    	vel[1][index] = cell_center[index].m_state.m_U[1];
	    	vel[2][index] = cell_center[index].m_state.m_U[2];
		vort3d[0][index] = getVorticityX(i,j,k);
		vort3d[1][index] = getVorticityY(i,j,k);
		vort3d[2][index] = getVorticityZ(i,j,k);
	    }
	    else
	    {
	    	pres[index] = 0.0;
		for (d = 0; d < 3; ++d)
		{
		    vel[d][index] = 0.0;
		    vort3d[d][index] = 0.0;
		}
	    }
	}
	FT_ParallelExchGridArrayBuffer(pres,front);
	FT_ParallelExchGridArrayBuffer(vel[0],front);
	FT_ParallelExchGridArrayBuffer(vel[1],front);
	FT_ParallelExchGridArrayBuffer(vel[2],front);
	FT_ParallelExchGridArrayBuffer(vort3d[0],front);
	FT_ParallelExchGridArrayBuffer(vort3d[1],front);
	FT_ParallelExchGridArrayBuffer(vort3d[2],front);
}	/* end copyMeshStates */


void Incompress_Solver_Smooth_3D_Cartesian::copyMeshStates_vd()
{
        int i,j,k,d,l,index;
        double **vel = field->vel;
//        double *pres = field->pres;
//        double *vort = field->vort;
//        double **vort3d = field->vort3d;
        // for vd
//        double *dens = field->dens;
//        double *dens_old = field->dens_old;
//        double *conc = field->conc;

        for (i = imin; i <= imax; ++i)
        for (j = jmin; j <= jmax; ++j)
        for (k = kmin; k <= kmax; ++k)
        {
            index = d_index3d(i,j,k,top_gmax);
            if (ifluid_comp(top_comp[index]))
            {
//                pres[index] = cell_center[index].m_state.m_P;
                vel[0][index] = cell_center[index].m_state.m_U[0];
                vel[1][index] = cell_center[index].m_state.m_U[1];
                vel[2][index] = cell_center[index].m_state.m_U[2];
//                vort3d[0][index] = getVorticityX(i,j,k);
//                vort3d[1][index] = getVorticityY(i,j,k);
//                vort3d[2][index] = getVorticityZ(i,j,k);
                //for vd
//                dens[index] = cell_center[index].m_state.m_rho;
//                dens_old[index] = cell_center[index].m_state.m_rho_old;
//                conc[index] = cell_center[index].m_state.m_c;
            }
            else
            {
//                pres[index] = 0;
                for (d = 0; d < 3; ++d)
                {
                    vel[d][index] = 0;
//                    vort3d[d][index] = 0;
                }
                //for vd
//                dens[index] = dens_old[index] = 0;
//                conc[index] = 0;
                assert (false);
            }
        }
//        FT_ParallelExchGridArrayBuffer(pres,front);
        FT_ParallelExchGridArrayBuffer(vel[0],front);
        FT_ParallelExchGridArrayBuffer(vel[1],front);
        FT_ParallelExchGridArrayBuffer(vel[2],front);
        // removal tag: HAOZ
        enforceReflectionState(vel);// on m_U
        for (l = 0; l < dim; ++l)
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_U[l] = vel[l][index];
        }
        // Update Reflection Boundary Condition TODO && FIXME: Solute_Reflect contribute divergence non-free situation.
        //Solute_Reflect(0, vel[0]);
        //Solute_Reflect(1, vel[1]);
        //Solute_Reflect(2, vel[2]);
//        FT_ParallelExchGridArrayBuffer(vort3d[0],front);
//        FT_ParallelExchGridArrayBuffer(vort3d[1],front);
//        FT_ParallelExchGridArrayBuffer(vort3d[2],front);
        // for vd
//        FT_ParallelExchGridArrayBuffer(dens,front);
//        FT_ParallelExchGridArrayBuffer(dens_old,front);
//        FT_ParallelExchGridArrayBuffer(conc,front);

        if (debugging("interpolate") && pp_mynode()==0 && front->step==20)
        {
            for (i = imin; i <= imax; ++i)
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index3d(i,j,14,top_gmax);
                printf("For (%d,%d,14) on proc #0, vel[index=%d] = (%12.8g,%12.8g,%12.8g) at ts = %d\n",
                       i,j,index,vel[0][index],vel[1][index],vel[2][index],front->step);
            }
        }

        if (debugging("interpolate") && pp_mynode()==1 && front->step==20)
        {
            for (i = imin; i <= imax; ++i)
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index3d(i,j,7,top_gmax);
                printf("For (%d,%d,7) on proc #1, vel[index=%d] = (%12.8g,%12.8g,%12.8g) at ts = %d\n",
                       i,j,index,vel[0][index],vel[1][index],vel[2][index],front->step);
            }
        }

        printf("**********************************************************\n");
        printf("**********************************************************\n");
        printf("**********************************************************\n");
        printf("**********************************************************\n");

        if (debugging("interpolate") && pp_mynode()==0 && front->step==20)
        {
            for (i = imin; i <= imax; ++i)
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index3d(i,j,7,top_gmax);
                printf("For (%d,%d,7) on proc #0, vel[index=%d] = (%12.8g,%12.8g,%12.8g) at ts = %d\n",
                       i,j,index,vel[0][index],vel[1][index],vel[2][index],front->step);
            }
        }

        if (debugging("interpolate") && pp_mynode()==1 && front->step==20)
        {
            for (i = imin; i <= imax; ++i)
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index3d(i,j,0,top_gmax);
                printf("For (%d,%d,0) on proc #1, vel[index=%d] = (%12.8g,%12.8g,%12.8g) at ts = %d\n",
                       i,j,index,vel[0][index],vel[1][index],vel[2][index],front->step);
            }
        }

} /* end copyMeshStates_vd */


void Incompress_Solver_Smooth_3D_Cartesian::
	compDiffWithSmoothProperty_1st_decoupled(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
	int i,j,k,l,nb,icoords[MAXD];
        L_STATE state;
        INTERFACE *intfc = front->interf;
	double coords[MAXD], crx_coords[MAXD];
	double coeff[6],mu[6],mu0,rho,corner[6],rhs,U_nb[6];
        double speed;
        double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	double (*getStateVel[3])(POINTER) =
			{getStateXvel,getStateYvel,getStateZvel};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;

        setIndexMap();

	max_speed = 0.0;

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	for (l = 0; l < dim; ++l)
	{
            PETSc solver;
            solver.Create(ilower, iupper-1, 7, 7);
	    solver.Reset_A();
	    solver.Reset_b();
	    solver.Reset_x();

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
            	I = ijk_to_I[i][j][k];
            	if (I == -1) continue;

            	index = d_index3d(i,j,k,top_gmax);
            	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

		icoords[0] = i;
		icoords[1] = j;
		icoords[2] = k;
		comp = top_comp[index];

            	I_nb[0] = ijk_to_I[i-1][j][k]; //west
            	I_nb[1] = ijk_to_I[i+1][j][k]; //east
            	I_nb[2] = ijk_to_I[i][j-1][k]; //south
            	I_nb[3] = ijk_to_I[i][j+1][k]; //north
            	I_nb[4] = ijk_to_I[i][j][k-1]; //lower
            	I_nb[5] = ijk_to_I[i][j][k+1]; //upper


            	mu0 = cell_center[index].m_state.m_mu;
            	rho = cell_center[index].m_state.m_rho;

            	for (nb = 0; nb < 6; nb++)
            	{
                    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                                comp,&intfc_state,&hs,crx_coords,m_t_new) &&
                                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		    {
			U_nb[nb] = getStateVel[l](intfc_state);
			if (wave_type(hs) == DIRICHLET_BOUNDARY ||
			    wave_type(hs) == NEUMANN_BOUNDARY)
			    mu[nb] = mu0;
			else
			    mu[nb] = 1.0/2*(mu0 +
				cell_center[index_nb[nb]].m_state.m_mu);
		    }
                    else
		    {
                    	U_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[l];
			mu[nb] = 1.0/2*(mu0 +
				cell_center[index_nb[nb]].m_state.m_mu);
		    }
            	}

            	coeff[0] = 0.5*m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            	coeff[1] = 0.5*m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            	coeff[2] = 0.5*m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            	coeff[3] = 0.5*m_dt/rho*mu[3]/(top_h[1]*top_h[1]);
            	coeff[4] = 0.5*m_dt/rho*mu[4]/(top_h[2]*top_h[2]);
            	coeff[5] = 0.5*m_dt/rho*mu[5]/(top_h[2]*top_h[2]);

            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, state);

            	solver.Set_A(I,I,1+coeff[0]+coeff[1]+coeff[2]+coeff[3]+
				coeff[4]+coeff[5]);
		rhs = (1-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*
		      		cell_center[index].m_state.m_U[l];

		for(nb = 0; nb < 6; nb++)
		{
		    if(I_nb[nb] != -1)
		    {
			solver.Set_A(I,I_nb[nb],-coeff[nb]);
			rhs += coeff[nb]*U_nb[nb];
		    }
		    else
			rhs += 2.0*coeff[nb]*U_nb[nb];
		}
		rhs += m_dt*state.m_U[l];
		rhs += m_dt*cell_center[index].m_state.f_surf[l];
		rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho;

		solver.Set_b(I, rhs);
            }

            solver.SetMaxIter(40000);
            solver.SetTol(1e-10);

	    start_clock("Befor Petsc solve");
            solver.Solve_GMRES();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

	    //if(rel_residual > 1)
	    //{
	    //	printf("\nThe solution diverges! The residual is %le. Solve again using GMRES!\n",rel_residual);
	    //	solver.Reset_x();
	    //	solver.Solve_GMRES();
	    //	solver.GetNumIterations(&num_iter);
            //	solver.GetFinalRelativeResidualNorm(&rel_residual);
	    //}
	    stop_clock("After Petsc solve");

            // get back the solution
            solver.Get_x(x);

            if (debugging("PETSc"))
                (void) printf("L_CARTESIAN::"
			"compDiffWithSmoothProperty_1st_decoupled: "
                        "num_iter = %d, rel_residual = %le. \n",
                        num_iter,rel_residual);

	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ijk_to_I[i][j][k];
                index = d_index3d(i,j,k,top_gmax);
                if (I >= 0)
                {
                    cell_center[index].m_state.m_U[l] = x[I-ilower];
                }
                else
                {
                    cell_center[index].m_state.m_U[l] = 0.0;
                }
            }

	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }
	for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
            speed = fabs(cell_center[index].m_state.m_U[0]) +
                    fabs(cell_center[index].m_state.m_U[1]) +
                    fabs(cell_center[index].m_state.m_U[2]);
            if (speed > max_speed)
                    max_speed = speed;
	}
        pp_global_max(&max_speed,1);

        FT_FreeThese(1,x);
}       /* end compDiffWithSmoothProperty3d_1st_decoupled */


void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmI(void)
{
        int i,j,k,index;

        for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_P += cell_center[index].m_state.m_phi;
	    array[index] = cell_center[index].m_state.m_P;
	}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_P = array[index];
	    cell_center[index].m_state.m_q = array[index];
	}
}        /* end computePressurePmI3d */


void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmII(void)
{
        int i,j,k,index;
        double mu0;

        for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
            mu0 = 0.5*cell_center[index].m_state.m_mu;
            cell_center[index].m_state.m_P = (cell_center[index].m_state.m_q +
				              cell_center[index].m_state.m_phi -
                        	              mu0*cell_center[index].m_state.div_U);
	    array[index] = cell_center[index].m_state.m_P;
	}
	scatMeshArray();
	for(k = 0; k <= top_gmax[2]; k++)
	for(j = 0; j <= top_gmax[1]; j++)
	for(i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    cell_center[index].m_state.m_P = array[index];
	    cell_center[index].m_state.m_q = array[index];
	}
}        /* end computePressurePmII3d */


void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmII_vd(void)
{
}


void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmIII(void)
{
        int i,j,k,index;
        double mu0;

	for(k = 0; k <= top_gmax[2]; k++)
	for(j = 0; j <= top_gmax[1]; j++)
	for(i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    mu0 = 0.5*cell_center[index].m_state.m_mu;
	    cell_center[index].m_state.m_P =
		cell_center[index].m_state.m_phi -
		accum_dt*mu0*cell_center[index].m_state.div_U;
	    cell_center[index].m_state.m_q = 0.0;
	}
}        /* end computePressurePmIII3d */


void Incompress_Solver_Smooth_3D_Cartesian::computePressure(void)
{
	switch (iFparams->num_scheme)
	{
	case BELL_COLELLA:
	    computePressurePmI();
	    break;
	case KIM_MOIN:
	    computePressurePmII();
	    break;
	case SIMPLE:
	case PEROT_BOTELLA:
	    computePressurePmIII();
	    break;
	}
	computeGradientQ();
}	/* end computePressure */

void Incompress_Solver_Smooth_3D_Cartesian::computePressure_vd(void)
{
}


void Incompress_Solver_Smooth_3D_Cartesian::computePressure_MAC_vd(void)
{
	switch (iFparams->num_scheme)
	{
	case BELL_COLELLA:
	    computePressurePmI();
//            computePressurePmII();
	    break;
	case KIM_MOIN:
	    computePressurePmII();
	    break;
	case SIMPLE:
	case PEROT_BOTELLA:
	    computePressurePmIII();
	    break;
	}
	computeGradientQ_MAC_vd();
} /* end computePressure_MAC_vd */


//set q_z = 0
void Incompress_Solver_Smooth_3D_Cartesian::computePressure_MAC_zeroW_vd(void)
{
        switch (iFparams->num_scheme)
        {
        case BELL_COLELLA:
            computePressurePmII();
            break;
        case KIM_MOIN:
            computePressurePmII();
            break;
        case SIMPLE:
        case PEROT_BOTELLA:
            computePressurePmIII();
            break;
        }
        computeGradientQ_MAC_zeroW_vd();
} /* end computePressure_MAC_zeroW_vd */


//set q_y = 0
void Incompress_Solver_Smooth_3D_Cartesian::computePressure_MAC_zeroV_vd(void)
{
        switch (iFparams->num_scheme)
        {
        case BELL_COLELLA:
            computePressurePmII();
            break;
        case KIM_MOIN:
            computePressurePmII();
            break;
        case SIMPLE:
        case PEROT_BOTELLA:
            computePressurePmIII();
            break;
        }
        computeGradientQ_MAC_zeroV_vd();
} /* end computePressure_MAC_zeroV_vd */


void Incompress_Solver_Smooth_3D_Cartesian::computeGradientQ(void)
{
	int i,j,k,l,index;
	double *grad_q;
	int icoords[MAXD];

	for (k = kmin; k <= kmax; ++k)
	for (j = jmin; j <= jmax; ++j)
	for (i = imin; i <= imax; ++i)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = cell_center[index].m_state.m_q;
	}
	scatMeshArray();

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    grad_q = cell_center[index].m_state.grad_q;
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    computeFieldPointGrad(icoords,array,grad_q);
	}
	for (l = 0; l < dim; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
		array[index] = cell_center[index].m_state.grad_q[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
		cell_center[index].m_state.grad_q[l] = array[index];
	    }
	}
} /* end computeGradientQ */


void Incompress_Solver_Smooth_3D_Cartesian::computeGradientQ_MAC_vd(void)
{
        bool bNoBoundary[6];
	int index_nb[6];
	int i,j,k,l,index;
	double rho,g;
	int icoords[MAXD];
        L_STATE gravity;
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        COMPONENT comp;

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
            index = d_index3d(i,j,k,top_gmax);
            comp = top_comp[index];

     	    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
            //4 directions
            bNoBoundary[0] = YES;
            bNoBoundary[1] = YES;
            bNoBoundary[2] = YES;
            bNoBoundary[3] = YES;
            //LOWER
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
                    comp,&intfc_state,&hs,crx_coords,m_t_int) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                bNoBoundary[4] = NO;
            else
                bNoBoundary[4] = YES;
            //UPPER
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
                    comp,&intfc_state,&hs,crx_coords,m_t_int) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                bNoBoundary[5] = NO;
            else
                bNoBoundary[5] = YES;

            if (!bNoBoundary[5]) //cells on UPPER bdry
            {
                //get gravity force terms
                getRectangleCenter(index, coords);
                computeSourceTerm(coords, gravity);

                g = gravity.m_U[2];
                rho = 0.5*(cell_center[index].m_state.m_rho + cell_center[index].m_state.m_rho_old);

                cell_center[index].m_state.grad_q[0] = (cell_center[index_nb[1]].m_state.m_q -
                                                        cell_center[index].m_state.m_q)/top_h[0];
                cell_center[index].m_state.grad_q[1] = (cell_center[index_nb[3]].m_state.m_q -
                                                        cell_center[index].m_state.m_q)/top_h[1];
                cell_center[index].m_state.grad_q[2] = rho*g;//p_z = rho*g on UPPER/LOWER bdry
            }
            else //other cells
            {
                cell_center[index].m_state.grad_q[0] = (cell_center[index_nb[1]].m_state.m_q -
                                                        cell_center[index].m_state.m_q)/top_h[0];
                cell_center[index].m_state.grad_q[1] = (cell_center[index_nb[3]].m_state.m_q -
                                                        cell_center[index].m_state.m_q)/top_h[1];
                cell_center[index].m_state.grad_q[2] = (cell_center[index_nb[5]].m_state.m_q -
                                                        cell_center[index].m_state.m_q)/top_h[2];
            }
	}
        //scatter states
        /*
         * removal tag: HAOZ
         * there is a vector field grad_q
         * should this be done in a fashion same as for velocity?????
         * */
	for (l = 0; l < dim; ++l)
	{
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
		array[index] = cell_center[index].m_state.grad_q[l];
	    }
	    scatMeshArray();
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index = d_index3d(i,j,k,top_gmax);
		cell_center[index].m_state.grad_q[l] = array[index];
        vecarray[l][index] = array[index];
	    }
	}
    enforceReflectionState(vecarray);// on grad_q
    for (l = 0; l < dim; l++)
    for (k = 0; k <= top_gmax[2]; k++)
    for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
    {
        index = d_index3d(i,j,k,top_gmax);
        cell_center[index].m_state.grad_q[l] = vecarray[l][index];
    }
} /* end computeGradientQ_MAC_vd */


//set q_z = 0
void Incompress_Solver_Smooth_3D_Cartesian::computeGradientQ_MAC_zeroW_vd(void)
{
        bool bNoBoundary[6];
        int index_nb[6];
        int i,j,k,l,index;
        double rho,g;
        int icoords[MAXD];
        L_STATE gravity;
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        COMPONENT comp;

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            comp = top_comp[index];

            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);

            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            // 4 directions
            bNoBoundary[0] = YES;
            bNoBoundary[1] = YES;
            bNoBoundary[2] = YES;
            bNoBoundary[3] = YES;
            // LOWER
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
                    comp,&intfc_state,&hs,crx_coords,m_t_int) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                bNoBoundary[4] = NO;
            else
                bNoBoundary[4] = YES;
            // UPPER
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
                    comp,&intfc_state,&hs,crx_coords,m_t_int) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                bNoBoundary[5] = NO;
            else
                bNoBoundary[5] = YES;

            if (!bNoBoundary[5]) //cells on UPPER bdry
            {
                //get gravity force terms
                getRectangleCenter(index, coords);
                computeSourceTerm(coords, gravity);

                g = gravity.m_U[2];
                rho = 0.5*(cell_center[index].m_state.m_rho + cell_center[index].m_state.m_rho_old);

                cell_center[index].m_state.grad_q[0] = (cell_center[index_nb[1]].m_state.m_q -
                                                        cell_center[index].m_state.m_q)/top_h[0];
                cell_center[index].m_state.grad_q[1] = (cell_center[index_nb[3]].m_state.m_q -
                                                        cell_center[index].m_state.m_q)/top_h[1];
                cell_center[index].m_state.grad_q[2] = rho*g; //p_z = rho*g on UPPER/LOWER bdry
            }
            else //other cells
            {
                cell_center[index].m_state.grad_q[0] = (cell_center[index_nb[1]].m_state.m_q -
                                                        cell_center[index].m_state.m_q)/top_h[0];
                cell_center[index].m_state.grad_q[1] = (cell_center[index_nb[3]].m_state.m_q -
                                                        cell_center[index].m_state.m_q)/top_h[1];
                cell_center[index].m_state.grad_q[2] = (cell_center[index_nb[5]].m_state.m_q -
                                                        cell_center[index].m_state.m_q)/top_h[2];
            }
        }
        //scatter states
        for (l = 0; l < dim; ++l)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.grad_q[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.grad_q[l] = array[index];
            }
        }
        //set q_z = 0
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.grad_q[2] = 0.0;
        }
} /* end computeGradientQ_MAC_zeroW_vd */


//set q_y = 0
void Incompress_Solver_Smooth_3D_Cartesian::computeGradientQ_MAC_zeroV_vd(void)
{
        bool bNoBoundary[6];
        int index_nb[6];
        int i,j,k,l,index;
        double rho,g;
        int icoords[MAXD];
        L_STATE gravity;
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        COMPONENT comp;

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            comp = top_comp[index];

            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);

            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            // 4 directions
            bNoBoundary[0] = YES;
            bNoBoundary[1] = YES;
            bNoBoundary[2] = YES;
            bNoBoundary[3] = YES;
            // LOWER
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
                    comp,&intfc_state,&hs,crx_coords,m_t_int) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                bNoBoundary[4] = NO;
            else
                bNoBoundary[4] = YES;
            // UPPER
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
                    comp,&intfc_state,&hs,crx_coords,m_t_int) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                bNoBoundary[5] = NO;
            else
                bNoBoundary[5] = YES;

            if (!bNoBoundary[5]) //cells on UPPER bdry
            {
                //get gravity force terms
                getRectangleCenter(index, coords);
                computeSourceTerm(coords, gravity);

                g = gravity.m_U[2];
                rho = 0.5*(cell_center[index].m_state.m_rho + cell_center[index].m_state.m_rho_old);

                cell_center[index].m_state.grad_q[0] = (cell_center[index_nb[1]].m_state.m_q -
                                                        cell_center[index].m_state.m_q)/top_h[0];
                cell_center[index].m_state.grad_q[1] = (cell_center[index_nb[3]].m_state.m_q -
                                                        cell_center[index].m_state.m_q)/top_h[1];
                cell_center[index].m_state.grad_q[2] = rho*g; //p_z = rho*g on UPPER/LOWER bdry
            }
            else //other cells
            {
                cell_center[index].m_state.grad_q[0] = (cell_center[index_nb[1]].m_state.m_q -
                                                        cell_center[index].m_state.m_q)/top_h[0];
                cell_center[index].m_state.grad_q[1] = (cell_center[index_nb[3]].m_state.m_q -
                                                        cell_center[index].m_state.m_q)/top_h[1];
                cell_center[index].m_state.grad_q[2] = (cell_center[index_nb[5]].m_state.m_q -
                                                        cell_center[index].m_state.m_q)/top_h[2];
            }
        }
        //scatter states
        for (l = 0; l < dim; ++l)
        {
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                array[index] = cell_center[index].m_state.grad_q[l];
            }
            scatMeshArray();
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
                cell_center[index].m_state.grad_q[l] = array[index];
            }
        }
        //set q_y = 0
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.grad_q[1] = 0.0;
        }
} /* end computeGradientQ_MAC_zeroV_vd */


#define		MAX_TRI_FOR_INTEGRAL		100
void Incompress_Solver_Smooth_3D_Cartesian::surfaceTension(
	double *d,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *force,
	double dd)
{
    	/*
    	double cellVolume = top_h[0]*top_h[1]*top_h[2];
	int i,j,k,num_tris;
	TRI *tri,*tri_list[MAX_TRI_FOR_INTEGRAL];
	double kappa_tmp,kappa,mag_nor,area,delta;
	double median[MAXD],nor[MAXD];
	POINT *p;

	TriAndFirstRing(hse,hs,&num_tris,tri_list);
	for (i = 0; i < num_tris; ++i)
	{
	    kappa = 0.0;
	    tri = tri_list[i];
	    for (j = 0; j < 3; ++j) median[j] = 0.0;
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		for (k = 0; k < 3; ++k)
		    median[k] += Coords(p)[k];
	    	GetFrontCurvature(p,Hyper_surf_element(tri),hs,
				&kappa_tmp,front);
		kappa += kappa_tmp;
		nor[j] = Tri_normal(tri)[j];
	    }
	    kappa /= 3.0;
	    mag_nor = mag_vector(nor,3);
	    area = 0.5*mag_nor;
	    for (j = 0; j < 3; ++j)
	    {
		nor[j] /= mag_nor;
		median[j] /= 3.0;
	    }
	    delta = smoothedDeltaFunction(coords,median);
	    if (delta == 0.0) continue;
	    for (j = 0; j < dim; ++j)
	    {
		force[j] += delta*sigma*area*kappa*nor[j];
	    }
	}
	for (i = 0; i < dim; i++)
	    force[j] /= cellVolume;
	*/
}	/* end surfaceTension3d */

void Incompress_Solver_Smooth_3D_Cartesian::setInitialCondition()
{
	int i;
	COMPONENT comp;
	double coords[MAXD];

	FT_MakeGridIntfc(front);
	setDomain();

        m_rho[0] = iFparams->rho1;
        m_rho[1] = iFparams->rho2;
        m_mu[0] = iFparams->mu1;
        m_mu[1] = iFparams->mu2;
	m_comp[0] = iFparams->m_comp1;
	m_comp[1] = iFparams->m_comp2;
	m_smoothing_radius = iFparams->smoothing_radius;
	m_sigma = iFparams->surf_tension;
	mu_min = rho_min = HUGE;
	for (i = 0; i < 2; ++i)
	{
	    if (ifluid_comp(m_comp[i]))
	    {
        	mu_min = std::min(mu_min,m_mu[i]);
        	rho_min = std::min(rho_min,m_rho[i]);
	    }
	}

	// Initialize state at cell_center
        for (i = 0; i < cell_center.size(); i++)
        {
            getRectangleCenter(i, coords);
	    cell_center[i].m_state.setZero();
	    comp = top_comp[i];
	    if (getInitialState != NULL)
	    	(*getInitialState)(comp,coords,cell_center[i].m_state,dim,
						iFparams);
        }
	computeGradientQ();
        copyMeshStates();
	setAdvectionDt();
}       /* end setInitialCondition */

double Incompress_Solver_Smooth_3D_Cartesian::computeFieldPointDiv(
        int *icoords,
        double **field)
{
        int index;
        COMPONENT comp;
        int i,j,k,nb;
	int index_nb[6];
        double div,u_nb[2],v_nb[2],w_nb[2];
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
      	double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};

	i = icoords[0];
	j = icoords[1];
	k = icoords[2];
	index = d_index3d(i,j,k,top_gmax);
        comp = top_comp[index];

	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,WEST,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    u_nb[0] = getStateXvel(intfc_state);
	else
	    u_nb[0] = (field[0][index] + field[0][index_nb[0]])/2.0;

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,EAST,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    u_nb[1] = getStateXvel(intfc_state);
	else
	    u_nb[1] = (field[0][index] + field[0][index_nb[1]])/2.0;

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,SOUTH,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    v_nb[0] = getStateYvel(intfc_state);
	else
	    v_nb[0] = (field[1][index] + field[1][index_nb[2]])/2.0;

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,NORTH,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    v_nb[1] = getStateYvel(intfc_state);
	else
	    v_nb[1] = (field[1][index] + field[1][index_nb[3]])/2.0;

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    w_nb[0] = getStateZvel(intfc_state);
	else
	    w_nb[0] = (field[2][index] + field[2][index_nb[4]])/2.0;

	if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
		comp,&intfc_state,&hs,crx_coords,m_t_new) &&
	        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    w_nb[1] = getStateZvel(intfc_state);
	else
	    w_nb[1] = (field[2][index] + field[2][index_nb[5]])/2.0;

        div = (u_nb[1]-u_nb[0])/top_h[0] + (v_nb[1]-v_nb[0])/top_h[1] + (w_nb[1]-w_nb[0])/top_h[2];
        return div;
} /* end computeFieldPointDiv */


double Incompress_Solver_Smooth_3D_Cartesian::computeFieldPointDiv_Neumann_vd(
        int *icoords,
        double **field)
{
        int index;
        COMPONENT comp;
        int i,j,k,nb;
        int index_nb[6];
        double div,u_nb[2],v_nb[2],w_nb[2];
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;

        i = icoords[0];
        j = icoords[1];
        k = icoords[2];
        index = d_index3d(i,j,k,top_gmax);
        comp = top_comp[index];

        index_nb[0] = d_index3d(i-1,j,k,top_gmax);
        index_nb[1] = d_index3d(i+1,j,k,top_gmax);
        index_nb[2] = d_index3d(i,j-1,k,top_gmax);
        index_nb[3] = d_index3d(i,j+1,k,top_gmax);
        index_nb[4] = d_index3d(i,j,k-1,top_gmax);
        index_nb[5] = d_index3d(i,j,k+1,top_gmax);

        //use homogeneous Neumann B.C.
        if (FT_StateStructAtGridCrossing_tmp(front,icoords,WEST,
                comp,&intfc_state,&hs,crx_coords,m_t_new) &&
                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
            u_nb[0] = field[0][index];
        else
            u_nb[0] = (field[0][index] + field[0][index_nb[0]])/2.0;

        if (FT_StateStructAtGridCrossing_tmp(front,icoords,EAST,
                comp,&intfc_state,&hs,crx_coords,m_t_new) &&
                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
            u_nb[1] = field[0][index];
        else
            u_nb[1] = (field[0][index] + field[0][index_nb[1]])/2.0;

        if (FT_StateStructAtGridCrossing_tmp(front,icoords,SOUTH,
                comp,&intfc_state,&hs,crx_coords,m_t_new) &&
                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
            v_nb[0] = field[1][index];
        else
            v_nb[0] = (field[1][index] + field[1][index_nb[2]])/2.0;

        if (FT_StateStructAtGridCrossing_tmp(front,icoords,NORTH,
                comp,&intfc_state,&hs,crx_coords,m_t_new) &&
                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
            v_nb[1] = field[1][index];
        else
            v_nb[1] = (field[1][index] + field[1][index_nb[3]])/2.0;

        if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
                comp,&intfc_state,&hs,crx_coords,m_t_new) &&
                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
            w_nb[0] = field[2][index];
        else
            w_nb[0] = (field[2][index] + field[2][index_nb[4]])/2.0;

        if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
                comp,&intfc_state,&hs,crx_coords,m_t_new) &&
                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
            w_nb[1] = field[2][index];
        else
            w_nb[1] = (field[2][index] + field[2][index_nb[5]])/2.0;

        div = (u_nb[1]-u_nb[0])/top_h[0] + (v_nb[1]-v_nb[0])/top_h[1] + (w_nb[1]-w_nb[0])/top_h[2];
        return div;
} /* end computeFieldPointDiv_Neumann_vd */


double Incompress_Solver_Smooth_3D_Cartesian::computeFieldPointDiv_MAC_vd(
        int *icoords,
        double **field)
{
    int bNoBoundary[6];//bNoBoundary change type from bool to int
        int index;
        COMPONENT comp;
        int i,j,k,nb;
	int index_nb[6];
        double div;
        double coords[MAXD],crx_coords[MAXD];
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

	i = icoords[0];
	j = icoords[1];
	k = icoords[2];
	index = d_index3d(i,j,k,top_gmax);
        comp = top_comp[index];

	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	index_nb[4] = d_index3d(i,j,k-1,top_gmax);

    /*
     * removal tag: HAOZ
     * for PERIODIC_BOUNDARY,
     * // if a computation cell is at Boundary, use its corresponding one due to the nature of PERIODIC
     *
     * for REFLECTION_BOUNDARY,
     * // if a computation cell is at Boundary, simply use this cell since the buffer size and REFLECTION nature
     *
     * for FT_StateStructAtGridCrossing_tmp,
     * this function can only examine DIRICHLET and NEUMANN B.C.
     *
     * for rect_boundary_type(interface, dir, side)
     * this macro can only examine PERIODIC and REFLECTION B.C.
     *
     * and there is no crossover
     *
     * without changing the code, when using REFLECTION_BOUNDARY,
     * the current function calls FT_StateStructAtGridCrossing_tmp
     * which means function choose to pick neighbor cells or
     * For Reflection B.C., there is no Boundary Cell.
     *
    */
        // 4 directions
        // This should be removed. removal tag: HAOZ
        for (nb = 0; nb < 6; nb++)
        {
            checkBoundaryCondition(dir[nb],icoords,&bNoBoundary[nb],m_t_new,comp);
        }
        // TODO && FIXME: remove later
        /*
        // LOWER
        if (FT_StateStructAtGridCrossing_tmp(front,icoords,LOWER,
                comp,&intfc_state,&hs,crx_coords,m_t_new) &&
                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
            bNoBoundary[4] = NO;
        else
            bNoBoundary[4] = YES;
        // UPPER
        if (FT_StateStructAtGridCrossing_tmp(front,icoords,UPPER,
                comp,&intfc_state,&hs,crx_coords,m_t_new) &&
                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
            bNoBoundary[5] = NO;
        else
            bNoBoundary[5] = YES;
        */
        div = 0;
        if (bNoBoundary[4]==2) //cells on LOWER bdry NEUMANN
        {
            div += (field[0][index] - field[0][index_nb[0]])/top_h[0];
            div += (field[1][index] - field[1][index_nb[2]])/top_h[1];
            div += (field[2][index] - 0)/top_h[2];
            /*
             * removal tag: HAOZ
             * when using REFLECTION, this block should never be called
             * print function as debugging line
             * PERIODIC_BOUNDARY debugging print will carry out here.
            */
        }
        else //other cells
        {
            /*
             * removal tag: HAOZ
             *
             *
             * */
            div += (field[0][index] - field[0][index_nb[0]])/top_h[0];
            div += (field[1][index] - field[1][index_nb[2]])/top_h[1];
            div += (field[2][index] - field[2][index_nb[4]])/top_h[2];
        }
        return div;
} /* end computeFieldPointDiv_MAC_vd */


void Incompress_Solver_Smooth_3D_Cartesian::computeFieldPointGradPhi(
        int *icoords,
        double *field,
        double *grad_field)
{
        int index;
        COMPONENT comp;
        int i,j,k,nb;
	int index_nb[6];
        double p_nbedge[6],p0;  //the p values on the cell edges and cell center
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

	i = icoords[0];
	j = icoords[1];
	k = icoords[2];

	index = d_index3d(i,j,k,top_gmax);
        comp = top_comp[index];
	p0 = field[index];

	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	for (nb = 0; nb < 6; nb++)
	{
	    if(FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
			comp,&intfc_state,&hs,crx_coords) &&
		        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
		    boundary_state(hs) == NULL)
		    p_nbedge[nb] = 0.0;
		else //Neumann B.C.
		    p_nbedge[nb] = p0;
	    }
	    else
		p_nbedge[nb] = 0.5*(p0 + field[index_nb[nb]]);
	}
	grad_field[0] = (p_nbedge[1] - p_nbedge[0])/top_h[0];
	grad_field[1] = (p_nbedge[3] - p_nbedge[2])/top_h[1];
	grad_field[2] = (p_nbedge[5] - p_nbedge[4])/top_h[2];
}


void Incompress_Solver_Smooth_3D_Cartesian::computeFieldPointGradRho(
        int *icoords,
        double *field,
        double *grad_field)
{
        int index;
        COMPONENT comp;
        int i,j,k,nb;
        int index_nb[6];
        double p_nbedge[6],p0;  //the p values on the cell faces and cell center
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

        i = icoords[0];
        j = icoords[1];
        k = icoords[2];

        index = d_index3d(i,j,k,top_gmax);
        comp = top_comp[index];
        p0 = field[index];

        index_nb[0] = d_index3d(i-1,j,k,top_gmax);
        index_nb[1] = d_index3d(i+1,j,k,top_gmax);
        index_nb[2] = d_index3d(i,j-1,k,top_gmax);
        index_nb[3] = d_index3d(i,j+1,k,top_gmax);
        index_nb[4] = d_index3d(i,j,k-1,top_gmax);
        index_nb[5] = d_index3d(i,j,k+1,top_gmax);

        for (nb = 0; nb < 6; nb++)
        {
            if(FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                        comp,&intfc_state,&hs,crx_coords) &&
                        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
            {
                if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    boundary_state(hs) == NULL)
                {
                    if (nb == 4)
                        p_nbedge[nb] = m_rho[1];
                    if (nb == 5)
                        p_nbedge[nb] = m_rho[0];
                }
                else // homogenous Neumann B.C. for rho
                    p_nbedge[nb] = p0;
            }
            else
                p_nbedge[nb] = 0.5*(p0 + field[index_nb[nb]]);
        }
        grad_field[0] = (p_nbedge[1] - p_nbedge[0])/top_h[0];
        grad_field[1] = (p_nbedge[3] - p_nbedge[2])/top_h[1];
        grad_field[2] = (p_nbedge[5] - p_nbedge[4])/top_h[2];
}


void Incompress_Solver_Smooth_3D_Cartesian::computeFieldPointGradRho_MAC_vd(
        int *icoords,
        double *field,
        double *grad_field)
{
        int index,index_nb[6];
        COMPONENT comp;
        int i,j,k,nb;
        double p_nb[6],p0; //the p values on the cell faces and cell center
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

        i = icoords[0];
        j = icoords[1];
        k = icoords[2];
        index = d_index3d(i,j,k,top_gmax);
        comp = top_comp[index];
        p0 = field[index];

        index_nb[0] = d_index3d(i-1,j,k,top_gmax);
        index_nb[1] = d_index3d(i+1,j,k,top_gmax);
        index_nb[2] = d_index3d(i,j-1,k,top_gmax);
        index_nb[3] = d_index3d(i,j+1,k,top_gmax);
        index_nb[4] = d_index3d(i,j,k-1,top_gmax);
        index_nb[5] = d_index3d(i,j,k+1,top_gmax);

        for (nb = 0; nb < 6; ++nb)
        {
            if(FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                        comp,&intfc_state,&hs,crx_coords) &&
                        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
            {
                if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    boundary_state(hs) == NULL)
                {
                    if (nb == 4)
                        p_nb[nb] = m_rho[1];
                    if (nb == 5)
                        p_nb[nb] = m_rho[0];
                }
                else // homogenous Neumann B.C. for rho
                    p_nb[nb] = p0;
            }
            else
                p_nb[nb] = field[index_nb[nb]];
        }
        grad_field[0] = (p_nb[1] - p0)/top_h[0];
        grad_field[1] = (p_nb[3] - p0)/top_h[1];
        grad_field[2] = (p_nb[5] - p0)/top_h[2];
}


void Incompress_Solver_Smooth_3D_Cartesian::computeFieldPointGrad_MAC_vd(
        int *icoords,
        double *field,
        double *grad_field)
{
        int index;
        COMPONENT comp;
        int i,j,k,nb;
        int index_nb[6];
        double p_nb[6],p0;
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

        i = icoords[0];
        j = icoords[1];
        k = icoords[2];
        index = d_index3d(i,j,k,top_gmax);
        comp = top_comp[index];
        p0 = field[index];

        index_nb[0] = d_index3d(i-1,j,k,top_gmax);
        index_nb[1] = d_index3d(i+1,j,k,top_gmax);
        index_nb[2] = d_index3d(i,j-1,k,top_gmax);
        index_nb[3] = d_index3d(i,j+1,k,top_gmax);
        index_nb[4] = d_index3d(i,j,k-1,top_gmax);
        index_nb[5] = d_index3d(i,j,k+1,top_gmax);

        /*
         * removal tag: HAOZ
         * REFLECTION B.C.
         * */
        for (nb = 0; nb < 6; nb++)
        {
            if(FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                        comp,&intfc_state,&hs,crx_coords) &&
                        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
            {
                if ((wave_type(hs) == DIRICHLET_BOUNDARY &&
                     boundary_state(hs) == NULL)
                     || wave_type(hs) == NEUMANN_BOUNDARY)
                    p_nb[nb] = p0;
                else
                    p_nb[nb] = field[index_nb[nb]];
            }
            else
                p_nb[nb] = field[index_nb[nb]];
        }
        grad_field[0] = (p_nb[1] - p0)/top_h[0];
        grad_field[1] = (p_nb[3] - p0)/top_h[1];
        grad_field[2] = (p_nb[5] - p0)/top_h[2];
} /* end computeFieldPointGrad_MAC_vd */


void Incompress_Solver_Smooth_3D_Cartesian::computeFieldPointGrad(
        int *icoords,
        double *field,
        double *grad_field)
{
        int index;
        COMPONENT comp;
        int i,j,k,nb;
	int index_nb[6];
        double p_nbedge[6],p0;  //the p values on the cell edges and cell center
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

	i = icoords[0];
	j = icoords[1];
	k = icoords[2];

	index = d_index3d(i,j,k,top_gmax);
        comp = top_comp[index];
	p0 = field[index];

	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

	for (nb = 0; nb < 6; nb++)
	{
	    if(FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
			comp,&intfc_state,&hs,crx_coords) &&
		        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
		    boundary_state(hs) == NULL)
		    p_nbedge[nb] = 0.0;
		else
		{
		    //p_nbedge[nb] = p0;

		    if (nb == 0)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[1]] - p0);
		    else if (nb == 1)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[0]]);
		    else if (nb == 2)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[3]] - p0);
		    else if (nb == 3)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[2]]);
		    else if (nb == 4)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[5]] - p0);
		    else if (nb == 5)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[4]]);
		}
	    }
	    else
		p_nbedge[nb] = (p0 + field[index_nb[nb]])/2.0;
	}
	grad_field[0] = (p_nbedge[1] - p_nbedge[0])/top_h[0];
	grad_field[1] = (p_nbedge[3] - p_nbedge[2])/top_h[1];
	grad_field[2] = (p_nbedge[5] - p_nbedge[4])/top_h[2];
}


void Incompress_Solver_Smooth_3D_Cartesian::compDiffWithSmoothProperty_2nd_decoupled(void)
{

    printf("\nEntering the Incompress_debug and the m_dt for diffusion step is %.16g\n",m_dt);

    COMPONENT comp;
    int index,index_nb[6],size;
    int I,I_nb[6];
    int i,j,k,l,nb,icoords[MAXD];
    L_STATE source_term;
    INTERFACE *intfc = front->interf;
    double coords[MAXD], crx_coords[MAXD];
    double coeff[6],mu[6],mu0,rho,corner[6],rhs;

    // U_nb contains states at neighbor cell or states on the boundary.
    double U_nb[6], U_nb_new[6], U_center;

    double speed;
    double *x;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    double (*getStateVel[3])(POINTER) =
    {getStateXvel,getStateYvel,getStateZvel};
    POINTER intfc_state;
    HYPER_SURF *hs;

    setIndexMap();


    max_speed = 0.0;

    size = iupper - ilower;
    FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

    for (l = 0; l < dim; ++l)
    {
	PETSc solver;

	solver.Create(ilower, iupper-1, 7, 7);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
		for (i = imin; i <= imax; i++)
		{
		    I = ijk_to_I[i][j][k];
		    if (I == -1) continue;

		    index = d_index3d(i,j,k,top_gmax);
		    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
		    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
		    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
		    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
		    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
		    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

		    icoords[0] = i;
		    icoords[1] = j;
		    icoords[2] = k;
		    comp = top_comp[index];

		    I_nb[0] = ijk_to_I[i-1][j][k]; //west
		    I_nb[1] = ijk_to_I[i+1][j][k]; //east
		    I_nb[2] = ijk_to_I[i][j-1][k]; //south
		    I_nb[3] = ijk_to_I[i][j+1][k]; //north
		    I_nb[4] = ijk_to_I[i][j][k-1]; //lower
		    I_nb[5] = ijk_to_I[i][j][k+1]; //upper


		    mu0 = cell_center[index].m_state.m_mu;
		    rho = cell_center[index].m_state.m_rho;
		    U_center =  cell_center[index].m_state.m_U[l];

		    for (nb = 0; nb < 6; nb++)
		    {
			if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
				comp,&intfc_state,&hs,crx_coords,m_t_old) &&
				wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
			{
			    // old boundary condition
			    U_nb[nb] = getStateVel[l](intfc_state);
			    // new boundary condition
			    FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
				    comp,&intfc_state,&hs,crx_coords,m_t_new);
			    U_nb_new[nb] = getStateVel[l](intfc_state);

			    if (wave_type(hs) == DIRICHLET_BOUNDARY ||
				    wave_type(hs) == NEUMANN_BOUNDARY)
				mu[nb] = mu0;
			    else
				mu[nb] = 1.0/2*(mu0 +
					cell_center[index_nb[nb]].m_state.m_mu);
			}
			else
			{
			    U_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[l];
			    mu[nb] = 1.0/2*(mu0 +
				    cell_center[index_nb[nb]].m_state.m_mu);
			}
		    }


		    getRectangleCenter(index, coords);
		    computeSourceTerm(coords, source_term);


		    rhs = 0;
		    double dh[3] = {top_h[0], top_h[1], top_h[2]};
		    for(nb = 0; nb<6; nb++)
		    {
			// use dh[1]*dh[2] as the face area
			if(nb<2)
			{
			    dh[0] = top_h[0];
			    dh[1] = top_h[1];
			    dh[2] = top_h[2];
			}
			else if(nb<4 && nb>=2)
			{
			    dh[0] = top_h[1];
			    dh[1] = top_h[2];
			    dh[2] = top_h[0];
			}
			else if(nb<6 && nb>=4)
			{
			    dh[0] = top_h[2];
			    dh[1] = top_h[0];
			    dh[2] = top_h[1];
			}

			if(I_nb[nb]>=0)  // interior
			{
			    // u^{*}
			    solver.Add_A(I, I,        (-1) * -0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]));
			    solver.Add_A(I, I_nb[nb], (-1) *  0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]));
			    // u^{n}
			    rhs += 0.5*m_dt/rho*mu[nb] *
				    (U_nb[nb]-U_center) /(dh[0]*dh[0]);
			}
			else		// boundary
			{
			    // u^{*}
			    solver.Add_A(I, I,       (-1) * -0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]) * 2);
			    rhs += 0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]) * 2 * U_nb_new[nb];
			    // u^{n}
			    rhs += 0.5*m_dt/rho*mu[nb] *
				    (U_nb[nb]-U_center) /(dh[0]*dh[0]) * 2;
			}
		    }

		    rhs += m_dt*source_term.m_U[l];
		    rhs += m_dt*cell_center[index].m_state.f_surf[l];
		    rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho;
		    //printf("gradq = %.16g\n",cell_center[index].m_state.grad_q[l]);

		    rhs -= m_dt * cell_center[index].m_state.m_adv[l]; //advection source term

		    solver.Add_A(I, I, 1.0);
		    rhs += U_center;

		    solver.Add_b(I, rhs);
		}

	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);
	solver.Solve_GMRES();

	// get back the solution
	solver.Get_x(x);

	PetscInt num_iter;
	double rel_residual;
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);

	if (debugging("PETSc"))
	    (void) printf("L_CARTESIAN::"
		    "compDiffWithSmoothProperty_2nd_decoupled: "
		    "num_iter = %d, rel_residual = %le. \n",
		    num_iter,rel_residual);

	for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
		for (i = imin; i <= imax; i++)
		{
		    I = ijk_to_I[i][j][k];
		    index = d_index3d(i,j,k,top_gmax);
		    if (I >= 0)
		    {
			cell_center[index].m_state.m_U[l] = x[I-ilower];
		    }
		    else
		    {
			cell_center[index].m_state.m_U[l] = 0.0;
		    }
		}

	for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
		for (i = imin; i <= imax; i++)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    array[index] = cell_center[index].m_state.m_U[l];
		}
	scatMeshArray();
	for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
		for (i = 0; i <= top_gmax[0]; i++)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    cell_center[index].m_state.m_U[l] = array[index];
		}
    }
    for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index3d(i,j,k,top_gmax);
		speed = fabs(cell_center[index].m_state.m_U[0]) +
			fabs(cell_center[index].m_state.m_U[1]) +
			fabs(cell_center[index].m_state.m_U[2]);
		if (speed > max_speed)
		    max_speed = speed;
	    }
    pp_global_max(&max_speed,1);

    FT_FreeThese(1,x);
}


void Incompress_Solver_Smooth_3D_Cartesian::compAdvectionTerm_decoupled(void)
{
    int index;
    int I;
    int i,j,k,icoords[MAXD];

    setIndexMap();


	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    I = ijk_to_I[i][j][k];
	    if (I == -1) continue;

	    index = d_index3d(i,j,k,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    double convectionTerm[3];
	    getAdvectionTerm_decoupled(icoords,convectionTerm);
	    cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	    cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	    cell_center[index].m_state.m_adv[2] = convectionTerm[2];
	}
}


void Incompress_Solver_Smooth_3D_Cartesian::compAdvectionTerm_coupled(void)
{
    int index;
    int I;
    int i,j,k,icoords[MAXD];

    setIndexMap();


	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    I = ijk_to_I[i][j][k];
	    if (I == -1) continue;

	    index = d_index3d(i,j,k,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    double convectionTerm[3];
	    getAdvectionTerm_coupled(icoords,convectionTerm);
	    cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	    cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	    cell_center[index].m_state.m_adv[2] = convectionTerm[2];
	}
}


void Incompress_Solver_Smooth_3D_Cartesian::compAdvectionTerm_decoupled_upgraded(void)
{
    int index;
    int I;
    int i,j,k,icoords[MAXD];

    setIndexMap();


	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    I = ijk_to_I[i][j][k];
	    if (I == -1) continue;

	    index = d_index3d(i,j,k,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    double convectionTerm[3];
	    getAdvectionTerm_decoupled_upgraded(icoords,convectionTerm);
	    cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	    cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	    cell_center[index].m_state.m_adv[2] = convectionTerm[2];
	}
}


void Incompress_Solver_Smooth_3D_Cartesian::compAdvectionTerm_coupled_upgraded(void)
{
    int index;
    int I;
    int i,j,k,icoords[MAXD];

    setIndexMap();


	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    I = ijk_to_I[i][j][k];
	    if (I == -1) continue;

	    index = d_index3d(i,j,k,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

	    double convectionTerm[3];
	    getAdvectionTerm_coupled_upgraded(icoords,convectionTerm);
	    cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	    cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	    cell_center[index].m_state.m_adv[2] = convectionTerm[2];
	}
}


void Incompress_Solver_Smooth_3D_Cartesian::getAdvectionTerm_decoupled(
	int *icoords,
	double convectionTerm[3])
{
    bool bNoBoundary;
    int ICoords[3];
    L_STATE sl, sr, state_west, state_east, state_south, state_north, state_lower, state_upper;

    double dx = top_h[0];
    double dy = top_h[1];
    double dz = top_h[2];

    convectionTerm[0] = 0;
    convectionTerm[1] = 0;
    convectionTerm[2] = 0;

    // WEST
    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep(ICoords,EAST,sl);
    }
    getFaceVelocity_middleStep(icoords,WEST,sr);
    getRiemannSolution(COORD_X,sl,sr,state_west);

    // EAST
    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep(ICoords,WEST,sr);
    }
    getFaceVelocity_middleStep(icoords,EAST,sl);
    getRiemannSolution(COORD_X,sl,sr,state_east);

    // SOUTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep(ICoords,NORTH,sl);
    }
    getFaceVelocity_middleStep(icoords,SOUTH,sr);
    getRiemannSolution(COORD_Y,sl,sr,state_south);

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep(ICoords,SOUTH,sr);
    }
    getFaceVelocity_middleStep(icoords,NORTH,sl);
    getRiemannSolution(COORD_Y,sl,sr,state_north);

    // LOWER
    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep(ICoords,UPPER,sl);
    }
    getFaceVelocity_middleStep(icoords,LOWER,sr);
    getRiemannSolution(COORD_Z,sl,sr,state_lower);

    // UPPER
    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep(ICoords,LOWER,sr);
    }
    getFaceVelocity_middleStep(icoords,UPPER,sl);
    getRiemannSolution(COORD_Z,sl,sr,state_upper);

    convectionTerm[0] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[0]-state_west.m_U[0])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[0]-state_south.m_U[0])/dy +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]) * (state_upper.m_U[0]-state_lower.m_U[0])/dz;
    convectionTerm[1] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[1]-state_west.m_U[1])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[1]-state_south.m_U[1])/dy +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]) * (state_upper.m_U[1]-state_lower.m_U[1])/dz;
    convectionTerm[2] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[2]-state_west.m_U[2])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[2]-state_south.m_U[2])/dy +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]) * (state_upper.m_U[2]-state_lower.m_U[2])/dz;
}


void Incompress_Solver_Smooth_3D_Cartesian::getAdvectionTerm_coupled(
	int *icoords,
	double convectionTerm[3])
{
    bool bNoBoundary;
    int ICoords[3];
    L_STATE sl, sr, state_west, state_east, state_south, state_north, state_lower, state_upper;

    double dx = top_h[0];
    double dy = top_h[1];
    double dz = top_h[2];

    convectionTerm[0] = 0;
    convectionTerm[1] = 0;
    convectionTerm[2] = 0;

    // WEST
    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled(ICoords,EAST,sl);
    }
    getFaceVelocity_middleStep_coupled(icoords,WEST,sr);
    getRiemannSolution(COORD_X,sl,sr,state_west);

    // EAST
    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled(ICoords,WEST,sr);
    }
    getFaceVelocity_middleStep_coupled(icoords,EAST,sl);
    getRiemannSolution(COORD_X,sl,sr,state_east);

    // SOUTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled(ICoords,NORTH,sl);
    }
    getFaceVelocity_middleStep_coupled(icoords,SOUTH,sr);
    getRiemannSolution(COORD_Y,sl,sr,state_south);

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled(ICoords,SOUTH,sr);
    }
    getFaceVelocity_middleStep_coupled(icoords,NORTH,sl);
    getRiemannSolution(COORD_Y,sl,sr,state_north);

    // LOWER
    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep_coupled(ICoords,UPPER,sl);
    }
    getFaceVelocity_middleStep_coupled(icoords,LOWER,sr);
    getRiemannSolution(COORD_Z,sl,sr,state_lower);

    // UPPER
    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,sr,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep_coupled(ICoords,LOWER,sr);
    }
    getFaceVelocity_middleStep_coupled(icoords,UPPER,sl);
    getRiemannSolution(COORD_Z,sl,sr,state_upper);

    convectionTerm[0] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[0]-state_west.m_U[0])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[0]-state_south.m_U[0])/dy +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]) * (state_upper.m_U[0]-state_lower.m_U[0])/dz;
    convectionTerm[1] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[1]-state_west.m_U[1])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[1]-state_south.m_U[1])/dy +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]) * (state_upper.m_U[1]-state_lower.m_U[1])/dz;
    convectionTerm[2] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[2]-state_west.m_U[2])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[2]-state_south.m_U[2])/dy +
	    1.0/2*(state_upper.m_U[2]+state_lower.m_U[2]) * (state_upper.m_U[2]-state_lower.m_U[2])/dz;
}


void Incompress_Solver_Smooth_3D_Cartesian::getAdvectionTerm_decoupled_upgraded(
	int *icoords,
	double convectionTerm[3])
{
    bool bNoBoundary;
    int ICoords[3];
    L_STATE sl, sr, state_west_bar, state_east_bar, state_south_bar, state_north_bar, state_lower_bar, state_upper_bar;

    L_STATE state_west_hat, state_east_hat, state_south_hat, state_north_hat, state_lower_hat, state_upper_hat;

    L_STATE state_west_hat_l, state_west_hat_r;
    L_STATE state_east_hat_l, state_east_hat_r;
    L_STATE state_south_hat_l, state_south_hat_r;
    L_STATE state_north_hat_l, state_north_hat_r;
    L_STATE state_lower_hat_l, state_lower_hat_r;
    L_STATE state_upper_hat_l, state_upper_hat_r;

    double transverseD[3];

    double dx = top_h[0];
    double dy = top_h[1];
    double dz = top_h[2];

    convectionTerm[0] = 0;
    convectionTerm[1] = 0;
    convectionTerm[2] = 0;


    ///////////////////////////// Get the state_hat on six faces first //////////////////////

    // WEST
    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,state_west_hat_l,m_t_int);
    if(!bNoBoundary)
    {
	state_west_hat = state_west_hat_l;
	//getFaceVelocity_middleStep_hat(icoords,WEST,state_west_hat_r);
	//getRiemannSolution(COORD_X,state_west_hat_l,state_west_hat_r,state_west_hat);
    }
    else
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,EAST,state_west_hat_l);
	getFaceVelocity_middleStep_hat(icoords,WEST,state_west_hat_r);
	getRiemannSolution(COORD_X,state_west_hat_l,state_west_hat_r,state_west_hat);
    }

    // EAST
    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,state_east_hat_r,m_t_int);
    if(!bNoBoundary)
    {
	state_east_hat = state_east_hat_r;
	//getFaceVelocity_middleStep_hat(icoords,EAST,state_east_hat_l);
	//getRiemannSolution(COORD_X,state_east_hat_l,state_east_hat_r,state_east_hat);
    }
    else
    {
	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,WEST,state_east_hat_r);
	getFaceVelocity_middleStep_hat(icoords,EAST,state_east_hat_l);
	getRiemannSolution(COORD_X,state_east_hat_l,state_east_hat_r,state_east_hat);
    }

    // SOUTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,state_south_hat_l,m_t_int);
    if(!bNoBoundary)
    {
	state_south_hat = state_south_hat_l;
	//getFaceVelocity_middleStep_hat(icoords,SOUTH,state_south_hat_r);
	//getRiemannSolution(COORD_Y,state_south_hat_l,state_south_hat_r,state_south_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,NORTH,state_south_hat_l);
	getFaceVelocity_middleStep_hat(icoords,SOUTH,state_south_hat_r);
	getRiemannSolution(COORD_Y,state_south_hat_l,state_south_hat_r,state_south_hat);
    }

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,state_north_hat_r,m_t_int);
    if(!bNoBoundary)
    {
	state_north_hat = state_north_hat_r;
	//getFaceVelocity_middleStep_hat(icoords,NORTH,state_north_hat_l);
	//getRiemannSolution(COORD_Y,state_north_hat_l,state_north_hat_r,state_north_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,SOUTH,state_north_hat_r);
	getFaceVelocity_middleStep_hat(icoords,NORTH,state_north_hat_l);
	getRiemannSolution(COORD_Y,state_north_hat_l,state_north_hat_r,state_north_hat);

    }

    // LOWER
    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,state_lower_hat_l,m_t_int);
    if(!bNoBoundary)
    {
	state_lower_hat = state_lower_hat_l;
	//getFaceVelocity_middleStep_hat(icoords,SOUTH,state_south_hat_r);
	//getRiemannSolution(COORD_Y,state_south_hat_l,state_south_hat_r,state_south_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep_hat(ICoords,UPPER,state_lower_hat_l);
	getFaceVelocity_middleStep_hat(icoords,LOWER,state_lower_hat_r);
	getRiemannSolution(COORD_Z,state_lower_hat_l,state_lower_hat_r,state_lower_hat);
    }

    // UPPER
    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,state_upper_hat_r,m_t_int);
    if(!bNoBoundary)
    {
	state_upper_hat = state_upper_hat_r;
	//getFaceVelocity_middleStep_hat(icoords,NORTH,state_north_hat_l);
	//getRiemannSolution(COORD_Y,state_north_hat_l,state_north_hat_r,state_north_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep_hat(ICoords,LOWER,state_upper_hat_r);
	getFaceVelocity_middleStep_hat(icoords,UPPER,state_upper_hat_l);
	getRiemannSolution(COORD_Z,state_upper_hat_l,state_upper_hat_r,state_upper_hat);

    }

    ///////////////////////////// get the state_bar on four edges ////////////////////////////////

    // WEST

    transverseD[0] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);


    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,sl,m_t_int);
    if(!bNoBoundary)
    {
	state_west_bar = sl;
	//getFaceVelocity_middleStep_bar(icoords,WEST,sr,transverseD,state_west_hat_r);
	//getRiemannSolution(COORD_X,sl,sr,state_west_bar);
    }
    else
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_bar(ICoords,EAST,sl,transverseD,state_west_hat_l);
	getFaceVelocity_middleStep_bar(icoords,WEST,sr,transverseD,state_west_hat_r);
	getRiemannSolution(COORD_X,sl,sr,state_west_bar);
    }

    // EAST

    transverseD[0] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);


    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,sr,m_t_int);
    if(!bNoBoundary)
    {
	state_east_bar = sr;
	//getFaceVelocity_middleStep_bar(icoords,EAST,sl,transverseD,state_east_hat_l);
	//getRiemannSolution(COORD_X,sl,sr,state_east_bar);

    }
    else
    {

	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_bar(ICoords,WEST,sr,transverseD,state_east_hat_r);
	getFaceVelocity_middleStep_bar(icoords,EAST,sl,transverseD,state_east_hat_l);
	getRiemannSolution(COORD_X,sl,sr,state_east_bar);
    }

    // SOUTH

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);



    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,sl,m_t_int);
    if(!bNoBoundary)
    {
	state_south_bar = sl;
	//getFaceVelocity_middleStep_bar(icoords,SOUTH,sr,transverseD,state_south_hat_r);
	//getRiemannSolution(COORD_Y,sl,sr,state_south_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_bar(ICoords,NORTH,sl,transverseD,state_south_hat_l);
	getFaceVelocity_middleStep_bar(icoords,SOUTH,sr,transverseD,state_south_hat_r);
	getRiemannSolution(COORD_Y,sl,sr,state_south_bar);
    }

    // NORTH

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);



    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,sr,m_t_int);
    if(!bNoBoundary)
    {
	state_north_bar = sr;
	//getFaceVelocity_middleStep_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	//getRiemannSolution(COORD_Y,sl,sr,state_north_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_bar(ICoords,SOUTH,sr,transverseD,state_north_hat_r);
	getFaceVelocity_middleStep_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	getRiemannSolution(COORD_Y,sl,sr,state_north_bar);

    }

    // LOWER

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);


    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,sl,m_t_int);
    if(!bNoBoundary)
    {
	state_lower_bar = sl;
	//getFaceVelocity_middleStep_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	//getRiemannSolution(COORD_Y,sl,sr,state_north_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep_bar(ICoords,UPPER,sl,transverseD,state_lower_hat_l);
	getFaceVelocity_middleStep_bar(icoords,LOWER,sr,transverseD,state_lower_hat_r);
	getRiemannSolution(COORD_Z,sl,sr,state_lower_bar);

    }

    // UPPER

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);


    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,sr,m_t_int);
    if(!bNoBoundary)
    {
	state_upper_bar = sr;
	//getFaceVelocity_middleStep_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	//getRiemannSolution(COORD_Y,sl,sr,state_north_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep_bar(ICoords,LOWER,sr,transverseD,state_upper_hat_r);
	getFaceVelocity_middleStep_bar(icoords,UPPER,sl,transverseD,state_upper_hat_l);
	getRiemannSolution(COORD_Z,sl,sr,state_upper_bar);

    }

    convectionTerm[0] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[0]-state_west_bar.m_U[0])/dx +
	    1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[0]-state_south_bar.m_U[0])/dy +
	    1.0/2*(state_upper_bar.m_U[2]+state_lower_bar.m_U[2]) * (state_upper_bar.m_U[0]-state_lower_bar.m_U[0])/dz;
    convectionTerm[1] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[1]-state_west_bar.m_U[1])/dx +
	    1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[1]-state_south_bar.m_U[1])/dy +
	    1.0/2*(state_upper_bar.m_U[2]+state_lower_bar.m_U[2]) * (state_upper_bar.m_U[1]-state_lower_bar.m_U[1])/dz;
    convectionTerm[2] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[2]-state_west_bar.m_U[2])/dx +
	    1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[2]-state_south_bar.m_U[2])/dy +
	    1.0/2*(state_upper_bar.m_U[2]+state_lower_bar.m_U[2]) * (state_upper_bar.m_U[2]-state_lower_bar.m_U[2])/dz;
}


void Incompress_Solver_Smooth_3D_Cartesian::getAdvectionTerm_coupled_upgraded(
	int *icoords,
	double convectionTerm[3])
{
    bool bNoBoundary;
    int ICoords[3];
    L_STATE sl, sr, state_west_bar, state_east_bar, state_south_bar, state_north_bar, state_lower_bar, state_upper_bar;

    L_STATE state_west_hat, state_east_hat, state_south_hat, state_north_hat, state_lower_hat, state_upper_hat;

    L_STATE state_west_hat_l, state_west_hat_r;
    L_STATE state_east_hat_l, state_east_hat_r;
    L_STATE state_south_hat_l, state_south_hat_r;
    L_STATE state_north_hat_l, state_north_hat_r;
    L_STATE state_lower_hat_l, state_lower_hat_r;
    L_STATE state_upper_hat_l, state_upper_hat_r;

    double transverseD[3];

    double dx = top_h[0];
    double dy = top_h[1];
    double dz = top_h[2];

    convectionTerm[0] = 0;
    convectionTerm[1] = 0;
    convectionTerm[2] = 0;


    ///////////////////////////// Get the state_hat on six faces first //////////////////////

    // WEST
    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,state_west_hat_l,m_t_int);
    if(!bNoBoundary)
    {
	state_west_hat = state_west_hat_l;
	//getFaceVelocity_middleStep_hat(icoords,WEST,state_west_hat_r);
	//getRiemannSolution(COORD_X,state_west_hat_l,state_west_hat_r,state_west_hat);
    }
    else
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,EAST,state_west_hat_l);
	getFaceVelocity_middleStep_hat(icoords,WEST,state_west_hat_r);
	getRiemannSolution(COORD_X,state_west_hat_l,state_west_hat_r,state_west_hat);
    }

    // EAST
    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,state_east_hat_r,m_t_int);
    if(!bNoBoundary)
    {
	state_east_hat = state_east_hat_r;
	//getFaceVelocity_middleStep_hat(icoords,EAST,state_east_hat_l);
	//getRiemannSolution(COORD_X,state_east_hat_l,state_east_hat_r,state_east_hat);
    }
    else
    {
	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,WEST,state_east_hat_r);
	getFaceVelocity_middleStep_hat(icoords,EAST,state_east_hat_l);
	getRiemannSolution(COORD_X,state_east_hat_l,state_east_hat_r,state_east_hat);
    }

    // SOUTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,state_south_hat_l,m_t_int);
    if(!bNoBoundary)
    {
	state_south_hat = state_south_hat_l;
	//getFaceVelocity_middleStep_hat(icoords,SOUTH,state_south_hat_r);
	//getRiemannSolution(COORD_Y,state_south_hat_l,state_south_hat_r,state_south_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,NORTH,state_south_hat_l);
	getFaceVelocity_middleStep_hat(icoords,SOUTH,state_south_hat_r);
	getRiemannSolution(COORD_Y,state_south_hat_l,state_south_hat_r,state_south_hat);
    }

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,state_north_hat_r,m_t_int);
    if(!bNoBoundary)
    {
	state_north_hat = state_north_hat_r;
	//getFaceVelocity_middleStep_hat(icoords,NORTH,state_north_hat_l);
	//getRiemannSolution(COORD_Y,state_north_hat_l,state_north_hat_r,state_north_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_hat(ICoords,SOUTH,state_north_hat_r);
	getFaceVelocity_middleStep_hat(icoords,NORTH,state_north_hat_l);
	getRiemannSolution(COORD_Y,state_north_hat_l,state_north_hat_r,state_north_hat);

    }

    // LOWER
    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,state_lower_hat_l,m_t_int);
    if(!bNoBoundary)
    {
	state_lower_hat = state_lower_hat_l;
	//getFaceVelocity_middleStep_hat(icoords,SOUTH,state_south_hat_r);
	//getRiemannSolution(COORD_Y,state_south_hat_l,state_south_hat_r,state_south_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep_hat(ICoords,UPPER,state_lower_hat_l);
	getFaceVelocity_middleStep_hat(icoords,LOWER,state_lower_hat_r);
	getRiemannSolution(COORD_Z,state_lower_hat_l,state_lower_hat_r,state_lower_hat);
    }

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,state_upper_hat_r,m_t_int);
    if(!bNoBoundary)
    {
	state_upper_hat = state_upper_hat_r;
	//getFaceVelocity_middleStep_hat(icoords,NORTH,state_north_hat_l);
	//getRiemannSolution(COORD_Y,state_north_hat_l,state_north_hat_r,state_north_hat);
    }
    else
    {
	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep_hat(ICoords,LOWER,state_upper_hat_r);
	getFaceVelocity_middleStep_hat(icoords,UPPER,state_upper_hat_l);
	getRiemannSolution(COORD_Z,state_upper_hat_l,state_upper_hat_r,state_upper_hat);

    }

    ///////////////////////////// get the state_bar on four edges ////////////////////////////////

    // WEST

    transverseD[0] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);


    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,sl,m_t_int);
    if(!bNoBoundary)
    {
	state_west_bar = sl;
	//getFaceVelocity_middleStep_bar(icoords,WEST,sr,transverseD,state_west_hat_r);
	//getRiemannSolution(COORD_X,sl,sr,state_west_bar);
    }
    else
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled_bar(ICoords,EAST,sl,transverseD,state_west_hat_l);
	getFaceVelocity_middleStep_coupled_bar(icoords,WEST,sr,transverseD,state_west_hat_r);
	getRiemannSolution(COORD_X,sl,sr,state_west_bar);
    }

    // EAST


    transverseD[0] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);


    bNoBoundary = getNeighborOrBoundaryState(icoords,EAST,sr,m_t_int);
    if(!bNoBoundary)
    {
	state_east_bar = sr;
	//getFaceVelocity_middleStep_bar(icoords,EAST,sl,transverseD,state_east_hat_l);
	//getRiemannSolution(COORD_X,sl,sr,state_east_bar);

    }
    else
    {

	ICoords[0] = icoords[0] + 1;
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled_bar(ICoords,WEST,sr,transverseD,state_east_hat_r);
	getFaceVelocity_middleStep_coupled_bar(icoords,EAST,sl,transverseD,state_east_hat_l);
	getRiemannSolution(COORD_X,sl,sr,state_east_bar);
    }

    // SOUTH
    //

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);



    bNoBoundary = getNeighborOrBoundaryState(icoords,SOUTH,sl,m_t_int);
    if(!bNoBoundary)
    {
	state_south_bar = sl;
	//getFaceVelocity_middleStep_bar(icoords,SOUTH,sr,transverseD,state_south_hat_r);
	//getRiemannSolution(COORD_Y,sl,sr,state_south_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] - 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled_bar(ICoords,NORTH,sl,transverseD,state_south_hat_l);
	getFaceVelocity_middleStep_coupled_bar(icoords,SOUTH,sr,transverseD,state_south_hat_r);
	getRiemannSolution(COORD_Y,sl,sr,state_south_bar);
    }

    // NORTH

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[0] - state_lower_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[1] - state_lower_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dz) * (state_lower_hat.m_U[2] + state_upper_hat.m_U[2]) * (state_upper_hat.m_U[2] - state_lower_hat.m_U[2]);



    bNoBoundary = getNeighborOrBoundaryState(icoords,NORTH,sr,m_t_int);
    if(!bNoBoundary)
    {
	state_north_bar = sr;
	//getFaceVelocity_middleStep_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	//getRiemannSolution(COORD_Y,sl,sr,state_north_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1] + 1;
	ICoords[2] = icoords[2];
	getFaceVelocity_middleStep_coupled_bar(ICoords,SOUTH,sr,transverseD,state_north_hat_r);
	getFaceVelocity_middleStep_coupled_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	getRiemannSolution(COORD_Y,sl,sr,state_north_bar);

    }

    // LOWER

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);


    bNoBoundary = getNeighborOrBoundaryState(icoords,LOWER,sl,m_t_int);
    if(!bNoBoundary)
    {
	state_lower_bar = sl;
	//getFaceVelocity_middleStep_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	//getRiemannSolution(COORD_Y,sl,sr,state_north_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] - 1;
	getFaceVelocity_middleStep_coupled_bar(ICoords,UPPER,sl,transverseD,state_lower_hat_l);
	getFaceVelocity_middleStep_coupled_bar(icoords,LOWER,sr,transverseD,state_lower_hat_r);
	getRiemannSolution(COORD_Z,sl,sr,state_lower_bar);

    }

    // UPPER

    transverseD[0] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[0] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);

    transverseD[1] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    transverseD[1] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);

    transverseD[2] =  1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[2] - state_west_hat.m_U[2]);
    transverseD[2] += 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[2] - state_south_hat.m_U[2]);


    bNoBoundary = getNeighborOrBoundaryState(icoords,UPPER,sr,m_t_int);
    if(!bNoBoundary)
    {
	state_upper_bar = sr;
	//getFaceVelocity_middleStep_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	//getRiemannSolution(COORD_Y,sl,sr,state_north_bar);
    }
    else
    {

	ICoords[0] = icoords[0];
	ICoords[1] = icoords[1];
	ICoords[2] = icoords[2] + 1;
	getFaceVelocity_middleStep_coupled_bar(ICoords,LOWER,sr,transverseD,state_upper_hat_r);
	getFaceVelocity_middleStep_coupled_bar(icoords,UPPER,sl,transverseD,state_upper_hat_l);
	getRiemannSolution(COORD_Z,sl,sr,state_upper_bar);

    }

    convectionTerm[0] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[0]-state_west_bar.m_U[0])/dx +
	    1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[0]-state_south_bar.m_U[0])/dy +
	    1.0/2*(state_upper_bar.m_U[2]+state_lower_bar.m_U[2]) * (state_upper_bar.m_U[0]-state_lower_bar.m_U[0])/dz;
    convectionTerm[1] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[1]-state_west_bar.m_U[1])/dx +
	    1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[1]-state_south_bar.m_U[1])/dy +
	    1.0/2*(state_upper_bar.m_U[2]+state_lower_bar.m_U[2]) * (state_upper_bar.m_U[1]-state_lower_bar.m_U[1])/dz;
    convectionTerm[2] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[2]-state_west_bar.m_U[2])/dx +
	    1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[2]-state_south_bar.m_U[2])/dy +
	    1.0/2*(state_upper_bar.m_U[2]+state_lower_bar.m_U[2]) * (state_upper_bar.m_U[2]-state_lower_bar.m_U[2])/dz;
}


void Incompress_Solver_Smooth_3D_Cartesian::getFaceVelocity_middleStep_hat(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_hat)
{
    L_STATE  state_orig;
    int index;
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    state_orig = cell_center[index].m_state;
    double sL;

    double dx = 0, dy = 0, dz = 0, slope_x_limited[3] = {0,0,0}, slope_y_limited[3] = {0,0,0}, slope_z_limited[3] = {0,0,0};



    switch(dir)
    {
    case WEST:
	getLimitedSlope(icoords,COORD_X,slope_x_limited);
	dx = top_h[0];
	if (state_orig.m_U[0] <= 0)
	    sL = 1.0;
	else
	    sL = 0.0;
	state_hat.m_U[0] = state_orig.m_U[0] + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[0];
	state_hat.m_U[1] = state_orig.m_U[1] + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[1];
	state_hat.m_U[2] = state_orig.m_U[2] + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[2];
	break;


    case EAST:
	getLimitedSlope(icoords,COORD_X,slope_x_limited);
	dx = top_h[0];
	if (state_orig.m_U[0] >= 0)
	    sL = 1.0;
	else
	    sL = 0.0;
	state_hat.m_U[0] = state_orig.m_U[0] + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[0];
	state_hat.m_U[1] = state_orig.m_U[1] + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[1];
	state_hat.m_U[2] = state_orig.m_U[2] + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[2];
	break;


    case SOUTH:
	getLimitedSlope(icoords,COORD_Y,slope_y_limited);
	dy = top_h[1];
	if (state_orig.m_U[1] <= 0)
	    sL = 1.0;
	else
	    sL = 0.0;
	state_hat.m_U[0] = state_orig.m_U[0] + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[0];
	state_hat.m_U[1] = state_orig.m_U[1] + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[1];
	state_hat.m_U[2] = state_orig.m_U[2] + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[2];
	break;

    case NORTH:
	getLimitedSlope(icoords,COORD_Y,slope_y_limited);
	dy = top_h[1];
	if (state_orig.m_U[1] >= 0)
	    sL = 1.0;
	else
	    sL = 0.0;
	state_hat.m_U[0] = state_orig.m_U[0] + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[0];
	state_hat.m_U[1] = state_orig.m_U[1] + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[1];
	state_hat.m_U[2] = state_orig.m_U[2] + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[2];
	break;

    case LOWER:
	getLimitedSlope(icoords,COORD_Z,slope_z_limited);
	dz = top_h[2];
	if (state_orig.m_U[2] <= 0)
	    sL = 1.0;
	else
	    sL = 0.0;
	state_hat.m_U[0] = state_orig.m_U[0] + sL*(-dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[0];
	state_hat.m_U[1] = state_orig.m_U[1] + sL*(-dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[1];
	state_hat.m_U[2] = state_orig.m_U[2] + sL*(-dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[2];
	break;

    case UPPER:
	getLimitedSlope(icoords,COORD_Z,slope_z_limited);
	dz = top_h[2];
	if (state_orig.m_U[2] >= 0)
	    sL = 1.0;
	else
	    sL = 0.0;
	state_hat.m_U[0] = state_orig.m_U[0] + sL*(dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[0];
	state_hat.m_U[1] = state_orig.m_U[1] + sL*(dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[1];
	state_hat.m_U[2] = state_orig.m_U[2] + sL*(dz/2.0 - m_dt/2.0 * state_orig.m_U[2])*slope_z_limited[2];
	break;
    default:
	assert(false);
    }
}


void Incompress_Solver_Smooth_3D_Cartesian::getFaceVelocity_middleStep_bar(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_bar,
	double transverseD[3],
	L_STATE state_hat)
{
    int index;
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double rho = cell_center[index].m_state.m_rho;

    double diffusion[3], gradP[3];
    gradP[0] = cell_center[index].m_state.grad_q[0];
    gradP[1] = cell_center[index].m_state.grad_q[1];
    gradP[2] = cell_center[index].m_state.grad_q[2];

    getDifffusion(icoords,diffusion);
    state_bar.m_U[0] = state_hat.m_U[0] + m_dt/2.0 * (-transverseD[0] + diffusion[0]/rho - gradP[0]/rho);
    state_bar.m_U[1] = state_hat.m_U[1] + m_dt/2.0 * (-transverseD[1] + diffusion[1]/rho - gradP[1]/rho);
    state_bar.m_U[2] = state_hat.m_U[2] + m_dt/2.0 * (-transverseD[2] + diffusion[2]/rho - gradP[2]/rho);

    double coords[3];
    L_STATE source_term;

    getRectangleCenter(index, coords);
    computeSourceTerm_Adv(coords, source_term);

    state_bar.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_bar.m_U[1] += m_dt/2 * source_term.m_U[1];
    state_bar.m_U[2] += m_dt/2 * source_term.m_U[2];
}


void Incompress_Solver_Smooth_3D_Cartesian::getFaceVelocity_middleStep_coupled_bar(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_bar,
	double transverseD[3],
	L_STATE state_hat)
{
    int index;
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    double rho = cell_center[index].m_state.m_rho;

    double diffusion[3], gradP[3];
    gradP[0] = cell_center[index].m_state.grad_q[0];
    gradP[1] = cell_center[index].m_state.grad_q[1];
    gradP[2] = cell_center[index].m_state.grad_q[2];

    getDiffusion_coupled(icoords,diffusion);
    state_bar.m_U[0] = state_hat.m_U[0] + m_dt/2.0 * (-transverseD[0] + diffusion[0]/rho - gradP[0]/rho);
    state_bar.m_U[1] = state_hat.m_U[1] + m_dt/2.0 * (-transverseD[1] + diffusion[1]/rho - gradP[1]/rho);
    state_bar.m_U[2] = state_hat.m_U[2] + m_dt/2.0 * (-transverseD[2] + diffusion[2]/rho - gradP[2]/rho);

    double coords[3];
    L_STATE source_term;

    getRectangleCenter(index, coords);
    computeSourceTerm_Adv(coords, source_term);

    state_bar.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_bar.m_U[1] += m_dt/2 * source_term.m_U[1];
    state_bar.m_U[2] += m_dt/2 * source_term.m_U[2];
}


void Incompress_Solver_Smooth_3D_Cartesian::getFaceVelocity_middleStep(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_face)
{
    int index;
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    L_STATE state;
    state = cell_center[index].m_state;
    double rho = cell_center[index].m_state.m_rho;


    double dx = 0, dy = 0,dz = 0, slope_x_limited[3] = {0,0,0}, slope_y_limited[3] = {0,0,0}, slope_z_limited[3] = {0,0,0};

    switch(dir)
    {
    case WEST:
	dx = -top_h[0];
	break;
    case EAST:
	dx =  top_h[0];
	break;
    case SOUTH:
	dy = -top_h[1];
	break;
    case NORTH:
	dy =  top_h[1];
	break;
    case LOWER:
	dz = -top_h[2];
	break;
    case UPPER:
	dz =  top_h[2];
	break;
    default:
	assert(false);
    }

    state_face.m_U[0] = state.m_U[0];
    state_face.m_U[1] = state.m_U[1];
    state_face.m_U[2] = state.m_U[2];

//    return;

    getLimitedSlope(icoords,COORD_X,slope_x_limited);
    getLimitedSlope(icoords,COORD_Y,slope_y_limited);
    getLimitedSlope(icoords,COORD_Z,slope_z_limited);

    // dx/2, dy/2, dz/2
    state_face.m_U[0] += dx/2 * slope_x_limited[0] + dy/2 * slope_y_limited[0] + dz/2 * slope_z_limited[0];
    state_face.m_U[1] += dx/2 * slope_x_limited[1] + dy/2 * slope_y_limited[1] + dz/2 * slope_z_limited[1];
    state_face.m_U[2] += dx/2 * slope_x_limited[2] + dy/2 * slope_y_limited[2] + dz/2 * slope_z_limited[2];

    //    return;
    // dt/2
    double diffusion[3], gradP[3];
    getDifffusion(icoords,diffusion);
    gradP[0] = cell_center[index].m_state.grad_q[0];
    gradP[1] = cell_center[index].m_state.grad_q[1];
    gradP[2] = cell_center[index].m_state.grad_q[2];

    state_face.m_U[0] += m_dt/2 *
	    (diffusion[0]/rho - state.m_U[0]*slope_x_limited[0] - state.m_U[1]*slope_y_limited[0] - state.m_U[2]*slope_z_limited[0] - gradP[0]/rho);
    state_face.m_U[1] += m_dt/2 *
	    (diffusion[1]/rho - state.m_U[0]*slope_x_limited[1] - state.m_U[1]*slope_y_limited[1] - state.m_U[2]*slope_z_limited[1] - gradP[1]/rho);
    state_face.m_U[2] += m_dt/2 *
    	    (diffusion[2]/rho - state.m_U[0]*slope_x_limited[2] - state.m_U[1]*slope_y_limited[2] - state.m_U[2]*slope_z_limited[2] - gradP[2]/rho);

    //    return;
    // rhs
    double coords[3];
    L_STATE source_term;

    getRectangleCenter(d_index3d(icoords[0],icoords[1],icoords[2],top_gmax), coords);
    computeSourceTerm_Adv(coords, source_term);

    state_face.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_face.m_U[1] += m_dt/2 * source_term.m_U[1];
    state_face.m_U[2] += m_dt/2 * source_term.m_U[2];
}

void Incompress_Solver_Smooth_3D_Cartesian::getFaceVelocity_middleStep_coupled(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_face)
{
    int index;
    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    L_STATE state;
    state = cell_center[index].m_state;
    double rho = cell_center[index].m_state.m_rho;


    double dx = 0, dy = 0,dz = 0, slope_x_limited[3] = {0,0,0}, slope_y_limited[3] = {0,0,0}, slope_z_limited[3] = {0,0,0};

    switch(dir)
    {
    case WEST:
	dx = -top_h[0];
	break;
    case EAST:
	dx =  top_h[0];
	break;
    case SOUTH:
	dy = -top_h[1];
	break;
    case NORTH:
	dy =  top_h[1];
	break;
    case LOWER:
	dz = -top_h[2];
	break;
    case UPPER:
	dz =  top_h[2];
	break;
    default:
	assert(false);
    }

    state_face.m_U[0] = state.m_U[0];
    state_face.m_U[1] = state.m_U[1];
    state_face.m_U[2] = state.m_U[2];

//    return;

    getLimitedSlope(icoords,COORD_X,slope_x_limited);
    getLimitedSlope(icoords,COORD_Y,slope_y_limited);
    getLimitedSlope(icoords,COORD_Z,slope_z_limited);

    // dx/2, dy/2, dz/2
    state_face.m_U[0] += dx/2 * slope_x_limited[0] + dy/2 * slope_y_limited[0] + dz/2 * slope_z_limited[0];
    state_face.m_U[1] += dx/2 * slope_x_limited[1] + dy/2 * slope_y_limited[1] + dz/2 * slope_z_limited[1];
    state_face.m_U[2] += dx/2 * slope_x_limited[2] + dy/2 * slope_y_limited[2] + dz/2 * slope_z_limited[2];

    //    return;
    // dt/2
    double diffusion[3], gradP[3];
    getDiffusion_coupled(icoords,diffusion);
    gradP[0] = cell_center[index].m_state.grad_q[0];
    gradP[1] = cell_center[index].m_state.grad_q[1];
    gradP[2] = cell_center[index].m_state.grad_q[2];

    state_face.m_U[0] += m_dt/2 *
	    (diffusion[0]/rho - state.m_U[0]*slope_x_limited[0] - state.m_U[1]*slope_y_limited[0] - state.m_U[2]*slope_z_limited[0] - gradP[0]/rho);
    state_face.m_U[1] += m_dt/2 *
	    (diffusion[1]/rho - state.m_U[0]*slope_x_limited[1] - state.m_U[1]*slope_y_limited[1] - state.m_U[2]*slope_z_limited[1] - gradP[1]/rho);
    state_face.m_U[2] += m_dt/2 *
    	    (diffusion[2]/rho - state.m_U[0]*slope_x_limited[2] - state.m_U[1]*slope_y_limited[2] - state.m_U[2]*slope_z_limited[2] - gradP[2]/rho);

    //    return;
    // rhs
    double coords[3];
    L_STATE source_term;

    getRectangleCenter(d_index3d(icoords[0],icoords[1],icoords[2],top_gmax), coords);
    computeSourceTerm_Adv(coords, source_term);

    state_face.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_face.m_U[1] += m_dt/2 * source_term.m_U[1];
    state_face.m_U[2] += m_dt/2 * source_term.m_U[2];
}



/**
* compute mu * (Uxx+Uyy+Uzz)
* @param icoords
* @param diffusion
* @param gradP
*/


void Incompress_Solver_Smooth_3D_Cartesian::getDifffusion(
	int *icoords,
	double diffusion[3])
{
    double Uxx[3], Uyy[3], Uzz[3];
    getDU2(icoords,COORD_X,Uxx);
    getDU2(icoords,COORD_Y,Uyy);
    getDU2(icoords,COORD_Z,Uzz);

    double mu = cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_state.m_mu;

    diffusion[0] = mu * (Uxx[0] + Uyy[0] + Uzz[0]);
    diffusion[1] = mu * (Uxx[1] + Uyy[1] + Uzz[1]);
    diffusion[2] = mu * (Uxx[2] + Uyy[2] + Uzz[2]);
}


void Incompress_Solver_Smooth_3D_Cartesian::getDiffusion_coupled(
	int *icoords,
	double diffusion[3])
{
    int index,index_nb[18];
    double mu[6],mu_edge[6],mu0;
    L_STATE Unb,corner_state;
    double U0_nb[18],U1_nb[18],U2_nb[18],U0_center,U1_center,U2_center;
    // U0[6] -- U0[17] U1[6] -- U1[17]  U2[6] -- U2[17] are corner values
    int nb;
    GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
    bool bNoBoundary[6];
    double dh[3],dh0[3],dh1[3];
    double coords[MAXD],corner_coords[MAXD];

    diffusion[0] = 0.0;
    diffusion[1] = 0.0;
    diffusion[2] = 0.0;

    int i = icoords[0];
    int j = icoords[1];
    int k = icoords[2];

    index = d_index3d(i,j,k,top_gmax);
    //6 neighbours of the center cell

    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
    index_nb[5] = d_index3d(i,j,k+1,top_gmax);

    //xy cut neighbours
    index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
    index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
    index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
    index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);

    //yz cut neighbours
    index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
    index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
    index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
    index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);

    //xz cut neighbours
    index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
    index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
    index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
    index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);

    mu0 = cell_center[index].m_state.m_mu;

    U0_center = cell_center[index].m_state.m_U[0];
    U1_center = cell_center[index].m_state.m_U[1];
    U2_center = cell_center[index].m_state.m_U[2];

    for (nb = 0; nb < 6; nb++)
    {
	bNoBoundary[nb] = getNeighborOrBoundaryState(icoords,dir[nb],Unb,m_t_old);
	U0_nb[nb] = Unb.m_U[0];
	U1_nb[nb] = Unb.m_U[1];
	U2_nb[nb] = Unb.m_U[2];
	if(!bNoBoundary[nb])
	{
	    mu[nb] = mu0;
	    mu_edge[nb] = mu0;
	}
	else
	{
	    mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
	    mu_edge[nb] = 0.5*(mu[nb] + mu0);
	}
    }

    // non-cross derivative terms

    dh[0] = top_h[0];
    dh[1] = top_h[1];
    dh[2] = top_h[2];

    if (bNoBoundary[0])
	dh0[0] = top_h[0];
    else
	dh0[0] = top_h[0]/2.0;
    if (bNoBoundary[1])
	dh1[0] = top_h[0];
    else
	dh1[0] = top_h[0]/2.0;

    if (bNoBoundary[2])
	dh0[1] = top_h[1];
    else
	dh0[1] = top_h[1]/2.0;
    if (bNoBoundary[3])
	dh1[1] = top_h[1];
    else
	dh1[1] = top_h[1]/2.0;

    if (bNoBoundary[4])
	dh0[2] = top_h[2];
    else
	dh0[2] = top_h[2]/2.0;
    if (bNoBoundary[5])
	dh1[2] = top_h[2];
    else
	dh1[2] = top_h[2]/2.0;

    diffusion[0] += 2.0*(mu_edge[1]*(U0_nb[1]-U0_center)/dh1[0] - mu_edge[0]*(U0_center-U0_nb[0])/dh0[0])/dh[0];// (2*mu*u_x)_x
    diffusion[1] +=     (mu_edge[1]*(U1_nb[1]-U1_center)/dh1[0] - mu_edge[0]*(U1_center-U1_nb[0])/dh0[0])/dh[0];// (mu*v_x)_x
    diffusion[2] +=     (mu_edge[1]*(U2_nb[1]-U2_center)/dh1[0] - mu_edge[0]*(U2_center-U2_nb[0])/dh0[0])/dh[0];// (mu*w_x)_x

    diffusion[0] +=     (mu_edge[3]*(U0_nb[3]-U0_center)/dh1[1] - mu_edge[2]*(U0_center-U0_nb[2])/dh0[1])/dh[1];// (mu*u_y)_y
    diffusion[1] += 2.0*(mu_edge[3]*(U1_nb[3]-U1_center)/dh1[1] - mu_edge[2]*(U1_center-U1_nb[2])/dh0[1])/dh[1];// (2*mu*v_y)_y
    diffusion[2] +=     (mu_edge[3]*(U2_nb[3]-U2_center)/dh1[1] - mu_edge[2]*(U2_center-U2_nb[2])/dh0[1])/dh[1];// (mu*w_y)_y

    diffusion[0] +=     (mu_edge[5]*(U0_nb[5]-U0_center)/dh1[2] - mu_edge[4]*(U0_center-U0_nb[4])/dh0[2])/dh[2];// (mu*u_z)_z
    diffusion[1] +=     (mu_edge[5]*(U1_nb[5]-U1_center)/dh1[2] - mu_edge[4]*(U1_center-U1_nb[4])/dh0[2])/dh[2];// (mu*v_z)_z
    diffusion[2] += 2.0*(mu_edge[5]*(U2_nb[5]-U2_center)/dh1[2] - mu_edge[4]*(U2_center-U2_nb[4])/dh0[2])/dh[2];// (2*mu*w_z)_z


    // get the coords in the cell center
    getRectangleCenter(index, coords);

    //cross derivative terms

    //traverse the corners on 3 cut planes

    //corner (i-1/2,j-1/2,k)
    /*
    if (!bNoBoundary[0] && bNoBoundary[2])
    {
	U0_nb[6] = U0_nb[0];
	U1_nb[6] = U1_nb[0];
	U2_nb[6] = U2_nb[0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[2])
    {
	U0_nb[6] = U0_nb[2];
	U1_nb[6] = U1_nb[2];
	U2_nb[6] = U2_nb[2];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[2])
    {
	U0_nb[6] = U0_nb[0];
	U1_nb[6] = U1_nb[0];
	U2_nb[6] = U2_nb[0];
    }
    else
    {
	U0_nb[6] = (U0_nb[0]+U0_nb[2]+cell_center[index_nb[6]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[6] = (U1_nb[0]+U1_nb[2]+cell_center[index_nb[6]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[6] = (U2_nb[0]+U2_nb[2]+cell_center[index_nb[6]].m_state.m_U[2]+U2_center)/4.0;
    } */

    if (bNoBoundary[0] && bNoBoundary[2])
    {
        U0_nb[6] = (U0_nb[0]+U0_nb[2]+cell_center[index_nb[6]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[6] = (U1_nb[0]+U1_nb[2]+cell_center[index_nb[6]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[6] = (U2_nb[0]+U2_nb[2]+cell_center[index_nb[6]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] - 0.5*top_h[0];
        corner_coords[1] = coords[1] - 0.5*top_h[1];
        corner_coords[2] = coords[2] - 0.0*top_h[2];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[6] = corner_state.m_U[0];
        U1_nb[6] = corner_state.m_U[1];
        U2_nb[6] = corner_state.m_U[2];
    }

    //corner (i+1/2,j-1/2,k)
    /*
    if (!bNoBoundary[1] && bNoBoundary[2])
    {
	U0_nb[7] = U0_nb[1];
	U1_nb[7] = U1_nb[1];
	U2_nb[7] = U2_nb[1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[2])
    {
	U0_nb[7] = U0_nb[2];
	U1_nb[7] = U1_nb[2];
	U2_nb[7] = U2_nb[2];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[2])
    {
	U0_nb[7] = U0_nb[1];
	U1_nb[7] = U1_nb[1];
	U2_nb[7] = U2_nb[1];
    }
    else
    {
	U0_nb[7] = (U0_nb[1]+U0_nb[2]+cell_center[index_nb[7]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[7] = (U1_nb[1]+U1_nb[2]+cell_center[index_nb[7]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[7] = (U2_nb[1]+U2_nb[2]+cell_center[index_nb[7]].m_state.m_U[2]+U2_center)/4.0;
    } */

    if (bNoBoundary[1] && bNoBoundary[2])
    {
        U0_nb[7] = (U0_nb[1]+U0_nb[2]+cell_center[index_nb[7]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[7] = (U1_nb[1]+U1_nb[2]+cell_center[index_nb[7]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[7] = (U2_nb[1]+U2_nb[2]+cell_center[index_nb[7]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] + 0.5*top_h[0];
        corner_coords[1] = coords[1] - 0.5*top_h[1];
        corner_coords[2] = coords[2] - 0.0*top_h[2];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[7] = corner_state.m_U[0];
        U1_nb[7] = corner_state.m_U[1];
        U2_nb[7] = corner_state.m_U[2];
    }

    //corner (i+1/2,j+1/2,k)
    /*
    if (!bNoBoundary[1] && bNoBoundary[3])
    {
	U0_nb[8] = U0_nb[1];
	U1_nb[8] = U1_nb[1];
	U2_nb[8] = U2_nb[1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[3])
    {
	U0_nb[8] = U0_nb[3];
	U1_nb[8] = U1_nb[3];
	U2_nb[8] = U2_nb[3];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[3])
    {
	U0_nb[8] = U0_nb[1];
	U1_nb[8] = U1_nb[1];
	U2_nb[8] = U2_nb[1];
    }
    else
    {
	U0_nb[8] = (U0_nb[1]+U0_nb[3]+cell_center[index_nb[8]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[8] = (U1_nb[1]+U1_nb[3]+cell_center[index_nb[8]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[8] = (U2_nb[1]+U2_nb[3]+cell_center[index_nb[8]].m_state.m_U[2]+U2_center)/4.0;
    } */

    if (bNoBoundary[1] && bNoBoundary[3])
    {
        U0_nb[8] = (U0_nb[1]+U0_nb[3]+cell_center[index_nb[8]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[8] = (U1_nb[1]+U1_nb[3]+cell_center[index_nb[8]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[8] = (U2_nb[1]+U2_nb[3]+cell_center[index_nb[8]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] + 0.5*top_h[0];
        corner_coords[1] = coords[1] + 0.5*top_h[1];
        corner_coords[2] = coords[2] - 0.0*top_h[2];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[8] = corner_state.m_U[0];
        U1_nb[8] = corner_state.m_U[1];
        U2_nb[8] = corner_state.m_U[2];
    }

    //corner (i-1/2,j+1/2,k)
    /*
    if (!bNoBoundary[0] && bNoBoundary[3])
    {
	U0_nb[9] = U0_nb[0];
	U1_nb[9] = U1_nb[0];
	U2_nb[9] = U2_nb[0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[3])
    {
	U0_nb[9] = U0_nb[3];
	U1_nb[9] = U1_nb[3];
	U2_nb[9] = U2_nb[3];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[3])
    {
	U0_nb[9] = U0_nb[0];
	U1_nb[9] = U1_nb[0];
	U2_nb[9] = U2_nb[0];
    }
    else
    {
	U0_nb[9] = (U0_nb[0]+U0_nb[3]+cell_center[index_nb[9]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[9] = (U1_nb[0]+U1_nb[3]+cell_center[index_nb[9]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[9] = (U2_nb[0]+U2_nb[3]+cell_center[index_nb[9]].m_state.m_U[2]+U2_center)/4.0;
    } */

    if (bNoBoundary[0] && bNoBoundary[3])
    {
        U0_nb[9] = (U0_nb[0]+U0_nb[3]+cell_center[index_nb[9]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[9] = (U1_nb[0]+U1_nb[3]+cell_center[index_nb[9]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[9] = (U2_nb[0]+U2_nb[3]+cell_center[index_nb[9]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] - 0.5*top_h[0];
        corner_coords[1] = coords[1] + 0.5*top_h[1];
        corner_coords[2] = coords[2] - 0.0*top_h[2];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[9] = corner_state.m_U[0];
        U1_nb[9] = corner_state.m_U[1];
        U2_nb[9] = corner_state.m_U[2];
    }

    //corner (i,j-1/2,k-1/2)
    /*
    if (!bNoBoundary[2] && bNoBoundary[4])
    {
	U0_nb[10] = U0_nb[2];
	U1_nb[10] = U1_nb[2];
	U2_nb[10] = U2_nb[2];
    }
    else if(bNoBoundary[2] && !bNoBoundary[4])
    {
	U0_nb[10] = U0_nb[4];
	U1_nb[10] = U1_nb[4];
	U2_nb[10] = U2_nb[4];
    }
    else if(!bNoBoundary[2] && !bNoBoundary[4])
    {
	U0_nb[10] = U0_nb[2];
	U1_nb[10] = U1_nb[2];
	U2_nb[10] = U2_nb[2];
    }
    else
    {
	U0_nb[10] = (U0_nb[2]+U0_nb[4]+cell_center[index_nb[10]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[10] = (U1_nb[2]+U1_nb[4]+cell_center[index_nb[10]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[10] = (U2_nb[2]+U2_nb[4]+cell_center[index_nb[10]].m_state.m_U[2]+U2_center)/4.0;
    } */

    if (bNoBoundary[2] && bNoBoundary[4])
    {
        U0_nb[10] = (U0_nb[2]+U0_nb[4]+cell_center[index_nb[10]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[10] = (U1_nb[2]+U1_nb[4]+cell_center[index_nb[10]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[10] = (U2_nb[2]+U2_nb[4]+cell_center[index_nb[10]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] - 0.0*top_h[0];
        corner_coords[1] = coords[1] - 0.5*top_h[1];
        corner_coords[2] = coords[2] - 0.5*top_h[2];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[10] = corner_state.m_U[0];
        U1_nb[10] = corner_state.m_U[1];
        U2_nb[10] = corner_state.m_U[2];
    }

    //corner (i,j+1/2,k-1/2)
    /*
    if (!bNoBoundary[3] && bNoBoundary[4])
    {
	U0_nb[11] = U0_nb[3];
	U1_nb[11] = U1_nb[3];
	U2_nb[11] = U2_nb[3];
    }
    else if(bNoBoundary[3] && !bNoBoundary[4])
    {
	U0_nb[11] = U0_nb[4];
	U1_nb[11] = U1_nb[4];
	U2_nb[11] = U2_nb[4];
    }
    else if(!bNoBoundary[3] && !bNoBoundary[4])
    {
	U0_nb[11] = U0_nb[3];
	U1_nb[11] = U1_nb[3];
	U2_nb[11] = U2_nb[3];
    }
    else
    {
	U0_nb[11] = (U0_nb[3]+U0_nb[4]+cell_center[index_nb[11]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[11] = (U1_nb[3]+U1_nb[4]+cell_center[index_nb[11]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[11] = (U2_nb[3]+U2_nb[4]+cell_center[index_nb[11]].m_state.m_U[2]+U2_center)/4.0;
    } */

    if (bNoBoundary[3] && bNoBoundary[4])
    {
        U0_nb[11] = (U0_nb[3]+U0_nb[4]+cell_center[index_nb[11]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[11] = (U1_nb[3]+U1_nb[4]+cell_center[index_nb[11]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[11] = (U2_nb[3]+U2_nb[4]+cell_center[index_nb[11]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] - 0.0*top_h[0];
        corner_coords[1] = coords[1] + 0.5*top_h[1];
        corner_coords[2] = coords[2] - 0.5*top_h[2];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[11] = corner_state.m_U[0];
        U1_nb[11] = corner_state.m_U[1];
        U2_nb[11] = corner_state.m_U[2];
    }

    //corner (i,j+1/2,k+1/2)
    /*
    if (!bNoBoundary[3] && bNoBoundary[5])
    {
	U0_nb[12] = U0_nb[3];
	U1_nb[12] = U1_nb[3];
	U2_nb[12] = U2_nb[3];
    }
    else if(bNoBoundary[3] && !bNoBoundary[5])
    {
	U0_nb[12] = U0_nb[5];
	U1_nb[12] = U1_nb[5];
	U2_nb[12] = U2_nb[5];
    }
    else if(!bNoBoundary[3] && !bNoBoundary[5])
    {
	U0_nb[12] = U0_nb[3];
	U1_nb[12] = U1_nb[3];
	U2_nb[12] = U2_nb[3];
    }
    else
    {
	U0_nb[12] = (U0_nb[3]+U0_nb[5]+cell_center[index_nb[12]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[12] = (U1_nb[3]+U1_nb[5]+cell_center[index_nb[12]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[12] = (U2_nb[3]+U2_nb[5]+cell_center[index_nb[12]].m_state.m_U[2]+U2_center)/4.0;
    } */

    if (bNoBoundary[3] && bNoBoundary[5])
    {
        U0_nb[12] = (U0_nb[3]+U0_nb[5]+cell_center[index_nb[12]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[12] = (U1_nb[3]+U1_nb[5]+cell_center[index_nb[12]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[12] = (U2_nb[3]+U2_nb[5]+cell_center[index_nb[12]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] - 0.0*top_h[0];
        corner_coords[1] = coords[1] + 0.5*top_h[1];
        corner_coords[2] = coords[2] + 0.5*top_h[2];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[12] = corner_state.m_U[0];
        U1_nb[12] = corner_state.m_U[1];
        U2_nb[12] = corner_state.m_U[2];
    }


    //corner (i,j-1/2,k+1/2)
    /*
    if (!bNoBoundary[2] && bNoBoundary[5])
    {
	U0_nb[13] = U0_nb[2];
	U1_nb[13] = U1_nb[2];
	U2_nb[13] = U2_nb[2];
    }
    else if(bNoBoundary[2] && !bNoBoundary[5])
    {
	U0_nb[13] = U0_nb[5];
	U1_nb[13] = U1_nb[5];
	U2_nb[13] = U2_nb[5];
    }
    else if(!bNoBoundary[2] && !bNoBoundary[5])
    {
	U0_nb[13] = U0_nb[2];
	U1_nb[13] = U1_nb[2];
	U2_nb[13] = U2_nb[2];
    }
    else
    {
	U0_nb[13] = (U0_nb[2]+U0_nb[5]+cell_center[index_nb[13]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[13] = (U1_nb[2]+U1_nb[5]+cell_center[index_nb[13]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[13] = (U2_nb[2]+U2_nb[5]+cell_center[index_nb[13]].m_state.m_U[2]+U2_center)/4.0;
    } */

    if (bNoBoundary[2] && bNoBoundary[5])
    {
        U0_nb[13] = (U0_nb[2]+U0_nb[5]+cell_center[index_nb[13]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[13] = (U1_nb[2]+U1_nb[5]+cell_center[index_nb[13]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[13] = (U2_nb[2]+U2_nb[5]+cell_center[index_nb[13]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] - 0.0*top_h[0];
        corner_coords[1] = coords[1] - 0.5*top_h[1];
        corner_coords[2] = coords[2] + 0.5*top_h[2];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[13] = corner_state.m_U[0];
        U1_nb[13] = corner_state.m_U[1];
        U2_nb[13] = corner_state.m_U[2];
    }

    //corner (i-1/2,j,k-1/2)
    /*
    if (!bNoBoundary[0] && bNoBoundary[4])
    {
	U0_nb[14] = U0_nb[0];
	U1_nb[14] = U1_nb[0];
	U2_nb[14] = U2_nb[0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[4])
    {
	U0_nb[14] = U0_nb[4];
	U1_nb[14] = U1_nb[4];
	U2_nb[14] = U2_nb[4];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[4])
    {
	U0_nb[14] = U0_nb[0];
	U1_nb[14] = U1_nb[0];
	U2_nb[14] = U2_nb[0];
    }
    else
    {
	U0_nb[14] = (U0_nb[0]+U0_nb[4]+cell_center[index_nb[14]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[14] = (U1_nb[0]+U1_nb[4]+cell_center[index_nb[14]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[14] = (U2_nb[0]+U2_nb[4]+cell_center[index_nb[14]].m_state.m_U[2]+U2_center)/4.0;
    } */

    if (bNoBoundary[0] && bNoBoundary[4])
    {
        U0_nb[14] = (U0_nb[0]+U0_nb[4]+cell_center[index_nb[14]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[14] = (U1_nb[0]+U1_nb[4]+cell_center[index_nb[14]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[14] = (U2_nb[0]+U2_nb[4]+cell_center[index_nb[14]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] - 0.5*top_h[0];
        corner_coords[1] = coords[1] - 0.0*top_h[1];
        corner_coords[2] = coords[2] - 0.5*top_h[2];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[14] = corner_state.m_U[0];
        U1_nb[14] = corner_state.m_U[1];
        U2_nb[14] = corner_state.m_U[2];
    }

    //corner (i+1/2,j,k-1/2)
    /*
    if (!bNoBoundary[1] && bNoBoundary[4])
    {
	U0_nb[15] = U0_nb[1];
	U1_nb[15] = U1_nb[1];
	U2_nb[15] = U2_nb[1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[4])
    {
	U0_nb[15] = U0_nb[4];
	U1_nb[15] = U1_nb[4];
	U2_nb[15] = U2_nb[4];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[4])
    {
	U0_nb[15] = U0_nb[1];
	U1_nb[15] = U1_nb[1];
	U2_nb[15] = U2_nb[1];
    }
    else
    {
	U0_nb[15] = (U0_nb[1]+U0_nb[4]+cell_center[index_nb[15]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[15] = (U1_nb[1]+U1_nb[4]+cell_center[index_nb[15]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[15] = (U2_nb[1]+U2_nb[4]+cell_center[index_nb[15]].m_state.m_U[2]+U2_center)/4.0;
    } */

    if (bNoBoundary[1] && bNoBoundary[4])
    {
        U0_nb[15] = (U0_nb[1]+U0_nb[4]+cell_center[index_nb[15]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[15] = (U1_nb[1]+U1_nb[4]+cell_center[index_nb[15]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[15] = (U2_nb[1]+U2_nb[4]+cell_center[index_nb[15]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] + 0.5*top_h[0];
        corner_coords[1] = coords[1] - 0.0*top_h[1];
        corner_coords[2] = coords[2] - 0.5*top_h[2];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[15] = corner_state.m_U[0];
        U1_nb[15] = corner_state.m_U[1];
        U2_nb[15] = corner_state.m_U[2];
    }

    //corner (i+1/2,j,k+1/2)
    /*
    if (!bNoBoundary[1] && bNoBoundary[5])
    {
	U0_nb[16] = U0_nb[1];
	U1_nb[16] = U1_nb[1];
	U2_nb[16] = U2_nb[1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[5])
    {
	U0_nb[16] = U0_nb[5];
	U1_nb[16] = U1_nb[5];
	U2_nb[16] = U2_nb[5];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[5])
    {
	U0_nb[16] = U0_nb[1];
	U1_nb[16] = U1_nb[1];
	U2_nb[16] = U2_nb[1];
    }
    else
    {
	U0_nb[16] = (U0_nb[1]+U0_nb[5]+cell_center[index_nb[16]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[16] = (U1_nb[1]+U1_nb[5]+cell_center[index_nb[16]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[16] = (U2_nb[1]+U2_nb[5]+cell_center[index_nb[16]].m_state.m_U[2]+U2_center)/4.0;
    } */

    if (bNoBoundary[1] && bNoBoundary[5])
    {
        U0_nb[16] = (U0_nb[1]+U0_nb[5]+cell_center[index_nb[16]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[16] = (U1_nb[1]+U1_nb[5]+cell_center[index_nb[16]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[16] = (U2_nb[1]+U2_nb[5]+cell_center[index_nb[16]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] + 0.5*top_h[0];
        corner_coords[1] = coords[1] - 0.0*top_h[1];
        corner_coords[2] = coords[2] + 0.5*top_h[2];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[16] = corner_state.m_U[0];
        U1_nb[16] = corner_state.m_U[1];
        U2_nb[16] = corner_state.m_U[2];
    }

    //corner (i-1/2,j,k+1/2)
    /*
    if (!bNoBoundary[0] && bNoBoundary[5])
    {
	U0_nb[17] = U0_nb[0];
	U1_nb[17] = U1_nb[0];
	U2_nb[17] = U2_nb[0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[5])
    {
	U0_nb[17] = U0_nb[5];
	U1_nb[17] = U1_nb[5];
	U2_nb[17] = U2_nb[5];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[5])
    {
	U0_nb[17] = U0_nb[0];
	U1_nb[17] = U1_nb[0];
	U2_nb[17] = U2_nb[0];
    }
    else
    {
	U0_nb[17] = (U0_nb[0]+U0_nb[5]+cell_center[index_nb[17]].m_state.m_U[0]+U0_center)/4.0;
	U1_nb[17] = (U1_nb[0]+U1_nb[5]+cell_center[index_nb[17]].m_state.m_U[1]+U1_center)/4.0;
	U2_nb[17] = (U2_nb[0]+U2_nb[5]+cell_center[index_nb[17]].m_state.m_U[2]+U2_center)/4.0;
    } */

    if (bNoBoundary[0] && bNoBoundary[5])
    {
        U0_nb[17] = (U0_nb[0]+U0_nb[5]+cell_center[index_nb[17]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[17] = (U1_nb[0]+U1_nb[5]+cell_center[index_nb[17]].m_state.m_U[1]+U1_center)/4.0;
        U2_nb[17] = (U2_nb[0]+U2_nb[5]+cell_center[index_nb[17]].m_state.m_U[2]+U2_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] - 0.5*top_h[0];
        corner_coords[1] = coords[1] - 0.0*top_h[1];
        corner_coords[2] = coords[2] + 0.5*top_h[2];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[17] = corner_state.m_U[0];
        U1_nb[17] = corner_state.m_U[1];
        U2_nb[17] = corner_state.m_U[2];
    }

    diffusion[0] += (mu_edge[2]*U1_nb[6]- mu_edge[2]*U1_nb[7]+ mu_edge[3]*U1_nb[8]- mu_edge[3]*U1_nb[9])/ (top_h[0]*top_h[1]);//(mu*v_x)_y
    diffusion[0] += (mu_edge[4]*U2_nb[14]-mu_edge[4]*U2_nb[15]+mu_edge[5]*U2_nb[16]-mu_edge[5]*U2_nb[17])/(top_h[0]*top_h[2]);//(mu*w_x)_z

    diffusion[1] += (mu_edge[0]*U0_nb[6]- mu_edge[1]*U0_nb[7]+ mu_edge[1]*U0_nb[8]- mu_edge[0]*U0_nb[9])/ (top_h[0]*top_h[1]);//(mu*u_y)_x
    diffusion[1] += (mu_edge[4]*U2_nb[10]-mu_edge[4]*U2_nb[11]+mu_edge[5]*U2_nb[12]-mu_edge[5]*U2_nb[13])/(top_h[1]*top_h[2]);//(mu*w_y)_z


    diffusion[2] += (mu_edge[0]*U0_nb[14]-mu_edge[1]*U0_nb[15]+mu_edge[1]*U0_nb[16]-mu_edge[0]*U0_nb[17])/(top_h[0]*top_h[2]);//(mu*u_z)_x
    diffusion[2] += (mu_edge[2]*U1_nb[10]-mu_edge[3]*U1_nb[11]+mu_edge[3]*U1_nb[12]-mu_edge[2]*U1_nb[13])/(top_h[1]*top_h[2]);//(mu*v_z)_y

}


/**
* calculate Uxx,Uyy,Uzz and Px,Py,Pz.
* @param dir
* @param icoords
* @param dU2
* @param dP
*/
void Incompress_Solver_Smooth_3D_Cartesian::getDU2(
	int *icoords,
	EBM_COORD xyz,
	double dU2[3])
{
    double dh0, dh1, dh;
    L_STATE U0, U1, U2;

    U1 = cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_state;

    bool bNoBoundary[2];

    if(xyz==COORD_X)
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,WEST,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[0];
	else
	    dh0 = top_h[0]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,EAST,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[0];
	else
	    dh1 = top_h[0]/2;

	dh = top_h[0];
    }
    else if(xyz==COORD_Y)	//
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,SOUTH,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[1];
	else
	    dh0 = top_h[1]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,NORTH,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[1];
	else
	    dh1 = top_h[1]/2;

	dh = top_h[1];
    }
    else
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,LOWER,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[2];
	else
	    dh0 = top_h[2]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,UPPER,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[2];
	else
	    dh1 = top_h[2]/2;

	dh = top_h[2];
    }

    dU2[0] = ((U2.m_U[0] - U1.m_U[0])/dh1 - (U1.m_U[0] - U0.m_U[0])/dh0) / dh;
    dU2[1] = ((U2.m_U[1] - U1.m_U[1])/dh1 - (U1.m_U[1] - U0.m_U[1])/dh0) / dh;
    dU2[2] = ((U2.m_U[2] - U1.m_U[2])/dh1 - (U1.m_U[2] - U0.m_U[2])/dh0) / dh;

}

// Minmod slope limiter
void Incompress_Solver_Smooth_3D_Cartesian::getLimitedSlope(
	int *icoords,
	EBM_COORD xyz,
	double slope[3])
{
    double dh0, dh1;
    L_STATE U0, U1, U2;

    U1 = cell_center[d_index3d(icoords[0],icoords[1],icoords[2],top_gmax)].m_state;

    bool bNoBoundary[2];

    if(xyz==COORD_X)
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,WEST,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[0];
	else
	    dh0 = top_h[0]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,EAST,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[0];
	else
	    dh1 = top_h[0]/2;
    }
    else if(xyz==COORD_Y)	//
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,SOUTH,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[1];
	else
	    dh0 = top_h[1]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,NORTH,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[1];
	else
	    dh1 = top_h[1]/2;
    }
    else  //xyz == COORD_Z
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,LOWER,U0,m_t_old);
	if(bNoBoundary[0])
	    dh0 = top_h[1];
	else
	    dh0 = top_h[1]/2;

	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,UPPER,U2,m_t_old);
	if(bNoBoundary[1])
	    dh1 = top_h[1];
	else
	    dh1 = top_h[1]/2;
    }
    slope[0] = EBM_minmod((U1.m_U[0]-U0.m_U[0])/dh0, (U2.m_U[0]-U1.m_U[0])/dh1);
    slope[1] = EBM_minmod((U1.m_U[1]-U0.m_U[1])/dh0, (U2.m_U[1]-U1.m_U[1])/dh1);
    slope[2] = EBM_minmod((U1.m_U[2]-U0.m_U[2])/dh0, (U2.m_U[2]-U1.m_U[2])/dh1);
}

// Van Leer Slope Limiter


void Incompress_Solver_Smooth_3D_Cartesian::getLimitedSlope_Vanleer(
	int *icoords,
	EBM_COORD xyz,
	double slope[3])
{
    double dx,dy,dz;
    L_STATE U0, U1, U2;
    double u_lim, v_lim, w_lim;
    double u_slope, v_slope, w_slope;

    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);

    U1 = cell_center[index].m_state;
    dx = top_h[0];
    dy = top_h[1];
    dz = top_h[2];

    bool bNoBoundary[2];


    if(xyz==COORD_X)
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,WEST,U0,m_t_old);
	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,EAST,U2,m_t_old);
	if(bNoBoundary[0] || bNoBoundary[1])
	{

	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dx, 2.0*(U1.m_U[0]-U0.m_U[0])/dx);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dx, 2.0*(U1.m_U[1]-U0.m_U[1])/dx);
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/dx, 2.0*(U1.m_U[2]-U0.m_U[2])/dx);

	    u_slope = (U2.m_U[0] - U0.m_U[0])/(2*dx);
	    v_slope = (U2.m_U[1] - U0.m_U[1])/(2*dx);
	    w_slope = (U2.m_U[2] - U0.m_U[2])/(2*dx);


	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);

	}
	else if ( (!bNoBoundary[0]) || bNoBoundary[1])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dx, 2.0*(U1.m_U[0]-U0.m_U[0])/(dx/2.0));
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dx, 2.0*(U1.m_U[1]-U0.m_U[1])/(dx/2.0));
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/dx, 2.0*(U1.m_U[2]-U0.m_U[2])/(dx/2.0));

	    u_slope = (U2.m_U[0] + U1.m_U[0] + 2.0*U0.m_U[0])/(2.0*dx);
	    v_slope = (U2.m_U[1] + U1.m_U[1] + 2.0*U0.m_U[1])/(2.0*dx);
	    w_slope = (U2.m_U[2] + U1.m_U[2] + 2.0*U0.m_U[2])/(2.0*dx);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);
	}
	else if ( (!bNoBoundary[1]) || bNoBoundary[0])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/(dx/2.0), 2.0*(U1.m_U[0]-U0.m_U[0])/dx);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/(dx/2.0), 2.0*(U1.m_U[1]-U0.m_U[1])/dx);
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/(dx/2.0), 2.0*(U1.m_U[2]-U0.m_U[2])/dx);

	    u_slope = (2.0*U2.m_U[0] - U1.m_U[0] - U0.m_U[0])/(2.0*dx);
	    v_slope = (2.0*U2.m_U[1] - U1.m_U[1] - U0.m_U[1])/(2.0*dx);
	    w_slope = (2.0*U2.m_U[2] - U1.m_U[2] - U0.m_U[2])/(2.0*dx);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);

	}
	else
	{
	    printf("\nThe number of cell in x direction is less than 1!!\n");
	    slope[0] = slope[1] = slope[2] = 0.0;
	}

    }
    else if (xyz == COORD_Y)	//
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,SOUTH,U0,m_t_old);
	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,NORTH,U2,m_t_old);

	if(bNoBoundary[0] || bNoBoundary[1])
	{

	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dy, 2.0*(U1.m_U[0]-U0.m_U[0])/dy);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dy, 2.0*(U1.m_U[1]-U0.m_U[1])/dy);
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/dy, 2.0*(U1.m_U[2]-U0.m_U[2])/dy);

	    u_slope = (U2.m_U[0] - U0.m_U[0])/(2*dy);
	    v_slope = (U2.m_U[1] - U0.m_U[1])/(2*dy);
	    w_slope = (U2.m_U[2] - U0.m_U[2])/(2*dy);


	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);

	}
	else if ( (!bNoBoundary[0]) || bNoBoundary[1])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dy, 2.0*(U1.m_U[0]-U0.m_U[0])/(dy/2.0));
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dy, 2.0*(U1.m_U[1]-U0.m_U[1])/(dy/2.0));
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/dy, 2.0*(U1.m_U[2]-U0.m_U[2])/(dy/2.0));

	    u_slope = (U2.m_U[0] + U1.m_U[0] + 2.0*U0.m_U[0])/(2.0*dy);
	    v_slope = (U2.m_U[1] + U1.m_U[1] + 2.0*U0.m_U[1])/(2.0*dy);
	    w_slope = (U2.m_U[2] + U1.m_U[2] + 2.0*U0.m_U[2])/(2.0*dy);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);
	}
	else if ( (!bNoBoundary[1]) || bNoBoundary[0])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/(dy/2.0), 2.0*(U1.m_U[0]-U0.m_U[0])/dy);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/(dy/2.0), 2.0*(U1.m_U[1]-U0.m_U[1])/dy);
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/(dy/2.0), 2.0*(U1.m_U[2]-U0.m_U[2])/dy);

	    u_slope = (2.0*U2.m_U[0] - U1.m_U[0] - U0.m_U[0])/(2.0*dy);
	    v_slope = (2.0*U2.m_U[1] - U1.m_U[1] - U0.m_U[1])/(2.0*dy);
	    w_slope = (2.0*U2.m_U[2] - U1.m_U[2] - U0.m_U[2])/(2.0*dy);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);

	}
	else
	{
	    printf("\nThe number of cell in y direction is less than 1!!\n");
	    slope[0] = slope[1] = slope[2] = 0.0;
	}
    }
    else // xyz == COORD_Z
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,LOWER,U0,m_t_old);
	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,UPPER,U2,m_t_old);

	if(bNoBoundary[0] || bNoBoundary[1])
	{

	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dz, 2.0*(U1.m_U[0]-U0.m_U[0])/dz);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dz, 2.0*(U1.m_U[1]-U0.m_U[1])/dz);
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/dz, 2.0*(U1.m_U[2]-U0.m_U[2])/dz);

	    u_slope = (U2.m_U[0] - U0.m_U[0])/(2*dz);
	    v_slope = (U2.m_U[1] - U0.m_U[1])/(2*dz);
	    w_slope = (U2.m_U[2] - U0.m_U[2])/(2*dz);


	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);

	}
	else if ( (!bNoBoundary[0]) || bNoBoundary[1])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dz, 2.0*(U1.m_U[0]-U0.m_U[0])/(dz/2.0));
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dz, 2.0*(U1.m_U[1]-U0.m_U[1])/(dz/2.0));
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/dz, 2.0*(U1.m_U[2]-U0.m_U[2])/(dz/2.0));

	    u_slope = (U2.m_U[0] + U1.m_U[0] + 2.0*U0.m_U[0])/(2.0*dz);
	    v_slope = (U2.m_U[1] + U1.m_U[1] + 2.0*U0.m_U[1])/(2.0*dz);
	    w_slope = (U2.m_U[2] + U1.m_U[2] + 2.0*U0.m_U[2])/(2.0*dz);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);
	}
	else if ( (!bNoBoundary[1]) || bNoBoundary[0])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/(dz/2.0), 2.0*(U1.m_U[0]-U0.m_U[0])/dz);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/(dz/2.0), 2.0*(U1.m_U[1]-U0.m_U[1])/dz);
	    w_lim = EBM_minmod(2.0*(U2.m_U[2]-U1.m_U[2])/(dz/2.0), 2.0*(U1.m_U[2]-U0.m_U[2])/dz);

	    u_slope = (2.0*U2.m_U[0] - U1.m_U[0] - U0.m_U[0])/(2.0*dz);
	    v_slope = (2.0*U2.m_U[1] - U1.m_U[1] - U0.m_U[1])/(2.0*dz);
	    w_slope = (2.0*U2.m_U[2] - U1.m_U[2] - U0.m_U[2])/(2.0*dz);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	    slope[2] = EBM_minmod(w_slope, w_lim);

	}
	else
	{
	    printf("\nThe number of cell in z direction is less than 1!!\n");
	    slope[0] = slope[1] = slope[2] = 0.0;
	}
    }
}


double Incompress_Solver_Smooth_3D_Cartesian::EBM_minmod(
	double x,
	double y)
{

    double sign = x*y;

    if(sign<0)
	return 0;
    else if(sign>=0)
    {
	if(fabs(x)<fabs(y))
	    return x;
	else
	    return y;
    }

    return 0;
}

/**
* get the state from the neighbor cell or boundary.
* @param icoords
* @param dir
* @param comp
* @param state
* @return true,  valid state from neighbor cell
* 	   false, valid state from boundary
*/
bool Incompress_Solver_Smooth_3D_Cartesian::getNeighborOrBoundaryState(
	int icoords[3],
	GRID_DIRECTION dir,
	L_STATE &state,
	double t)
{
    double crx_coords[MAXD];
    static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};
    POINTER intfc_state;
    HYPER_SURF *hs;

    int index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
    int comp = cell_center[index].comp;

    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir,
	    comp,&intfc_state,&hs,crx_coords,t) &&
	    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
    {
	state.m_U[0] = getStateVel[0](intfc_state);
	state.m_U[1] = getStateVel[1](intfc_state);
	state.m_U[2] = getStateVel[2](intfc_state);
	return false;
    }
    else
    {
	int index_nb;
	switch(dir)
	{
	case WEST:
	    index_nb = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
	    break;
	case EAST:
	    index_nb = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
	    break;
	case SOUTH:
	    index_nb = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
	    break;
	case NORTH:
	    index_nb = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
	    break;
	case LOWER:
	    index_nb = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
	    break;
	case UPPER:
	    index_nb = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);
	    break;
	default:
	    assert(false);
	}
	state = cell_center[index_nb].m_state;
	return true;
    }
}


/**
*
* @param state_left
* @param state_right
* @param ans
*/
void Incompress_Solver_Smooth_3D_Cartesian::getRiemannSolution(
	EBM_COORD xyz,
	L_STATE &state_left,
	L_STATE &state_right,
	L_STATE &ans)
{
    L_STATE sl, sr;


    // rotate state
    if(xyz==COORD_X)
    {
	sl.m_U[0] = state_left.m_U[0];
	sl.m_U[1] = state_left.m_U[1];
	sl.m_U[2] = state_left.m_U[2];

	sr.m_U[0] = state_right.m_U[0];
	sr.m_U[1] = state_right.m_U[1];
	sr.m_U[2] = state_right.m_U[2];
    }
    else if(xyz==COORD_Y)
    {
	sl.m_U[0] = state_left.m_U[1];
	sl.m_U[1] = state_left.m_U[0];
	sl.m_U[2] = state_left.m_U[2];

	sr.m_U[0] = state_right.m_U[1];
	sr.m_U[1] = state_right.m_U[0];
	sr.m_U[2] = state_right.m_U[2];
    }
    else
    {
	sl.m_U[0] = state_left.m_U[2];
	sl.m_U[1] = state_left.m_U[1];
	sl.m_U[2] = state_left.m_U[0];

	sr.m_U[0] = state_right.m_U[2];
	sr.m_U[1] = state_right.m_U[1];
	sr.m_U[2] = state_right.m_U[0];
    }

    // calculate the Riemann solution
    double uL = sl.m_U[0];
    double uR = sr.m_U[0];

    // BCG, JCP 85, 257-283 (1989)
    // ut + uux = 0
    if(uL>=0 && (uL+uR)>=0)
	ans.m_U[0] = uL;
    else if(uL<0 && uR>0)
	ans.m_U[0] = 0;
    else
	ans.m_U[0] = uR;

    // vt + uvx = 0
    if(ans.m_U[0]>0)
	ans.m_U[1] = sl.m_U[1];
    else if(ans.m_U[0]<0)
	ans.m_U[1] = sr.m_U[1];
    else
	ans.m_U[1] = 1.0/2*(sl.m_U[1]+sr.m_U[1]);

    if(ans.m_U[0]>0)
	ans.m_U[2] = sl.m_U[2];
    else if(ans.m_U[0]<0)
	ans.m_U[2] = sr.m_U[2];
    else
	ans.m_U[2] = 1.0/2*(sl.m_U[2]+sr.m_U[2]);

    // rotate state
    if(xyz==COORD_X)
	; // do nothing
    else if(xyz==COORD_Y)
	std::swap(ans.m_U[0],ans.m_U[1]);
    else
	std::swap(ans.m_U[0],ans.m_U[2]);
}


//rho*nu_turbulent for momemtum eqn and nu_turbulent/Sc_turbulent for mass diffusion
void Incompress_Solver_Smooth_3D_Cartesian::computeSubgridModel_vd(void)
{
    int i,j,k,l;
    int index, index1, index2, index_nb[18];
    int index_nbb, coord;
    int size;
    L_STATE state;
    double *u, *v, *w;
    double *rho, *c, rho_face, c_face;
    double ux,uy,uz,vx,vy,vz,wx,wy,wz;
    double *s, *s11, *s12, *s13, *s22, *s23, *s33;
    double *vel_u, *vel_v, *vel_w;
    double *cx, *cy, *cz;
    double sum_p,sum_pu,sum_pv,sum_pw,sum_pc; //suppose k=1 for sum_pc
    double sum_puu,sum_pvv,sum_pww,sum_puv,sum_pvw,sum_puw,sum_puc,sum_pvc,sum_pwc;
    double sum_s11,sum_s12,sum_s13,sum_s22,sum_s23,sum_s33;
    double sum_pss11,sum_pss12,sum_pss13,sum_pss22,sum_pss23,sum_pss33,sum_s;
    double sum_pscx,sum_pscy,sum_pscz,sum_cx,sum_cy,sum_cz;
    double ma11, ma12, ma13, ma22, ma23, ma33, la11, la12, la13, la22, la23, la33;
    double my1, my2, my3, ly1, ly2, ly3; //suppose k=1
    double cs,max_cs[MAXD+1],cc,max_cc[MAXD+1];
    double max_turbulent_kinematic_viscosity[MAXD+1],max_turbulent_diffusity[MAXD+1];
    double avg_turbulent_kinematic_viscosity[MAXD+1],avg_turbulent_diffusity[MAXD+1];
    int    count;
    int    ii,jj,kk,iii,jjj,kkk;
    int    NB = 2;
    int    NBC = pow(NB,3);
    double gdeno, gnume;
    double ***ldeno, ***lnume;
    double ***cdeno, ***cnume;
    double delta2, tdelta2;
    delta2 = sqr(pow((top_h[0]*top_h[1]*top_h[2]),(1.0/3.0)));
    tdelta2 = sqr(pow((NB*top_h[0]*NB*top_h[1]*NB*top_h[2]),(1.0/3.0)));

    size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
    FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&w,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&rho,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&c,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&cx,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&cy,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&cz,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&s,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&s11,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&s12,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&s13,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&s22,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&s23,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&s33,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&vel_u,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&vel_v,size,sizeof(double));
    FT_VectorMemoryAlloc((POINTER*)&vel_w,size,sizeof(double));
    FT_TriArrayMemoryAlloc((POINTER*)&ldeno,(imax-imin+1)/NB,(jmax-jmin+1)/NB,(kmax-kmin+1)/NB,sizeof(double));
    FT_TriArrayMemoryAlloc((POINTER*)&lnume,(imax-imin+1)/NB,(jmax-jmin+1)/NB,(kmax-kmin+1)/NB,sizeof(double));
    FT_TriArrayMemoryAlloc((POINTER*)&cdeno,(imax-imin+1)/NB,(jmax-jmin+1)/NB,(kmax-kmin+1)/NB,sizeof(double));
    FT_TriArrayMemoryAlloc((POINTER*)&cnume,(imax-imin+1)/NB,(jmax-jmin+1)/NB,(kmax-kmin+1)/NB,sizeof(double));

    const int  nn = pp_numnodes();
    int        myid = pp_mynode();
    int   *ppgmax = front->pp_grid->gmax;
    int   ppx = myid % ppgmax[0];
    int   ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
    int   ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

    for (k = 0; k <= top_gmax[2]; k++)
    for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
    {
        index = d_index3d(i,j,k,top_gmax);

        //LOWER bdry
        if (ppz==0 && k==kmin-1)
        {
            index_nbb = d_index3d(i,j,kmin,top_gmax);
            u[index] = -cell_center[index_nbb].m_state.m_U[0];
            v[index] = -cell_center[index_nbb].m_state.m_U[1];
            //w[4] = 0
            w[index] = 0;
            rho[index] = cell_center[index_nbb].m_state.m_rho;
            c[index] = cell_center[index_nbb].m_state.m_c;
            continue;
        }

        //UPPER bdry
        if (ppz==ppgmax[2]-1 && k==kmax+1)
        {
            index_nbb = d_index3d(i,j,kmax,top_gmax);
            u[index] = -cell_center[index_nbb].m_state.m_U[0];
            v[index] = -cell_center[index_nbb].m_state.m_U[1];
            rho[index] = cell_center[index_nbb].m_state.m_rho;
            c[index] = cell_center[index_nbb].m_state.m_c;
            //w[5] = -w[4]
            index_nbb = d_index3d(i,j,kmax-1,top_gmax);
            w[index] = -cell_center[index_nbb].m_state.m_U[2];
            continue;
        }

        u[index] = cell_center[index].m_state.m_U[0];
        v[index] = cell_center[index].m_state.m_U[1];
        w[index] = cell_center[index].m_state.m_U[2];
        rho[index] = cell_center[index].m_state.m_rho;
        c[index] = cell_center[index].m_state.m_c;
    }

    //calc turbulent viscosities and diffusities on 3 cell-faces/cell center in each cell
    for (coord = 0; coord < MAXD+1; coord++)
    {
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);

            //6 neighbours of the center cell
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);

            //xy cut neighbours
            index_nb[6] = d_index3d(i-1,j-1,k,top_gmax);
            index_nb[7] = d_index3d(i+1,j-1,k,top_gmax);
            index_nb[8] = d_index3d(i+1,j+1,k,top_gmax);
            index_nb[9] = d_index3d(i-1,j+1,k,top_gmax);

            //yz cut neighbours
            index_nb[10] = d_index3d(i,j-1,k-1,top_gmax);
            index_nb[11] = d_index3d(i,j+1,k-1,top_gmax);
            index_nb[12] = d_index3d(i,j+1,k+1,top_gmax);
            index_nb[13] = d_index3d(i,j-1,k+1,top_gmax);

            //xz cut neighbours
            index_nb[14] = d_index3d(i-1,j,k-1,top_gmax);
            index_nb[15] = d_index3d(i+1,j,k-1,top_gmax);
            index_nb[16] = d_index3d(i+1,j,k+1,top_gmax);
            index_nb[17] = d_index3d(i-1,j,k+1,top_gmax);

            //u-face
            if (coord==0) {
                ux = (u[index_nb[1]]-u[index_nb[0]]) / (2.0*top_h[0]);
                uy = (u[index_nb[3]]-u[index_nb[2]]) / (2.0*top_h[1]);
                uz = (u[index_nb[5]]-u[index_nb[4]]) / (2.0*top_h[2]);

                vx = (v[index_nb[1]]+v[index_nb[7]]-v[index_nb[2]]-v[index]) / (2.0*top_h[0]);
                vy = (v[index_nb[1]]+v[index]-v[index_nb[2]]-v[index_nb[7]]) / (2.0*top_h[1]);
                index1 = d_index3d(i+1,j-1,k+1,top_gmax);
                index2 = d_index3d(i+1,j-1,k-1,top_gmax);
                vz = (v[index_nb[5]]+v[index_nb[16]]+v[index_nb[13]]+v[index1]
                     -v[index_nb[4]]-v[index_nb[15]]-v[index_nb[10]]-v[index2]) / (8.0*top_h[2]);

                wx = (w[index_nb[1]]+w[index_nb[15]]-w[index_nb[4]]-w[index]) / (2.0*top_h[0]);
                index1 = d_index3d(i+1,j+1,k-1,top_gmax);
                index2 = d_index3d(i+1,j-1,k-1,top_gmax);
                wy = (w[index_nb[3]]+w[index_nb[11]]+w[index_nb[8]]+w[index1]
                     -w[index_nb[2]]-w[index_nb[10]]-w[index_nb[7]]-w[index2]) / (8.0*top_h[1]);
                wz = (w[index_nb[1]]+w[index]-w[index_nb[4]]-w[index_nb[15]]) / (2.0*top_h[2]);

                cx[index] = (c[index_nb[1]]-c[index]) / top_h[0];
                cy[index] = (c[index_nb[3]]+c[index_nb[8]]-c[index_nb[2]]-c[index_nb[7]]) / (4.0*top_h[1]);
                cz[index] = (c[index_nb[5]]+c[index_nb[16]]-c[index_nb[4]]-c[index_nb[15]]) / (4.0*top_h[2]);

                vel_u[index] = u[index];
                vel_v[index] = (v[index_nb[1]]+v[index_nb[7]]+v[index_nb[2]]+v[index])/4;
                vel_w[index] = (w[index_nb[1]]+w[index_nb[15]]+w[index_nb[4]]+w[index])/4;
            }

            //v-face
            else if (coord==1) {
                ux = (u[index_nb[3]]+u[index]-u[index_nb[0]]-u[index_nb[9]]) / (2.0*top_h[0]);
                uy = (u[index_nb[3]]+u[index_nb[9]]-u[index_nb[0]]-u[index]) / (2.0*top_h[1]);
                index1 = d_index3d(i-1,j+1,k+1,top_gmax);
                index2 = d_index3d(i-1,j+1,k-1,top_gmax);
                uz = (u[index_nb[5]]+u[index_nb[17]]+u[index_nb[12]]+u[index1]
                     -u[index_nb[4]]-u[index_nb[14]]-u[index_nb[11]]-u[index2]) / (8.0*top_h[2]);

                vx = (v[index_nb[1]]-v[index_nb[0]]) / (2.0*top_h[0]);
                vy = (v[index_nb[3]]-v[index_nb[2]]) / (2.0*top_h[1]);
                vz = (v[index_nb[5]]-v[index_nb[4]]) / (2.0*top_h[2]);

                index1 = d_index3d(i+1,j+1,k-1,top_gmax);
                index2 = d_index3d(i-1,j+1,k-1,top_gmax);
                wx = (w[index_nb[1]]+w[index_nb[15]]+w[index_nb[8]]+w[index1]
                     -w[index_nb[0]]-w[index_nb[14]]-w[index_nb[9]]-w[index2]) / (8.0*top_h[0]);
                wy = (w[index_nb[3]]+w[index_nb[11]]-w[index_nb[4]]-w[index]) / (2.0*top_h[1]);
                wz = (w[index_nb[3]]+w[index]-w[index_nb[4]]-w[index_nb[11]]) / (2.0*top_h[2]);

                cx[index] = (c[index_nb[1]]+c[index_nb[8]]-c[index_nb[0]]-c[index_nb[9]]) / (4.0*top_h[0]);
                cy[index] = (c[index_nb[3]]-c[index]) / top_h[1];
                cz[index] = (c[index_nb[5]]+c[index_nb[12]]-c[index_nb[4]]-c[index_nb[11]]) / (4.0*top_h[2]);

                vel_u[index] = (u[index_nb[0]]+u[index_nb[9]]+u[index_nb[3]]+u[index])/4;
                vel_v[index] = v[index];
                vel_w[index] = (w[index_nb[4]]+w[index_nb[11]]+w[index_nb[3]]+w[index])/4;
            }

            //w-face
            else if (coord==2) {
                ux = (u[index_nb[5]]+u[index]-u[index_nb[0]]-u[index_nb[17]]) / (2.0*top_h[0]);
                index1 = d_index3d(i-1,j+1,k+1,top_gmax);
                index2 = d_index3d(i-1,j-1,k+1,top_gmax);
                uy = (u[index_nb[3]]+u[index_nb[9]]+u[index_nb[12]]+u[index1]
                     -u[index_nb[2]]-u[index_nb[6]]-u[index_nb[13]]-u[index2]) / (8.0*top_h[1]);
                uz = (u[index_nb[5]]+u[index_nb[17]]-u[index_nb[0]]-u[index]) / (2.0*top_h[2]);

                index1 = d_index3d(i+1,j-1,k+1,top_gmax);
                index2 = d_index3d(i-1,j-1,k+1,top_gmax);
                vx = (v[index_nb[1]]+v[index_nb[7]]+v[index_nb[16]]+v[index1]
                     -v[index_nb[0]]-v[index_nb[6]]-v[index_nb[17]]-v[index2]) / (8.0*top_h[0]);
                vy = (v[index_nb[5]]+v[index]-v[index_nb[2]]-v[index_nb[13]]) / (2.0*top_h[1]);
                vz = (v[index_nb[5]]+v[index_nb[13]]-v[index_nb[2]]-v[index]) / (2.0*top_h[2]);

                wx = (w[index_nb[1]]-w[index_nb[0]]) / (2.0*top_h[0]);
                wy = (w[index_nb[3]]-w[index_nb[2]]) / (2.0*top_h[1]);
                wz = (w[index_nb[5]]-w[index_nb[4]]) / (2.0*top_h[2]);

                cx[index] = (c[index_nb[1]]+c[index_nb[16]]-c[index_nb[0]]-c[index_nb[17]]) / (4.0*top_h[0]);
                cy[index] = (c[index_nb[3]]+c[index_nb[12]]-c[index_nb[2]]-c[index_nb[13]]) / (4.0*top_h[1]);
                cz[index] = (c[index_nb[5]]-c[index]) / top_h[2];

                vel_u[index] = (u[index_nb[0]]+u[index_nb[17]]+u[index_nb[5]]+u[index])/4;
                vel_v[index] = (v[index_nb[5]]+v[index_nb[13]]+v[index_nb[2]]+v[index])/4;
                vel_w[index] = w[index];
            }

            //cell center
            else if (coord==3) {
                ux = (u[index]-u[index_nb[0]]) / top_h[0];
                uy = (u[index_nb[3]]+u[index_nb[9]]-u[index_nb[2]]-u[index_nb[6]]) / (2.0*top_h[1]);
                uz = (u[index_nb[5]]+u[index_nb[17]]-u[index_nb[4]]-u[index_nb[14]]) / (2.0*top_h[2]);

                vx = (v[index_nb[1]]+v[index_nb[7]]-v[index_nb[0]]-v[index_nb[6]]) / (2.0*top_h[0]);
                vy = (v[index]-v[index_nb[2]]) / top_h[1];
                vz = (v[index_nb[5]]+v[index_nb[13]]-v[index_nb[4]]-v[index_nb[10]]) / (2.0*top_h[2]);

                wx = (w[index_nb[1]]+w[index_nb[15]]-w[index_nb[0]]-w[index_nb[14]]) / (2.0*top_h[0]);
                wy = (w[index_nb[3]]+w[index_nb[11]]-w[index_nb[2]]-w[index_nb[10]]) / (2.0*top_h[1]);
                wz = (w[index]-w[index_nb[4]]) / top_h[2];

                cx[index] = (c[index_nb[1]]-c[index_nb[0]]) / (2.0*top_h[0]);
                cy[index] = (c[index_nb[3]]-c[index_nb[2]]) / (2.0*top_h[1]);
                cz[index] = (c[index_nb[5]]-c[index_nb[4]]) / (2.0*top_h[2]);

                vel_u[index] = (u[index_nb[0]]+u[index])/2;
                vel_v[index] = (v[index_nb[2]]+v[index])/2;
                vel_w[index] = (w[index_nb[4]]+w[index])/2;
            }

            //S_ij
            s11[index] = ux;
            s12[index] = (uy+vx)/2.0;
            s13[index] = (uz+wx)/2.0;
            s22[index] = vy;
            s23[index] = (vz+wy)/2.0;
            s33[index] = wz;
            s[index] = sqrt(2*( sqr(s11[index]) + sqr(s12[index])
                         + sqr(s13[index]) + sqr(s22[index])
                         + sqr(s23[index]) + sqr(s33[index])
                         + sqr(s12[index]) + sqr(s13[index])
                         + sqr(s23[index]) ));

            //Sa_ij
            s11[index] -= (ux+vy+wz)/3.0;
            s22[index] -= (ux+vy+wz)/3.0;
            s33[index] -= (ux+vy+wz)/3.0;
        }

        for (k = 0; k <= ((kmax-kmin+1)/NB)-1; k++)
        for (j = 0; j <= ((jmax-jmin+1)/NB)-1; j++)
        for (i = 0; i <= ((imax-imin+1)/NB)-1; i++)
        {
            kk = (NB*k)+kmin;
            jj = (NB*j)+jmin;
            ii = (NB*i)+imin;

            sum_p = sum_pu = sum_pv = sum_pw = sum_pc = 0.0;
            sum_puu = sum_pvv = sum_pww = 0.0;
            sum_puv = sum_pvw = sum_puw = 0.0;
            sum_puc = sum_pvc = sum_pwc = 0.0;
            sum_s11 = sum_s12 = sum_s13 = sum_s22 = sum_s23 = sum_s33 = 0.0;
            sum_pss11 = sum_pss12 = sum_pss13 = sum_pss22 = sum_pss23 = sum_pss33 = 0.0;
            sum_pscx = sum_pscy = sum_pscz = sum_cx = sum_cy = sum_cz = 0.0;
            sum_s = 0.0;
            for (kkk = kk; kkk < kk+NB; kkk++)
            for (jjj = jj; jjj < jj+NB; jjj++)
            for (iii = ii; iii < ii+NB; iii++)
            {
                index = d_index3d(iii,jjj,kkk,top_gmax);
                if (coord==0)		index_nbb = d_index3d(iii+1,jjj,kkk,top_gmax);
                else if (coord==1)	index_nbb = d_index3d(iii,jjj+1,kkk,top_gmax);
                else if (coord==2)	index_nbb = d_index3d(iii,jjj,kkk+1,top_gmax);

                if (coord < 3) { //3 cell faces
                    rho_face = 2.0/(1.0/rho[index] + 1.0/rho[index_nbb]);
                    c_face = (c[index] + c[index_nbb])/2;
                }
                else { //cell center
                    rho_face = rho[index];
                    c_face = c[index];
                }

                sum_p += rho_face;
                sum_pu += rho_face*vel_u[index];
                sum_pv += rho_face*vel_v[index];
                sum_pw += rho_face*vel_w[index];
                sum_pc += rho_face*c_face;
                sum_puu += rho_face*vel_u[index]*vel_u[index];
                sum_pvv += rho_face*vel_v[index]*vel_v[index];
                sum_pww += rho_face*vel_w[index]*vel_w[index];
                sum_puv += rho_face*vel_u[index]*vel_v[index];
                sum_puw += rho_face*vel_u[index]*vel_w[index];
                sum_pvw += rho_face*vel_v[index]*vel_w[index];
                sum_puc += rho_face*vel_u[index]*c_face;
                sum_pvc += rho_face*vel_v[index]*c_face;
                sum_pwc += rho_face*vel_w[index]*c_face;
                sum_s11 += s11[index];
                sum_s12 += s12[index];
                sum_s13 += s13[index];
                sum_s22 += s22[index];
                sum_s23 += s23[index];
                sum_s33 += s33[index];
                sum_pss11 += rho_face*s[index]*s11[index];
                sum_pss12 += rho_face*s[index]*s12[index];
                sum_pss13 += rho_face*s[index]*s13[index];
                sum_pss22 += rho_face*s[index]*s22[index];
                sum_pss23 += rho_face*s[index]*s23[index];
                sum_pss33 += rho_face*s[index]*s33[index];
                sum_pscx += rho_face*s[index]*cx[index];
                sum_pscy += rho_face*s[index]*cy[index];
                sum_pscz += rho_face*s[index]*cz[index];
                sum_cx += cx[index];
                sum_cy += cy[index];
                sum_cz += cz[index];
                sum_s += s[index];
            }
            //TODO: ma_ij differs from Eq.(37) of Tingguang Ma's PhD thesis in sum_s/(NBC), which is an approx.
            ma11 = (2.0*delta2*(sum_pss11/(NBC)))
                    - (2.0*tdelta2*(sum_p/(NBC))*(sum_s/(NBC))*(sum_s11/(NBC)));
            ma12 = (2.0*delta2*(sum_pss12/(NBC)))
                    - (2.0*tdelta2*(sum_p/(NBC))*(sum_s/(NBC))*(sum_s12/(NBC)));
            ma13 = (2.0*delta2*(sum_pss13/(NBC)))
                    - (2.0*tdelta2*(sum_p/(NBC))*(sum_s/(NBC))*(sum_s13/(NBC)));
            ma22 = (2.0*delta2*(sum_pss22/(NBC)))
                    - (2.0*tdelta2*(sum_p/(NBC))*(sum_s/(NBC))*(sum_s22/(NBC)));
            ma23 = (2.0*delta2*(sum_pss23/(NBC)))
                    - (2.0*tdelta2*(sum_p/(NBC))*(sum_s/(NBC))*(sum_s23/(NBC)));
            ma33 = (2.0*delta2*(sum_pss33/(NBC)))
                    - (2.0*tdelta2*(sum_p/(NBC))*(sum_s/(NBC))*(sum_s33/(NBC)));
            la11 = (sum_puu/(NBC))-((sum_pu/(NBC))*(sum_pu/(NBC))/(sum_p/(NBC)));
            la12 = (sum_puv/(NBC))-((sum_pu/(NBC))*(sum_pv/(NBC))/(sum_p/(NBC)));
            la13 = (sum_puw/(NBC))-((sum_pu/(NBC))*(sum_pw/(NBC))/(sum_p/(NBC)));
            la22 = (sum_pvv/(NBC))-((sum_pv/(NBC))*(sum_pv/(NBC))/(sum_p/(NBC)));
            la23 = (sum_pvw/(NBC))-((sum_pv/(NBC))*(sum_pw/(NBC))/(sum_p/(NBC)));
            la33 = (sum_pww/(NBC))-((sum_pw/(NBC))*(sum_pw/(NBC))/(sum_p/(NBC)));

//            printf("trace of L_ij = %20.16g\n", la11+la22+la33);

            //suppose k=1
            //TODO: my_1i differs from Eq.(48) of Tingguang Ma's PhD thesis in sum_s/(NBC) and sum_cxi/(NBC), which are approx.
            my1 = (delta2*(sum_pscx/(NBC)))
                    - (tdelta2*(sum_p/(NBC))*(sum_s/(NBC))*(sum_cx/(NBC)));
            my2 = (delta2*(sum_pscy/(NBC)))
                    - (tdelta2*(sum_p/(NBC))*(sum_s/(NBC))*(sum_cy/(NBC)));
            my3 = (delta2*(sum_pscz/(NBC)))
                    - (tdelta2*(sum_p/(NBC))*(sum_s/(NBC))*(sum_cz/(NBC)));
            ly1 = (sum_puc/(NBC))-((sum_pu/(NBC))*(sum_pc/(NBC))/(sum_p/(NBC)));
            ly2 = (sum_pvc/(NBC))-((sum_pv/(NBC))*(sum_pc/(NBC))/(sum_p/(NBC)));
            ly3 = (sum_pwc/(NBC))-((sum_pw/(NBC))*(sum_pc/(NBC))/(sum_p/(NBC)));

            ldeno[i][j][k] = ((ma11*ma11) + (ma12*ma12)
                            + (ma13*ma13) + (ma22*ma22)
                            + (ma23*ma23));
            lnume[i][j][k] = ((la11*ma11) + (la12*ma12)
                            + (la13*ma13) + (la22*ma22)
                            + (la23*ma23));
            cnume[i][j][k] = ((my1*my1) + (my2*my2) + (my3*my3));
            cdeno[i][j][k] = ((ly1*my1) + (ly2*my2) + (ly3*my3));
        }

        max_cs[coord] = 0;
        max_cc[coord] = 0;
        max_turbulent_kinematic_viscosity[coord] = 0;
        max_turbulent_diffusity[coord] = 0;
        avg_turbulent_kinematic_viscosity[coord] = 0;
        avg_turbulent_diffusity[coord] = 0;
	count = 0;

        for (k = 0; k <= ((kmax-kmin+1)/NB)-1; k++)
        for (j = 0; j <= ((jmax-jmin+1)/NB)-1; j++)
        for (i = 0; i <= ((imax-imin+1)/NB)-1; i++)
        {
            kk = (NB*k)+kmin;
            jj = (NB*j)+jmin;
            ii = (NB*i)+imin;
            cs = lnume[i][j][k]/ldeno[i][j][k];
            cc = cdeno[i][j][k]/cnume[i][j][k];

            if (cs < 0 || fabs(ldeno[i][j][k]) < 1e-16)	cs = 0;
            if (fabs(cnume[i][j][k]) < 1e-16)		cc = 0;
            if (cs > max_cs[coord])                     max_cs[coord] = cs;
            if (cc > max_cc[coord])                     max_cc[coord] = cc;


            for (kkk = kk; kkk < kk+NB; kkk++)
            for (jjj = jj; jjj < jj+NB; jjj++)
            for (iii = ii; iii < ii+NB; iii++)
            {
                index = d_index3d(iii,jjj,kkk,top_gmax);
                if (coord==0)           index_nbb = d_index3d(iii+1,jjj,kkk,top_gmax);
                else if (coord==1)      index_nbb = d_index3d(iii,jjj+1,kkk,top_gmax);
                else if (coord==2)	index_nbb = d_index3d(iii,jjj,kkk+1,top_gmax);

                if (coord < 3) { //3 cell faces
                    rho_face = 2.0/(1.0/rho[index] + 1.0/rho[index_nbb]);
                }
                else { //cell center
                    rho_face = rho[index];
                }

                //(rho*nu_t)
                cell_center[index].m_state.m_mu_turbulent[coord] = rho_face*cs*delta2*s[index];
                //(nu_t/Sc_t)
                if (cc*s[index] < 0)	cc = 0;
                cell_center[index].m_state.m_Dcoef_turbulent[coord] = cc*delta2*s[index];

                if (cs*delta2*s[index] > max_turbulent_kinematic_viscosity[coord])
                    max_turbulent_kinematic_viscosity[coord] = cs*delta2*s[index];

                if (cc*delta2*s[index] > max_turbulent_diffusity[coord])
                    max_turbulent_diffusity[coord] = cc*delta2*s[index];

		// calculate the average nu_t and Diff_t in the mixing layer
                double z_coord = cell_center[index].m_coords[2];
		if (z_coord >= zmin_vf && z_coord <= zmax_vf)
		{
		    count++;
		    avg_turbulent_kinematic_viscosity[coord] += cs*delta2*s[index];
		    avg_turbulent_diffusity[coord] += cc*delta2*s[index];
		}
            }
        }

        if (pp_numnodes() > 1)
	{
            pp_global_max(&max_cs[coord],1);
            pp_global_max(&max_cc[coord],1);
            pp_global_max(&max_turbulent_kinematic_viscosity[coord],1);
            pp_global_max(&max_turbulent_diffusity[coord],1);

	    pp_global_isum(&count,1);
            pp_global_sum(&avg_turbulent_kinematic_viscosity[coord],1);
            pp_global_sum(&avg_turbulent_diffusity[coord],1);
	}

	if (count > 0)
	{
	    avg_turbulent_kinematic_viscosity[coord] /= count;
	    avg_turbulent_diffusity[coord] /= count;
	}

        if (coord == 0)
	{
	    printf("\nIn computeSubgridModel_vd(): \n");
            printf("zmin_vf = %lf, zmax_vf = %lf, count = %d\n", zmin_vf ,zmax_vf, count);
	}
        printf("max cs[%d] = %20.16g\n", coord, max_cs[coord]);
        printf("max cc[%d] = %20.16g\n", coord, max_cc[coord]);
        printf("max nu_t[%d] = %20.16g\n", coord, max_turbulent_kinematic_viscosity[coord]);
        printf("max Diff_t[%d] = %20.16g\n", coord, max_turbulent_diffusity[coord]);
        printf("avg nu_t[%d] = %20.16g\n", coord, avg_turbulent_kinematic_viscosity[coord]);
        printf("avg Diff_t[%d] = %20.16g\n", coord, avg_turbulent_diffusity[coord]);
    }


    //print output file
    char filename[256];
    //sprintf(filename,"%s/sim_hc_sgs.txt",out_name);
    sprintf(filename,"sim_hc_sgs.txt");
    FILE *outfile;
    double A, g, H, tau;
    A = fabs(m_rho[1]-m_rho[0])/(m_rho[1]+m_rho[0]);
    g = fabs(iFparams->gravity[2]);
    H = 32;
    tau = sqrt(A*g/H)*front->time;
    if (pp_mynode() == 0)
    {
        if (front->step == 0)
        {
            outfile = fopen(filename,"w");
            fprintf(outfile,"  ts     time         tau      max_turb_nu  avg_turb_nu       nu       max_turb_D   avg_turb_D       D\n");
        }
        else
        {
            outfile = fopen(filename,"a");
        }
        fprintf(outfile,"%4d %e %e %e %e %e %e %e %e\n", front->step,front->time,tau,max_turbulent_kinematic_viscosity[3],avg_turbulent_kinematic_viscosity[3],(m_mu[0]+m_mu[1])/(m_rho[0]+m_rho[1]),max_turbulent_diffusity[3],avg_turbulent_diffusity[3],m_Dcoef[0]);

        fclose(outfile);
    }


    //scatter states
    for (coord = 0; coord < MAXD+1; coord++)
    {
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_mu_turbulent[coord];
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_mu_turbulent[coord] = array[index];
        }

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            array[index] = cell_center[index].m_state.m_Dcoef_turbulent[coord];
        }
        scatMeshArray();
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            cell_center[index].m_state.m_Dcoef_turbulent[coord] = array[index];
        }
    }

    FT_FreeThese(4,ldeno,lnume,cdeno,cnume);
    FT_FreeThese(8,u,v,w,rho,c,cx,cy,cz);
    FT_FreeThese(3,vel_u,vel_v,vel_w);
    FT_FreeThese(7,s,s11,s12,s13,s22,s23,s33);
}       /* end computeSubgridModel_vd */

/*
 *
 *      Smeeton Youngs Experiment 105
 *      with Reflecting Boundary Condition
 *      with Contact Angle
 *      with Meniscus
 *      NOT IN USE FOR NOW.
 *
 * */


//FUNC REFLECTION BOUNDARY CONDITION REVISED! Free-Slip Wall Condition
void Incompress_Solver_Smooth_3D_Cartesian::ReflectBC(
        int dir,
        int side,
        double *solute)
{
    int i, j, k;
    int InIndex, BIndex, OnIndex;
    //printf("HZ in func %s lbuf = [%d %d %d] ubuf = [%d %d %d]\n", __func__,lbuf[0], lbuf[1], lbuf[2], ubuf[0],ubuf[1], ubuf[2]);

    // This is where the code takes care of Reflecting Boundary Condition.
    if (side == 0)
    {
        switch (dir)
        {
            case 0: // X Lower
                 for (j = 0; j <= jmax+ubuf[1]; j++)
                     for (k = 0; k <= kmax+ubuf[2]; k++)
                     {
                         for (i = imin; i < imin+lbuf[0]-1; i++)
                         {
                              InIndex = d_index3d(i,j,k,top_gmax);//This is Interior State index
                              BIndex = d_index3d(lbuf[0]-2-i+imin,j,k,top_gmax);//This is Boundary State index
                              solute[BIndex] = -solute[InIndex];
                         }
                         OnIndex = d_index3d(imin-1,j,k,top_gmax);
                         solute[OnIndex] = 0.0; // On wall normal velocity component vanish
                     }
                 break;
            case 1: // Y Lower
                 for (i = 0; i <= imax+ubuf[0]; i++)
                     for (k = 0; k <= kmax+ubuf[2]; k++)
                     {
                         for (j = jmin; j < jmin+lbuf[1]-1; j++)
                         {
                              InIndex = d_index3d(i,j,k,top_gmax);//This is Interior State index
                              BIndex = d_index3d(i,lbuf[1]-2-j+jmin,k,top_gmax);//This is Boundary State index
                              solute[BIndex] = -solute[InIndex];
                         }
                         OnIndex = d_index3d(i,jmin-1,k,top_gmax);
                         solute[OnIndex] = 0.0; // On wall normal velocity component vanish
                     }
                 break;
            case 2: // Z Lower
                 for (i = 0; i <= imax+ubuf[0]; i++)
                     for (j = 0; j <= jmax+ubuf[1]; j++)
                     {
                         for (k = kmin; k < kmin+lbuf[2]-1; k++)
                         {
                              InIndex = d_index3d(i,j,k,top_gmax);//This is Interior State index
                              BIndex = d_index3d(i,j,lbuf[2]-2-k+kmin,top_gmax);//This is Boundary State index
                              solute[BIndex] = -solute[InIndex];
                         }
                         OnIndex = d_index3d(i,j,kmin-1,top_gmax);
                         solute[OnIndex] = 0.0; // On wall normal velocity component vanish
                     }
                 break;
        }
    }
    else
    {
        switch (dir)
        {
            case 0: // X Upper
                 for (j = 0; j <= jmax+ubuf[1]; j++)
                     for (k = 0; k <= kmax+ubuf[2]; k++)
                     {
                         for (i = imax-ubuf[0]+1; i <= imax-1; i++)// X Upper
                         {
                             InIndex = d_index3d(i,j,k,top_gmax);//This is Interior State index
                             BIndex = d_index3d(2*imax-i,j,k,top_gmax);//This is Boundary State index
                             solute[BIndex] = -solute[InIndex];

                         }
                         OnIndex = d_index3d(imax,j,k,top_gmax);
                         solute[OnIndex] = 0.0; // On wall normal velocity component vanish
                     }
                 break;
            case 1: // Y Upper
                 for(i = 0; i <= imax+ubuf[0]; i++)
                     for (k = 0; k <= kmax+ubuf[2]; k++)
                     {
                         for (j = jmax-ubuf[1]+1; j <= jmax-1; j++)//Y Upper
                         {
                             InIndex = d_index3d(i,j,k,top_gmax);//This is Interior State index
                             BIndex = d_index3d(i,2*jmax-j,k,top_gmax);//This is Boundary State index
                             solute[BIndex] = -solute[InIndex];
                         }
                         OnIndex = d_index3d(i,jmax,k,top_gmax);
                         solute[OnIndex] = 0.0; // On wall normal velocity component vanish
                     }
                 break;
            case 2: // Z Upper
                 for (i = 0; i <= imax+ubuf[0]; i++)
                     for (j = 0; j <= jmax+ubuf[1]; j++)
                     {
                         for (k = kmax-ubuf[2]+1; k <= kmax-1; k++)//Z Upper
                         {
                             InIndex = d_index3d(i,j,k,top_gmax);//This is Interior State index
                             BIndex = d_index3d(i,j,2*kmax-k,top_gmax);//This is Boundary State index
                             solute[BIndex] = -solute[InIndex];
                         }
                         OnIndex = d_index3d(i,j,kmax,top_gmax);
                         solute[OnIndex] = 0.0; // On wall normal velocity component vanish
                     }
                 break;
        }
    }
}

void Incompress_Solver_Smooth_3D_Cartesian::Solute_Reflect(
        int dir,
        double *solute)
{
     int side, i;
     INTERFACE *intfc = front->grid_intfc;

     for (side = 0; side < 2; side++)
     {
         if (rect_boundary_type(intfc,dir,side) == REFLECTION_BOUNDARY)
         {
              ReflectBC(dir, side, solute);
         }
     }
}

//Reflection Boundary Condition Treatment
void Incompress_Solver_Smooth_3D_Cartesian::enforceReflectionState(double **vel)
{
     int    dir;

     for (dir = 0; dir < dim; dir++)
     {
         Solute_Reflect(dir,vel[dir]);
    }
}

// Convert GRID_DIRECTION into dir and side.
// A helper function
extern void convertGridDirectionToDirSide(GRID_DIRECTION direct, int *dir, int *side)
{
     if (direct == WEST)
     {
         *dir = 0;
         *side = 0;
     }
     if (direct == EAST)
     {
          *dir = 0;
          *side = 1;
     }
     if (direct == SOUTH)
     {
         *dir = 1;
         *side = 0;
     }
     if (direct == NORTH)
     {
          *dir = 1;
          *side = 1;
     }
     if (direct == LOWER)
     {
         *dir = 2;
         *side = 0;
     }
     if (direct == UPPER)
     {
          *dir = 2;
          *side = 1;
     }
}
/**
 * Doxygen compatible:
 *
 * this function check
 * PERIODIC_BOUNDARY/INTERIOR_STATE  = 0
 * DIRICHLET_BOUNDARY = 1
 * NEUMANN_BOUNDARY = 2
 * REFLECTION_BOUNDARY = 3
 *
 * logic setup: if statement (!0), means a boundary occurred.
 * if (>=2) DERIVATIVE BOUNDARY(NEUMANN, REFLECTION)
 * if (1) SIMPLE BOUNDARY(DIRICHLET)
 *
 *
 * if (!0)
 * {
 *      if (1)
 *          // DIRICHLET TREATMENT
 *      if (2)
 *          // NEUMANN TREATMENT
 *      if (3)
 *          // REFLECTION TREATMENT
 * }
 * */
void Incompress_Solver_Smooth_3D_Cartesian::checkBoundaryCondition(GRID_DIRECTION dir, int *icoords, int *bNoBoundary, double m_t_new, COMPONENT comp)
{
        POINTER intfc_state;
        HYPER_SURF *hs;
        INTERFACE *intfc = front->interf;
        int dirr, side;
        double crx_coords[MAXD];

        convertGridDirectionToDirSide(dir, &dirr, &side);
        if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir,
                comp,&intfc_state,&hs,crx_coords,m_t_new) &&
                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
        {
            //examine NEUMANN Boundary or DIRICHLET Boundary
            if (wave_type(hs) == NEUMANN_BOUNDARY)
            {
                *bNoBoundary = 2;
            }
            if (wave_type(hs) == DIRICHLET_BOUNDARY)
            {
                *bNoBoundary = 1;
                printf("DIRICHLET BOUNDARY was NOT implemented in func %s\n", __func__);
                clean_up(ERROR);
            }
        }
        else if (rect_boundary_type(intfc,dirr, side) == REFLECTION_BOUNDARY && FT_Reflect(icoords, dirr, side))
            *bNoBoundary = 3;
        else // This is for either PERIODIC BC or simply INTERIOR CELL
            *bNoBoundary = 0;
}

// REFLECTION BOUNDARY CONDITION check if icoords is on the WALL.
bool Incompress_Solver_Smooth_3D_Cartesian::FT_Reflect(int *icoords, int dir, int side)
{
    int wallIndex[3][2] = {{imin,imax},{jmin,jmax},{kmin,kmax}};

    if (icoords[dir] == wallIndex[dir][side])
        return true;
    return false;
}
