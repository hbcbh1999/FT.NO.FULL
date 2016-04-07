/*******************************************************************
 * 	               iFcartsn2d.cpp	
 *******************************************************************/
#include "iFluid.h"
#include "solver.h"

void Incompress_Solver_Smooth_2D_Cartesian::computeNewDensity_vd(int flag)
{
        int i, j, k, index;
        COMPONENT comp;
        double density;

        max_density = -1;

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index2d(i,j,top_gmax);
            comp = top_comp[index];
            if (!ifluid_comp(comp))
            {
                cell_center[index].m_state.m_rho = 0.0;
                cell_center[index].m_state.m_rho_old = 0.0;
                continue;
            }

            if (!flag)
            {
                cell_center[index].m_state.m_rho_old = cell_center[index].m_state.m_rho;
                cell_center[index].m_state.m_rho -= m_dt*cell_center[index].m_state.m_rho_adv;
            }
            else
                cell_center[index].m_state.m_rho = cell_center[index].m_state.m_rho_old - 
                                                   m_dt*cell_center[index].m_state.m_rho_adv;

            density = fabs(cell_center[index].m_state.m_rho);
            if (density > max_density)
                max_density = density;
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            array[index] = cell_center[index].m_state.m_rho_old;
        }
            scatMeshArray();
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            cell_center[index].m_state.m_rho_old = array[index];
        }
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            array[index] = cell_center[index].m_state.m_rho;
        }
        scatMeshArray();
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            cell_center[index].m_state.m_rho = array[index];
        }
        pp_global_max(&max_density,1);
}       /* end computeNewDensity_vd */


void Incompress_Solver_Smooth_2D_Cartesian::compDiffWithSmoothProperty_velocity_vd(void)
{
        COMPONENT comp;
        int index,index_nb[8],size;
        int I,I_nb[8];
        double coords[MAXD],corner_coords[MAXD],crx_coords[MAXD];
        double coeff0[4],coeff1[4],coeff_temp1,coeff_temp2;
        double mu[4],mu_edge[4],mu0,rho_mid,rhs;
        double U0_nb[8],U1_nb[8],U0_center, U1_center;
        double U0_nb_new[4],U1_nb_new[4];
        int flag[4];  //denote whether this is dirichlet or neumann boundary
        L_STATE state,corner_state,corner_state_new;
        int i,j,k,nb,icoords[MAXD];
        INTERFACE *intfc = front->interf;
        double speed;
        double *x;
        GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
        double (*getStateVel[2])(POINTER) = {getStateXvel,getStateYvel};
        POINTER intfc_state;
        HYPER_SURF *hs;
        PetscInt num_iter;
        double rel_residual;

        max_speed = 0.0;
        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,2*size,sizeof(double));

        PETSc solver;

        solver.Create(2*ilower, 2*iupper-1, 15, 15);
                // two buffer zones, 5 u and 9 v for the first  equation
                //                   5 v and 9 u for the second equation

        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ij_to_I[i][j];
            if (I == -1) continue; 

            index  = d_index2d(i,j,top_gmax);
            index_nb[0] = d_index2d(i-1,j,top_gmax);
            index_nb[1] = d_index2d(i+1,j,top_gmax);
            index_nb[2] = d_index2d(i,j-1,top_gmax);
            index_nb[3] = d_index2d(i,j+1,top_gmax);
            index_nb[4] = d_index2d(i-1,j-1,top_gmax);
            index_nb[5] = d_index2d(i+1,j-1,top_gmax);
            index_nb[6] = d_index2d(i+1,j+1,top_gmax);
            index_nb[7] = d_index2d(i-1,j+1,top_gmax);

            /*
             *       |          |
             *   7   |    3     |   6
             *-------|----------|-------
             *       |          |
             *   0   |          |   1
             *       |          |
             *-------|----------|--------
             *   4   |    2     |   5
             *       |          |
             */

            I_nb[0] = ij_to_I[i-1][j];
            I_nb[1] = ij_to_I[i+1][j];
            I_nb[2] = ij_to_I[i][j-1];
            I_nb[3] = ij_to_I[i][j+1];
            I_nb[4] = ij_to_I[i-1][j-1];
            I_nb[5] = ij_to_I[i+1][j-1];
            I_nb[6] = ij_to_I[i+1][j+1];
            I_nb[7] = ij_to_I[i-1][j+1];

            icoords[0] = i;
            icoords[1] = j;
            comp = top_comp[index];

            mu0 = cell_center[index].m_state.m_mu;
            rho_mid = 0.5*(cell_center[index].m_state.m_rho
                      + cell_center[index].m_state.m_rho_old);
            U0_center = cell_center[index].m_state.m_U[0];
            U1_center = cell_center[index].m_state.m_U[1];

            for (nb = 0; nb < 4; nb++)
            {
                if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                {
                    flag[nb] = 1;
                    // old boundary condition
                    U0_nb[nb] = getStateVel[0](intfc_state);
                    U1_nb[nb] = getStateVel[1](intfc_state);
                    // new boundary condition
                    FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords,m_t_new);
                    U0_nb_new[nb] = getStateVel[0](intfc_state);
                    U1_nb_new[nb] = getStateVel[1](intfc_state);

                    if (wave_type(hs) == DIRICHLET_BOUNDARY ||
                        wave_type(hs) == NEUMANN_BOUNDARY)
                    {
                        mu_edge[nb] = mu0;
                        mu[nb] = mu0;
                    }
                    else
                    {
                        mu_edge[nb] = 0.5*(mu0 +
                                cell_center[index_nb[nb]].m_state.m_mu);
                        mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
                    }
                }
                else
                {
                    flag[nb] = 0;
                    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
                    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
                    mu_edge[nb] = 0.5*(mu0 +
                                cell_center[index_nb[nb]].m_state.m_mu);
                    mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
                }
            }

            for (nb = 4; nb < 8; nb++) // corner cell value, interior 
            {
                if (I_nb[nb] != -1)
                {
                    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
                    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
                }
            }

            // source term
            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);

            //Setting the coefficients and matrix for the U0 in first equation
            coeff0[0] = 2.0*m_dt/3.0/rho_mid * mu_edge[0]/(top_h[0]*top_h[0]);  // left
            coeff0[1] = 2.0*m_dt/3.0/rho_mid * mu_edge[1]/(top_h[0]*top_h[0]);  // right
            coeff0[2] = 0.5*m_dt/rho_mid * mu_edge[2]/(top_h[1]*top_h[1]);      // down 
            coeff0[3] = 0.5*m_dt/rho_mid * mu_edge[3]/(top_h[1]*top_h[1]);      // up


            solver.Add_A(I*2,I*2,1.0);
            rhs = U0_center;

            for (nb = 0; nb < 4; nb++)
            {
                if (I_nb[nb] != -1)
                {
                    solver.Add_A(I*2,I_nb[nb]*2,-coeff0[nb]);
                    rhs += coeff0[nb]*U0_nb[nb];
                }
                else 
                {
                    coeff0[nb] = 2.0*coeff0[nb];
                    rhs += coeff0[nb]*(U0_nb[nb] + U0_nb_new[nb]);
                }

                solver.Add_A(I*2,I*2,coeff0[nb]);
                rhs -= coeff0[nb]*U0_center;
            }

            //set the coefficients and matrix for U1 in first equation
            //traverse the four corners

            //corner 4

            coeff_temp1 = mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho_mid;
            coeff_temp2 = -2.0*mu_edge[0]*m_dt/(top_h[0]*top_h[1])/3.0/rho_mid;
            /*
            if (flag[0] == 1 && flag[2] == 0)
            {
                rhs += coeff_temp1*U1_nb[0];
                rhs += coeff_temp2*U1_nb[0];
            }
            else if(flag[0] == 0 && flag[2] ==1)
            {
                rhs += coeff_temp1*U1_nb[2];
                rhs += coeff_temp2*U1_nb[2];
            }
            else if(flag[0] ==1 && flag[2] == 1)
            {
                rhs += coeff_temp1*U1_nb[0];
                rhs += coeff_temp2*U1_nb[0];
            }
            else {
                rhs += coeff_temp1*(U1_nb[0]+U1_nb[2]+U1_nb[4]+U1_center)/8.0;
                rhs += coeff_temp2*(U1_nb[0]+U1_nb[2]+U1_nb[4]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,       -coeff_temp1/8.0);
                solver.Add_A(I*2,I_nb[0]*2+1, -coeff_temp1/8.0);
                solver.Add_A(I*2,I_nb[2]*2+1, -coeff_temp1/8.0);
                solver.Add_A(I*2,I_nb[4]*2+1, -coeff_temp1/8.0);

                solver.Add_A(I*2,I*2+1,       -coeff_temp2/8.0);
                solver.Add_A(I*2,I_nb[0]*2+1, -coeff_temp2/8.0);
                solver.Add_A(I*2,I_nb[2]*2+1, -coeff_temp2/8.0);
                solver.Add_A(I*2,I_nb[4]*2+1, -coeff_temp2/8.0);
            } */
            if (flag[0] ==0 && flag[2] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U1_nb[0]+U1_nb[2]+U1_nb[4]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,       -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2,I_nb[0]*2+1, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2,I_nb[2]*2+1, -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2,I_nb[4]*2+1, -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                getExactSolution_vd(corner_coords,m_t_old,corner_state);
                getExactSolution_vd(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*(coeff_temp1+coeff_temp2)*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }

            //corner 5

            coeff_temp1 = -mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho_mid;
            coeff_temp2 = 2.0*mu_edge[1]*m_dt/(top_h[0]*top_h[1])/3.0/rho_mid;
            /*
            if (flag[2] == 1 && flag[1] == 0)
            {
                rhs += coeff_temp1*U1_nb[2];
                rhs += coeff_temp2*U1_nb[2];
            }
            else if(flag[2] == 0 && flag[1] == 1)
            {
                rhs += coeff_temp1*U1_nb[1];
                rhs += coeff_temp2*U1_nb[1];
            }
            else if(flag[2] == 1 && flag[1] == 1)
            {
                rhs += coeff_temp1*U1_nb[2];
                rhs += coeff_temp2*U1_nb[2];
            }
            else {
                rhs += coeff_temp1*(U1_nb[2]+U1_nb[1]+U1_nb[5]+U1_center)/8.0;
                rhs += coeff_temp2*(U1_nb[2]+U1_nb[1]+U1_nb[5]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,        -coeff_temp1/8.0);
                solver.Add_A(I*2,I_nb[2]*2+1,  -coeff_temp1/8.0);
                solver.Add_A(I*2,I_nb[1]*2+1,  -coeff_temp1/8.0);
                solver.Add_A(I*2,I_nb[5]*2+1,  -coeff_temp1/8.0);

                solver.Add_A(I*2,I*2+1,        -coeff_temp2/8.0);
                solver.Add_A(I*2,I_nb[2]*2+1,  -coeff_temp2/8.0);
                solver.Add_A(I*2,I_nb[1]*2+1,  -coeff_temp2/8.0);
                solver.Add_A(I*2,I_nb[5]*2+1,  -coeff_temp2/8.0);
            } */
            if (flag[1] ==0 && flag[2] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U1_nb[2]+U1_nb[1]+U1_nb[5]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2,I_nb[2]*2+1,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2,I_nb[1]*2+1,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2,I_nb[5]*2+1,  -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                getExactSolution_vd(corner_coords,m_t_old,corner_state);
                getExactSolution_vd(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*(coeff_temp1+coeff_temp2)*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }


            //corner 6

            coeff_temp1 = mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho_mid;
            coeff_temp2 = -2.0*mu_edge[1]*m_dt/(top_h[0]*top_h[1])/3.0/rho_mid;
            /*
            if (flag[1] == 1 && flag[3] == 0)
            {
                rhs += coeff_temp1*U1_nb[1];
                rhs += coeff_temp2*U1_nb[1];
            }
            else if(flag[1] == 0 && flag[3] == 1)
            {
                rhs += coeff_temp1*U1_nb[3];
                rhs += coeff_temp2*U1_nb[3];
            }
            else if(flag[1] == 1 && flag[3] == 1)
            {
                rhs += coeff_temp1*U1_nb[1];
                rhs += coeff_temp2*U1_nb[1];
            }
            else {
                rhs += coeff_temp1*(U1_nb[1]+U1_nb[3]+U1_nb[6]+U1_center)/8.0;
                rhs += coeff_temp2*(U1_nb[1]+U1_nb[3]+U1_nb[6]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,        -coeff_temp1/8.0);
                solver.Add_A(I*2,I_nb[1]*2+1,  -coeff_temp1/8.0);
                solver.Add_A(I*2,I_nb[3]*2+1,  -coeff_temp1/8.0);
                solver.Add_A(I*2,I_nb[6]*2+1,  -coeff_temp1/8.0);

                solver.Add_A(I*2,I*2+1,        -coeff_temp2/8.0);
                solver.Add_A(I*2,I_nb[1]*2+1,  -coeff_temp2/8.0);
                solver.Add_A(I*2,I_nb[3]*2+1,  -coeff_temp2/8.0);
                solver.Add_A(I*2,I_nb[6]*2+1,  -coeff_temp2/8.0);
            } */
            if (flag[1] ==0 && flag[3] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U1_nb[1]+U1_nb[3]+U1_nb[6]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2,I_nb[1]*2+1,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2,I_nb[3]*2+1,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2,I_nb[6]*2+1,  -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                getExactSolution_vd(corner_coords,m_t_old,corner_state);
                getExactSolution_vd(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*(coeff_temp1+coeff_temp2)*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }


            //corner 7

            coeff_temp1 = -mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho_mid;
            coeff_temp2 = 2.0*mu_edge[0]*m_dt/(top_h[0]*top_h[1])/3.0/rho_mid;
            /*
            if (flag[3] == 1 && flag[0] == 0)
            {
                rhs += coeff_temp1*U1_nb[3];
                rhs += coeff_temp2*U1_nb[3];
            }
            else if(flag[3] == 0 && flag[0] == 1)
            {
                rhs += coeff_temp1*U1_nb[0];
                rhs += coeff_temp2*U1_nb[0];
            }
            else if(flag[3] == 1 && flag[0] == 1)
            {
                rhs += coeff_temp1*U1_nb[3];
                rhs += coeff_temp2*U1_nb[3];
            }
            else {
                rhs += coeff_temp1*(U1_nb[3]+U1_nb[0]+U1_nb[7]+U1_center)/8.0;
                rhs += coeff_temp2*(U1_nb[3]+U1_nb[0]+U1_nb[7]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,        -coeff_temp1/8.0);
                solver.Add_A(I*2,I_nb[3]*2+1,  -coeff_temp1/8.0);
                solver.Add_A(I*2,I_nb[0]*2+1,  -coeff_temp1/8.0);
                solver.Add_A(I*2,I_nb[7]*2+1,  -coeff_temp1/8.0);

                solver.Add_A(I*2,I*2+1,        -coeff_temp2/8.0);
                solver.Add_A(I*2,I_nb[3]*2+1,  -coeff_temp2/8.0);
                solver.Add_A(I*2,I_nb[0]*2+1,  -coeff_temp2/8.0);
                solver.Add_A(I*2,I_nb[7]*2+1,  -coeff_temp2/8.0);
            } */
            if (flag[3] ==0 && flag[0] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U1_nb[3]+U1_nb[0]+U1_nb[7]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2,I_nb[3]*2+1,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2,I_nb[0]*2+1,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2,I_nb[7]*2+1,  -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                getExactSolution_vd(corner_coords,m_t_old,corner_state);
                getExactSolution_vd(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*(coeff_temp1+coeff_temp2)*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }

            rhs += m_dt*state.m_U[0]/rho_mid;                            // source term
            rhs += m_dt*cell_center[index].m_state.f_surf[0];            // surface tension
            rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho_mid;    // grad_p
            rhs -= m_dt*cell_center[index].m_state.m_adv[0];             // advection term

            solver.Add_b(I*2, rhs);


            /****************************************************************/
            //Setting the coefficients for U1 in the second equation

            coeff1[0] = 0.5*m_dt/rho_mid * mu_edge[0]/(top_h[0]*top_h[0]);      // left
            coeff1[1] = 0.5*m_dt/rho_mid * mu_edge[1]/(top_h[0]*top_h[0]);      // right
            coeff1[2] = 2.0*m_dt/3.0/rho_mid * mu_edge[2]/(top_h[1]*top_h[1]);  // down 
            coeff1[3] = 2.0*m_dt/3.0/rho_mid * mu_edge[3]/(top_h[1]*top_h[1]);  // up


            solver.Add_A(I*2+1,I*2+1,1.0);
            rhs = U1_center;

            for (nb = 0; nb < 4; nb++)
            {
                if (I_nb[nb] != -1)
                {
                    solver.Add_A(I*2+1,I_nb[nb]*2+1,-coeff1[nb]);
                    rhs += coeff1[nb]*U1_nb[nb];
                }
                else
                {
                    coeff1[nb] = 2.0*coeff1[nb];
                    rhs += coeff1[nb]*(U1_nb[nb] + U1_nb_new[nb]);
                }

                solver.Add_A(I*2+1,I*2+1,coeff1[nb]);
                rhs -= coeff1[nb]*U1_center;
            }

            //set the coefficients and matrix for U0 in second equation
            //traverse the four corners

            //corner 4

            coeff_temp1 = mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho_mid;
            coeff_temp2 = -2.0*mu_edge[2]*m_dt/(top_h[0]*top_h[1])/3.0/rho_mid;
            /*
            if (flag[0] == 1 && flag[2] == 0)
            {
                rhs += coeff_temp1*U0_nb[0];
                rhs += coeff_temp2*U0_nb[0];
            }
            else if(flag[0] == 0 && flag[2] ==1)
            {
                rhs += coeff_temp1*U0_nb[2];
                rhs += coeff_temp2*U0_nb[2];
            }
            else if(flag[0] ==1 && flag[2] == 1)
            {
                rhs += coeff_temp1*U0_nb[0];
                rhs += coeff_temp2*U0_nb[0];
            }
            else {
                rhs += coeff_temp1*(U0_nb[0]+U0_nb[2]+U0_nb[4]+U0_center)/8.0;
                rhs += coeff_temp2*(U0_nb[0]+U0_nb[2]+U0_nb[4]+U0_center)/8.0;

                solver.Add_A(I*2+1,I*2,       -coeff_temp1/8.0);
                solver.Add_A(I*2+1,I_nb[0]*2, -coeff_temp1/8.0);
                solver.Add_A(I*2+1,I_nb[2]*2, -coeff_temp1/8.0);
                solver.Add_A(I*2+1,I_nb[4]*2, -coeff_temp1/8.0);

                solver.Add_A(I*2+1,I*2,       -coeff_temp2/8.0);
                solver.Add_A(I*2+1,I_nb[0]*2, -coeff_temp2/8.0);
                solver.Add_A(I*2+1,I_nb[2]*2, -coeff_temp2/8.0);
                solver.Add_A(I*2+1,I_nb[4]*2, -coeff_temp2/8.0);
            } */
            if (flag[0] ==0 && flag[2] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U0_nb[0]+U0_nb[2]+U0_nb[4]+U0_center)/8.0;

                solver.Add_A(I*2+1,I*2,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2+1,I_nb[0]*2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2+1,I_nb[2]*2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2+1,I_nb[4]*2,  -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                getExactSolution_vd(corner_coords,m_t_old,corner_state);
                getExactSolution_vd(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*(coeff_temp1+coeff_temp2)*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }


            //corner 5

            coeff_temp1 = -mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho_mid;
            coeff_temp2 = 2.0*mu_edge[2]*m_dt/(top_h[0]*top_h[1])/3.0/rho_mid;
            /*
            if(flag[2] == 1 && flag[1] == 0)
            {
                rhs += coeff_temp1*U0_nb[2];
                rhs += coeff_temp2*U0_nb[2];
            }
            else if(flag[2] == 0 && flag[1] ==1)
            {
                rhs += coeff_temp1*U0_nb[1];
                rhs += coeff_temp2*U0_nb[1];
            }
            else if(flag[2] ==1 && flag[1] == 1)
            {
                rhs += coeff_temp1*U0_nb[2];
                rhs += coeff_temp2*U0_nb[2];
            }
            else {
                rhs += coeff_temp1*(U0_nb[2]+U0_nb[1]+U0_nb[5]+U0_center)/8.0;
                rhs += coeff_temp2*(U0_nb[2]+U0_nb[1]+U0_nb[5]+U0_center)/8.0;

                solver.Add_A(I*2+1,I*2,       -coeff_temp1/8.0);
                solver.Add_A(I*2+1,I_nb[2]*2, -coeff_temp1/8.0);
                solver.Add_A(I*2+1,I_nb[1]*2, -coeff_temp1/8.0);
                solver.Add_A(I*2+1,I_nb[5]*2, -coeff_temp1/8.0);

                solver.Add_A(I*2+1,I*2,       -coeff_temp2/8.0);
                solver.Add_A(I*2+1,I_nb[2]*2, -coeff_temp2/8.0);
                solver.Add_A(I*2+1,I_nb[1]*2, -coeff_temp2/8.0);
                solver.Add_A(I*2+1,I_nb[5]*2, -coeff_temp2/8.0);
            } */
            if (flag[1] ==0 && flag[2] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U0_nb[2]+U0_nb[1]+U0_nb[5]+U0_center)/8.0;

                solver.Add_A(I*2+1,I*2,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2+1,I_nb[2]*2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2+1,I_nb[1]*2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2+1,I_nb[5]*2,  -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                getExactSolution_vd(corner_coords,m_t_old,corner_state);
                getExactSolution_vd(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*(coeff_temp1+coeff_temp2)*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }


            //corner 6

            coeff_temp1 = mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho_mid;
            coeff_temp2 = -2.0*mu_edge[3]*m_dt/(top_h[0]*top_h[1])/3.0/rho_mid;
            /*
            if(flag[1] == 1 && flag[3] == 0)
            {
                rhs += coeff_temp1*U0_nb[1];
                rhs += coeff_temp2*U0_nb[1];
            }
            else if(flag[1] == 0 && flag[3] ==1)
            {
                rhs += coeff_temp1*U0_nb[3];
                rhs += coeff_temp2*U0_nb[3];
            }
            else if(flag[1] ==1 && flag[3] == 1)
            {
                rhs += coeff_temp1*U0_nb[1];
                rhs += coeff_temp2*U0_nb[1];
            }
            else {
                rhs += coeff_temp1*(U0_nb[1]+U0_nb[3]+U0_nb[6]+U0_center)/8.0;
                rhs += coeff_temp2*(U0_nb[1]+U0_nb[3]+U0_nb[6]+U0_center)/8.0;

                solver.Add_A(I*2+1,I*2,       -coeff_temp1/8.0);
                solver.Add_A(I*2+1,I_nb[1]*2, -coeff_temp1/8.0);
                solver.Add_A(I*2+1,I_nb[3]*2, -coeff_temp1/8.0);
                solver.Add_A(I*2+1,I_nb[6]*2, -coeff_temp1/8.0);

                solver.Add_A(I*2+1,I*2,       -coeff_temp2/8.0);
                solver.Add_A(I*2+1,I_nb[1]*2, -coeff_temp2/8.0);
                solver.Add_A(I*2+1,I_nb[3]*2, -coeff_temp2/8.0);
                solver.Add_A(I*2+1,I_nb[6]*2, -coeff_temp2/8.0);
            } */
            if (flag[1] ==0 && flag[3] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U0_nb[1]+U0_nb[3]+U0_nb[6]+U0_center)/8.0;

                solver.Add_A(I*2+1,I*2,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2+1,I_nb[1]*2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2+1,I_nb[3]*2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2+1,I_nb[6]*2,  -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                getExactSolution_vd(corner_coords,m_t_old,corner_state);
                getExactSolution_vd(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*(coeff_temp1+coeff_temp2)*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }


            //corner 7

            coeff_temp1 = -mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho_mid;
            coeff_temp2 = 2.0*mu_edge[3]*m_dt/(top_h[0]*top_h[1])/3.0/rho_mid;
            /*
            if(flag[3] == 1 && flag[0] == 0)
            {
                rhs += coeff_temp1*U0_nb[3];
                rhs += coeff_temp2*U0_nb[3];
            }
            else if(flag[3] == 0 && flag[0] ==1)
            {
                rhs += coeff_temp1*U0_nb[0];
                rhs += coeff_temp2*U0_nb[0];
            }
            else if(flag[3] ==1 && flag[0] == 1)
            {
                rhs += coeff_temp1*U0_nb[3];
                rhs += coeff_temp2*U0_nb[3];
            }
            else {
                rhs += coeff_temp1*(U0_nb[3]+U0_nb[0]+U0_nb[7]+U0_center)/8.0;
                rhs += coeff_temp2*(U0_nb[3]+U0_nb[0]+U0_nb[7]+U0_center)/8.0;

                solver.Add_A(I*2+1,I*2,       -coeff_temp1/8.0);
                solver.Add_A(I*2+1,I_nb[3]*2, -coeff_temp1/8.0);
                solver.Add_A(I*2+1,I_nb[0]*2, -coeff_temp1/8.0);
                solver.Add_A(I*2+1,I_nb[7]*2, -coeff_temp1/8.0);

                solver.Add_A(I*2+1,I*2,       -coeff_temp2/8.0);
                solver.Add_A(I*2+1,I_nb[3]*2, -coeff_temp2/8.0);
                solver.Add_A(I*2+1,I_nb[0]*2, -coeff_temp2/8.0);
                solver.Add_A(I*2+1,I_nb[7]*2, -coeff_temp2/8.0);
            } */
            if (flag[3] ==0 && flag[0] == 0)
            {
                rhs += (coeff_temp1+coeff_temp2)*(U0_nb[3]+U0_nb[0]+U0_nb[7]+U0_center)/8.0;

                solver.Add_A(I*2+1,I*2,        -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2+1,I_nb[3]*2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2+1,I_nb[0]*2,  -(coeff_temp1+coeff_temp2)/8.0);
                solver.Add_A(I*2+1,I_nb[7]*2,  -(coeff_temp1+coeff_temp2)/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                getExactSolution_vd(corner_coords,m_t_old,corner_state);
                getExactSolution_vd(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*(coeff_temp1+coeff_temp2)*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }


            rhs += m_dt*state.m_U[1]/rho_mid;
            rhs += m_dt*cell_center[index].m_state.f_surf[1];
            rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho_mid;
            rhs -= m_dt*cell_center[index].m_state.m_adv[1];

            solver.Add_b(I*2+1, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

        start_clock("Before Petsc solver");

        solver.Solve_GMRES();
        solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);

        //if (rel_residual > 1)
        //{
        //    solver.Reset_x();
        //    solver.Solve_GMRES();
        //    solver.GetNumIterations(&num_iter);
        //    solver.GetFinalRelativeResidualNorm(&rel_residual);
        //}

        stop_clock("After Petsc solver");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
        {
            printf("\nIncompress_Solver_Smooth_2D_Cartesian::"
                        "compDiffWithSmoothProperty_velocity_vd: "
                        "num_iter = %d, rel_residual = %le. \n",
                        num_iter,rel_residual);
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ij_to_I[i][j];
            index = d_index2d(i,j,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*2-ilower*2];
                cell_center[index].m_state.m_U[1] = x[I*2+1-ilower*2];
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]);
                if (speed > max_speed)
                    max_speed = speed;
            }
            else
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]);
                if (speed > max_speed)
                    max_speed = speed;
            }
        }
        for (k = 0; k < 2; ++k) //scatter the velocity
        {
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                array[index] = cell_center[index].m_state.m_U[k];
            }
            scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                cell_center[index].m_state.m_U[k] = array[index];
            }
        }
        pp_global_max(&max_speed,1);
        FT_FreeThese(1,x);
}       /* end compDiffWithSmoothProperty2D_velocity_vd */


void Incompress_Solver_Smooth_2D_Cartesian::compDiffWithSmoothProperty_velocity_decoupled_vd(void)
{   
    COMPONENT comp;
    int index,index_nb[4],size;
    int I,I_nb[4];
    int i,j,l,nb,icoords[MAXD];
    L_STATE source_term;
    INTERFACE *intfc = front->interf;
    double coords[MAXD],crx_coords[MAXD];
    double coeff[4],mu[4],mu0,rho,corner[4],rhs;
    
    // U_nb contains states at neighbor cell or states on the boundary.   
    double U_nb[4],U_nb_new[4], U_center;
    
    double speed;
    double *x;
    GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
    double (*getStateVel[2])(POINTER) = {getStateXvel,getStateYvel};
    POINTER intfc_state;
    HYPER_SURF *hs;

    max_speed = 0.0;

    setIndexMap();

    size = iupper - ilower;
    FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

    for (l = 0; l < dim; ++l)
    {
        PETSc solver;
        solver.Create(ilower, iupper-1, 5, 5);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();

        for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I  = ij_to_I[i][j];
                if (I == -1) continue;

                index  = d_index2d(i,j,top_gmax);
                index_nb[0] = d_index2d(i-1,j,top_gmax);
                index_nb[1] = d_index2d(i+1,j,top_gmax);
                index_nb[2] = d_index2d(i,j-1,top_gmax);
                index_nb[3] = d_index2d(i,j+1,top_gmax);
                icoords[0] = i;
                icoords[1] = j;
                comp = top_comp[index];

                I_nb[0] = ij_to_I[i-1][j]; // left or west
                I_nb[1] = ij_to_I[i+1][j]; // right or east
                I_nb[2] = ij_to_I[i][j-1]; // down or south
                I_nb[3] = ij_to_I[i][j+1]; // up or north


                mu0 = cell_center[index].m_state.m_mu;
                rho = 0.5*(cell_center[index].m_state.m_rho +
                    cell_center[index].m_state.m_rho_old);
                U_center =  cell_center[index].m_state.m_U[l];

                for (nb = 0; nb < 4; nb++)
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
                double dh[2] = {top_h[0], top_h[1]};
                for(nb = 0; nb < 4; nb++)
                {
                    // use dh[1] as the edge length
                    if(nb >= 2)
                        std::swap(dh[0],dh[1]);

                    if(I_nb[nb] >= 0)
                    {
                        // u^{*}
                        solver.Add_A(I, I,        (-1) * -0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]));
                        solver.Add_A(I, I_nb[nb], (-1) *  0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]));
                        // u^{n}
                        rhs += 0.5*m_dt/rho*mu[nb] *
                                (U_nb[nb]-U_center) /(dh[0]*dh[0]);
                    }
                    else                // boundary
                    {
                        // u^{*}
                        solver.Add_A(I, I,       (-1) * -0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]) * 2);
                        rhs += 0.5*m_dt/rho*mu[nb] /(dh[0]*dh[0]) * 2 * U_nb_new[nb];
                        // u^{n}
                        rhs += 0.5*m_dt/rho*mu[nb] *
                                (U_nb[nb]-U_center) /(dh[0]*dh[0]) * 2;
                    }
                }

                // interior
                rhs += m_dt*source_term.m_U[l];
                rhs += m_dt*cell_center[index].m_state.f_surf[l];

               // rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho *dh[0]*dh[1];
                rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho;
                rhs -= m_dt * cell_center[index].m_state.m_adv[l];

                solver.Add_A(I, I, 1.0);
                rhs += U_center;
                solver.Add_b(I, rhs);
            }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-10);
        solver.Solve_GMRES();

        // get back the solution
        solver.Get_x(x);

        int num_iter;
        double rel_residual;
        solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);

        if (debugging("PETSc"))
            (void) printf("Incompress_Solver_Smooth_2D_Cartesian::"
                    "compDiffWithSmoothProperty_2nd_decoupled: "
                    "num_iter = %d, rel_residual = %le. \n",
                    num_iter,rel_residual);

        for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ij_to_I[i][j];
                index = d_index2d(i,j,top_gmax);
                if (I >= 0)
                {
                    cell_center[index].m_state.m_U[l] = x[I-ilower];
                }
                else
                {
                    cell_center[index].m_state.m_U[l] = 0.0;
                }
            }

        for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
        scatMeshArray();
        for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
    }

    for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index2d(i,j,top_gmax);
            speed = fabs(cell_center[index].m_state.m_U[0]) +
                    fabs(cell_center[index].m_state.m_U[1]);
            if (speed > max_speed)
                max_speed = speed;
        }
    pp_global_max(&max_speed,1);

    FT_FreeThese(1,x);
} /* end compDiffWithSmoothProperty_velocity_decoupled_vd */


void Incompress_Solver_Smooth_2D_Cartesian::computeNewConcentration_vd(void)
{
        COMPONENT comp;
        int index,index_nb[4],size;
        int I,I_nb[4];
        double coords[MAXD],crx_coords[MAXD];
        double coeff0[4],coeff1[4];
        double Dcoef[4],Dcoef_edge[4],Dcoef0,rhs;
        double rho[4], rho_edge[4], rho_old[4], rho_old_edge[4], rho_mid, rho0, rho0_old;
        double c_nb[4],c_center;
        int flag[4];  //denote whether this is dirichlet or neumann boundary
        L_STATE state;
        int i,j,k,nb,icoords[MAXD];
        INTERFACE *intfc = front->interf;
        double concentration;
        double *x;
        GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
        //double (*getStateConc)(POINTER) = getStateConc;
        POINTER intfc_state;
        HYPER_SURF *hs;
        PetscInt num_iter;
        double rel_residual;

        max_concentration = -1;
        min_concentration = HUGE;
        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

        PETSc solver;
        solver.Create(ilower, iupper-1, 5, 0);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ij_to_I[i][j];
            if (I == -1) continue;         // we do nothing with cells outside the boundary 

            index  = d_index2d(i,j,top_gmax);
            index_nb[0] = d_index2d(i-1,j,top_gmax);
            index_nb[1] = d_index2d(i+1,j,top_gmax);
            index_nb[2] = d_index2d(i,j-1,top_gmax);
            index_nb[3] = d_index2d(i,j+1,top_gmax);

            /*
             *       |          |
             *   7   |    3     |   6
             *-------|----------|-------
             *       |          |
             *   0   |          |   1
             *       |          |
             *-------|----------|--------
             *   4   |    2     |   5
             *       |          |
             */

            I_nb[0] = ij_to_I[i-1][j];
            I_nb[1] = ij_to_I[i+1][j];
            I_nb[2] = ij_to_I[i][j-1];
            I_nb[3] = ij_to_I[i][j+1];

            icoords[0] = i;
            icoords[1] = j;
            comp = top_comp[index];

            Dcoef0   = cell_center[index].m_state.m_Dcoef;
            rho0     = cell_center[index].m_state.m_rho;
            rho0_old = cell_center[index].m_state.m_rho_old;
            rho_mid  = 0.5*(rho0 + rho0_old);
            c_center = cell_center[index].m_state.m_c;

            for (nb = 0; nb < 4; nb++)
            {
                if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                {
                    flag[nb] = 1;
                   // c_nb[nb] = getStateConc(intfc_state);
                    c_nb[nb] = m_c[0];

                    if (wave_type(hs) == DIRICHLET_BOUNDARY ||
                        wave_type(hs) == NEUMANN_BOUNDARY)
                    {
                        Dcoef_edge[nb] = Dcoef[nb] = Dcoef0;
                        //rho_old[nb] = rho_old_edge[nb] = rho0_old;
                        //rho[nb] = rho_edge[nb] = rho0;
                        rho_old[nb] = rho_old_edge[nb] = m_rho[0];
                        rho[nb] = rho_edge[nb] = m_rho[0];
                    }
                    else
                    {
                        Dcoef[nb]        = cell_center[index_nb[nb]].m_state.m_Dcoef;
                        Dcoef_edge[nb]   = 0.5*(Dcoef0 + cell_center[index_nb[nb]].m_state.m_Dcoef);
                        rho_old[nb]      = cell_center[index_nb[nb]].m_state.m_rho_old;
                        rho_old_edge[nb] = 0.5*(rho0_old + cell_center[index_nb[nb]].m_state.m_rho_old);
                        rho[nb]          = cell_center[index_nb[nb]].m_state.m_rho;
                        rho_edge[nb]     = 0.5*(rho0 + cell_center[index_nb[nb]].m_state.m_rho);
                    }
                }
                else
                {
                    flag[nb]         = 0;
                    c_nb[nb]         = cell_center[index_nb[nb]].m_state.m_c;
                    Dcoef[nb]        = cell_center[index_nb[nb]].m_state.m_Dcoef;
                    Dcoef_edge[nb]   = 0.5*(Dcoef0 + cell_center[index_nb[nb]].m_state.m_Dcoef);
                    rho_old[nb]      = cell_center[index_nb[nb]].m_state.m_rho_old;
                    rho_old_edge[nb] = 0.5*(rho0_old + cell_center[index_nb[nb]].m_state.m_rho_old);
                    rho[nb]          = cell_center[index_nb[nb]].m_state.m_rho;
                    rho_edge[nb]     = 0.5*(rho0 + cell_center[index_nb[nb]].m_state.m_rho);
                }
            }

            //Setting the coefficients and matrix for c in concentration eqn.
            coeff0[0] = m_dt/2.0/rho_mid * rho_edge[0] * Dcoef_edge[0]/(top_h[0]*top_h[0]);  // left
            coeff0[1] = m_dt/2.0/rho_mid * rho_edge[1] * Dcoef_edge[1]/(top_h[0]*top_h[0]);  // right
            coeff0[2] = m_dt/2.0/rho_mid * rho_edge[2] * Dcoef_edge[2]/(top_h[1]*top_h[1]);  // down 
            coeff0[3] = m_dt/2.0/rho_mid * rho_edge[3] * Dcoef_edge[3]/(top_h[1]*top_h[1]);  // up

            coeff1[0] = m_dt/2.0/rho_mid * rho_old_edge[0] * Dcoef_edge[0]/(top_h[0]*top_h[0]);  // left
            coeff1[1] = m_dt/2.0/rho_mid * rho_old_edge[1] * Dcoef_edge[1]/(top_h[0]*top_h[0]);  // right
            coeff1[2] = m_dt/2.0/rho_mid * rho_old_edge[2] * Dcoef_edge[2]/(top_h[1]*top_h[1]);  // down 
            coeff1[3] = m_dt/2.0/rho_mid * rho_old_edge[3] * Dcoef_edge[3]/(top_h[1]*top_h[1]);  // up

            solver.Add_A(I,I,1.0);
            rhs = c_center;

            for (nb = 0; nb < 4; nb++)
            {
                if (I_nb[nb] != -1)                  // the neighbor cell is interior cell
                {
                    solver.Add_A(I,I_nb[nb],-coeff0[nb]);
                    rhs += coeff1[nb]*c_nb[nb];
                }
                else                                 // the neighbor cell on the bdry 
                {
                    coeff0[nb] = 2.0*coeff0[nb];     // compensate for the half space step on/near boundary, LHS
                    coeff1[nb] = 2.0*coeff1[nb];     // compensate for the half space step on/near boundary, RHS
                    rhs += (coeff0[nb]+coeff1[nb])*c_nb[nb]; // Dirichlet B.C., c_nb[nb]^(n+1)=c_nb[nb]^n, treated as rhs term 
                }

                solver.Add_A(I,I,coeff0[nb]);
                rhs -= coeff1[nb]*c_center;
            }

            rhs -= m_dt*cell_center[index].m_state.m_c_adv;    // advection term

            solver.Add_b(I, rhs);
        }
    
        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

        start_clock("Before Petsc solver");

        solver.Solve_GMRES();
        solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);

        //if (rel_residual > 1)
        //{
        //    solver.Reset_x();
        //    solver.Solve_GMRES();
        //    solver.GetNumIterations(&num_iter);
        //    solver.GetFinalRelativeResidualNorm(&rel_residual);
        //}

        stop_clock("After Petsc solver");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
        {
            printf("\nIncompress_Solver_Smooth_2D_Cartesian::"
                        "computeNewConcentration_vd: "
                        "num_iter = %d, rel_residual = %le. \n",
                        num_iter,rel_residual);
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ij_to_I[i][j];
            index = d_index2d(i,j,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_c = x[I-ilower];
                concentration = fabs(cell_center[index].m_state.m_c);
                if (concentration > max_concentration)
                    max_concentration = concentration;
                if (concentration < min_concentration)
                    min_concentration = concentration;
            }
            else
            {
                cell_center[index].m_state.m_c = 0.0;
                concentration = fabs(cell_center[index].m_state.m_c);
                if (concentration > max_concentration)
                    max_concentration = concentration;
                if (concentration < min_concentration)
                    min_concentration = concentration;
            }
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            array[index] = cell_center[index].m_state.m_c;
        }
        scatMeshArray();
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            cell_center[index].m_state.m_c = array[index];
        }
        pp_global_max(&max_concentration,1);
        pp_global_min(&min_concentration,1);
        FT_FreeThese(1,x);
}       /* end computeNewConcentration_vd */


void Incompress_Solver_Smooth_2D_Cartesian::computeProjection_vd(void)
{
        int index;
        int i,j,l,icoords[MAXD];
        double P_max,P_min;
        double **vel = iFparams->field->vel;
        double sum_div;
        double value;
        double diffusion[4];

        sum_div = 0.0;
        max_value = -1;

        for (l = 0; l < dim; ++l)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            vel[l][index] = cell_center[index].m_state.m_U[l]; //U*
        }

        /* Compute velocity divergence */
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            index = d_index2d(i,j,top_gmax);
            source_tmp[index] = computeFieldPointDiv(icoords,vel); // div(U*)
            getDivU_coupled_vd(icoords,diffusion,2);  // divergence constraint,i.e. div[u^(n+1)]
            source[index] = source_tmp[index] - diffusion[2];
            diff_coeff[index] = cell_center[index].m_state.m_rho;
            diff_coeff_old[index] = cell_center[index].m_state.m_rho_old;
        }
        FT_ParallelExchGridArrayBuffer(source,front);
        FT_ParallelExchGridArrayBuffer(source_tmp,front);
        FT_ParallelExchGridArrayBuffer(diff_coeff,front);
        FT_ParallelExchGridArrayBuffer(diff_coeff_old,front);
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index2d(i,j,top_gmax);
            cell_center[index].m_state.div_U = source_tmp[index]; // div(U*)
            source[index] /= accum_dt;
        }

        if(debugging("step_size"))
        {
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index2d(i,j,top_gmax);
                value = fabs(source_tmp[index]);
                sum_div = sum_div + source_tmp[index];
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            printf("\nThe summation of divergence of U^{star} is %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            printf("\nThe max value of divergence of U^{star} is %.16g\n",max_value);
            max_value = 0.0;
        }
        poisson_solver2d_vd(front,ilower,iupper,ij_to_I,source,diff_coeff,diff_coeff_old,
                        array,&P_max,&P_min);

        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index2d(i,j,top_gmax);
            cell_center[index].m_state.m_phi = array[index];
        }
}       /* end computeProjection_vd */


void Incompress_Solver_Smooth_2D_Cartesian::computeNewVelocity_vd(void)
{
        int i, j, k, index0, index;
        double grad_phi[2], rho;
        COMPONENT comp;
        double speed;
        int icoords[MAXD];

        int l;
        double **vel = iFparams->field->vel;
        double sum_div, value;

        max_speed = 0.0;

        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            array[index] = cell_center[index].m_state.m_phi;
        }
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index2d(i,j,top_gmax);
            comp = top_comp[index];
            if (!ifluid_comp(comp))
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
                continue;
            }
            rho = 0.5*(cell_center[index].m_state.m_rho + cell_center[index].m_state.m_rho_old);
            icoords[0] = i;
            icoords[1] = j;
            computeFieldPointGrad(icoords,array,grad_phi);
            cell_center[index].m_state.m_U[0] -= accum_dt/rho*grad_phi[0];
            cell_center[index].m_state.m_U[1] -= accum_dt/rho*grad_phi[1];
            speed = fabs(cell_center[index].m_state.m_U[0]) +
                    fabs(cell_center[index].m_state.m_U[1]);
            if (speed > max_speed)
                max_speed = speed;
        }
        for (k = 0; k < 2; ++k)
        {
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                array[index] = cell_center[index].m_state.m_U[k];
            }
            scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                cell_center[index].m_state.m_U[k] = array[index];
            }
        }
        pp_global_max(&max_speed,1);

        if(debugging("step_size"))
        {
            sum_div = 0.0;
            max_value = 0.0;

            for (l = 0; l < dim; ++l)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                vel[l][index] = cell_center[index].m_state.m_U[l];
            }

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                icoords[0] = i;
                icoords[1] = j;
                index = d_index2d(i,j,top_gmax);
                source[index] = computeFieldPointDiv(icoords,vel);
            }
            FT_ParallelExchGridArrayBuffer(source,front);

            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                cell_center[index].m_state.div_U = source[index];
            }

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index2d(i,j,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                sum_div = sum_div + cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U^{n+1} is %.16g\n",
                                        sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U^{n+1} is %.16g\n",
                                        max_value);
            max_value = 0.0;
        } 

} /* end computeNewVelocity_vd */


// for initial condition: 
//              setInitialCondition_vd();  
// this function should be called before solve_vd()
// for the source term of the momentum equation:        
//              computeSourceTerm();
void Incompress_Solver_Smooth_2D_Cartesian::solve_vd(double dt)
{

        printf("\nEntering Incompress Solve_vd! The dt for this time step is %.16g\n",dt);

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
        max_density = 0.0;

        start_clock("solve_vd");
        setDomain_vd();

        setComponent();
        if (debugging("trace"))
            printf("Passed setComponent()\n");
        setGlobalIndex();
        if (debugging("trace"))
            printf("Passed setGlobalIndex()\n");

        start_clock("setSmoothedProperties_vd");
        setSmoothedProperties_vd();
        stop_clock("setSmoothedProperties_vd");
        if (debugging("trace"))
            printf("Passed setSmoothedProperties_vd()\n");

        // 1) solve for intermediate U, density and concentration 
        // solve for estimated adv terms using lagged source term
        start_clock("compAdvectionTerm_coupled_vd(0)");
        compAdvectionTerm_coupled_vd(0);
        stop_clock("compAdvectionTerm_coupled_vd(0)");
        if (debugging("step_size"))
            printf("max_speed after computeAdvection_coupled_vd(0): %20.14f\n",
                                max_speed);

        // solve for estimated density explicitly
        start_clock("compNewDensity_vd(0)");
        computeNewDensity_vd(0);
        stop_clock("compNewDensity_vd(0)");

        // solve for accurate adv terms, i.e. corrector step of MAC
        start_clock("compAdvectionTerm_coupled_vd(1)");
        compAdvectionTerm_coupled_vd(1);
        stop_clock("compAdvectionTerm_coupled_vd(1)");
        if (debugging("step_size"))
            printf("max_speed after computeAdvection_coupled_vd(1): %20.14f\n",
                                max_speed);

        // solve for accurate density explicitly
        start_clock("compNewDensity_vd(1)");
        computeNewDensity_vd(1);
        stop_clock("compNewDensity_vd(1)");

        // solve for concentration implicitly
        start_clock("compNewConcentration_vd");
        computeNewConcentration_vd();
        stop_clock("compNewConcentration_vd");

        // solve for U^{\star} implicitly
        start_clock("compDiffWithSmoothProperty_velocity_vd");
        compDiffWithSmoothProperty_velocity_vd();
        stop_clock("compDiffWithSmoothProperty_velocity_vd");

        start_clock("compSGS");
        //compSGS();    //Subgrid model by Hyunkyun Lim
        stop_clock("compSGS");

        if (debugging("step_size"))
            printf("max_speed after compDiffWithSmoothProperty_velocity_vd(): %20.14f\n",
                                max_speed);

        // 2) projection step
        accum_dt += m_dt;
        if (accum_dt >= min_dt)
        {
            start_clock("computeProjection_vd");
            computeProjection_vd();
            stop_clock("computeProjection_vd");

            start_clock("computePressure");
            // Could PmII and PmIII be modified to be applied to vd?
            computePressure();
            stop_clock("computePressure");

            start_clock("computeNewVelocity_vd");
            computeNewVelocity_vd();
            stop_clock("computeNewVelocity_vd");
            accum_dt = 0.0;
        }

        if (debugging("sample_velocity"))
        {
            sampleVelocity();
        }

        if (debugging("step_size"))
            printf("max_speed after computeNewVelocity_vd(): %20.14f\n",
                                max_speed);

        start_clock("copyMeshStates_vd");
        copyMeshStates_vd();
        stop_clock("copyMeshStates_vd");

        setAdvectionDt();
        stop_clock("solve_vd");
}       /* end solve_vd */


void Incompress_Solver_Smooth_2D_Cartesian::copyMeshStates(void)
{
        int i,j,k,d,index;
        double **vel = field->vel;
        double *pres = field->pres;
        double *vort = field->vort;
        double **vort3d = field->vort3d;

        for (i = imin; i <= imax; ++i)
        for (j = jmin; j <= jmax; ++j)
        {
            index  = d_index2d(i,j,top_gmax);
            if (ifluid_comp(top_comp[index]))
            {
                pres[index] = cell_center[index].m_state.m_P;
                vel[0][index] = cell_center[index].m_state.m_U[0];
                vel[1][index] = cell_center[index].m_state.m_U[1];
                vort[index] = getVorticity(i,j);
            }
            else
            {
                pres[index] = 0.0;
                vel[0][index] = 0.0;
                vel[1][index] = 0.0;
                vort[index] = 0.0;
            }
        }

        FT_ParallelExchGridArrayBuffer(pres,front);
        FT_ParallelExchGridArrayBuffer(vort,front);
        FT_ParallelExchGridArrayBuffer(vel[0],front);
        FT_ParallelExchGridArrayBuffer(vel[1],front);
}       /* end copyMeshStates */


void Incompress_Solver_Smooth_2D_Cartesian::copyMeshStates_vd(void)
{
        int i,j,k,d,index;
        double **vel = field->vel;
        double *pres = field->pres;
        double *vort = field->vort;
        double **vort3d = field->vort3d;
        // for vd
        double *conc = field->conc;
        double *dens = field->dens;
        double *dens_old = field->dens_old;

        for (i = imin; i <= imax; ++i)
        for (j = jmin; j <= jmax; ++j)
        {
            index  = d_index2d(i,j,top_gmax);
            if (ifluid_comp(top_comp[index]))
            {
                pres[index] = cell_center[index].m_state.m_P;
                vel[0][index] = cell_center[index].m_state.m_U[0];
                vel[1][index] = cell_center[index].m_state.m_U[1];
                vort[index] = getVorticity(i,j);
                // for vd
                dens[index] = cell_center[index].m_state.m_rho;
                dens_old[index] = cell_center[index].m_state.m_rho_old;
                conc[index] = cell_center[index].m_state.m_c;
            }
            else
            {
                pres[index] = 0.0;
                vel[0][index] = 0.0;
                vel[1][index] = 0.0;
                vort[index] = 0.0;
                // for vd
                dens[index] = dens_old[index] = 0.0;
                conc[index] = 0.0;
            }
        }

        FT_ParallelExchGridArrayBuffer(pres,front);
        FT_ParallelExchGridArrayBuffer(vort,front);
        FT_ParallelExchGridArrayBuffer(vel[0],front);
        FT_ParallelExchGridArrayBuffer(vel[1],front);
        // for vd
        FT_ParallelExchGridArrayBuffer(dens,front);
        FT_ParallelExchGridArrayBuffer(dens_old,front);
        FT_ParallelExchGridArrayBuffer(conc,front);
}       /* end copyMeshStates_vd */


void Incompress_Solver_Smooth_2D_Cartesian::setInitialCondition_vd(LEVEL_FUNC_PACK *level_func_pack)
{
        int i;
        COMPONENT comp;
        double coords[MAXD];

        boolean status;
        int j,l,index,sign;
        double t[MAXD],*force;
        double center[MAXD],point[MAXD],nor[MAXD],phi,H,D;
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF *hs;
        int range = (int)(m_smoothing_radius+1);
        double rho;

        FT_MakeGridIntfc(front);
        setDomain_vd();
        setComponent();

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
            cell_center[i].m_state.setZero();
            comp = top_comp[i];
            if (getInitialState != NULL)
                (*getInitialState)(comp,coords,cell_center[i].m_state,dim,iFparams); // for U[k]
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            comp  = cell_center[index].comp;
            rho = cell_center[index].m_state.m_rho;
            if (!ifluid_comp(comp)) continue;

            getRectangleCenter(index, center);
            status = FT_FindNearestIntfcPointInRange(front,comp,center,point,
                                t,&hse,&hs,range);

            force = cell_center[index].m_state.f_surf;
            for (l = 0; l < dim; ++l) force[l] = 0.0;

            if (status == YES &&
                ifluid_comp(positive_component(hs)) &&
                ifluid_comp(negative_component(hs)))
            {
                sign = (comp == m_comp[0]) ? -1 : 1;
                D = smoothedDeltaFunction(center,point);
                H = smoothedStepFunction(center,point,sign);
                cell_center[index].m_state.m_c = m_c[0]  +
                                        (m_c[1]-m_c[0])*H; // ?????????????????
                cell_center[index].m_state.m_rho = m_rho[0] +
                                        (m_rho[1]-m_rho[0])*H;
                cell_center[index].m_state.m_rho_old = 
                               cell_center[index].m_state.m_rho;
                /*
                if (m_sigma != 0.0 && D != 0.0)
                {
                    surfaceTension_Fedkiw(D,hse,hs,force,m_sigma,rho,t);
                } */
            }
            else
            {
                switch (comp)
                {
                case LIQUID_COMP1:
                    cell_center[index].m_state.m_c = m_c[0];
                    cell_center[index].m_state.m_rho = m_rho[0];
                    cell_center[index].m_state.m_rho_old = 
                                   cell_center[index].m_state.m_rho;
                    break;
                case LIQUID_COMP2:
                    cell_center[index].m_state.m_c = m_c[1];
                    cell_center[index].m_state.m_rho = m_rho[1];
                    cell_center[index].m_state.m_rho_old = 
                                   cell_center[index].m_state.m_rho;
                    break;
                }
            }
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            array[index] = cell_center[index].m_state.m_c;
        }
        scatMeshArray();
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            cell_center[index].m_state.m_c = array[index];
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            array[index] = cell_center[index].m_state.m_rho;
        }
        scatMeshArray();
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            cell_center[index].m_state.m_rho = array[index];
            cell_center[index].m_state.m_rho_old = array[index];
        }

        for (l = 0; l < dim; l++) 
        {
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index2d(i,j,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index = d_index2d(i,j,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }

        computeGradientQ();
        copyMeshStates_vd();
        setAdvectionDt();
}       /* end setInitialCondition_vd */


void Incompress_Solver_Smooth_2D_Cartesian::compAdvectionTerm_coupled_vd(int flag)
{
    int I;
    int i,j,icoords[MAXD];
    int index,l;

    setIndexMap();

    if (!flag)
    {   
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ij_to_I[i][j];
            if (I == -1) continue;

            index  = d_index2d(i,j,top_gmax);
            icoords[0] = i;
            icoords[1] = j;
            //get U_bar, rho_bar, c_bar, and div(U_bar)
            getStatesBar_coupled_vd(icoords,m_t_int);
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            for (l=0; l<=3; l++)
            {
                cell_center[index].m_state.m_u_bar_tmp[l] =
                  cell_center[index].m_state.m_u_bar[l];
                cell_center[index].m_state.m_v_bar_tmp[l] =
                  cell_center[index].m_state.m_v_bar[l];
            }
        }

        //get phi in MAC step using Source^n
        computeMacPhi_vd(flag);
       
        //update U_bar using phi
        computeNewUbar_vd(m_t_int,flag);

        double div_U, value, sum_div;
        sum_div = 0.0;
        max_value = -1;
        if(debugging("step_size"))
        {
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index2d(i,j,top_gmax);
                div_U = (cell_center[index].m_state.m_u_bar[1]-cell_center[index].m_state.m_u_bar[0])/top_h[0]
                       +(cell_center[index].m_state.m_v_bar[3]-cell_center[index].m_state.m_v_bar[2])/top_h[1];
                value = fabs(div_U);
                sum_div = sum_div + div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            printf("\nThe summation of divergence of U^bar_0 is %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            printf("\nThe max value of divergence of U^bar_0 is %.16g\n",max_value);
            max_value = 0.0;
        }

        //get advection terms   
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            cell_center[index].m_state.m_rho_adv =
                (cell_center[index].m_state.m_rho_bar[1]*cell_center[index].m_state.m_u_bar[1] -
                 cell_center[index].m_state.m_rho_bar[0]*cell_center[index].m_state.m_u_bar[0])/top_h[0] +
                (cell_center[index].m_state.m_rho_bar[3]*cell_center[index].m_state.m_v_bar[3] -
                 cell_center[index].m_state.m_rho_bar[2]*cell_center[index].m_state.m_v_bar[2])/top_h[1];
        }
    }
    else // flag == 1
    {
        //get phi in MAC step using Source^(n+1/2)
        computeMacPhi_vd(flag);

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            for (l=0; l<=3; l++)
            {
                cell_center[index].m_state.m_u_bar[l] =
                  cell_center[index].m_state.m_u_bar_tmp[l];
		cell_center[index].m_state.m_v_bar[l] =
                  cell_center[index].m_state.m_v_bar_tmp[l];
            }
        }

        //update U_bar using phi
        computeNewUbar_vd(m_t_int,flag);
        
        double div_U, value, sum_div;
        sum_div = 0.0;
        max_value = -1;
        if(debugging("step_size"))
        {
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index2d(i,j,top_gmax);
                div_U = (cell_center[index].m_state.m_u_bar[1]-cell_center[index].m_state.m_u_bar[0])/top_h[0]
                       +(cell_center[index].m_state.m_v_bar[3]-cell_center[index].m_state.m_v_bar[2])/top_h[1];
                value = fabs(div_U);
                sum_div = sum_div + div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            printf("\nThe summation of divergence of U^bar_1 is %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            printf("\nThe max value of divergence of U^bar_1 is %.16g\n",max_value);
            max_value = 0.0;
        }

        //get advection terms   
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            cell_center[index].m_state.m_adv[0] =
                1.0/2*(cell_center[index].m_state.m_u_bar[1] + cell_center[index].m_state.m_u_bar[0]) *
                      (cell_center[index].m_state.m_u_bar[1] - cell_center[index].m_state.m_u_bar[0])/top_h[0] +
                1.0/2*(cell_center[index].m_state.m_v_bar[3] + cell_center[index].m_state.m_v_bar[2]) *
                      (cell_center[index].m_state.m_u_bar[3] - cell_center[index].m_state.m_u_bar[2])/top_h[1];
            cell_center[index].m_state.m_adv[1] =
                1.0/2*(cell_center[index].m_state.m_u_bar[1] + cell_center[index].m_state.m_u_bar[0]) *
                      (cell_center[index].m_state.m_v_bar[1] - cell_center[index].m_state.m_v_bar[0])/top_h[0] +
                1.0/2*(cell_center[index].m_state.m_v_bar[3] + cell_center[index].m_state.m_v_bar[2]) *
                      (cell_center[index].m_state.m_v_bar[3] - cell_center[index].m_state.m_v_bar[2])/top_h[1];
            // rho^(n+1/2) or rho^bar???
            cell_center[index].m_state.m_rho_adv =
                (cell_center[index].m_state.m_rho_bar[1]*cell_center[index].m_state.m_u_bar[1] -
                 cell_center[index].m_state.m_rho_bar[0]*cell_center[index].m_state.m_u_bar[0])/top_h[0] +
                (cell_center[index].m_state.m_rho_bar[3]*cell_center[index].m_state.m_v_bar[3] -
                 cell_center[index].m_state.m_rho_bar[2]*cell_center[index].m_state.m_v_bar[2])/top_h[1];
            cell_center[index].m_state.m_c_adv =
                1.0/2*(cell_center[index].m_state.m_u_bar[1] + cell_center[index].m_state.m_u_bar[0]) *
                      (cell_center[index].m_state.m_c_bar[1] - cell_center[index].m_state.m_c_bar[0])/top_h[0] +
                1.0/2*(cell_center[index].m_state.m_v_bar[3] + cell_center[index].m_state.m_v_bar[2]) *
                      (cell_center[index].m_state.m_c_bar[3] - cell_center[index].m_state.m_c_bar[2])/top_h[1];
        }
    }

        /*
        // convection terms of NS eqn's
        convectionTerm[0] =
                1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[0]-state_west_bar.m_U[0])/dx +
                1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[0]-state_south_bar.m_U[0])/dy;

        convectionTerm[1] =
                1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[1]-state_west_bar.m_U[1])/dx +
                1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[1]-state_south_bar.m_U[1])/dy;

        // convection term of continuity eqn: div(rho{\times}U)
        convectionTerm[2] =
                (state_east_bar.m_rho*state_east_bar.m_U[0] - state_west_bar.m_rho*state_west_bar.m_U[0])/dx +
                (state_north_bar.m_rho*state_north_bar.m_U[1] - state_south_bar.m_rho*state_south_bar.m_U[1])/dy;

        // convection term of concentration eqn: U{\cdot}grad_c
        convectionTerm[3] =
                1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_c-state_west_bar.m_c)/dx +
                1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_c-state_south_bar.m_c)/dy;
       */
} /* end compAdvectionTerm_coupled_vd */


void Incompress_Solver_Smooth_2D_Cartesian::computeMacPhi_vd(int flag)
{
        int index;
        int i,j,icoords[MAXD];
        double P_max,P_min;
        double sum_div;
        double value;
        double diffusion[4];

        sum_div = 0.0;
        max_value = -1;

        /* Compute velocity divergence */
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            index = d_index2d(i,j,top_gmax);
            //get the divergence constraint: div(U^{n}) or div(U^{n+1/2})
            getDivU_coupled_vd(icoords,diffusion,flag);
            source[index] = cell_center[index].m_state.div_U -
                            diffusion[2];
            if(!flag)
                diff_coeff[index] = cell_center[index].m_state.m_rho;
            else
                diff_coeff[index] = cell_center[index].m_state.m_rho_old;
        }
        FT_ParallelExchGridArrayBuffer(source,front);
        FT_ParallelExchGridArrayBuffer(diff_coeff,front);

        if(debugging("step_size") && !flag)
        {
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index2d(i,j,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                sum_div = sum_div + cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            printf("\nThe summation of divergence of U^bar is %.16g\n",sum_div);
            pp_global_max(&max_value,1);
            printf("\nThe max value of divergence of U^bar is %.16g\n",max_value);
            max_value = 0.0;
        }
        poisson_solver2d_MacPhi_vd(front,ilower,iupper,ij_to_I,source,diff_coeff,
                        array,&P_max,&P_min);

        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index2d(i,j,top_gmax);
            cell_center[index].m_state.m_phi = array[index];
        }
}       /* end computeMacPhi_vd */


void Incompress_Solver_Smooth_2D_Cartesian::computeNewUbar_vd(double t, int flag)
{
        int i, j, k, index;
        double grad_phi[2], phi, rho, rho_nb[4];
        COMPONENT comp;
        int icoords[MAXD];
        double crx_coords[MAXD];
        int l, nb;

        POINTER intfc_state;
        HYPER_SURF *hs;
        int ind[4],index_nb[8];

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index2d(i,j,top_gmax);
            index_nb[0] = d_index2d(i-1,j,top_gmax);
            index_nb[1] = d_index2d(i+1,j,top_gmax);
            index_nb[2] = d_index2d(i,j-1,top_gmax);
            index_nb[3] = d_index2d(i,j+1,top_gmax);
            index_nb[4] = d_index2d(i-1,j-1,top_gmax);
            index_nb[5] = d_index2d(i+1,j-1,top_gmax);
            index_nb[6] = d_index2d(i+1,j+1,top_gmax);
            index_nb[7] = d_index2d(i-1,j+1,top_gmax);
            comp = top_comp[index];
            //comp = cell_center[index].comp;
            phi = cell_center[index].m_state.m_phi;

            /*
             *       |          |
             *   7   |    3     |   6
             *-------|----------|-------
             *       |          |
             *   0   |          |   1
             *       |          |
             *-------|----------|--------
             *   4   |    2     |   5
             *       |          |
             */

            if (!ifluid_comp(comp))
            {
                for (l=0; l<=3; l++)
                {
                    cell_center[index].m_state.m_u_bar[l] = 0.0;
                    cell_center[index].m_state.m_v_bar[l] = 0.0;
                }
                continue;
            }

            icoords[0] = i;
            icoords[1] = j;
            // WEST, periodic B.C., actually no boundary
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,WEST,
                    comp,&intfc_state,&hs,crx_coords,t) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                ind[0] = 0;
            else
                ind[0] = 1;
            // EAST, periodic B.C., actually no boundary
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,EAST,
                    comp,&intfc_state,&hs,crx_coords,t) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                ind[1] = 0;
            else
                ind[1] = 1;
            // SOUTH, Neumann B.C.
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,SOUTH,
                    comp,&intfc_state,&hs,crx_coords,t) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                ind[2] = 0;
            else
                ind[2] = 1;
            // NORTH, Neumann B.C.
            if (FT_StateStructAtGridCrossing_tmp(front,icoords,NORTH,
                    comp,&intfc_state,&hs,crx_coords,t) &&
                    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                ind[3] = 0;
            else
                ind[3] = 1;

            if (!flag)
            {
                rho = cell_center[index].m_state.m_rho;
                for (nb=0; nb<4; nb++)
                    rho_nb[nb] = cell_center[index_nb[nb]].m_state.m_rho;
            }
            else // flag == 1
            {
                rho = cell_center[index].m_state.m_rho_old;
                for (nb=0; nb<4; nb++)
                    rho_nb[nb] = cell_center[index_nb[nb]].m_state.m_rho_old;
            }

            // homogeneous Neumann B.C. in SOUTH & NORTH for phi
            // cells on SOUTH boundary
            if (ind[2]==0)
            {
                cell_center[index].m_state.m_u_bar[0] -= (phi - cell_center[index_nb[0]].m_state.m_phi)/top_h[0]*2.0/(rho_nb[0]+rho);
                cell_center[index].m_state.m_v_bar[0] -= (cell_center[index_nb[3]].m_state.m_phi +
                                                          cell_center[index_nb[7]].m_state.m_phi -
                                                          cell_center[index_nb[0]].m_state.m_phi -
                                                          phi)/4/top_h[1]*2.0/(rho_nb[0]+rho);
                cell_center[index].m_state.m_u_bar[1] -= (cell_center[index_nb[1]].m_state.m_phi - phi)/top_h[0]*2.0/(rho_nb[1]+rho);
                cell_center[index].m_state.m_v_bar[1] -= (cell_center[index_nb[3]].m_state.m_phi +
                                                          cell_center[index_nb[6]].m_state.m_phi -
                                                          cell_center[index_nb[1]].m_state.m_phi -
                                                          phi)/4/top_h[1]*2.0/(rho_nb[1]+rho);
                cell_center[index].m_state.m_u_bar[3] -= (cell_center[index_nb[1]].m_state.m_phi +
                                                          cell_center[index_nb[6]].m_state.m_phi -
                                                          cell_center[index_nb[0]].m_state.m_phi -
                                                          cell_center[index_nb[7]].m_state.m_phi)/4/top_h[0]*2.0/(rho_nb[3]+rho);
                cell_center[index].m_state.m_v_bar[3] -= (cell_center[index_nb[3]].m_state.m_phi - phi)/top_h[1]*2.0/(rho_nb[3]+rho);
            }

            // cells on NORTH boundary
            if (ind[3]==0)
            {
                cell_center[index].m_state.m_u_bar[0] -= (phi - cell_center[index_nb[0]].m_state.m_phi)/top_h[0]*2.0/(rho_nb[0]+rho);
                cell_center[index].m_state.m_v_bar[0] -= -(cell_center[index_nb[2]].m_state.m_phi +
                                                          cell_center[index_nb[4]].m_state.m_phi -
                                                          cell_center[index_nb[0]].m_state.m_phi -
                                                          phi)/4/top_h[1]*2.0/(rho_nb[0]+rho);
                cell_center[index].m_state.m_u_bar[1] -= (cell_center[index_nb[1]].m_state.m_phi - phi)/top_h[0]*2.0/(rho_nb[1]+rho);
                cell_center[index].m_state.m_v_bar[1] -= -(cell_center[index_nb[2]].m_state.m_phi +
                                                          cell_center[index_nb[5]].m_state.m_phi -
                                                          cell_center[index_nb[1]].m_state.m_phi -
                                                          phi)/4/top_h[1]*2.0/(rho_nb[1]+rho);
                cell_center[index].m_state.m_u_bar[2] -= (cell_center[index_nb[1]].m_state.m_phi +
                                                          cell_center[index_nb[5]].m_state.m_phi -
                                                          cell_center[index_nb[0]].m_state.m_phi -
                                                          cell_center[index_nb[4]].m_state.m_phi)/4/top_h[0]*2.0/(rho_nb[2]+rho);
                cell_center[index].m_state.m_v_bar[2] -= (phi - cell_center[index_nb[2]].m_state.m_phi)/top_h[1]*2.0/(rho_nb[2]+rho);
            }

            // interior cells
            if (ind[2]==1 && ind[3]==1)
            {
                cell_center[index].m_state.m_u_bar[0] -= (phi - cell_center[index_nb[0]].m_state.m_phi)/top_h[0]*2.0/(rho_nb[0]+rho);
                cell_center[index].m_state.m_v_bar[0] -= (cell_center[index_nb[3]].m_state.m_phi +
                                                          cell_center[index_nb[7]].m_state.m_phi -
                                                          cell_center[index_nb[4]].m_state.m_phi -
                                                          cell_center[index_nb[2]].m_state.m_phi)/4/top_h[1]*2.0/(rho_nb[0]+rho);
                cell_center[index].m_state.m_u_bar[1] -= (cell_center[index_nb[1]].m_state.m_phi - phi)/top_h[0]*2.0/(rho_nb[1]+rho);
                cell_center[index].m_state.m_v_bar[1] -= (cell_center[index_nb[3]].m_state.m_phi +
                                                          cell_center[index_nb[6]].m_state.m_phi -
                                                          cell_center[index_nb[2]].m_state.m_phi -
                                                          cell_center[index_nb[5]].m_state.m_phi)/4/top_h[1]*2.0/(rho_nb[1]+rho);
                cell_center[index].m_state.m_u_bar[2] -= (cell_center[index_nb[1]].m_state.m_phi +
                                                          cell_center[index_nb[5]].m_state.m_phi -
                                                          cell_center[index_nb[0]].m_state.m_phi -
                                                          cell_center[index_nb[4]].m_state.m_phi)/4/top_h[0]*2.0/(rho_nb[2]+rho);
                cell_center[index].m_state.m_v_bar[2] -= (phi - cell_center[index_nb[2]].m_state.m_phi)/top_h[1]*2.0/(rho_nb[2]+rho);
                cell_center[index].m_state.m_u_bar[3] -= (cell_center[index_nb[1]].m_state.m_phi +
                                                          cell_center[index_nb[6]].m_state.m_phi -
                                                          cell_center[index_nb[0]].m_state.m_phi -
                                                          cell_center[index_nb[7]].m_state.m_phi)/4/top_h[0]*2.0/(rho_nb[3]+rho);
                cell_center[index].m_state.m_v_bar[3] -= (cell_center[index_nb[3]].m_state.m_phi - phi)/top_h[1]*2.0/(rho_nb[3]+rho);
            }
        }
} /* end computeNewUbar_vd */


void Incompress_Solver_Smooth_2D_Cartesian::getStatesBar_coupled_vd(
        int *icoords,
        double m_t_int)
{
    bool bNoBoundary;
    int ICoords[2];
    int index;
    L_STATE sl, sr, state_west_bar, state_east_bar, state_south_bar, state_north_bar;
    L_STATE state_west_hat, state_east_hat, state_south_hat, state_north_hat;
    L_STATE state_west_hat_l, state_west_hat_r;
    L_STATE state_east_hat_l, state_east_hat_r;
    L_STATE state_south_hat_l, state_south_hat_r;
    L_STATE state_north_hat_l, state_north_hat_r;

    double transverseD[4];
    double dx = top_h[0];
    double dy = top_h[1];

    //////////////////////// Get the state_hat on four edges first //////////////////////////

    // WEST
    bNoBoundary = getNeighborOrBoundaryState_vd(icoords,WEST,state_west_hat_l,m_t_int);
    if(!bNoBoundary)
    {
        state_west_hat = state_west_hat_l;
    }
    else
    {
        ICoords[0] = icoords[0] - 1;
        ICoords[1] = icoords[1];
        getFaceState_middleStep_hat_vd(ICoords,EAST,state_west_hat_l);
        getFaceState_middleStep_hat_vd(icoords,WEST,state_west_hat_r);
        getRiemannSolution_vd(COORD_X,state_west_hat_l,state_west_hat_r,state_west_hat);
    }

    // EAST
    bNoBoundary = getNeighborOrBoundaryState_vd(icoords,EAST,state_east_hat_r,m_t_int);
    if(!bNoBoundary)
    {
        state_east_hat = state_east_hat_r;
    }
    else
    {
        ICoords[0] = icoords[0] + 1;
        ICoords[1] = icoords[1];
        getFaceState_middleStep_hat_vd(ICoords,WEST,state_east_hat_r);
        getFaceState_middleStep_hat_vd(icoords,EAST,state_east_hat_l);
        getRiemannSolution_vd(COORD_X,state_east_hat_l,state_east_hat_r,state_east_hat);
    }

    // SOUTH
    bNoBoundary = getNeighborOrBoundaryState_vd(icoords,SOUTH,state_south_hat_l,m_t_int);
    if(!bNoBoundary)
    {
        state_south_hat = state_south_hat_l;
    }
    else
    {
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1] - 1;
        getFaceState_middleStep_hat_vd(ICoords,NORTH,state_south_hat_l);
        getFaceState_middleStep_hat_vd(icoords,SOUTH,state_south_hat_r);
        getRiemannSolution_vd(COORD_Y,state_south_hat_l,state_south_hat_r,state_south_hat);
    }

    // NORTH
    bNoBoundary = getNeighborOrBoundaryState_vd(icoords,NORTH,state_north_hat_r,m_t_int);
    if(!bNoBoundary)
    {
        state_north_hat = state_north_hat_r;
    }
    else
    {
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1] + 1;
        getFaceState_middleStep_hat_vd(ICoords,SOUTH,state_north_hat_r);
        getFaceState_middleStep_hat_vd(icoords,NORTH,state_north_hat_l);
        getRiemannSolution_vd(COORD_Y,state_north_hat_l,state_north_hat_r,state_north_hat);
    }

    ///////////////////////////// get the state_bar on four edges ////////////////////////////////

    // WEST

    transverseD[0] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]); // transverse m_U[0]
    transverseD[1] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]); // transverse m_U[1]
    transverseD[2] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_rho  - state_south_hat.m_rho) ; // transverse m_rho
    transverseD[3] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_c    - state_south_hat.m_c)   ; // transverse m_c

    bNoBoundary = getNeighborOrBoundaryState_vd(icoords,WEST,sl,m_t_int);
    if(!bNoBoundary)
    {
        state_west_bar = sl;
    }
    else
    {
        ICoords[0] = icoords[0] - 1;
        ICoords[1] = icoords[1];
        getFaceState_middleStep_bar_vd(ICoords,EAST,sl,transverseD,state_west_hat_l);
        getFaceState_middleStep_bar_vd(icoords,WEST,sr,transverseD,state_west_hat_r);
        getRiemannSolution_vd(COORD_X,sl,sr,state_west_bar);
    }

    // EAST

    transverseD[0] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]); // transverse m_U[0]
    transverseD[1] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]); // transverse m_U[1]
    transverseD[2] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_rho  - state_south_hat.m_rho) ; // transverse m_rho
    transverseD[3] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_c    - state_south_hat.m_c)   ; // transverse m_c


    bNoBoundary = getNeighborOrBoundaryState_vd(icoords,EAST,sr,m_t_int);
    if(!bNoBoundary)
    {
        state_east_bar = sr;
    }
    else
    {
        ICoords[0] = icoords[0] + 1;
        ICoords[1] = icoords[1];
        getFaceState_middleStep_bar_vd(ICoords,WEST,sr,transverseD,state_east_hat_r);
        getFaceState_middleStep_bar_vd(icoords,EAST,sl,transverseD,state_east_hat_l);
        getRiemannSolution_vd(COORD_X,sl,sr,state_east_bar);
    }

    // SOUTH

    transverseD[0] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]); // transverse m_U[0]
    transverseD[1] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]); // transverse m_U[1]
    transverseD[2] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_rho  - state_west_hat.m_rho) ; // transverse m_rho
    transverseD[3] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_c    - state_west_hat.m_c)   ; // transverse m_c


    bNoBoundary = getNeighborOrBoundaryState_vd(icoords,SOUTH,sl,m_t_int);
    if(!bNoBoundary)
    {
        state_south_bar = sl;
    }
    else
    {
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1] - 1;
        getFaceState_middleStep_bar_vd(ICoords,NORTH,sl,transverseD,state_south_hat_l);
        getFaceState_middleStep_bar_vd(icoords,SOUTH,sr,transverseD,state_south_hat_r);
        getRiemannSolution_vd(COORD_Y,sl,sr,state_south_bar);
    }

    // NORTH

    transverseD[0] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]); // transverse m_U[0]
    transverseD[1] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]); // transverse m_U[1]
    transverseD[2] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_rho  - state_west_hat.m_rho) ; // transverse m_rho    
    transverseD[3] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_c    - state_west_hat.m_c)   ; // transverse m_c    

    bNoBoundary = getNeighborOrBoundaryState_vd(icoords,NORTH,sr,m_t_int);
    if(!bNoBoundary)
    {
        state_north_bar = sr;
    }
    else
    {
        ICoords[0] = icoords[0];
        ICoords[1] = icoords[1] + 1;
        getFaceState_middleStep_bar_vd(ICoords,SOUTH,sr,transverseD,state_north_hat_r);
        getFaceState_middleStep_bar_vd(icoords,NORTH,sl,transverseD,state_north_hat_l);
        getRiemannSolution_vd(COORD_Y,sl,sr,state_north_bar);
    }

    index  = d_index2d(icoords[0],icoords[1],top_gmax);
    cell_center[index].m_state.m_u_bar[0] = state_west_bar.m_U[0];
    cell_center[index].m_state.m_u_bar[1] = state_east_bar.m_U[0];
    cell_center[index].m_state.m_u_bar[2] = state_south_bar.m_U[0];
    cell_center[index].m_state.m_u_bar[3] = state_north_bar.m_U[0];

    cell_center[index].m_state.m_v_bar[0] = state_west_bar.m_U[1];
    cell_center[index].m_state.m_v_bar[1] = state_east_bar.m_U[1];
    cell_center[index].m_state.m_v_bar[2] = state_south_bar.m_U[1];
    cell_center[index].m_state.m_v_bar[3] = state_north_bar.m_U[1];

    cell_center[index].m_state.m_rho_bar[0] = state_west_bar.m_rho;
    cell_center[index].m_state.m_rho_bar[1] = state_east_bar.m_rho;
    cell_center[index].m_state.m_rho_bar[2] = state_south_bar.m_rho;
    cell_center[index].m_state.m_rho_bar[3] = state_north_bar.m_rho;

    cell_center[index].m_state.m_c_bar[0] = state_west_bar.m_c;
    cell_center[index].m_state.m_c_bar[1] = state_east_bar.m_c;
    cell_center[index].m_state.m_c_bar[2] = state_south_bar.m_c;
    cell_center[index].m_state.m_c_bar[3] = state_north_bar.m_c;

    // get divergence of U_bar
    cell_center[index].m_state.div_U = (state_east_bar.m_U[0]-state_west_bar.m_U[0])/top_h[0]
                       + (state_north_bar.m_U[1]-state_south_bar.m_U[1])/top_h[1];
} /* end getStatesBar_coupled_vd */


void Incompress_Solver_Smooth_2D_Cartesian::getFaceState_middleStep_hat_vd(
        int *icoords,
        GRID_DIRECTION dir,
        L_STATE &state_hat)
{
    L_STATE  state_orig;
    int index;
    index = d_index2d(icoords[0],icoords[1],top_gmax);
    state_orig = cell_center[index].m_state;
    double sL;

    double dx = 0, dy = 0, slope_x_limited[4] = {0,0,0,0}, slope_y_limited[4] = {0,0,0,0};


    //    return the slope limiter, needed to be modified;

    //getLimitedSlope(icoords,COORD_X,slope_x_limited);
    //getLimitedSlope(icoords,COORD_Y,slope_y_limited);

    switch(dir)
    {
    case WEST:
        getLimitedSlope_vd(icoords,COORD_X,slope_x_limited);
        dx = top_h[0];
        if (state_orig.m_U[0] <= 0)
            sL = 1.0;
        else
            sL = 0.0;
        state_hat.m_U[0] = state_orig.m_U[0] + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[0];
        state_hat.m_U[1] = state_orig.m_U[1] + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[1];
        state_hat.m_rho  = state_orig.m_rho  + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_rho)*slope_x_limited[2];
        state_hat.m_c    = state_orig.m_c    + sL*(-dx/2.0 - m_dt/2.0 * state_orig.m_c)*slope_x_limited[3]; // non-conservative form of concentration eqn!
        break;


    case EAST:
        getLimitedSlope_vd(icoords,COORD_X,slope_x_limited);
        dx = top_h[0];
        if (state_orig.m_U[0] >= 0)
            sL = 1.0;
        else
            sL = 0.0;
        state_hat.m_U[0] = state_orig.m_U[0] + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[0];
        state_hat.m_U[1] = state_orig.m_U[1] + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_U[0])*slope_x_limited[1];
        state_hat.m_rho  = state_orig.m_rho  + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_rho)*slope_x_limited[2];
        state_hat.m_c    = state_orig.m_c    + sL*(dx/2.0 - m_dt/2.0 * state_orig.m_c)*slope_x_limited[3];
        break;


    case SOUTH:
        getLimitedSlope_vd(icoords,COORD_Y,slope_y_limited);
        dy = top_h[1];
        if (state_orig.m_U[1] <= 0)
            sL = 1.0;
        else
            sL = 0.0;
        state_hat.m_U[0] = state_orig.m_U[0] + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[0];
        state_hat.m_U[1] = state_orig.m_U[1] + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[1];
        state_hat.m_rho  = state_orig.m_rho  + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_rho)*slope_y_limited[2];
        state_hat.m_c    = state_orig.m_c    + sL*(-dy/2.0 - m_dt/2.0 * state_orig.m_c)*slope_y_limited[3];
        break;

    case NORTH:
        getLimitedSlope_vd(icoords,COORD_Y,slope_y_limited);
        dy = top_h[1];
        if (state_orig.m_U[1] >= 0)
            sL = 1.0;
        else
            sL = 0.0;
        state_hat.m_U[0] = state_orig.m_U[0] + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[0];
        state_hat.m_U[1] = state_orig.m_U[1] + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_U[1])*slope_y_limited[1];
        state_hat.m_rho  = state_orig.m_rho  + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_rho)*slope_y_limited[2];
        state_hat.m_c    = state_orig.m_c    + sL*(dy/2.0 - m_dt/2.0 * state_orig.m_c)*slope_y_limited[3];
        break;
    default:
        assert(false);
    }
} /* end getFaceState_middleStep_hat_vd */


void Incompress_Solver_Smooth_2D_Cartesian::getFaceState_middleStep_bar_vd(
        int *icoords,
        GRID_DIRECTION dir,
        L_STATE &state_bar,
        double transverseD[4],
        L_STATE state_hat)
{
    int index;
    index = d_index2d(icoords[0],icoords[1],top_gmax);
    double rho = cell_center[index].m_state.m_rho;

    double diffusion[4], gradP[2];
    gradP[0] = cell_center[index].m_state.grad_q[0];
    gradP[1] = cell_center[index].m_state.grad_q[1];

    getDiffusion_coupled_vd(icoords,diffusion);  // get viscous terms of NS eqn's 
    getDivU_coupled_vd(icoords,diffusion,0);     // get the divergence constraint
    getDiffusionC_coupled_vd(icoords,diffusion); // get the diffusion term of concentration eqn

    state_bar.m_U[0] = state_hat.m_U[0] + m_dt/2.0 * (-transverseD[0] + diffusion[0]/rho);
    state_bar.m_U[1] = state_hat.m_U[1] + m_dt/2.0 * (-transverseD[1] + diffusion[1]/rho);
    state_bar.m_rho  = state_hat.m_rho  + m_dt/2.0 * (-transverseD[2] - diffusion[2]*rho);
    state_bar.m_c    = state_hat.m_c    + m_dt/2.0 * (-transverseD[3] + diffusion[3]/rho);

    double coords[2];
    L_STATE source_term;

    getRectangleCenter(index, coords);
    computeSourceTerm(coords, source_term);

    state_bar.m_U[0] += m_dt/2 * source_term.m_U[0]/rho;
    state_bar.m_U[1] += m_dt/2 * source_term.m_U[1]/rho;
} /* end getFaceState_middleStep_bar_vd */


void Incompress_Solver_Smooth_2D_Cartesian::getDiffusion_coupled_vd(
        int *icoords,
        double diffusion[4])
{
    int index,index_nb[8];
    double mu[4],mu_edge[4],mu0;
    L_STATE Unb,corner_state;
    double U0_nb[8],U1_nb[8],U0_center, U1_center;
    int nb;
    GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
    bool bNoBoundary[4];
    double dh[2],dh0[2],dh1[2];
    double coords[MAXD],corner_coords[MAXD];

    diffusion[0] = 0.0;
    diffusion[1] = 0.0;

    int i = icoords[0];
    int j = icoords[1];
    index = d_index2d(i,j,top_gmax);
    index_nb[0] = d_index2d(i-1,j,top_gmax);
    index_nb[1] = d_index2d(i+1,j,top_gmax);
    index_nb[2] = d_index2d(i,j-1,top_gmax);
    index_nb[3] = d_index2d(i,j+1,top_gmax);

    index_nb[4] = d_index2d(i-1,j-1,top_gmax);
    index_nb[5] = d_index2d(i+1,j-1,top_gmax);
    index_nb[6] = d_index2d(i+1,j+1,top_gmax);
    index_nb[7] = d_index2d(i-1,j+1,top_gmax);

    mu0 = cell_center[index].m_state.m_mu;
    U0_center = cell_center[index].m_state.m_U[0];
    U1_center = cell_center[index].m_state.m_U[1];

    for (nb = 0; nb < 4; nb++)
    {
        bNoBoundary[nb] = getNeighborOrBoundaryState_vd(icoords,dir[nb],Unb,m_t_old);
        U0_nb[nb] = Unb.m_U[0];
        U1_nb[nb] = Unb.m_U[1];
        if(!bNoBoundary[nb])
        {
            mu[nb] = mu0;
            mu_edge[nb] = mu0;
        }
        else
        {
            mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
            mu_edge[nb] = 0.5*(cell_center[index_nb[nb]].m_state.m_mu + mu0);
        }
    }


    // non-cross derivative terms

    dh[0] = top_h[0];
    dh[1] = top_h[1];
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

    diffusion[0] += 4.0/3.0*(mu_edge[1]*(U0_nb[1]-U0_center)/dh1[0] - mu_edge[0]*(U0_center-U0_nb[0])/dh0[0])/dh[0];  // (4/3*mu*u_x)_x
    diffusion[1] +=         (mu_edge[1]*(U1_nb[1]-U1_center)/dh1[0] - mu_edge[0]*(U1_center-U1_nb[0])/dh0[0])/dh[0];  // (mu*v_x)_x

    diffusion[0] +=         (mu_edge[3]*(U0_nb[3]-U0_center)/dh1[1] - mu_edge[2]*(U0_center-U0_nb[2])/dh0[1])/dh[1];  // (mu*u_y)_y
    diffusion[1] += 4.0/3.0*(mu_edge[3]*(U1_nb[3]-U1_center)/dh1[1] - mu_edge[2]*(U1_center-U1_nb[2])/dh0[1])/dh[1];  // (4/3*mu*v_y)_y

    // get the coords in the cell center
    getRectangleCenter(index, coords);

    //cross derivative terms

    //traverse the four corners to get corner values

    //corner (i-1/2,j-1/2)
    /*
    if (!bNoBoundary[0] && bNoBoundary[2])
    {
        U0_nb[4] = U0_nb[0];
        U1_nb[4] = U1_nb[0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[2])
    {
        U0_nb[4] = U0_nb[2];
        U1_nb[4] = U1_nb[2];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[2])
    {
        U0_nb[4] = U0_nb[0];
        U1_nb[4] = U1_nb[0];
    }
    else
    {
        U0_nb[4] = (U0_nb[0]+U0_nb[2]+cell_center[index_nb[4]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[4] = (U1_nb[0]+U1_nb[2]+cell_center[index_nb[4]].m_state.m_U[1]+U1_center)/4.0;
    } */
    if (bNoBoundary[0] && bNoBoundary[2])
    {
        U0_nb[4] = (U0_nb[0]+U0_nb[2]+cell_center[index_nb[4]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[4] = (U1_nb[0]+U1_nb[2]+cell_center[index_nb[4]].m_state.m_U[1]+U1_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] - 0.5*top_h[0];
        corner_coords[1] = coords[1] - 0.5*top_h[1];
        getExactSolution_vd(corner_coords,m_t_old,corner_state);
        U0_nb[4] = corner_state.m_U[0];
        U1_nb[4] = corner_state.m_U[1];
    }

    //corner (i+1/2,j-1/2)
    /*
    if (!bNoBoundary[1] && bNoBoundary[2])
    {
        U0_nb[5] = U0_nb[1];
        U1_nb[5] = U1_nb[1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[2])
    {
        U0_nb[5] = U0_nb[2];
        U1_nb[5] = U1_nb[2];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[2])
    {
        U0_nb[5] = U0_nb[1];
        U1_nb[5] = U1_nb[1];
    }
    else
    {
        U0_nb[5] = (U0_nb[1]+U0_nb[2]+cell_center[index_nb[5]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[5] = (U1_nb[1]+U1_nb[2]+cell_center[index_nb[5]].m_state.m_U[1]+U1_center)/4.0;
    } */
    if (bNoBoundary[1] && bNoBoundary[2])
    {
        U0_nb[5] = (U0_nb[1]+U0_nb[2]+cell_center[index_nb[5]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[5] = (U1_nb[1]+U1_nb[2]+cell_center[index_nb[5]].m_state.m_U[1]+U1_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] + 0.5*top_h[0];
        corner_coords[1] = coords[1] - 0.5*top_h[1];
        getExactSolution_vd(corner_coords,m_t_old,corner_state);
        U0_nb[5] = corner_state.m_U[0];
        U1_nb[5] = corner_state.m_U[1];
    }

    //corner (i+1/2,j+1/2)
    /*
    if (!bNoBoundary[1] && bNoBoundary[3])
    {
        U0_nb[6] = U0_nb[1];
        U1_nb[6] = U1_nb[1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[3])
    {
        U0_nb[6] = U0_nb[3];
        U1_nb[6] = U1_nb[3];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[3])
    {
        U0_nb[6] = U0_nb[1];
        U1_nb[6] = U1_nb[1];
    }
    else
    {
        U0_nb[6] = (U0_nb[1]+U0_nb[3]+cell_center[index_nb[6]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[6] = (U1_nb[1]+U1_nb[3]+cell_center[index_nb[6]].m_state.m_U[1]+U1_center)/4.0;
    } */
    if (bNoBoundary[1] && bNoBoundary[3])
    {
        U0_nb[6] = (U0_nb[1]+U0_nb[3]+cell_center[index_nb[6]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[6] = (U1_nb[1]+U1_nb[3]+cell_center[index_nb[6]].m_state.m_U[1]+U1_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] + 0.5*top_h[0];
        corner_coords[1] = coords[1] + 0.5*top_h[1];
        getExactSolution_vd(corner_coords,m_t_old,corner_state);
        U0_nb[6] = corner_state.m_U[0];
        U1_nb[6] = corner_state.m_U[1];
    }

    //corner (i-1/2,j+1/2)
    /*
    if (!bNoBoundary[0] && bNoBoundary[3])
    {
        U0_nb[7] = U0_nb[0];
        U1_nb[7] = U1_nb[0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[3])
    {
        U0_nb[7] = U0_nb[3];
        U1_nb[7] = U1_nb[3];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[3])
    {
        U0_nb[7] = U0_nb[0];
        U1_nb[7] = U1_nb[0];
    }
    else
    {
        U0_nb[7] = (U0_nb[0]+U0_nb[3]+cell_center[index_nb[7]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[7] = (U1_nb[0]+U1_nb[3]+cell_center[index_nb[7]].m_state.m_U[1]+U1_center)/4.0;
    } */
    if (bNoBoundary[0] && bNoBoundary[3])
    {
        U0_nb[7] = (U0_nb[0]+U0_nb[3]+cell_center[index_nb[7]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[7] = (U1_nb[0]+U1_nb[3]+cell_center[index_nb[7]].m_state.m_U[1]+U1_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] - 0.5*top_h[0];
        corner_coords[1] = coords[1] + 0.5*top_h[1];
        getExactSolution_vd(corner_coords,m_t_old,corner_state);
        U0_nb[7] = corner_state.m_U[0];
        U1_nb[7] = corner_state.m_U[1];
    }

    diffusion[0] +=         (mu_edge[2]*U1_nb[4]-mu_edge[2]*U1_nb[5]+mu_edge[3]*U1_nb[6]-mu_edge[3]*U1_nb[7])/(top_h[0]*top_h[1]);  // (mu*v_x)_y
    diffusion[0] -= 2.0/3.0*(mu_edge[0]*U1_nb[4]-mu_edge[1]*U1_nb[5]+mu_edge[1]*U1_nb[6]-mu_edge[0]*U1_nb[7])/(top_h[0]*top_h[1]);  // -(2/3*mu*v_y)_x
    diffusion[1] +=         (mu_edge[0]*U0_nb[4]-mu_edge[1]*U0_nb[5]+mu_edge[1]*U0_nb[6]-mu_edge[0]*U0_nb[7])/(top_h[0]*top_h[1]);  // (mu*u_y)_x
    diffusion[1] -= 2.0/3.0*(mu_edge[2]*U0_nb[4]-mu_edge[2]*U0_nb[5]+mu_edge[3]*U0_nb[6]-mu_edge[3]*U0_nb[7])/(top_h[0]*top_h[1]);  // -(2/3*mu*u_x)_y
} /* end getDiffusion_coupled_vd */


void Incompress_Solver_Smooth_2D_Cartesian::getDivU_coupled_vd(
        int *icoords,
        double diffusion[4],
        int flag)
{
    int index,index_nb[4];
    double Dcoef[4],Dcoef_edge[4],Dcoef0;
    double rho_edge[4];
    L_STATE rhonb;
    double rho_nb[4],rho_center;
    int nb;
    GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
    bool bNoBoundary[4];
    double dh[2],dh0[2],dh1[2];

    diffusion[2] = 0.0;

    int i = icoords[0];
    int j = icoords[1];
    index = d_index2d(i,j,top_gmax);
    index_nb[0] = d_index2d(i-1,j,top_gmax);
    index_nb[1] = d_index2d(i+1,j,top_gmax);
    index_nb[2] = d_index2d(i,j-1,top_gmax);
    index_nb[3] = d_index2d(i,j+1,top_gmax);

    Dcoef0 = cell_center[index].m_state.m_Dcoef;
    if (flag == 0 || flag == 2)
        rho_center = cell_center[index].m_state.m_rho;
    if (flag == 1)
        rho_center = 0.5*(cell_center[index].m_state.m_rho + 
                     cell_center[index].m_state.m_rho_old);

    for (nb = 0; nb < 4; nb++)
    {
        if (flag == 0)
            bNoBoundary[nb] = getNeighborOrBoundaryState_vd(icoords,dir[nb],rhonb,m_t_old);
        if (flag == 1)
            bNoBoundary[nb] = getNeighborOrBoundaryState_vd(icoords,dir[nb],rhonb,m_t_int);
        if (flag == 2)
            bNoBoundary[nb] = getNeighborOrBoundaryState_vd(icoords,dir[nb],rhonb,m_t_new);
        if (flag == 0 || flag == 2)
            rho_nb[nb] = rhonb.m_rho;
        if (flag == 1)
            rho_nb[nb] = 0.5*(rhonb.m_rho + rhonb.m_rho_old);
            
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
            rho_edge[nb] = 0.5*(rho_nb[nb] + rho_center);
        }
    }

    dh[0] = top_h[0];
    dh[1] = top_h[1];
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

    // -div(Dcoef/rho*Grad_rho)
    diffusion[2] -= (Dcoef_edge[1]/rho_edge[1]*(rho_nb[1]-rho_center)/dh1[0] - Dcoef_edge[0]/rho_edge[0]*(rho_center-rho_nb[0])/dh0[0])/dh[0];
                                                  // -(Dcoef/rho*rho_x)_x
    diffusion[2] -= (Dcoef_edge[3]/rho_edge[3]*(rho_nb[3]-rho_center)/dh1[1] - Dcoef_edge[2]/rho_edge[2]*(rho_center-rho_nb[2])/dh0[1])/dh[1];
                                                  // -(Dcoef/rho*rho_y)_y
} /* end getDivU_coupled_vd */


void Incompress_Solver_Smooth_2D_Cartesian::getDiffusionC_coupled_vd(
        int *icoords,
        double diffusion[4])
{
    int index,index_nb[4];
    double Dcoef[4],Dcoef_edge[4],Dcoef0;
    double rho_edge[4], c_edge[4];
    L_STATE rhocnb;
    double rho_nb[4], c_nb[4],rho_center, c_center;
    int nb;
    GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
    bool bNoBoundary[4];
    double dh[2],dh0[2],dh1[2];

    diffusion[3] = 0.0;

    int i = icoords[0];
    int j = icoords[1];
    index = d_index2d(i,j,top_gmax);
    index_nb[0] = d_index2d(i-1,j,top_gmax);
    index_nb[1] = d_index2d(i+1,j,top_gmax);
    index_nb[2] = d_index2d(i,j-1,top_gmax);
    index_nb[3] = d_index2d(i,j+1,top_gmax);

    Dcoef0 = cell_center[index].m_state.m_Dcoef;
    c_center = cell_center[index].m_state.m_c;
    rho_center = cell_center[index].m_state.m_rho;

    for (nb = 0; nb < 4; nb++)
    {
        bNoBoundary[nb] = getNeighborOrBoundaryState_vd(icoords,dir[nb],rhocnb,m_t_old);
        c_nb[nb] = rhocnb.m_c;
        rho_nb[nb] = rhocnb.m_rho;

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
            rho_edge[nb] = 0.5*(rho_nb[nb] + rho_center);
        }
    }

    dh[0] = top_h[0];
    dh[1] = top_h[1];
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

    // div(Dcoef*rho*Grad_c)
    diffusion[3] += (Dcoef_edge[1]*rho_edge[1]*(c_nb[1]-c_center)/dh1[0] - Dcoef_edge[0]*rho_edge[0]*(c_center-c_nb[0])/dh0[0])/dh[0];  
                                                     // (Dcoef*rho*c_x)_x
    diffusion[3] += (Dcoef_edge[3]*rho_edge[3]*(c_nb[3]-c_center)/dh1[1] - Dcoef_edge[2]*rho_edge[2]*(c_center-c_nb[2])/dh0[1])/dh[1];
                                                     // (Dcoef*rho*c_y)_y
} /* end getDiffusionC_coupled_vd */


void Incompress_Solver_Smooth_2D_Cartesian::getLimitedSlope_vd(
        int *icoords,
        EBM_COORD xyz,
        double slope[4])
{
    double dh0, dh1;
    L_STATE U0, U1, U2;

    U1 = cell_center[d_index2d(icoords[0],icoords[1],top_gmax)].m_state;

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
    else        //
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
    slope[0] = EBM_minmod((U1.m_U[0]-U0.m_U[0])/dh0, (U2.m_U[0]-U1.m_U[0])/dh1);
    slope[1] = EBM_minmod((U1.m_U[1]-U0.m_U[1])/dh0, (U2.m_U[1]-U1.m_U[1])/dh1);
    slope[2] = EBM_minmod((U1.m_rho -U0.m_rho) /dh0, (U2.m_rho -U1.m_rho) /dh1);
    slope[3] = EBM_minmod((U1.m_c   -U0.m_c)   /dh0, (U2.m_c   -U1.m_c)   /dh1);
} /* end getLimitedSlope_vd */


bool Incompress_Solver_Smooth_2D_Cartesian::getNeighborOrBoundaryState_vd(
        int icoords[2],
        GRID_DIRECTION dir,
        L_STATE &state,
        double t)
{
    double crx_coords[MAXD];
    static double (*getStateVel[2])(POINTER) = {getStateXvel,getStateYvel};
    //static double (*getStateDens)(POINTER) = getStateDens;
    //static double (*getStateConc)(POINTER) = getStateConc;

    POINTER intfc_state;
    HYPER_SURF *hs;

    int index = d_index2d(icoords[0],icoords[1],top_gmax);
    int comp = cell_center[index].comp;

    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir,
            comp,&intfc_state,&hs,crx_coords,t) &&
            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
    {
        state.m_U[0] = getStateVel[0](intfc_state);
        state.m_U[1] = getStateVel[1](intfc_state);
       // state.m_rho  = getStateDens(intfc_state);
       // state.m_c    = getStateConc(intfc_state);
        state.m_rho  = m_rho[0];
        state.m_c    = m_c[0];
        return false;
    }
    else
    {
        int index_nb;
        switch(dir)
        {
        case WEST:
            index_nb = d_index2d(icoords[0]-1,icoords[1],top_gmax);
            break;
        case EAST:
            index_nb = d_index2d(icoords[0]+1,icoords[1],top_gmax);
            break;
        case SOUTH:
            index_nb = d_index2d(icoords[0],icoords[1]-1,top_gmax);
            break;
        case NORTH:
            index_nb = d_index2d(icoords[0],icoords[1]+1,top_gmax);
            break;
        default:
            assert(false);
        }
        state = cell_center[index_nb].m_state;
        return true;
    }
} /* end getNeighborOrBoundaryState_debug_vd */


/**
*
* @param state_left
* @param state_right
* @param ans
*/
void Incompress_Solver_Smooth_2D_Cartesian::getRiemannSolution_vd(
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
        sl.m_rho  = state_left.m_rho;
        sl.m_c    = state_left.m_c;

        sr.m_U[0] = state_right.m_U[0];
        sr.m_U[1] = state_right.m_U[1];
        sr.m_rho  = state_right.m_rho;
        sr.m_c    = state_right.m_c;
    }
    else
    {
        sl.m_U[0] = state_left.m_U[1];
        sl.m_U[1] = state_left.m_U[0];
        sl.m_rho  = state_left.m_rho;
        sl.m_c    = state_left.m_c;

        sr.m_U[0] = state_right.m_U[1];
        sr.m_U[1] = state_right.m_U[0];
        sr.m_rho  = state_right.m_rho;
        sr.m_c    = state_right.m_c;
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
    // rho_t + u * rho_x = 0
    // c_t + u * c_x = 0
    if(ans.m_U[0]>0)
    {
        ans.m_U[1] = sl.m_U[1];
        ans.m_rho  = sl.m_rho;
        ans.m_c    = sl.m_c;
    }
    else if(ans.m_U[0]<0)
    {
        ans.m_U[1] = sr.m_U[1];
        ans.m_rho  = sr.m_rho;
        ans.m_c    = sr.m_c;
    }
    else
    {
        ans.m_U[1] = 1.0/2*(sl.m_U[1]+sr.m_U[1]);
        ans.m_rho  = 1.0/2*(sl.m_rho + sr.m_rho);
        ans.m_c    = 1.0/2*(sl.m_c + sr.m_c);
    }

    // rotate state
    if(xyz==COORD_X)
        ; // do nothing
    else
        std::swap(ans.m_U[0],ans.m_U[1]);
} /* end getRiemannSolution_vd */


void Incompress_Solver_Smooth_2D_Cartesian::computeAdvection(void)
{
	int i,j,k,index,index00,index01,index10,index11,size;
	L_STATE state;
	COMPONENT comp;
	double speed;
	double *u, *v;
	double u0,v0,u00,u01,u10,u11,v00,v01,v10,v11;
	double crx_coords[MAXD];
	int icoords[MAXD];
	POINTER intfc_state;
	HYPER_SURF *hs;

	max_speed = 0.0;
	size = (top_gmax[0]+1)*(top_gmax[1]+1);
	FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    u[index] = cell_center[index].m_state.m_U[0];
	    v[index] = cell_center[index].m_state.m_U[1];
	}

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{	
	    icoords[0] = i;
	    icoords[1] = j;
	    index  = d_index2d(i,j,top_gmax);
	    comp = top_comp[index];
	    if (comp == SOLID_COMP)
	    {
	    	cell_center[index].m_state.m_U[0] = 0.0; 
	    	cell_center[index].m_state.m_U[1] = 0.0; 
		continue;
	    }
	    u0 = u[index];
	    v0 = v[index];
	    index00 = d_index2d(i-1,j,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,WEST,comp,
			        &intfc_state,&hs,crx_coords,m_t_old) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
	        u00 = getStateXvel(intfc_state);
	        v00 = getStateYvel(intfc_state);
	    }
	    else
	    {
		u00 = u[index00];
		v00 = v[index00];
	    }
	    index01 = d_index2d(i+1,j,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,EAST,comp,
			        &intfc_state,&hs,crx_coords,m_t_old) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
	    	u01 = getStateXvel(intfc_state);
	    	v01 = getStateYvel(intfc_state);
	    }
	    else
	    {
		u01 = u[index01];
		v01 = v[index01];
	    }
	    index10 = d_index2d(i,j-1,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,SOUTH,comp,
			        &intfc_state,&hs,crx_coords,m_t_old) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
	    	u10 = getStateXvel(intfc_state);
	    	v10 = getStateYvel(intfc_state);
	    }
	    else
	    {
		u10 = u[index10];
		v10 = v[index10];
	    }
	    index11 = d_index2d(i,j+1,top_gmax);
	    if (FT_StateStructAtGridCrossing_tmp(front,icoords,NORTH,comp,
			        &intfc_state,&hs,crx_coords,m_t_old) &&
		                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		u11 = getStateXvel(intfc_state);
		v11 = getStateYvel(intfc_state);
	    }
	    else
	    {
		u11 = u[index11];
		v11 = v[index11];
	    }

	    cell_center[index].m_state.m_U[0] += -m_dt*(
				burger_flux(u00,u0,u01)/top_h[0] +
				linear_flux(v0,u10,u0,u11)/top_h[1]);

	    cell_center[index].m_state.m_U[1] += - m_dt*(
				linear_flux(u0,v00,v0,v01)/top_h[0] +
			 	burger_flux(v10,v0,v11)/top_h[1]);

	    speed = fabs(cell_center[index].m_state.m_U[0]) +
		    fabs(cell_center[index].m_state.m_U[1]);
	    if (speed > max_speed)
		max_speed = speed;
	}
	for (k = 0; k < 2; ++k)
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	array[index] = cell_center[index].m_state.m_U[k];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	cell_center[index].m_state.m_U[k] = array[index];
	    }
	}
	pp_global_max(&max_speed,1);
	FT_FreeThese(2,u,v);
}	/* end computeAdvection2d */


void Incompress_Solver_Smooth_2D_Cartesian::
	compDiffWithSmoothProperty_1st_coupled(void)
{
        COMPONENT comp;
	int index,index_nb[8],size;
	int I,I_nb[8];
	double coords[MAXD],crx_coords[MAXD];
	double coeff[8],mu[4],mu_edge[4],mu0,rho,rhs,U0_nb[8],U1_nb[8];
	int flag[4]; //denote whether this is dirichlet or neumann boundary
	L_STATE state;
	int i,j,k,nb,icoords[MAXD];
	INTERFACE *intfc = front->interf;
	double speed;
	double *x;
	GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
	double (*getStateVel[2])(POINTER) = {getStateXvel,getStateYvel};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;

	max_speed = 0.0;
	setIndexMap();

	size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,2*size,sizeof(double));

	PETSc solver;

	solver.Create(2*ilower, 2*iupper-1, 9, 9);
		// two buffer zones, 5 u and 4 v for the first  equation
	        //                   5 v and 4 u for the second equation

	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    I  = ij_to_I[i][j];
	    if (I == -1) continue;

	    index  = d_index2d(i,j,top_gmax);
	    index_nb[0] = d_index2d(i-1,j,top_gmax);
	    index_nb[1] = d_index2d(i+1,j,top_gmax);
	    index_nb[2] = d_index2d(i,j-1,top_gmax);
	    index_nb[3] = d_index2d(i,j+1,top_gmax);
	    index_nb[4] = d_index2d(i-1,j-1,top_gmax);
	    index_nb[5] = d_index2d(i+1,j-1,top_gmax);
	    index_nb[6] = d_index2d(i+1,j+1,top_gmax);
	    index_nb[7] = d_index2d(i-1,j+1,top_gmax);

	    /*
	     *       |		|
	     *   7   |    3	|   6
	     *-------|----------|-------
	     *	     |		|
	     * 	 0   | 		|   1
	     *       |		|
	     *-------|----------|--------
	     *   4   |    2	|   5
	     *       |		|
	     */

	    I_nb[0] = ij_to_I[i-1][j];
	    I_nb[1] = ij_to_I[i+1][j];
	    I_nb[2] = ij_to_I[i][j-1];
	    I_nb[3] = ij_to_I[i][j+1];
	    I_nb[4] = ij_to_I[i-1][j-1];
	    I_nb[5] = ij_to_I[i+1][j-1];
	    I_nb[6] = ij_to_I[i+1][j+1];
	    I_nb[7] = ij_to_I[i-1][j+1];

	    icoords[0] = i;
	    icoords[1] = j;
	    comp = top_comp[index];

	    mu0 = cell_center[index].m_state.m_mu;
	    rho = cell_center[index].m_state.m_rho;

	    for (nb = 0; nb < 4; nb++)
	    {
		if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
			    comp,&intfc_state,&hs,crx_coords,m_t_old) &&
		            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		{
		    flag[nb] = 1;
		    U0_nb[nb] = getStateVel[0](intfc_state);
		    U1_nb[nb] = getStateVel[1](intfc_state);
		    
		    if (wave_type(hs) == DIRICHLET_BOUNDARY ||
		        wave_type(hs) == NEUMANN_BOUNDARY)
		    {
			mu_edge[nb] = mu0;
			mu[nb] = mu0;
		    }
		    else
		    {
			mu_edge[nb] = 0.5*(mu0 + 
				cell_center[index_nb[nb]].m_state.m_mu);
			mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		    }
		}
		else
		{
		    flag[nb] = 0;
		    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
		    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
		    mu_edge[nb] = 0.5*(mu0 + 
				cell_center[index_nb[nb]].m_state.m_mu);
		    mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
		}
	    }

	    // Neighbour cell 4
	    if (flag[0] == 1 && flag[2] == 1)
	    {
		U0_nb[4] = 0.5*(U0_nb[0] + U0_nb[2]);
		U1_nb[4] = 0.5*(U1_nb[0] + U1_nb[2]);
	    }
	    else if(flag[0] == 1 && flag[2] == 0)
	    {
		U0_nb[4] = U0_nb[0];
		U1_nb[4] = U1_nb[0];
	    }
	    else if(flag[0] == 0 && flag[2] == 1)
	    {
		U0_nb[4] = U0_nb[2];
		U1_nb[4] = U1_nb[2];
	    }
	    else  
	    {
		U0_nb[4] = cell_center[index_nb[4]].m_state.m_U[0];
		U1_nb[4] = cell_center[index_nb[4]].m_state.m_U[1];
	    }
	    // Neighbour cell 5
	    if (flag[1] == 1 && flag[2] == 1)
	    {
		U0_nb[5] = 0.5*(U0_nb[2] + U0_nb[1]);
		U1_nb[5] = 0.5*(U1_nb[2] + U1_nb[1]);
	    }
	    else if(flag[1] == 1 && flag[2] == 0)
	    {
		U0_nb[5] = U0_nb[1];
		U1_nb[5] = U1_nb[1];
	    }
	    else if(flag[1] == 0 && flag[2] == 1)
	    {
		U0_nb[5] = U0_nb[2];
		U1_nb[5] = U1_nb[2];
	    }
	    else
	    {
		U0_nb[5] = cell_center[index_nb[5]].m_state.m_U[0];
		U1_nb[5] = cell_center[index_nb[5]].m_state.m_U[1];
	    }
	    // Neighbour cell 6 
	    if (flag[1] == 1 && flag[3] == 1)
	    {
		U0_nb[6] = 0.5*(U0_nb[1] + U0_nb[3]);
		U1_nb[6] = 0.5*(U1_nb[1] + U1_nb[3]);
	    }
	    else if(flag[1] == 1 && flag[3] == 0)
	    {
		U0_nb[6] = U0_nb[1];
		U1_nb[6] = U1_nb[1];
	    }
	    else if(flag[1] == 0 && flag[3] == 1)
	    {
		U0_nb[6] = U0_nb[3];
		U1_nb[6] = U1_nb[3];
	    }
	    else
	    {
		U0_nb[6] = cell_center[index_nb[6]].m_state.m_U[0];
		U1_nb[6] = cell_center[index_nb[6]].m_state.m_U[1];
	    }
	    // Neighbour cell 7
	    if (flag[3] == 1 && flag[0] == 1)
	    {
		U0_nb[7] = 0.5*(U0_nb[3] + U0_nb[0]);
		U1_nb[7] = 0.5*(U1_nb[3] + U1_nb[0]);
	    }
	    else if(flag[3] == 1 && flag[0] == 0)
	    {
		U0_nb[7] = U0_nb[3];
		U1_nb[7] = U1_nb[3];
	    }
	    else if(flag[3] == 0 && flag[0] == 1)
	    {
		U0_nb[7] = U0_nb[0];
		U1_nb[7] = U1_nb[0];
	    }
	    else
	    {
		U0_nb[7] = cell_center[index_nb[7]].m_state.m_U[0];
		U1_nb[7] = cell_center[index_nb[7]].m_state.m_U[1];
	    }


	    // source term
	    getRectangleCenter(index, coords);
	    computeSourceTerm(coords, state);

	    //Setting the coefficients for the first equation
	    coeff[0] = m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);	// down
	    coeff[1] = m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);	// right
	    coeff[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);	// up	
	    coeff[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);	// left

	    coeff[4] = 0.5*m_dt/rho * mu[2]/(4*top_h[0]*top_h[1]);
	    coeff[5] = -0.5*m_dt/rho * mu[2]/(4*top_h[0]*top_h[1]);
	    coeff[6] = 0.5*m_dt/rho * mu[3]/(4*top_h[0]*top_h[1]);
	    coeff[7] = -0.5*m_dt/rho * mu[3]/(4*top_h[0]*top_h[1]);

	    solver.Set_A(I*2,I*2,1+coeff[0]+coeff[1]+coeff[2]+coeff[3]);
	    rhs = (1-coeff[0]-coeff[1]-coeff[2]-coeff[3])*
				cell_center[index].m_state.m_U[0];
	    

	    for (nb = 0; nb < 4; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*2,I_nb[nb]*2,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }
	    for (nb = 4; nb < 8; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*2,I_nb[nb]*2+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }

	    rhs += m_dt*state.m_U[0];
	    rhs += m_dt*cell_center[index].m_state.f_surf[0];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;

	    solver.Set_b(I*2, rhs);


	    /****************************************************************/
	    //Setting the coefficients for the second equation

	    coeff[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);	// down
	    coeff[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);	// right
	    coeff[2] = m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);	// up	
	    coeff[3] = m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);	// left

	    coeff[4] = 0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[1]);
	    coeff[5] = -0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[1]);
	    coeff[6] = 0.5*m_dt/rho * mu[1]/(4*top_h[0]*top_h[1]);
	    coeff[7] = -0.5*m_dt/rho * mu[0]/(4*top_h[0]*top_h[1]);


	    solver.Set_A(I*2+1,I*2+1,1+coeff[0]+coeff[1]+coeff[2]+coeff[3]);
	    rhs = (1-coeff[0]-coeff[1]-coeff[2]-coeff[3])*
				cell_center[index].m_state.m_U[1];

	    for (nb = 0; nb < 4; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*2+1,I_nb[nb]*2+1,-coeff[nb]);
		    rhs += coeff[nb]*U1_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U1_nb[nb];
	    }
	    for (nb = 4; nb < 8; nb++)
	    {
		if (I_nb[nb] != -1)
		{
		    solver.Set_A(I*2+1,I_nb[nb]*2,-coeff[nb]);
		    rhs += coeff[nb]*U0_nb[nb];
		}
		else
		    rhs += 2.0*coeff[nb]*U0_nb[nb];
	    }

	    rhs += m_dt*state.m_U[1];
	    rhs += m_dt*cell_center[index].m_state.f_surf[1];
	    rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;

	    solver.Set_b(I*2+1, rhs);
	}

	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);

	start_clock("Before Petsc solver");
	
	solver.Solve_GMRES();
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);

	//if (rel_residual > 1)
	//{
	//    solver.Reset_x();
	//    solver.Solve_GMRES();
	//    solver.GetNumIterations(&num_iter);
	//    solver.GetFinalRelativeResidualNorm(&rel_residual);
	//}

	stop_clock("After Petsc solver");

	// get back the solution
	solver.Get_x(x);

	if (debugging("PETSc"))
	{
	    printf("\nIncompress_Solver_Smooth_2D_Cartesian::"
			"compDiffWithSmoothProperty_1st_coupled: "
	       		"num_iter = %d, rel_residual = %le. \n", 
			num_iter,rel_residual); 
	}

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    I = ij_to_I[i][j];
	    index = d_index2d(i,j,top_gmax);
	    if (I >= 0)
	    {
		cell_center[index].m_state.m_U[0] = x[I*2-ilower*2];
		cell_center[index].m_state.m_U[1] = x[I*2+1-ilower*2];
		speed = fabs(cell_center[index].m_state.m_U[0]) +
		    	fabs(cell_center[index].m_state.m_U[1]);
		if (speed > max_speed)
		    max_speed = speed;
	    }
	    else
	    {
		cell_center[index].m_state.m_U[0] = 0.0;
		cell_center[index].m_state.m_U[1] = 0.0;
		speed = fabs(cell_center[index].m_state.m_U[0]) +
		    	fabs(cell_center[index].m_state.m_U[1]);
		if (speed > max_speed)
		    max_speed = speed;
	    }
	}
	for (k = 0; k < 2; ++k) //scatter the velocity
	{
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	array[index] = cell_center[index].m_state.m_U[k];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	cell_center[index].m_state.m_U[k] = array[index];
	    }
	}
	pp_global_max(&max_speed,1);
	FT_FreeThese(1,x);
}	/* end compDiffWithSmoothProperty2D */


void Incompress_Solver_Smooth_2D_Cartesian::compDiffWithSmoothProperty_2nd_coupled(void) 
{
        COMPONENT comp;
        int index,index_nb[8],size;
        int I,I_nb[8];
        double coords[MAXD],corner_coords[MAXD],crx_coords[MAXD];
        double coeff0[4],coeff1[4],coeff_temp;
        double mu[4],mu_edge[4],mu0,rho,rhs;
        double U0_nb[8],U1_nb[8],U0_center, U1_center;
        double U0_nb_new[4],U1_nb_new[4];
        int flag[4]; //denote whether this is dirichlet or neumann boundary       
        L_STATE state,corner_state,corner_state_new;
        int i,j,k,nb,icoords[MAXD];
        INTERFACE *intfc = front->interf;
        double speed;
        double *x;
        GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
        double (*getStateVel[2])(POINTER) = {getStateXvel,getStateYvel};          
        POINTER intfc_state;
        HYPER_SURF *hs;
        PetscInt num_iter;
        double rel_residual;

        max_speed = 0.0;
        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,2*size,sizeof(double));

        PETSc solver;

        solver.Create(2*ilower, 2*iupper-1, 15, 15);
                // two buffer zones, 5 u and 9 v for the first  equation
                //                   5 v and 9 u for the second equation
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ij_to_I[i][j];
            if (I == -1) continue;

            index  = d_index2d(i,j,top_gmax);
            index_nb[0] = d_index2d(i-1,j,top_gmax);
            index_nb[1] = d_index2d(i+1,j,top_gmax);
            index_nb[2] = d_index2d(i,j-1,top_gmax);
            index_nb[3] = d_index2d(i,j+1,top_gmax);
            index_nb[4] = d_index2d(i-1,j-1,top_gmax);
            index_nb[5] = d_index2d(i+1,j-1,top_gmax);
            index_nb[6] = d_index2d(i+1,j+1,top_gmax);
            index_nb[7] = d_index2d(i-1,j+1,top_gmax);

            /*
             *       |          |
             *   7   |    3     |   6
             *-------|----------|-------
             *       |          |
             *   0   |          |   1
             *       |          |
             *-------|----------|--------
             *   4   |    2     |   5
             *       |          |
             */

            I_nb[0] = ij_to_I[i-1][j];
            I_nb[1] = ij_to_I[i+1][j];
            I_nb[2] = ij_to_I[i][j-1];
            I_nb[3] = ij_to_I[i][j+1];
            I_nb[4] = ij_to_I[i-1][j-1];
            I_nb[5] = ij_to_I[i+1][j-1];
            I_nb[6] = ij_to_I[i+1][j+1];
            I_nb[7] = ij_to_I[i-1][j+1];

            icoords[0] = i;
            icoords[1] = j;
            comp = top_comp[index];

            mu0 = cell_center[index].m_state.m_mu;
            rho = cell_center[index].m_state.m_rho;
            U0_center = cell_center[index].m_state.m_U[0];
            U1_center = cell_center[index].m_state.m_U[1];

            for (nb = 0; nb < 4; nb++)
            {
                if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords,m_t_old) &&
                            wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                {
                    flag[nb] = 1;
                    // old boundary condition
                    U0_nb[nb] = getStateVel[0](intfc_state);
                    U1_nb[nb] = getStateVel[1](intfc_state);
                    // new boundary condition
                    FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
                            comp,&intfc_state,&hs,crx_coords,m_t_new);
                    U0_nb_new[nb] = getStateVel[0](intfc_state);
                    U1_nb_new[nb] = getStateVel[1](intfc_state);

                    if (wave_type(hs) == DIRICHLET_BOUNDARY ||
                        wave_type(hs) == NEUMANN_BOUNDARY)
                    {
                        mu_edge[nb] = mu0;
                        mu[nb] = mu0;
                    }
                    else
                    {
                        mu_edge[nb] = 0.5*(mu0 +
                                cell_center[index_nb[nb]].m_state.m_mu);
                        mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
                    }
                }
                else
                {
                    flag[nb] = 0;
                    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
                    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
                    mu_edge[nb] = 0.5*(mu0 +
                                cell_center[index_nb[nb]].m_state.m_mu);
                    mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
                }
            }

            for (nb = 4; nb < 8; nb++) // corner cell value, interior
            {
                if (I_nb[nb] != -1)
                {
                    U0_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[0];
                    U1_nb[nb] = cell_center[index_nb[nb]].m_state.m_U[1];
                }
            }

            // source term
            getRectangleCenter(index, coords);
            computeSourceTerm(coords, state);

            //Setting the coefficients and matrix for the U0 in first equation
            coeff0[0] = m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);      // down
            coeff0[1] = m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);      // right
            coeff0[2] = 0.5*m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);  // up   
            coeff0[3] = 0.5*m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);  // left


            solver.Add_A(I*2,I*2,1.0);
            rhs = U0_center;

            for (nb = 0; nb < 4; nb++)
            {
                if (flag[nb] == 0)
                {
                    solver.Add_A(I*2,I_nb[nb]*2,-coeff0[nb]);
                    rhs += coeff0[nb]*U0_nb[nb];
                }
                else if (flag[nb] == 1)
                {
                    coeff0[nb] = 2.0*coeff0[nb];
                    rhs += coeff0[nb]*(U0_nb_new[nb] + U0_nb[nb]);
                }

                solver.Add_A(I*2,I*2,coeff0[nb]);
                rhs -= coeff0[nb]*U0_center;
            }

            //set the coefficients and matrix for U1 in first equation
            //traverse the four corners
            //corner 4

            coeff_temp = mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
            if (flag[0] == 1 && flag[2] == 0)
                rhs += 0.5*coeff_temp*(U1_nb_new[0]+U1_nb[0]);
            else if(flag[0] == 0 && flag[2] ==1)
                rhs += 0.5*coeff_temp*(U1_nb_new[2]+U1_nb[2]);
            else if(flag[0] ==1 && flag[2] == 1)
                rhs += 0.5*coeff_temp*(U1_nb_new[0]+U1_nb[0]);
            else (flag[0] == 0 && flag[2] ==0)
            {
                rhs += coeff_temp*(U1_nb[0]+U1_nb[2]+U1_nb[4]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,       -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[0]*2+1, -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[2]*2+1, -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[4]*2+1, -coeff_temp/8.0);
            } */
            if (flag[0] == 0 && flag[2] ==0)
            {
                rhs += coeff_temp*(U1_nb[0]+U1_nb[2]+U1_nb[4]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,       -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[0]*2+1, -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[2]*2+1, -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[4]*2+1, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }

            //corner 5

            coeff_temp = -mu_edge[2]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
            if (flag[2] == 1 && flag[1] == 0)
                rhs += 0.5*coeff_temp*(U1_nb_new[2]+U1_nb[2]);
            else if(flag[2] == 0 && flag[1] == 1)
                rhs += 0.5*coeff_temp*(U1_nb_new[1]+U1_nb[1]);
            else if(flag[2] == 1 && flag[1] == 1)
                rhs += 0.5*coeff_temp*(U1_nb_new[2]+U1_nb[2]);
            else {
                rhs += coeff_temp*(U1_nb[2]+U1_nb[1]+U1_nb[5]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,        -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[2]*2+1,  -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[1]*2+1,  -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[5]*2+1,  -coeff_temp/8.0);
            } */
            if (flag[2] == 0 && flag[1] ==0)
            {
                rhs += coeff_temp*(U1_nb[2]+U1_nb[1]+U1_nb[5]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,       -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[2]*2+1, -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[1]*2+1, -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[5]*2+1, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }

            //corner 6

            coeff_temp = mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
            if (flag[1] == 1 && flag[3] == 0)
                rhs += 0.5*coeff_temp*(U1_nb_new[1]+U1_nb[1]);
            else if(flag[1] == 0 && flag[3] == 1)
                rhs += 0.5*coeff_temp*(U1_nb_new[3]+U1_nb[3]);
            else if(flag[1] == 1 && flag[3] == 1)
                rhs += 0.5*coeff_temp*(U1_nb_new[1]+U1_nb[1]);
            else {
                rhs += coeff_temp*(U1_nb[1]+U1_nb[3]+U1_nb[6]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,        -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[1]*2+1,  -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[3]*2+1,  -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[6]*2+1,  -coeff_temp/8.0);
            }*/
            if (flag[1] == 0 && flag[3] == 0)
            {
                rhs += coeff_temp*(U1_nb[1]+U1_nb[3]+U1_nb[6]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,        -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[1]*2+1,  -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[3]*2+1,  -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[6]*2+1,  -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }

            //corner 7

            coeff_temp = -mu_edge[3]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
            if (flag[3] == 1 && flag[0] == 0)
                rhs += 0.5*coeff_temp*(U1_nb_new[3]+U1_nb[3]);
            else if(flag[3] == 0 && flag[0] == 1)
                rhs += 0.5*coeff_temp*(U1_nb_new[0]+U1_nb[0]);
            else if(flag[3] == 1 && flag[0] == 1)
                rhs += 0.5*coeff_temp*(U1_nb_new[3]+U1_nb[3]);
            else {
                rhs += coeff_temp*(U1_nb[3]+U1_nb[0]+U1_nb[7]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,        -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[3]*2+1,  -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[0]*2+1,  -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[7]*2+1,  -coeff_temp/8.0);
            } */
            if (flag[3] == 0 && flag[0] == 0)
            {
                rhs += coeff_temp*(U1_nb[3]+U1_nb[0]+U1_nb[7]+U1_center)/8.0;

                solver.Add_A(I*2,I*2+1,        -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[3]*2+1,  -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[0]*2+1,  -coeff_temp/8.0);
                solver.Add_A(I*2,I_nb[7]*2+1,  -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[1]+
                    corner_state_new.m_U[1]);
            }


            rhs += m_dt*state.m_U[0]/rho;
            rhs += m_dt*cell_center[index].m_state.f_surf[0];
            rhs -= m_dt*cell_center[index].m_state.grad_q[0]/rho;
            rhs -= m_dt*cell_center[index].m_state.m_adv[0];

            solver.Add_b(I*2, rhs);


            /****************************************************************/
            //Setting the coefficients for U1 in the second equation

            coeff1[0] = 0.5*m_dt/rho * mu_edge[0]/(top_h[0]*top_h[0]);  // down
            coeff1[1] = 0.5*m_dt/rho * mu_edge[1]/(top_h[0]*top_h[0]);  // right
            coeff1[2] = m_dt/rho * mu_edge[2]/(top_h[1]*top_h[1]);      // up   
            coeff1[3] = m_dt/rho * mu_edge[3]/(top_h[1]*top_h[1]);      // left


            solver.Add_A(I*2+1,I*2+1,1.0);
            rhs = U1_center;


            for (nb = 0; nb < 4; nb++)
            {
                if (flag[nb] == 0)
                {
                    solver.Add_A(I*2+1,I_nb[nb]*2+1,-coeff1[nb]);
                    rhs += coeff1[nb]*U1_nb[nb];
                }
                else if (flag[nb] == 1)
                {
                    coeff1[nb] = 2.0*coeff1[nb];
                    rhs += coeff1[nb]*(U1_nb_new[nb] + U1_nb[nb]);
                }

                solver.Add_A(I*2+1,I*2+1,coeff1[nb]);
                rhs -= coeff1[nb]*U1_center;
            }

            //set the coefficients and matrix for U0 in second equation
            //traverse the four corners

            //corner 4

            coeff_temp = mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
            if (flag[0] == 1 && flag[2] == 0)
                rhs += 0.5*coeff_temp*(U0_nb_new[0]+U0_nb[0]);
            else if(flag[0] == 0 && flag[2] ==1)
                rhs += 0.5*coeff_temp*(U0_nb_new[2]+U0_nb[2]);
            else if(flag[0] ==1 && flag[2] == 1)
                rhs += 0.5*coeff_temp*(U0_nb_new[0]+U0_nb[0]);
            else {
                rhs += coeff_temp*(U0_nb[0]+U0_nb[2]+U0_nb[4]+U0_center)/8.0;

                solver.Add_A(I*2+1,I*2,       -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[0]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[2]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[4]*2, -coeff_temp/8.0);
            } */
            if (flag[0] ==0 && flag[2] == 0)
            {
                rhs += coeff_temp*(U0_nb[0]+U0_nb[2]+U0_nb[4]+U0_center)/8.0;

                solver.Add_A(I*2+1,I*2,       -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[0]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[2]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[4]*2, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }


            //corner 5

            coeff_temp = -mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
            if (flag[2] == 1 && flag[1] == 0)
                rhs += 0.5*coeff_temp*(U0_nb_new[2]+U0_nb[2]);
            else if(flag[2] == 0 && flag[1] ==1)
                rhs += 0.5*coeff_temp*(U0_nb_new[1]+U0_nb[1]);
            else if(flag[2] ==1 && flag[1] == 1)
                rhs += 0.5*coeff_temp*(U0_nb_new[2]+U0_nb[2]);
            else {
                rhs += coeff_temp*(U0_nb[2]+U0_nb[1]+U0_nb[5]+U0_center)/8.0;

                solver.Add_A(I*2+1,I*2,       -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[2]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[1]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[5]*2, -coeff_temp/8.0);
            } */
            if (flag[2] ==0 && flag[1] == 0)
            {
                rhs += coeff_temp*(U0_nb[2]+U0_nb[1]+U0_nb[5]+U0_center)/8.0;

                solver.Add_A(I*2+1,I*2,       -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[2]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[1]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[5]*2, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] - 0.5*top_h[1];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }

            //corner 6

            coeff_temp = mu_edge[1]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
            if (flag[1] == 1 && flag[3] == 0)
                rhs += 0.5*coeff_temp*(U0_nb_new[1]+U0_nb[1]);
            else if(flag[1] == 0 && flag[3] ==1)
                rhs += 0.5*coeff_temp*(U0_nb_new[3]+U0_nb[3]);
            else if(flag[1] ==1 && flag[3] == 1)
                rhs += 0.5*coeff_temp*(U0_nb_new[1]+U0_nb[1]);
            else {
                rhs += coeff_temp*(U0_nb[1]+U0_nb[3]+U0_nb[6]+U0_center)/8.0;
                                                    985,1-8       22%
                solver.Add_A(I*2+1,I*2,       -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[1]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[3]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[6]*2, -coeff_temp/8.0);
            } */
            if (flag[1] == 0 && flag[3] == 0)
            {
                rhs += coeff_temp*(U0_nb[1]+U0_nb[3]+U0_nb[6]+U0_center)/8.0;

                solver.Add_A(I*2+1,I*2,       -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[1]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[3]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[6]*2, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] + 0.5*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }

            //corner 7

            coeff_temp = -mu_edge[0]*m_dt/(top_h[0]*top_h[1])/rho;
            /*
            if (flag[3] == 1 && flag[0] == 0)
                rhs += 0.5*coeff_temp*(U0_nb_new[3]+U0_nb[3]);
            else if(flag[3] == 0 && flag[0] ==1)
                rhs += 0.5*coeff_temp*(U0_nb_new[0]+U0_nb[0]);
            else if(flag[3] ==1 && flag[0] == 1)
                rhs += 0.5*coeff_temp*(U0_nb_new[3]+U0_nb[3]);
            else {
                rhs += coeff_temp*(U0_nb[3]+U0_nb[0]+U0_nb[7]+U0_center)/8.0;
                                                    1022,1-8      23%
                solver.Add_A(I*2+1,I*2,       -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[3]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[0]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[7]*2, -coeff_temp/8.0);
            } */
            if (flag[3] == 0 && flag[0] == 0)
            {
                rhs += coeff_temp*(U0_nb[3]+U0_nb[0]+U0_nb[7]+U0_center)/8.0;

                solver.Add_A(I*2+1,I*2,       -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[3]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[0]*2, -coeff_temp/8.0);
                solver.Add_A(I*2+1,I_nb[7]*2, -coeff_temp/8.0);
            }
            else
            {
                corner_coords[0] = coords[0] - 0.5*top_h[0];
                corner_coords[1] = coords[1] + 0.5*top_h[1];
                getExactSolution(corner_coords,m_t_old,corner_state);
                getExactSolution(corner_coords,m_t_new,corner_state_new);
                rhs += 0.5*coeff_temp*(corner_state.m_U[0]+
                    corner_state_new.m_U[0]);
            }


            rhs += m_dt*state.m_U[1]/rho;
            rhs += m_dt*cell_center[index].m_state.f_surf[1];
            rhs -= m_dt*cell_center[index].m_state.grad_q[1]/rho;
            rhs -= m_dt*cell_center[index].m_state.m_adv[1];

            solver.Add_b(I*2+1, rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-10);

        start_clock("Before Petsc solver");

        solver.Solve_GMRES();
        solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);

        //if (rel_residual > 1)
        //{
        //    solver.Reset_x();
        //    solver.Solve_GMRES();
        //    solver.GetNumIterations(&num_iter);
        //    solver.GetFinalRelativeResidualNorm(&rel_residual);
        //}

        stop_clock("After Petsc solver");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
        {
            printf("\nIncompress_Solver_Smooth_2D_Cartesian::"
                        "compDiffWithSmoothProperty_2nd_coupled: "
                        "num_iter = %d, rel_residual = %le. \n",
                        num_iter,rel_residual);
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ij_to_I[i][j];
            index = d_index2d(i,j,top_gmax);
            if (I >= 0)
            {
                cell_center[index].m_state.m_U[0] = x[I*2-ilower*2];
                cell_center[index].m_state.m_U[1] = x[I*2+1-ilower*2];
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]);
                if (speed > max_speed)
                    max_speed = speed;

                printf("\t#### (i,j)={%3d,%3d}\n",i,j);
                printf("\t#### U_star_couple={%12.8f,%12.8f}\n",cell_center[index].m_state.m_U[0],
                      cell_center[index].m_state.m_U[1]);
            }
            else
            {
                cell_center[index].m_state.m_U[0] = 0.0;
                cell_center[index].m_state.m_U[1] = 0.0;
                speed = fabs(cell_center[index].m_state.m_U[0]) +
                        fabs(cell_center[index].m_state.m_U[1]);
                if (speed > max_speed)
                    max_speed = speed;

                printf("\t#### (i,j)={%3d,%3d}\n",i,j);
                printf("\t#### U_star_couple={%12.8f,%12.8f}\n",cell_center[index].m_state.m_U[0],
                      cell_center[index].m_state.m_U[1]);
            }
        }
        for (k = 0; k < 2; ++k) //scatter the velocity
        {
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                array[index] = cell_center[index].m_state.m_U[k];
            }
            scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                cell_center[index].m_state.m_U[k] = array[index];
            }
        }
        pp_global_max(&max_speed,1);
        FT_FreeThese(1,x);
} /* end compDiffWithSmoothProperty2D */


void Incompress_Solver_Smooth_2D_Cartesian::computeProjection_new(void)
{
	int index;
	int i,j,l,icoords[MAXD];
	double P_max,P_min;
	double **vel = iFparams->field->vel;
	double sum_div;
	double value;

	sum_div = 0.0;
	max_value = 0.0;

	for (l = 0; l < dim; ++l)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    vel[l][index] = cell_center[index].m_state.m_U[l];
	}

	/* Compute velocity (U^star) divergence */
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    index  = d_index(icoords,top_gmax,dim);
	    source[index] = computeFieldPointDiv(icoords,vel);
	    diff_coeff[index] = cell_center[index].m_state.m_rho;
	}
	FT_ParallelExchGridArrayBuffer(source,front);
	FT_ParallelExchGridArrayBuffer(diff_coeff,front);
	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    cell_center[index].m_state.div_U = source[index];
	    source[index] /= accum_dt;
	    array[index] = cell_center[index].m_state.m_phi;
	}

	if(debugging("step_size"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index2d(i,j,top_gmax);
	        value = fabs(cell_center[index].m_state.div_U);
		sum_div = sum_div + cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
	    }
	    pp_global_sum(&sum_div,1);
	    (void) printf("\nThe summation of divergence of U^{star} is %.16g\n",
					sum_div);
	    pp_global_max(&max_value,1);
	    (void) printf("\nThe max value of divergence of U^{star} is %.16g\n",
					max_value);
	    max_value = 0.0;
	}
	poisson_solver2d(front,ilower,iupper,ij_to_I,source,diff_coeff,
			array,&P_max,&P_min);

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    cell_center[index].m_state.m_phi = array[index];
	}
}	/* end computeProjection_new */

void Incompress_Solver_Smooth_2D_Cartesian::computeProjection(void)
{
	int index, index_nb[4], size;
	double rhs, coeff[4], rho[4], rho0;
	int I,I_nb[4];
	int i,j,l,icoords[MAXD];
	INTERFACE *intfc = front->interf;
	double P_max,P_min;
	int icrds_Pmax[MAXD],icrds_Pmin[MAXD];
	COMPONENT comp;
	double aII,rho_nb[4];
	double coords[MAXD],crx_coords[MAXD];
	double **vel = iFparams->field->vel;
	int num_nb;
	GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
	boolean use_neumann_solver = YES;
	max_value = 0.0;
	double value;
	double sum_div;
	sum_div = 0.0;
	PetscInt num_iter = 0;
	double rel_residual = 0.0;

	PETSc solver;
	solver.Create(ilower, iupper-1, 5, 5);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();
	size = iupper - ilower;

	setIndexMap();

	for (l = 0; l < dim; ++l)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    vel[l][index] = cell_center[index].m_state.m_U[l];
	}

	/* Compute velocity divergence */
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    index  = d_index(icoords,top_gmax,dim);
	    array[index] = computeFieldPointDiv(icoords,vel);
	}
	scatMeshArray();
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    cell_center[index].m_state.div_U = array[index];	    
	}

	if(debugging("step_size"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index2d(i,j,top_gmax);
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
	
	    rho0 = cell_center[index].m_state.m_rho;
	    num_nb = 0;
	    for (l = 0; l < 4; ++l)
	    {
		if (I_nb[l] == -1)
		    index_nb[l] = index;
		else num_nb++;
	    	rho[l] = 0.5*(rho0 + cell_center[index_nb[l]].m_state.m_rho);
	    	coeff[l] = 1.0/rho[l]/(top_h[l/2]*top_h[l/2]); 
	    }

	    rhs = cell_center[index].m_state.div_U/accum_dt;

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
            //if (aII != 0.0)
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
		printf("\nUnsolvable psudo-pressure, using original pressure!\n");
                solver.Set_A(I,I,1.0);
                rhs = cell_center[index].m_state.m_P;
            }
            solver.Set_b(I,rhs);
	}
	use_neumann_solver = pp_min_status(use_neumann_solver);
	
	solver.SetMaxIter(40000);
	solver.SetTol(1e-14);

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
	    (void) printf("Incompress_Solver_Smooth_2D_Cartesian::"
			"computeProjection: "
	       		"num_iter = %d, rel_residual = %le \n", 
			num_iter, rel_residual);
	
	P_max = -HUGE;		P_min = HUGE;
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    I = ij_to_I[i][j];
	    array[index] = x[I-ilower];
	}
	scatMeshArray();
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    cell_center[index].m_state.m_phi = array[index];
	}

	if(debugging("step_size"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index2d(i,j,top_gmax);
		value = fabs(cell_center[index].m_state.m_phi);
		if (value > max_value)
		    max_value = value;
	    }
	    pp_global_max(&max_value,1);
	    printf("\nThe max value of phi is %.16g\n",max_value);
	}

	FT_FreeThese(1,x);
}	/* end computeProjection */

void Incompress_Solver_Smooth_2D_Cartesian::computeNewVelocity(void)
{
	int i, j, k, index0, index;
	double grad_phi[2], rho;
	COMPONENT comp;
	double speed;
	int icoords[MAXD];
        /*
        int l;
        double **vel = iFparams->field->vel;
        double sum_div, value;
        */

	max_speed = 0.0;

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    array[index] = cell_center[index].m_state.m_phi;
	}
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    comp = top_comp[index];
	    if (!ifluid_comp(comp))
	    {
		cell_center[index].m_state.m_U[0] = 0.0;
		cell_center[index].m_state.m_U[1] = 0.0;
		continue;
	    }
	    rho = cell_center[index].m_state.m_rho;
	    icoords[0] = i;
	    icoords[1] = j;
	    computeFieldPointGradPhi(icoords,array,grad_phi);
	    cell_center[index].m_state.m_U[0] -= accum_dt/rho*grad_phi[0];
	    cell_center[index].m_state.m_U[1] -= accum_dt/rho*grad_phi[1];
	    speed = fabs(cell_center[index].m_state.m_U[0]) +
		    fabs(cell_center[index].m_state.m_U[1]);
	    if (speed > max_speed)
		max_speed = speed;
	}
	for (k = 0; k < 2; ++k)
	{
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	array[index] = cell_center[index].m_state.m_U[k];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	cell_center[index].m_state.m_U[k] = array[index];
	    }
	}
	pp_global_max(&max_speed,1);

        /*
        if(debugging("step_size"))
        {
            sum_div = 0.0;
            max_value = 0.0;

            for (l = 0; l < dim; ++l)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                vel[l][index] = cell_center[index].m_state.m_U[l];
            }

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                icoords[0] = i;
                icoords[1] = j;
                index = d_index2d(i,j,top_gmax);
                source[index] = computeFieldPointDiv(icoords,vel);
            }
            FT_ParallelExchGridArrayBuffer(source,front);

            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                cell_center[index].m_state.div_U = source[index];
            }

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index2d(i,j,top_gmax);
                value = fabs(cell_center[index].m_state.div_U);
                sum_div = sum_div + cell_center[index].m_state.div_U;
                if(value > max_value)
                    max_value = value;
            }
            pp_global_sum(&sum_div,1);
            (void) printf("\nThe summation of divergence of U^{n+1} is %.16g\n",
                                        sum_div);
            pp_global_max(&max_value,1);
            (void) printf("\nThe max value of divergence of U^{n+1} is %.16g\n",
                                        max_value);
            max_value = 0.0;
        }
        */
} /* end computeNewVelocity */


void Incompress_Solver_Smooth_2D_Cartesian::computeSourceTerm_Adv(double *coords, L_STATE &state) 
{
	int i;
	for (i = 0; i < dim; ++i)
	    state.m_U[i] = iFparams->gravity[i];

	state.m_P = HUGE_VAL;
} 


void Incompress_Solver_Smooth_2D_Cartesian::computeSourceTerm(double *coords, L_STATE &state) 
{
	int i;
	for (i = 0; i < dim; ++i)
	    state.m_U[i] = iFparams->gravity[i];

	state.m_P = HUGE_VAL;
}
void Incompress_Solver_Smooth_2D_Cartesian::computeSourceTerm(double *coords, double t, L_STATE &state) 
{
	computeSourceTerm(coords, state);
}


// for initial condition: 
// 		setInitialCondition();	
// this function should be called before solve()
// for the source term of the momentum equation: 	
// 		computeSourceTerm();
void Incompress_Solver_Smooth_2D_Cartesian::solve(double dt)
{

        printf("\nEntering Incompress Solve! The dt for this time step is %.16g\n",dt);

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
	if (debugging("step_size"))
	    printf("max_speed after computeAdvection(): %20.14f\n",
				max_speed);
	
	start_clock("compDiffWithSmoothProperty");
	compDiffWithSmoothProperty_1st_decoupled();
	//compDiffWithSmoothProperty();
	stop_clock("compDiffWithSmoothProperty");
	if (debugging("sample_velocity"))
	{
	    sampleVelocity();
	}

        start_clock("compSGS");
        //compSGS();	//Subgrid model by Hyunkyun Lim
        stop_clock("compSGS");

	if (debugging("step_size"))
	    printf("max_speed after compDiffWithSmoothProperty(): %20.14f\n",
				max_speed);

	// 2) projection step
	accum_dt += m_dt;
	if (accum_dt >= min_dt)
	{
	    start_clock("computeProjection");
	    //computeProjection_new();
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
	{
	    sampleVelocity();
	}

	if (debugging("step_size"))
	    printf("max_speed after computeNewVelocity(): %20.14f\n",
				max_speed);

	start_clock("copyMeshStates");
	copyMeshStates();
	stop_clock("copyMeshStates");

	setAdvectionDt();
	stop_clock("solve");
}	/* end solve */


double Incompress_Solver_Smooth_2D_Cartesian::getVorticity(int i, int j)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dy;
	double vorticity;

	dx = top_h[0];
	dy = top_h[1];
	index0 = d_index2d(i,j,top_gmax);
	index00 = d_index2d(i-1,j,top_gmax);
	index01 = d_index2d(i+1,j,top_gmax);
	index10 = d_index2d(i,j-1,top_gmax);
	index11 = d_index2d(i,j+1,top_gmax);
	v00 = -cell_center[index00].m_state.m_U[1];
	v01 =  cell_center[index01].m_state.m_U[1];
	v10 =  cell_center[index10].m_state.m_U[0];
	v11 = -cell_center[index11].m_state.m_U[0];

	vorticity = (v00 + v01)/2.0/dy + (v10 + v11)/2.0/dx;
	return vorticity;
}	/* end getVorticity */


void Incompress_Solver_Smooth_2D_Cartesian::
	compDiffWithSmoothProperty_1st_decoupled(void)
{
	COMPONENT comp;
        int index,index_nb[4],size;
        int I,I_nb[4];
        int i,j,l,nb,icoords[MAXD];
        L_STATE state;
        INTERFACE *intfc = front->interf;
        double coords[MAXD],crx_coords[MAXD];
	double coeff[4],mu[4],mu0,rho,corner[4],rhs,U_nb[4];
        double speed;
        double *x;
	GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
	double (*getStateVel[2])(POINTER) = {getStateXvel,getStateYvel};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;

	max_speed = 0.0;

        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	for (l = 0; l < dim; ++l)
	{
            PETSc solver;
            solver.Create(ilower, iupper-1, 5, 5);
	    solver.Reset_A();
	    solver.Reset_b();
	    solver.Reset_x();

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
            	I  = ij_to_I[i][j];
            	if (I == -1) continue;

            	index  = d_index2d(i,j,top_gmax);
            	index_nb[0] = d_index2d(i-1,j,top_gmax);
            	index_nb[1] = d_index2d(i+1,j,top_gmax);
            	index_nb[2] = d_index2d(i,j-1,top_gmax);
            	index_nb[3] = d_index2d(i,j+1,top_gmax);
		icoords[0] = i;
		icoords[1] = j;
		comp = top_comp[index];

            	I_nb[0] = ij_to_I[i-1][j]; // left or west
            	I_nb[1] = ij_to_I[i+1][j]; // right or east
            	I_nb[2] = ij_to_I[i][j-1]; // down or south
            	I_nb[3] = ij_to_I[i][j+1]; // up or north


		mu0 = cell_center[index].m_state.m_mu;
		rho = cell_center[index].m_state.m_rho;

            	for (nb = 0; nb < 4; nb++)
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

            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, state);

        	//first equation  decoupled, some terms may be lost
            	solver.Set_A(I,I,1+coeff[0]+coeff[1]+coeff[2]+coeff[3]);
            	rhs = (1-coeff[0]-coeff[1]-coeff[2]-coeff[3])
				*cell_center[index].m_state.m_U[l];

            	for (nb = 0; nb < 4; nb++)
            	{
		    if (I_nb[nb] != -1)
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
	    //	solver.Reset_x();
	    //	solver.Solve_GMRES();
	    //	solver.GetNumIterations(&num_iter);
            //	solver.GetFinalRelativeResidualNorm(&rel_residual);
	    //}
	    stop_clock("After Petsc solve");


            // get back the solution
            solver.Get_x(x);

            if (debugging("PETSc"))
                (void) printf("Incompress_Solver_Smooth_2D_Cartesian::"
			"compDiffWithSmoothProperty_1st_decoupled: "
                        "num_iter = %d, rel_residual = %le. \n",
                        num_iter,rel_residual);

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ij_to_I[i][j];
                index = d_index2d(i,j,top_gmax);
                if (I >= 0)
                {
                    cell_center[index].m_state.m_U[l] = x[I-ilower];
                }
                else
                {
                    cell_center[index].m_state.m_U[l] = 0.0;
                }
            }

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                array[index] = cell_center[index].m_state.m_U[l];
            }
            scatMeshArray();
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
                cell_center[index].m_state.m_U[l] = array[index];
            }
        }
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
            speed = fabs(cell_center[index].m_state.m_U[0]) +
                    fabs(cell_center[index].m_state.m_U[1]);
            if (speed > max_speed)
                max_speed = speed;
	}
        pp_global_max(&max_speed,1);

        FT_FreeThese(1,x);
}       /* end compDiffWithSmoothProperty2d_1st_decoupled */


void Incompress_Solver_Smooth_2D_Cartesian::computePressurePmI(void)
{
        int i,j,index;

	//scatMeshArray();
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
            index = d_index2d(i,j,top_gmax);
            cell_center[index].m_state.m_P += cell_center[index].m_state.m_phi;
	    array[index] = cell_center[index].m_state.m_P;
	}
	scatMeshArray();
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    cell_center[index].m_state.m_P = array[index];
	    cell_center[index].m_state.m_q = array[index];
	}

}        /* end computePressurePmI2d */


void Incompress_Solver_Smooth_2D_Cartesian::computePressurePmII(void)
{
        int i,j,index;
        double mu0;
	int icoords[MAXD];

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
            index = d_index2d(i,j,top_gmax);
            mu0 = 0.5*cell_center[index].m_state.m_mu;
	    cell_center[index].m_state.m_P = cell_center[index].m_state.m_q +
			cell_center[index].m_state.m_phi -
                        mu0*cell_center[index].m_state.div_U;
	    array[index] = cell_center[index].m_state.m_P;
	}
	scatMeshArray();
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    cell_center[index].m_state.m_P = array[index];
	    cell_center[index].m_state.m_q = array[index];
	}
}        /* end computePressurePmII2d */

void Incompress_Solver_Smooth_2D_Cartesian::computePressurePmIII(void)
{
        int i,j,index;
        double mu0;

	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
            index = d_index2d(i,j,top_gmax);
            mu0 = 0.5*cell_center[index].m_state.m_mu;
            cell_center[index].m_state.m_P = 
				cell_center[index].m_state.m_phi -
                        	accum_dt*mu0*cell_center[index].m_state.div_U;
	    cell_center[index].m_state.m_q = 0.0;
	}
}        /* end computePressurePmIII2d */


void Incompress_Solver_Smooth_2D_Cartesian::computePressure(void)
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


void Incompress_Solver_Smooth_2D_Cartesian::computeGradientQ(void)
{
	int i,j,l,index;
	double *grad_q;
	int icoords[MAXD];

	for (j = jmin; j <= jmax; ++j)
	for (i = imin; i <= imax; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = cell_center[index].m_state.m_q;
	}

	scatMeshArray();

	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    grad_q = cell_center[index].m_state.grad_q;
	    icoords[0] = i;
	    icoords[1] = j;
	    computeFieldPointGrad(icoords,array,grad_q);
	}
	for (l = 0; l < dim; ++l)
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index = d_index2d(i,j,top_gmax);
		array[index] = cell_center[index].m_state.grad_q[l];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	index = d_index2d(i,j,top_gmax);
		cell_center[index].m_state.grad_q[l] = array[index];
	    }
	}
}	/* end computeGradientQ2d */


void Incompress_Solver_Smooth_2D_Cartesian::surfaceTension(
	double *center,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *force,
        double sigma)
{
    	/*
    	double coords[3];
	int i,k,nb;
	BOND *pb,*bonds[100];
	BOND *b =  Bond_of_hse(hse);
	CURVE *c = Curve_of_hs(hs);
	double kappa,nor[MAXD];
	double kappa0,nor0[MAXD];
	double kappa1,nor1[MAXD];
	double len,delta;
	HYPER_SURF_ELEMENT *phse;
	double p[MAXD];
	

	BondAndNeighbors(hse,hs,&nb,bonds,3);

	for (i = 0; i < dim; ++i) force[i] = 0.0;
	for (k = 0; k < nb; ++k)
	{
	    pb = bonds[k];
	    for (i = 0; i < dim; ++i) 
		p[i] = 0.5*(Coords(pb->start)[i] + Coords(pb->end)[i]);
	    delta = smoothedDeltaFunction(coords,p);
	    if (delta == 0.0) continue;

	    len = bond_length(pb);
	    phse = Hyper_surf_element(pb);
	    GetFrontNormal(pb->start,phse,hs,nor0,front);
	    GetFrontNormal(pb->end,phse,hs,nor1,front);
	    for (i = 0; i < dim; ++i) nor[i] = 0.5*(nor0[i] + nor1[i]);
	    GetFrontCurvature(pb->start,phse,hs,&kappa0,front);
	    GetFrontCurvature(pb->end,phse,hs,&kappa1,front);
	    kappa = 0.5*(kappa0 + kappa1);
	    for (i = 0; i < dim; ++i) 
	    {
		force[i] += delta*sigma*len*kappa*nor[i];
	    }
	}
	for (i = 0; i < dim; i++)
	    force[i] /= cellVolume;
	*/
}	/* end surfaceTension2d */

void Incompress_Solver_Smooth_2D_Cartesian::setInitialCondition()
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


double Incompress_Solver_Smooth_2D_Cartesian::computeFieldPointDiv(
        int *icoords,
        double **field)
{
        int index;
        COMPONENT comp;
        int i, j, nb;
	int index_nb[4];
        double div,u_nb[2],v_nb[2];
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};      
      	double (*getStateVel[2])(POINTER) = {getStateXvel,getStateYvel};

	i = icoords[0];
	j = icoords[1];
	index = d_index2d(i,j,top_gmax);
        comp = top_comp[index];

	index_nb[0] = d_index2d(i-1,j,top_gmax);
	index_nb[1] = d_index2d(i+1,j,top_gmax);
	index_nb[2] = d_index2d(i,j-1,top_gmax);
	index_nb[3] = d_index2d(i,j+1,top_gmax);

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


        div = (u_nb[1] - u_nb[0])/top_h[0] + (v_nb[1] - v_nb[0])/top_h[1];
        return div;
}       /* end computeFieldPointDiv */


void Incompress_Solver_Smooth_2D_Cartesian::computeFieldPointGradPhi(
        int *icoords,
        double *field,
        double *grad_field)
{
        int index;
        COMPONENT comp;
        int i,j,nb;
	int index_nb[4];
        double p_nbedge[4],p0;  //the p values on the cell edges and cell center
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};      

	i = icoords[0];
	j = icoords[1];
	index = d_index2d(i,j,top_gmax);
        comp = top_comp[index];
	p0 = field[index];

	index_nb[0] = d_index2d(i-1,j,top_gmax);
	index_nb[1] = d_index2d(i+1,j,top_gmax);
	index_nb[2] = d_index2d(i,j-1,top_gmax);
	index_nb[3] = d_index2d(i,j+1,top_gmax);

	for (nb = 0; nb < 4; nb++)
	{
	    if(FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
			comp,&intfc_state,&hs,crx_coords) &&
		        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
		    boundary_state(hs) == NULL) //flow throught B.C.
		    p_nbedge[nb] = 0.0;
		else
		{  
		    p_nbedge[nb] = p0;
		    
		    //linear extrapolation to the boundary
		    /*
		    if (nb == 0)
			p_nbedge[nb] = p0 - 0.5 * (field[index_nb[1]] - p0);
		    else if (nb == 1)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[0]]);
		    else if (nb == 2)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[3]] - p0);
		    else if (nb == 3)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[2]]);
		    */

		}
	    }
	    else
		p_nbedge[nb] = (p0 + field[index_nb[nb]])/2.0;
	}
	grad_field[0] = (p_nbedge[1] - p_nbedge[0])/top_h[0];
	grad_field[1] = (p_nbedge[3] - p_nbedge[2])/top_h[1];
}


void Incompress_Solver_Smooth_2D_Cartesian::computeFieldPointGrad(
        int *icoords,
        double *field,
        double *grad_field)
{
        int index;
        COMPONENT comp;
        int i,j,nb;
	int index_nb[4];
        double p_nbedge[4],p0;  //the p values on the cell edges and cell center
        double coords[MAXD],crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};      

	i = icoords[0];
	j = icoords[1];
	index = d_index2d(i,j,top_gmax);
        comp = top_comp[index];
	p0 = field[index];

	index_nb[0] = d_index2d(i-1,j,top_gmax);
	index_nb[1] = d_index2d(i+1,j,top_gmax);
	index_nb[2] = d_index2d(i,j-1,top_gmax);
	index_nb[3] = d_index2d(i,j+1,top_gmax);

	for (nb = 0; nb < 4; nb++)
	{
	    if(FT_StateStructAtGridCrossing_tmp(front,icoords,dir[nb],
			comp,&intfc_state,&hs,crx_coords) &&
		        wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
	    {
		if (wave_type(hs) == DIRICHLET_BOUNDARY &&
		    boundary_state(hs) == NULL) //flow throught B.C.
		    p_nbedge[nb] = 0.0;
		else
		{  
		    //p_nbedge[nb] = p0;
		    
		    //linear extrapolation to the boundary
		    if (nb == 0)
			p_nbedge[nb] = p0 - 0.5 * (field[index_nb[1]] - p0);
		    else if (nb == 1)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[0]]);
		    else if (nb == 2)
			p_nbedge[nb] = p0 - 0.5*(field[index_nb[3]] - p0);
		    else if (nb == 3)
			p_nbedge[nb] = p0 + 0.5*(p0 - field[index_nb[2]]);
		    

		}
	    }
	    else
		p_nbedge[nb] = (p0 + field[index_nb[nb]])/2.0;
	}
	grad_field[0] = (p_nbedge[1] - p_nbedge[0])/top_h[0];
	grad_field[1] = (p_nbedge[3] - p_nbedge[2])/top_h[1];
}


void Incompress_Solver_Smooth_2D_Cartesian::compDiffWithSmoothProperty_2nd_decoupled(void)
{
    COMPONENT comp;
    int index,index_nb[4],size;
    int I,I_nb[4];
    int i,j,l,nb,icoords[MAXD];
    L_STATE source_term;
    INTERFACE *intfc = front->interf;
    double coords[MAXD],crx_coords[MAXD];
    double coeff[4],mu[4],mu0,rho,corner[4],rhs;

    // U_nb contains states at neighbor cell or states on the boundary.
    double U_nb[4],U_nb_new[4], U_center;

    double speed;
    double *x;
    GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
    double (*getStateVel[2])(POINTER) = {getStateXvel,getStateYvel};
    POINTER intfc_state;
    HYPER_SURF *hs;

    max_speed = 0.0;

    setIndexMap();

    size = iupper - ilower;
    FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

    for (l = 0; l < dim; ++l)
    {
	PETSc solver;
	solver.Create(ilower, iupper-1, 5, 5);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		I  = ij_to_I[i][j];
		if (I == -1) continue;

		index  = d_index2d(i,j,top_gmax);
		index_nb[0] = d_index2d(i-1,j,top_gmax);
		index_nb[1] = d_index2d(i+1,j,top_gmax);
		index_nb[2] = d_index2d(i,j-1,top_gmax);
		index_nb[3] = d_index2d(i,j+1,top_gmax);
		icoords[0] = i;
		icoords[1] = j;
		comp = top_comp[index];

		I_nb[0] = ij_to_I[i-1][j]; // left or west
		I_nb[1] = ij_to_I[i+1][j]; // right or east
		I_nb[2] = ij_to_I[i][j-1]; // down or south
		I_nb[3] = ij_to_I[i][j+1]; // up or north


		mu0 = cell_center[index].m_state.m_mu;
		rho = cell_center[index].m_state.m_rho;
		U_center = cell_center[index].m_state.m_U[l];

		for (nb = 0; nb < 4; nb++)
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
		double dh[2] = {top_h[0], top_h[1]};
		for(nb = 0; nb < 4; nb++)
		{
		    // use dh[1] as the edge length
		    if(nb >= 2)
			std::swap(dh[0],dh[1]);

		    if(I_nb[nb] >= 0)
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

		// interior
		rhs += m_dt*source_term.m_U[l];
		/*rhs += m_dt*cell_center[index].m_state.f_surf[l]; */

		// rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho *dh[0]*dh[1];
		rhs -= m_dt*cell_center[index].m_state.grad_q[l]/rho;
		rhs -= m_dt * cell_center[index].m_state.m_adv[l];

		solver.Add_A(I, I, 1.0);
		rhs += U_center;
		solver.Add_b(I, rhs);
	    }

	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);
	solver.Solve_GMRES();

	// get back the solution
	solver.Get_x(x);

	int num_iter;
	double rel_residual;
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);

	if (debugging("PETSc"))
	    (void) printf("Incompress_Solver_Smooth_2D_Cartesian::"
		    "compDiffWithSmoothProperty_2nd_decoupled: "
		    "num_iter = %d, rel_residual = %le. \n",
		    num_iter,rel_residual);

	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		I = ij_to_I[i][j];
		index = d_index2d(i,j,top_gmax);
		if (I >= 0)
		{
		    cell_center[index].m_state.m_U[l] = x[I-ilower];
		}
		else
		{
		    cell_center[index].m_state.m_U[l] = 0.0;
		}
	    }

	for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index  = d_index2d(i,j,top_gmax);
		array[index] = cell_center[index].m_state.m_U[l];
	    }
	scatMeshArray();
	for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
		index  = d_index2d(i,j,top_gmax);
		cell_center[index].m_state.m_U[l] = array[index];
	    }
    }

    for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    speed = fabs(cell_center[index].m_state.m_U[0]) +
		    fabs(cell_center[index].m_state.m_U[1]);
	    if (speed > max_speed)
		max_speed = speed;
	}
    pp_global_max(&max_speed,1);

    FT_FreeThese(1,x);
} /* end compDiffWithSmoothProperty_2nd_decoupled */


void Incompress_Solver_Smooth_2D_Cartesian::compAdvectionTerm_decoupled(void)
{
    int I;
    int i,j,icoords[MAXD];
    int index;
 
    setIndexMap();

	
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    I  = ij_to_I[i][j];
	    if (I == -1) continue;

	    index  = d_index2d(i,j,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;

	    double convectionTerm[2];
	    getAdvectionTerm_decoupled(icoords,convectionTerm);
	    cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	    cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	}
}


void Incompress_Solver_Smooth_2D_Cartesian::compAdvectionTerm_coupled(void)
{
    int I;
    int i,j,icoords[MAXD];
    int index;
 
    setIndexMap();

	
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    I  = ij_to_I[i][j];
	    if (I == -1) continue;

	    index  = d_index2d(i,j,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;

	    double convectionTerm[2];
	    getAdvectionTerm_coupled(icoords,convectionTerm);
	    cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	    cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	}
} /* end compAdvectionTerm_coupled */


void Incompress_Solver_Smooth_2D_Cartesian::compAdvectionTerm_decoupled_upgraded(void)
{
    int I;
    int i,j,icoords[MAXD];
    int index;
 
    setIndexMap();

	
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    I  = ij_to_I[i][j];
	    if (I == -1) continue;

	    index  = d_index2d(i,j,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;

	    double convectionTerm[2];
	    getAdvectionTerm_decoupled_upgraded(icoords,convectionTerm);
	    cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	    cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	}
} /* end compAdvectionTerm_decoupled_upgraded */


void Incompress_Solver_Smooth_2D_Cartesian::compAdvectionTerm_coupled_upgraded(void)
{
    int I;
    int i,j,icoords[MAXD];
    int index;
 
    setIndexMap();

	
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    I  = ij_to_I[i][j];
	    if (I == -1) continue;

	    index  = d_index2d(i,j,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;

	    double convectionTerm[2];
	    getAdvectionTerm_coupled_upgraded(icoords,convectionTerm);
	    cell_center[index].m_state.m_adv[0] = convectionTerm[0];
	    cell_center[index].m_state.m_adv[1] = convectionTerm[1];
	}
} /* end compAdvectionTerm_coupled_upgraded */


void Incompress_Solver_Smooth_2D_Cartesian::getAdvectionTerm_decoupled(
	int *icoords,
	double convectionTerm[2])
{
    bool bNoBoundary;
    int ICoords[2];
    L_STATE sl, sr, state_west, state_east, state_south, state_north;

    double dx = top_h[0];
    double dy = top_h[1];

    convectionTerm[0] = 0;
    convectionTerm[1] = 0;

    // WEST
    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
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
	getFaceVelocity_middleStep(ICoords,SOUTH,sr);
    }
    getFaceVelocity_middleStep(icoords,NORTH,sl);
    getRiemannSolution(COORD_Y,sl,sr,state_north);

    convectionTerm[0] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[0]-state_west.m_U[0])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[0]-state_south.m_U[0])/dy;
    convectionTerm[1] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[1]-state_west.m_U[1])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[1]-state_south.m_U[1])/dy;
} /* end getAdvectionTerm_decoupled */


void Incompress_Solver_Smooth_2D_Cartesian::getAdvectionTerm_coupled(
	int *icoords,
	double convectionTerm[2])
{
    bool bNoBoundary;
    int ICoords[2];
    L_STATE sl, sr, state_west, state_east, state_south, state_north;

    double dx = top_h[0];
    double dy = top_h[1];

    convectionTerm[0] = 0;
    convectionTerm[1] = 0;

    // WEST
    bNoBoundary = getNeighborOrBoundaryState(icoords,WEST,sl,m_t_int);
    if(bNoBoundary)
    {
	ICoords[0] = icoords[0] - 1;
	ICoords[1] = icoords[1];
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
	getFaceVelocity_middleStep_coupled(ICoords,SOUTH,sr);
    }
    getFaceVelocity_middleStep_coupled(icoords,NORTH,sl);
    getRiemannSolution(COORD_Y,sl,sr,state_north);

    convectionTerm[0] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[0]-state_west.m_U[0])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[0]-state_south.m_U[0])/dy;
    convectionTerm[1] =
	    1.0/2*(state_east.m_U[0]+state_west.m_U[0]) * (state_east.m_U[1]-state_west.m_U[1])/dx +
	    1.0/2*(state_north.m_U[1]+state_south.m_U[1]) * (state_north.m_U[1]-state_south.m_U[1])/dy;
} /* end getAdvectionTerm_coupled */


void Incompress_Solver_Smooth_2D_Cartesian::getAdvectionTerm_decoupled_upgraded(
	int *icoords,
	double convectionTerm[2])
{
    bool bNoBoundary;
    int ICoords[2];
    L_STATE sl, sr, state_west_bar, state_east_bar, state_south_bar, state_north_bar;

    L_STATE state_west_hat, state_east_hat, state_south_hat, state_north_hat;

    L_STATE state_west_hat_l, state_west_hat_r;
    L_STATE state_east_hat_l, state_east_hat_r;
    L_STATE state_south_hat_l, state_south_hat_r;
    L_STATE state_north_hat_l, state_north_hat_r;

    double transverseD[2];

    double dx = top_h[0];
    double dy = top_h[1];

    convectionTerm[0] = 0;
    convectionTerm[1] = 0;

    //////////////////////// Get the state_hat on four edges first //////////////////////////

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
	getFaceVelocity_middleStep_hat(ICoords,SOUTH,state_north_hat_r);
	getFaceVelocity_middleStep_hat(icoords,NORTH,state_north_hat_l);
	getRiemannSolution(COORD_Y,state_north_hat_l,state_north_hat_r,state_north_hat);

    }

    ///////////////////////////// get the state_bar on four edges ////////////////////////////////

    // WEST

    transverseD[0] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);
    transverseD[1] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);

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
	getFaceVelocity_middleStep_bar(ICoords,EAST,sl,transverseD,state_west_hat_l);
	getFaceVelocity_middleStep_bar(icoords,WEST,sr,transverseD,state_west_hat_r);
	getRiemannSolution(COORD_X,sl,sr,state_west_bar);
    }

    // EAST
    
    transverseD[0] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);
    transverseD[1] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);
    
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
	getFaceVelocity_middleStep_bar(ICoords,WEST,sr,transverseD,state_east_hat_r);
	getFaceVelocity_middleStep_bar(icoords,EAST,sl,transverseD,state_east_hat_l);
	getRiemannSolution(COORD_X,sl,sr,state_east_bar);
    }

    // SOUTH
    //
    transverseD[0] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[1] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    
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
	getFaceVelocity_middleStep_bar(ICoords,NORTH,sl,transverseD,state_south_hat_l);
	getFaceVelocity_middleStep_bar(icoords,SOUTH,sr,transverseD,state_south_hat_r);
	getRiemannSolution(COORD_Y,sl,sr,state_south_bar);
    }

    // NORTH

    transverseD[0] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[1] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    
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
	getFaceVelocity_middleStep_bar(ICoords,SOUTH,sr,transverseD,state_north_hat_r);
	getFaceVelocity_middleStep_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	getRiemannSolution(COORD_Y,sl,sr,state_north_bar);
    }

    convectionTerm[0] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[0]-state_west_bar.m_U[0])/dx +
	    1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[0]-state_south_bar.m_U[0])/dy;
    convectionTerm[1] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[1]-state_west_bar.m_U[1])/dx +
	    1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[1]-state_south_bar.m_U[1])/dy;
} /* end getAdvectionTerm_decoupled_upgraded */


void Incompress_Solver_Smooth_2D_Cartesian::getAdvectionTerm_coupled_upgraded(
	int *icoords,
	double convectionTerm[2])
{
    bool bNoBoundary;
    int ICoords[2];
    L_STATE sl, sr, state_west_bar, state_east_bar, state_south_bar, state_north_bar;

    L_STATE state_west_hat, state_east_hat, state_south_hat, state_north_hat;

    L_STATE state_west_hat_l, state_west_hat_r;
    L_STATE state_east_hat_l, state_east_hat_r;
    L_STATE state_south_hat_l, state_south_hat_r;
    L_STATE state_north_hat_l, state_north_hat_r;

    double transverseD[2];

    double dx = top_h[0];
    double dy = top_h[1];

    convectionTerm[0] = 0;
    convectionTerm[1] = 0;

    //////////////////////// Get the state_hat on four edges first //////////////////////////

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
	getFaceVelocity_middleStep_hat(ICoords,SOUTH,state_north_hat_r);
	getFaceVelocity_middleStep_hat(icoords,NORTH,state_north_hat_l);
	getRiemannSolution(COORD_Y,state_north_hat_l,state_north_hat_r,state_north_hat);
    }

    ///////////////////////////// get the state_bar on four edges ////////////////////////////////

    // WEST

    transverseD[0] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);
    transverseD[1] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);

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
	getFaceVelocity_middleStep_coupled_bar(ICoords,EAST,sl,transverseD,state_west_hat_l);
	getFaceVelocity_middleStep_coupled_bar(icoords,WEST,sr,transverseD,state_west_hat_r);
	getRiemannSolution(COORD_X,sl,sr,state_west_bar);
    }

    // EAST
    
    transverseD[0] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[0] - state_south_hat.m_U[0]);
    transverseD[1] = 1/(2.0*dy) * (state_south_hat.m_U[1] + state_north_hat.m_U[1]) * (state_north_hat.m_U[1] - state_south_hat.m_U[1]);
    
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
	getFaceVelocity_middleStep_coupled_bar(ICoords,WEST,sr,transverseD,state_east_hat_r);
	getFaceVelocity_middleStep_coupled_bar(icoords,EAST,sl,transverseD,state_east_hat_l);
	getRiemannSolution(COORD_X,sl,sr,state_east_bar);
    }

    // SOUTH
    //
    transverseD[0] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[1] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    
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
	getFaceVelocity_middleStep_coupled_bar(ICoords,NORTH,sl,transverseD,state_south_hat_l);
	getFaceVelocity_middleStep_coupled_bar(icoords,SOUTH,sr,transverseD,state_south_hat_r);
	getRiemannSolution(COORD_Y,sl,sr,state_south_bar);
    }

    // NORTH

    transverseD[0] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[0] - state_west_hat.m_U[0]);
    transverseD[1] = 1/(2.0*dx) * (state_west_hat.m_U[0] + state_east_hat.m_U[0]) * (state_east_hat.m_U[1] - state_west_hat.m_U[1]);
    
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
	getFaceVelocity_middleStep_coupled_bar(ICoords,SOUTH,sr,transverseD,state_north_hat_r);
	getFaceVelocity_middleStep_coupled_bar(icoords,NORTH,sl,transverseD,state_north_hat_l);
	getRiemannSolution(COORD_Y,sl,sr,state_north_bar);
    }

    convectionTerm[0] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[0]-state_west_bar.m_U[0])/dx +
	    1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[0]-state_south_bar.m_U[0])/dy;
    convectionTerm[1] =
	    1.0/2*(state_east_bar.m_U[0]+state_west_bar.m_U[0]) * (state_east_bar.m_U[1]-state_west_bar.m_U[1])/dx +
            1.0/2*(state_north_bar.m_U[1]+state_south_bar.m_U[1]) * (state_north_bar.m_U[1]-state_south_bar.m_U[1])/dy;
} /* end getAdvectionTerm_coupled_upgraded */


void Incompress_Solver_Smooth_2D_Cartesian::getFaceVelocity_middleStep(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_face)
{
    L_STATE state;
    int index;
    index = d_index2d(icoords[0],icoords[1],top_gmax);
    state = cell_center[index].m_state;
    double rho = cell_center[index].m_state.m_rho;


    double dx = 0, dy = 0, slope_x_limited[2] = {0,0}, slope_y_limited[2] = {0,0};

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
    default:
	assert(false);
    }

    state_face.m_U[0] = state.m_U[0];
    state_face.m_U[1] = state.m_U[1];

    //    return;

    getLimitedSlope(icoords,COORD_X,slope_x_limited);
    getLimitedSlope(icoords,COORD_Y,slope_y_limited);

    // dx/2, dy/2
    state_face.m_U[0] += dx/2 * slope_x_limited[0] + dy/2 * slope_y_limited[0];
    state_face.m_U[1] += dx/2 * slope_x_limited[1] + dy/2 * slope_y_limited[1];

    //    return;
    // dt/2
    double diffusion[2];
    double gradP[2];
    gradP[0] = cell_center[index].m_state.grad_q[0];
    gradP[1] = cell_center[index].m_state.grad_q[1];

    getDifffusion(icoords,diffusion);
    state_face.m_U[0] += m_dt/2 *
	    (diffusion[0]/rho - state.m_U[0]*slope_x_limited[0] - state.m_U[1]*slope_y_limited[0] - gradP[0]/rho);
    state_face.m_U[1] += m_dt/2 *
	    (diffusion[1]/rho - state.m_U[0]*slope_x_limited[1] - state.m_U[1]*slope_y_limited[1] - gradP[1]/rho);

    //    return;
    // rhs
    double coords[2];
    L_STATE source_term;

    getRectangleCenter(index, coords);
    computeSourceTerm_Adv(coords, source_term);

    state_face.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_face.m_U[1] += m_dt/2 * source_term.m_U[1];
}


void Incompress_Solver_Smooth_2D_Cartesian::getFaceVelocity_middleStep_coupled(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_face)
{
    L_STATE state;
    int index;
    index = d_index2d(icoords[0],icoords[1],top_gmax);
    state = cell_center[index].m_state;
    double rho = cell_center[index].m_state.m_rho;


    double dx = 0, dy = 0, slope_x_limited[2] = {0,0}, slope_y_limited[2] = {0,0};

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
    default:
	assert(false);
    }

    state_face.m_U[0] = state.m_U[0];
    state_face.m_U[1] = state.m_U[1];

    //    return;

    getLimitedSlope(icoords,COORD_X,slope_x_limited);
    getLimitedSlope(icoords,COORD_Y,slope_y_limited);

    // dx/2, dy/2
    state_face.m_U[0] += dx/2 * slope_x_limited[0] + dy/2 * slope_y_limited[0];
    state_face.m_U[1] += dx/2 * slope_x_limited[1] + dy/2 * slope_y_limited[1];

    //    return;
    // dt/2
    double diffusion[2];
    double gradP[2];
    gradP[0] = cell_center[index].m_state.grad_q[0];
    gradP[1] = cell_center[index].m_state.grad_q[1];

    getDiffusion_coupled(icoords,diffusion);
    state_face.m_U[0] += m_dt/2 *
	    (diffusion[0]/rho - state.m_U[0]*slope_x_limited[0] - state.m_U[1]*slope_y_limited[0] - gradP[0]/rho);
    state_face.m_U[1] += m_dt/2 *
	    (diffusion[1]/rho - state.m_U[0]*slope_x_limited[1] - state.m_U[1]*slope_y_limited[1] - gradP[1]/rho);

    //    return;
    // rhs
    double coords[2];
    L_STATE source_term;

    getRectangleCenter(index, coords);
    computeSourceTerm_Adv(coords, source_term);

    state_face.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_face.m_U[1] += m_dt/2 * source_term.m_U[1];
}


void Incompress_Solver_Smooth_2D_Cartesian::getFaceVelocity_middleStep_hat(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_hat)
{
    L_STATE  state_orig;
    int index;
    index = d_index2d(icoords[0],icoords[1],top_gmax);
    state_orig = cell_center[index].m_state;
    double sL;

    double dx = 0, dy = 0, slope_x_limited[2] = {0,0}, slope_y_limited[2] = {0,0};


    //    return the slope limiter, needed to be modified;

    //getLimitedSlope(icoords,COORD_X,slope_x_limited);
    //getLimitedSlope(icoords,COORD_Y,slope_y_limited);

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
	break;
    default:
	assert(false);
    }
}


void Incompress_Solver_Smooth_2D_Cartesian::getFaceVelocity_middleStep_bar(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_bar,
	double transverseD[2],
	L_STATE state_hat)
{
    int index;
    index = d_index2d(icoords[0],icoords[1],top_gmax);
    double rho = cell_center[index].m_state.m_rho;

    double diffusion[2], gradP[2];
    gradP[0] = cell_center[index].m_state.grad_q[0];
    gradP[1] = cell_center[index].m_state.grad_q[1];


    getDifffusion(icoords,diffusion);
    state_bar.m_U[0] = state_hat.m_U[0] + m_dt/2.0 * (-transverseD[0] + diffusion[0]/rho - gradP[0]/rho);
    state_bar.m_U[1] = state_hat.m_U[1] + m_dt/2.0 * (-transverseD[1] + diffusion[1]/rho - gradP[1]/rho);

    double coords[2];
    L_STATE source_term;

    getRectangleCenter(index, coords);
    computeSourceTerm_Adv(coords, source_term);

    state_bar.m_U[0] += m_dt/2 * source_term.m_U[0];
    state_bar.m_U[1] += m_dt/2 * source_term.m_U[1];
}

void Incompress_Solver_Smooth_2D_Cartesian::getFaceVelocity_middleStep_coupled_bar(
	int *icoords,
	GRID_DIRECTION dir,
	L_STATE &state_bar,
	double transverseD[2],
	L_STATE state_hat)
{
    int index;
    index = d_index2d(icoords[0],icoords[1],top_gmax);
    double rho = cell_center[index].m_state.m_rho;

    double diffusion[2], gradP[2];
    gradP[0] = cell_center[index].m_state.grad_q[0];
    gradP[1] = cell_center[index].m_state.grad_q[1];


    getDiffusion_coupled(icoords,diffusion);
    state_bar.m_U[0] = state_hat.m_U[0] + m_dt/2.0 * (-transverseD[0] + diffusion[0]/rho - gradP[0]/rho);
    state_bar.m_U[1] = state_hat.m_U[1] + m_dt/2.0 * (-transverseD[1] + diffusion[1]/rho - gradP[1]/rho);

    double coords[2];
    L_STATE source_term;

    getRectangleCenter(index, coords);
    computeSourceTerm_Adv(coords, source_term);

    state_bar.m_U[0] += m_dt/2 * source_term.m_U[0]/rho;
    state_bar.m_U[1] += m_dt/2 * source_term.m_U[1]/rho;
}


/**
* compute mu * (Uxx+Uyy)
* @param icoords
* @param diffusion
* @param gradP
*/
void Incompress_Solver_Smooth_2D_Cartesian::getDiffusion_coupled(
	int *icoords,
	double diffusion[2])
{
    int index,index_nb[8];
    double mu[4],mu_edge[4],mu0;
    L_STATE Unb,corner_state;
    double U0_nb[8],U1_nb[8],U0_center, U1_center;
    int nb;
    GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
    bool bNoBoundary[4];
    double dh[2],dh0[2],dh1[2];
    double coords[MAXD],corner_coords[MAXD];

    diffusion[0] = 0.0;
    diffusion[1] = 0.0;

    int i = icoords[0];
    int j = icoords[1];
    index = d_index2d(i,j,top_gmax);
    index_nb[0] = d_index2d(i-1,j,top_gmax);
    index_nb[1] = d_index2d(i+1,j,top_gmax);
    index_nb[2] = d_index2d(i,j-1,top_gmax);
    index_nb[3] = d_index2d(i,j+1,top_gmax);

    index_nb[4] = d_index2d(i-1,j-1,top_gmax);
    index_nb[5] = d_index2d(i+1,j-1,top_gmax);
    index_nb[6] = d_index2d(i+1,j+1,top_gmax);
    index_nb[7] = d_index2d(i-1,j+1,top_gmax);

    mu0 = cell_center[index].m_state.m_mu;

    U0_center = cell_center[index].m_state.m_U[0];
    U1_center = cell_center[index].m_state.m_U[1];
   

    for (nb = 0; nb < 4; nb++)
    {
	bNoBoundary[nb] = getNeighborOrBoundaryState(icoords,dir[nb],Unb,m_t_old);
	U0_nb[nb] = Unb.m_U[0];
	U1_nb[nb] = Unb.m_U[1];
	if(!bNoBoundary[nb])
	{
	    mu[nb] = mu0;
	    mu_edge[nb] = mu0;
	}
	else
	{
	    mu[nb] = cell_center[index_nb[nb]].m_state.m_mu;
	    mu_edge[nb] = 0.5*(cell_center[index_nb[nb]].m_state.m_mu + mu0);
	}
    }


    // non-cross derivative terms

    dh[0] = top_h[0];
    dh[1] = top_h[1];
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

    diffusion[0] += 2.0*(mu_edge[1]*(U0_nb[1]-U0_center)/dh1[0] - mu_edge[0]*(U0_center-U0_nb[0])/dh0[0])/dh[0];    // first diffusion coefficient (2*mu*u_x)_x
    diffusion[1] +=     (mu_edge[1]*(U1_nb[1]-U1_center)/dh1[0] - mu_edge[0]*(U1_center-U1_nb[0])/dh0[0])/dh[0];

    diffusion[0] +=     (mu_edge[3]*(U0_nb[3]-U0_center)/dh1[1] - mu_edge[2]*(U0_center-U0_nb[2])/dh0[1])/dh[1];
    diffusion[1] += 2.0*(mu_edge[3]*(U1_nb[3]-U1_center)/dh1[1] - mu_edge[2]*(U1_center-U1_nb[2])/dh0[1])/dh[1];

    // get the coords in the cell center
    getRectangleCenter(index, coords);

    //cross derivative terms

    //traverse the four corners to get corner values

    //corner (i-1/2,j-1/2)
    /*
    if (!bNoBoundary[0] && bNoBoundary[2])
    {
        U0_nb[4] = U0_nb[0];
        U1_nb[4] = U1_nb[0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[2])
    {
        U0_nb[4] = U0_nb[2];
        U1_nb[4] = U1_nb[2];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[2])
    {
        U0_nb[4] = U0_nb[0];
        U1_nb[4] = U1_nb[0];
    }
    else
    {
        U0_nb[4] = (U0_nb[0]+U0_nb[2]+cell_center[index_nb[4]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[4] = (U1_nb[0]+U1_nb[2]+cell_center[index_nb[4]].m_state.m_U[1]+U1_center)/4.0;
    } */
    if (bNoBoundary[0] && bNoBoundary[2])
    {
        U0_nb[4] = (U0_nb[0]+U0_nb[2]+cell_center[index_nb[4]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[4] = (U1_nb[0]+U1_nb[2]+cell_center[index_nb[4]].m_state.m_U[1]+U1_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] - 0.5*top_h[0];
        corner_coords[1] = coords[1] - 0.5*top_h[1];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[4] = corner_state.m_U[0];
        U1_nb[4] = corner_state.m_U[1];
    }
    //corner (i+1/2,j-1/2)
    /*
    if (!bNoBoundary[1] && bNoBoundary[2])
    {
        U0_nb[5] = U0_nb[1];
        U1_nb[5] = U1_nb[1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[2])
    {
        U0_nb[5] = U0_nb[2];
        U1_nb[5] = U1_nb[2];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[2])
    {
        U0_nb[5] = U0_nb[1];
        U1_nb[5] = U1_nb[1];
    }
    else
    {
        U0_nb[5] = (U0_nb[1]+U0_nb[2]+cell_center[index_nb[5]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[5] = (U1_nb[1]+U1_nb[2]+cell_center[index_nb[5]].m_state.m_U[1]+U1_center)/4.0;
    } */
    if (bNoBoundary[1] && bNoBoundary[2])
    {
        U0_nb[5] = (U0_nb[1]+U0_nb[2]+cell_center[index_nb[5]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[5] = (U1_nb[1]+U1_nb[2]+cell_center[index_nb[5]].m_state.m_U[1]+U1_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] + 0.5*top_h[0];
        corner_coords[1] = coords[1] - 0.5*top_h[1];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[5] = corner_state.m_U[0];
        U1_nb[5] = corner_state.m_U[1];
    }

    //corner (i+1/2,j+1/2)
    /*
    if (!bNoBoundary[1] && bNoBoundary[3])
    {
        U0_nb[6] = U0_nb[1];
        U1_nb[6] = U1_nb[1];
    }
    else if(bNoBoundary[1] && !bNoBoundary[3])
    {
        U0_nb[6] = U0_nb[3];
        U1_nb[6] = U1_nb[3];
    }
    else if(!bNoBoundary[1] && !bNoBoundary[3])
    {
        U0_nb[6] = U0_nb[1];
        U1_nb[6] = U1_nb[1];
    }
    else
    {
        U0_nb[6] = (U0_nb[1]+U0_nb[3]+cell_center[index_nb[6]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[6] = (U1_nb[1]+U1_nb[3]+cell_center[index_nb[6]].m_state.m_U[1]+U1_center)/4.0;
    } */
    if (bNoBoundary[1] && bNoBoundary[3])
    {
        U0_nb[6] = (U0_nb[1]+U0_nb[3]+cell_center[index_nb[6]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[6] = (U1_nb[1]+U1_nb[3]+cell_center[index_nb[6]].m_state.m_U[1]+U1_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] + 0.5*top_h[0];
        corner_coords[1] = coords[1] + 0.5*top_h[1];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[6] = corner_state.m_U[0];
        U1_nb[6] = corner_state.m_U[1];
    }

    //corner (i-1/2,j+1/2)
    /*
    if (!bNoBoundary[0] && bNoBoundary[3])
    {
        U0_nb[7] = U0_nb[0];
        U1_nb[7] = U1_nb[0];
    }
    else if(bNoBoundary[0] && !bNoBoundary[3])
    {
        U0_nb[7] = U0_nb[3];
        U1_nb[7] = U1_nb[3];
    }
    else if(!bNoBoundary[0] && !bNoBoundary[3])
    {
        U0_nb[7] = U0_nb[0];
        U1_nb[7] = U1_nb[0];
    }
    else
    {
        U0_nb[7] = (U0_nb[0]+U0_nb[3]+cell_center[index_nb[7]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[7] = (U1_nb[0]+U1_nb[3]+cell_center[index_nb[7]].m_state.m_U[1]+U1_center)/4.0;
    } */
    if (bNoBoundary[0] && bNoBoundary[3])
    {
        U0_nb[7] = (U0_nb[0]+U0_nb[3]+cell_center[index_nb[7]].m_state.m_U[0]+U0_center)/4.0;
        U1_nb[7] = (U1_nb[0]+U1_nb[3]+cell_center[index_nb[7]].m_state.m_U[1]+U1_center)/4.0;
    }
    else
    {
        corner_coords[0] = coords[0] - 0.5*top_h[0];
        corner_coords[1] = coords[1] + 0.5*top_h[1];
        getExactSolution(corner_coords,m_t_old,corner_state);
        U0_nb[7] = corner_state.m_U[0];
        U1_nb[7] = corner_state.m_U[1];
    }

    diffusion[0] += (mu_edge[2]*U1_nb[4]-mu_edge[2]*U1_nb[5]+mu_edge[3]*U1_nb[6]-mu_edge[3]*U1_nb[7])/(top_h[0]*top_h[1]); // (mu*v_x)_y
    diffusion[1] += (mu_edge[0]*U0_nb[4]-mu_edge[1]*U0_nb[5]+mu_edge[1]*U0_nb[6]-mu_edge[0]*U0_nb[7])/(top_h[0]*top_h[1]); // (mu*u_y)_x
} /* end getDiffusion_coupled */


void Incompress_Solver_Smooth_2D_Cartesian::getDifffusion(
	int *icoords,
	double diffusion[2])
{
    double Uxx[2], Uyy[2];
    getDU2(icoords,COORD_X,Uxx);
    getDU2(icoords,COORD_Y,Uyy);

    double mu = cell_center[d_index2d(icoords[0],icoords[1],top_gmax)].m_state.m_mu;

    diffusion[0] = mu * (Uxx[0] + Uyy[0]);
    diffusion[1] = mu * (Uxx[1] + Uyy[1]);
}

/**
* calculate Uxx or Uyy and Px or Py.
* @param dir
* @param icoords
* @param dU2
* @param dP
*/
void Incompress_Solver_Smooth_2D_Cartesian::getDU2(
	int *icoords,
	EBM_COORD xyz,
	double dU2[2])
{
    double dh0, dh1, dh;
    L_STATE U0, U1, U2;

    U1 = cell_center[d_index2d(icoords[0],icoords[1],top_gmax)].m_state;

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
    else	//
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

    dU2[0] = ((U2.m_U[0] - U1.m_U[0])/dh1 - (U1.m_U[0] - U0.m_U[0])/dh0) / dh;
    dU2[1] = ((U2.m_U[1] - U1.m_U[1])/dh1 - (U1.m_U[1] - U0.m_U[1])/dh0) / dh;
}

void Incompress_Solver_Smooth_2D_Cartesian::getLimitedSlope(
	int *icoords,
	EBM_COORD xyz,
	double slope[2])
{
    double dh0, dh1;
    L_STATE U0, U1, U2;

    U1 = cell_center[d_index2d(icoords[0],icoords[1],top_gmax)].m_state;

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
    else	//
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
    slope[0] = EBM_minmod((U1.m_U[0]-U0.m_U[0])/dh0, (U2.m_U[0]-U1.m_U[0])/dh1);
    slope[1] = EBM_minmod((U1.m_U[1]-U0.m_U[1])/dh0, (U2.m_U[1]-U1.m_U[1])/dh1);
}



void Incompress_Solver_Smooth_2D_Cartesian::getLimitedSlope_Vanleer(
	int *icoords,
	EBM_COORD xyz,
	double slope[2])
{
    double dx,dy;
    L_STATE U0, U1, U2;
    double u_lim,  v_lim ;
    double u_slope, v_slope;

    int index = d_index2d(icoords[0],icoords[1],top_gmax);
    U1 = cell_center[index].m_state;
    dx = top_h[0];
    dy = top_h[1];

    bool bNoBoundary[2];

    if(xyz==COORD_X)
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,WEST,U0,m_t_old);
	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,EAST,U2,m_t_old);
	if(bNoBoundary[0] || bNoBoundary[1])
	{

	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dx, 2.0*(U1.m_U[0]-U0.m_U[0])/dx);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dx, 2.0*(U1.m_U[1]-U0.m_U[1])/dx);

	    u_slope = (U2.m_U[0] - U0.m_U[0])/(2*dx);
	    v_slope = (U2.m_U[1] - U0.m_U[1])/(2*dx);
	    

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);

	}
	else if ( (!bNoBoundary[0]) || bNoBoundary[1])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dx, 2.0*(U1.m_U[0]-U0.m_U[0])/(dx/2.0));
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dx, 2.0*(U1.m_U[1]-U0.m_U[1])/(dx/2.0));

	    u_slope = (U2.m_U[0] + U1.m_U[0] + 2.0*U0.m_U[0])/(2.0*dx);
	    v_slope = (U2.m_U[1] + U1.m_U[1] + 2.0*U0.m_U[1])/(2.0*dx);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	}
	else if ( (!bNoBoundary[1]) || bNoBoundary[0])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/(dx/2.0), 2.0*(U1.m_U[0]-U0.m_U[0])/dx);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/(dx/2.0), 2.0*(U1.m_U[1]-U0.m_U[1])/dx);

	    u_slope = (2.0*U2.m_U[0] - U1.m_U[0] - U0.m_U[0])/(2.0*dx);
	    v_slope = (2.0*U2.m_U[1] - U1.m_U[1] - U0.m_U[1])/(2.0*dx);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);

	}
	else
	{
	    printf("\nThe number of cell in x direction is less than 1!!\n");
	    slope[0] = slope[1] = 0.0;
	}

    }
    else	//
    {
	bNoBoundary[0] = getNeighborOrBoundaryState(icoords,SOUTH,U0,m_t_old);
	bNoBoundary[1] = getNeighborOrBoundaryState(icoords,NORTH,U2,m_t_old);

	if(bNoBoundary[0] || bNoBoundary[1])
	{

	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dy, 2.0*(U1.m_U[0]-U0.m_U[0])/dy);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dy, 2.0*(U1.m_U[1]-U0.m_U[1])/dy);

	    u_slope = (U2.m_U[0] - U0.m_U[0])/(2*dy);
	    v_slope = (U2.m_U[1] - U0.m_U[1])/(2*dy);
	    

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);

	}
	else if ( (!bNoBoundary[0]) || bNoBoundary[1])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/dy, 2.0*(U1.m_U[0]-U0.m_U[0])/(dy/2.0));
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/dy, 2.0*(U1.m_U[1]-U0.m_U[1])/(dy/2.0));

	    u_slope = (U2.m_U[0] + U1.m_U[0] + 2.0*U0.m_U[0])/(2.0*dy);
	    v_slope = (U2.m_U[1] + U1.m_U[1] + 2.0*U0.m_U[1])/(2.0*dy);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);
	}
	else if ( (!bNoBoundary[1]) || bNoBoundary[0])
	{
	    u_lim = EBM_minmod(2.0*(U2.m_U[0]-U1.m_U[0])/(dy/2.0), 2.0*(U1.m_U[0]-U0.m_U[0])/dy);
	    v_lim = EBM_minmod(2.0*(U2.m_U[1]-U1.m_U[1])/(dy/2.0), 2.0*(U1.m_U[1]-U0.m_U[1])/dy);

	    u_slope = (2.0*U2.m_U[0] - U1.m_U[0] - U0.m_U[0])/(2.0*dy);
	    v_slope = (2.0*U2.m_U[1] - U1.m_U[1] - U0.m_U[1])/(2.0*dy);

	    slope[0] = EBM_minmod(u_slope, u_lim);
	    slope[1] = EBM_minmod(v_slope, v_lim);

	}
	else
	{
	    printf("\nThe number of cell in y direction is less than 1!!\n");
	    slope[0] = slope[1] = 0.0;
	}
    }
 }


double Incompress_Solver_Smooth_2D_Cartesian::EBM_minmod(
	double x,
	double y)
{
    //	return x;
    //	static double tol = 1e-5;
    //	static double tol = -1;
    double sign = x*y;

    //	if(fabs(sign)>tol)
    //	{
    if(sign<0)
	return 0;
    else if(sign>=0)
    {
	if(fabs(x)<fabs(y))
	    return x;
	else
	    return y;
    }
    //	}
    //	else
    //	{
    //		if(fabs(x)<fabs(y))
    //			return y;
    //		else
    //			return x;
    //	}

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
bool Incompress_Solver_Smooth_2D_Cartesian::getNeighborOrBoundaryState(
	int icoords[2],
	GRID_DIRECTION dir,
	L_STATE &state,
	double t)
{
    double crx_coords[MAXD];
    static double (*getStateVel[2])(POINTER) = {getStateXvel,getStateYvel};
    POINTER intfc_state;
    HYPER_SURF *hs;

    int index = d_index2d(icoords[0],icoords[1],top_gmax);
    int comp = cell_center[index].comp;

    if (FT_StateStructAtGridCrossing_tmp(front,icoords,dir,
	    comp,&intfc_state,&hs,crx_coords,t) &&
	    wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
    {
	state.m_U[0] = getStateVel[0](intfc_state);
	state.m_U[1] = getStateVel[1](intfc_state);
	return false;
    }
    else
    {
	int index_nb;
	switch(dir)
	{
	case WEST:
	    index_nb = d_index2d(icoords[0]-1,icoords[1],top_gmax);
	    break;
	case EAST:
	    index_nb = d_index2d(icoords[0]+1,icoords[1],top_gmax);
	    break;
	case SOUTH:
	    index_nb = d_index2d(icoords[0],icoords[1]-1,top_gmax);
	    break;
	case NORTH:
	    index_nb = d_index2d(icoords[0],icoords[1]+1,top_gmax);
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
void Incompress_Solver_Smooth_2D_Cartesian::getRiemannSolution(
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

	sr.m_U[0] = state_right.m_U[0];
	sr.m_U[1] = state_right.m_U[1];
    }
    else
    {
	sl.m_U[0] = state_left.m_U[1];
	sl.m_U[1] = state_left.m_U[0];

	sr.m_U[0] = state_right.m_U[1];
	sr.m_U[1] = state_right.m_U[0];
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

    // rotate state
    if(xyz==COORD_X)
	; // do nothing
    else
	std::swap(ans.m_U[0],ans.m_U[1]);
}


void Incompress_Solver_Smooth_2D_Cartesian::computeSubgridModel(void)
{
        int i,j,k,index,index0,index1,index2,index3,index4,size;  
        L_STATE state;
        double *u, *v;
        double u0x,ulx,urx,vlx,vrx;
        double u0y,uly,ury,vly,vry;
        double ux,uy,vx,vy;
        double s, *s11, *s12, *s22;
        double *ss11, *ss12, *ss22;
        double *tau00, *tau01, *tau10, *tau11;
        double *vel_u, *vel_v, *vel_uu, *vel_uv, *vel_vv;
        double sum_vel_u,sum_vel_v,sum_vel_uu,sum_vel_uv,sum_vel_vv;
        double sum_s11,sum_s12,sum_s22,sum_ss11,sum_ss12,sum_ss22,sum_s;
        double *ma11, *ma12, *la11, *la12, *la22;
        double *cs, *cs_ave, *deno, *nume, *co_coords_y;
        double coords[2];
        int    *r, num_r;
        int    ii,jj,iii,jjj;
        const int nn = pp_numnodes();
        num_r = (int)(((top_U[1]-top_L[1])/top_h[1])+1);

        size = (top_gmax[0]+1)*(top_gmax[1]+1);
        FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s22,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss22,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau00,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau01,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau10,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_u,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_v,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_uu,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_uv,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_vv,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&co_coords_y,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ma11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ma12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&la11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&la12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&la22,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&r,size,sizeof(int));
        FT_VectorMemoryAlloc((POINTER*)&cs,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&cs_ave,num_r,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&deno,num_r,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&nume,num_r,sizeof(double));

        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index2d(i,j,top_gmax);
            u[index] = cell_center[index].m_state.m_U[0];
            v[index] = cell_center[index].m_state.m_U[1];
            getRectangleCenter(index, coords);
            co_coords_y[index] = coords[1] + (top_h[1]/2.0);
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin-2; i <= imax+2; i++)
        {
            index0  = d_index2d(i,j,top_gmax);
            index1  = d_index2d(i-1,j,top_gmax);
            ulx = u[index1];
            uly = v[index1];
            index2  = d_index2d(i+1,j,top_gmax);
            urx = u[index2];
            ury = v[index2];
            index3  = d_index2d(i,j-1,top_gmax);
            vlx = u[index3];
            vly = v[index3];
            index4  = d_index2d(i,j+1,top_gmax);
            vrx = u[index4];
            vry = v[index4];

            ux = (urx - ulx) / (2.0*top_h[0]);
            uy = (ury - uly) / (2.0*top_h[1]);
            vx = (vrx - vlx) / (2.0*top_h[0]);
            vy = (vry - vly) / (2.0*top_h[1]);
            s11[index0] = ux;
            s12[index0] = (uy + vx)/2;
            s22[index0] = vy;
            s = sqrt(2*( (s11[index0]*s11[index0]) 
                         + (2*(s12[index0]*s12[index0]))
                         + (s22[index0]*s22[index0])));
            ss11[index0] = s*s11[index0];
            ss12[index0] = s*s12[index0];
            ss22[index0] = s*s22[index0];
            vel_u[index0] = u[index0];
            vel_v[index0] = v[index0];  
            vel_uu[index0] = u[index0]*u[index0]; 
            vel_uv[index0] = u[index0]*v[index0];  
            vel_vv[index0] = v[index0]*v[index0];      
        }

        for (j = jmin; j <= (jmax/2); j++)
        {
            jj = (2*j)-1;
            for (i = imin-1; i <= (imax/2)+1; i++)
            {
                ii = (2*i)-4;
                index = d_index2d(ii,jj,top_gmax);
                sum_vel_u = sum_vel_v = 0.0;
                sum_vel_uu = sum_vel_uv = sum_vel_vv = 0.0;
                sum_s11 = sum_s12 = sum_s22 = 0.0;
                sum_ss11 = sum_ss12 = sum_ss22 = sum_s = 0.0;
                for(jjj = jj; jjj < jj+2; jjj++)
                for(iii = ii; iii < ii+2; iii++)
                {
                    index0  = d_index2d(iii,jjj,top_gmax);
                    sum_vel_u += vel_u[index0];
                    sum_vel_v += vel_v[index0];
                    sum_vel_uu += vel_uu[index0];
                    sum_vel_uv += vel_uv[index0];
                    sum_vel_vv += vel_vv[index0];
                    sum_s11 += s11[index0];
                    sum_s12 += s12[index0];
                    sum_s22 += s22[index0];
                    sum_ss11 += ss11[index0];
                    sum_ss12 += ss12[index0];
                    sum_ss22 += ss22[index0];
                    sum_s += sqrt(2*( (s11[index0]*s11[index0]) 
                                  + (2*(s12[index0]*s12[index0])) 
                                  + (s22[index0]*s22[index0])));
                } 
                ma11[index] = (2.0*top_h[1]*top_h[1]*(sum_ss11/4.0))
                        - (2.0*4*top_h[1]*top_h[1]*(sum_s/4.0)*(sum_s11/4.0));
                ma12[index] = (2.0*top_h[1]*top_h[1]*(sum_ss12/4.0))
                        - (2.0*4*top_h[1]*top_h[1]*(sum_s/4.0)*(sum_s12/4.0));
                la11[index] = (sum_vel_uu/4.0)-((sum_vel_u/4.0)*
			(sum_vel_u/4.0));
                la12[index] = (sum_vel_uv/4.0)-((sum_vel_u/4.0)*
			(sum_vel_v/4.0));
                la12[index] = (sum_vel_vv/4.0)-((sum_vel_v/4.0)*
			(sum_vel_v/4.0));
                r[index] = (int)(co_coords_y[index]/(2*top_h[1]));
            }
        }

        for (k = 0; k < num_r; k++)
        {
            deno[k] = 0.0;
            nume[k] = 0.0;
        }

        for (k = 0; k < num_r; k++)
        for (j = jmin; j <= (jmax/2); j++)
        {
            jj = (2*j)-1;
            for (i = imin-1; i <= (imax/2)+1; i++)
            {
                ii = (2*i)-4;
                index0 = d_index2d(ii,jj,top_gmax);
                if(k == r[index0])
                {
                    deno[k] += (ma11[index0]*ma11[index0]) + 
				(ma12[index0]*ma12[index0]);
                    nume[k] += (((la11[index0]/2.0)-(la22[index0]/2.0))*
				ma11[index0]) + (la12[index0]*ma12[index0]);
                }
            }
        }

        pp_gsync();
        
        if (nn > 1)
        {
           for (k = 0; k < num_r; k++)
           {
              pp_global_sum(&deno[k],1L);
              pp_global_sum(&nume[k],1L);
           }
        }

        for (k = 0; k < num_r; k++)
        {
            if(deno[k] < 10e-16)
                cs_ave[k] = 0.0;
            else
                cs_ave[k] = nume[k]/deno[k];
        }

        for (j = jmin; j <= (jmax/2); j++)
        {
            jj = (2*j)-1;
            for (i = imin-1; i <= (imax/2)+1; i++)
            {
                ii = (2*i)-4;
                index = d_index2d(ii,jj,top_gmax);
                for(jjj = jj; jjj < jj+2; jjj++)
                for(iii = ii; iii < ii+2; iii++)
                {
                    index0 = d_index2d(iii,jjj,top_gmax);
                    cs[index0] = cs_ave[r[index]];
                }
            }
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin-1; i <= imax+1; i++)
        {
            index0  = d_index2d(i,j,top_gmax);
            s = sqrt(2*( (s11[index0]*s11[index0]) 
                          + (2*(s12[index0]*s12[index0]))
                          + (s22[index0]*s22[index0])));
            tau00[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*
                              s*((s11[index0]/2.0)-(s22[index0]/2.0));
            tau01[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*s*(s12[index0]);
            tau10[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*s*(s12[index0]);
            tau11[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*
                              s*((s22[index0]/2.0)-(s11[index0]/2.0));
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index0 = d_index2d(i,j,top_gmax);
            index1 = d_index2d(i-1,j,top_gmax);
            index2 = d_index2d(i+1,j,top_gmax);
            index3 = d_index2d(i,j-1,top_gmax);
            index4 = d_index2d(i,j+1,top_gmax);

            cell_center[index0].m_state.m_U[0] += -m_dt*(
                              ((tau00[index2]-tau00[index1])/(2.0*top_h[0])) + 
                                ((tau01[index4]-tau01[index3])/(2.0*top_h[1])));
            cell_center[index0].m_state.m_U[1] += -m_dt*(
                              ((tau10[index2]-tau10[index1])/(2.0*top_h[0])) + 
                              ((tau11[index4]-tau11[index3])/(2.0*top_h[1])));
        }
        FT_FreeThese(2,u,v);
        FT_FreeThese(4,tau00,tau01,tau10,tau11);
        FT_FreeThese(6,s11,s12,s22,ss11,ss12,ss22);
        FT_FreeThese(5,vel_u,vel_v,vel_uu,vel_uv,vel_vv);
        FT_FreeThese(11,co_coords_y,ma11,ma12,la11,la12,la22,r,cs,cs_ave,
					deno,nume);
}       /* end 2D_Cartesian::computeSubgridModel */
