/**
* iFdebug2d.cpp
*
*  Created on: Jun 23, 2010
*      Author: shuqiangwang & Yijie Zhou
*/
#include "iFluid_debug.h"
#include "ifluid_basic.h"
#include "solver.h"

/*****************************************************************************
*			Incompress_Solver_Smooth_2D_Cartesian_Debug
*
* 2D testing of the algorithm using the exact solution
*****************************************************************************/
Incompress_Solver_Smooth_2D_Cartesian_Debug::Incompress_Solver_Smooth_2D_Cartesian_Debug(Front &front)
: Incompress_Solver_Smooth_2D_Cartesian(front)
{
    m_bStartDebugging = false;
}

/**
* see also:
* Incompress_Solver_Smooth_2D_Cartesian::setInitialCondition
*/
void Incompress_Solver_Smooth_2D_Cartesian_Debug::setInitialCondition(void)
{
    // Incompress_Solver_Smooth_2D_Cartesian::setInitialCondition();
    printf("\nEnter 2D_Debug::setInitialCondition()\n");
    int i,j,l;
    int index;
    COMPONENT comp;
    double coords[MAXD];
        int reflect[MAXD];
        reflect[0] = reflect[1] = YES;

    FT_MakeGridIntfc(front);
    setDomain();

    m_rho[0] = 1;
    m_rho[1] = 1;
    m_mu[0] = 1;
    m_mu[1] = 1;
    m_comp[0] = iFparams->m_comp1;
    m_comp[1] = iFparams->m_comp2;
    m_smoothing_radius = iFparams->smoothing_radius;
    m_sigma = 0;
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

    double CFL = 0.1;
    front->max_time = 0.5;
    double dh;
    dh = std::min(front->rect_grid->h[0], front->rect_grid->h[1]);
    double dt = CFL * dh;
    printf("\nBefore adjusting, the dt = %20.16g\n",dt);
    int nStep = int (front->max_time/dt) + 1;
    dt = front->max_time/ nStep;
    printf("\nAfter adjusting, the dt = %20.16g\n",dt);
    front->dt = dt;
    for (i = 0; i < cell_center.size(); i++)
    {
	getRectangleCenter(i, coords);
	cell_center[i].m_state.setZero();
	getExactSolution(coords,front->time,cell_center[i].m_state);
	cell_center[i].m_state.m_q = cell_center[i].m_state.m_P;
        printf( "\t#### i = %3d\n", i);
        printf("\t#### m_U = {%12.8f, %12.8f}\n", cell_center[i].m_state.m_U[0],cell_center[i].m_state.m_U[1]);
        printf("\t#### m_P = {%12.8f}\n", cell_center[i].m_state.m_P);
        printf("\t#### m_q = {%12.8f}\n", cell_center[i].m_state.m_q);
    }
    front->dt = 0.0;

    for (l = 0; l < dim; l++)
    {
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = cell_center[index].m_state.m_U[l];
	}
	scatMeshArray(reflect);
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    cell_center[index].m_state.m_U[l] = array[index];
	}
    }

    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
	index = d_index2d(i,j,top_gmax);
	array[index] = cell_center[index].m_state.m_P;
    }
    scatMeshArray(reflect);
    for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
    {
	index = d_index2d(i,j,top_gmax);
	cell_center[index].m_state.m_P = array[index];
	cell_center[index].m_state.m_q = array[index];
    }

    computeGradientQ();
    copyMeshStates();
    setAdvectionDt(); //Not used in debug mode
    printf("\nAfter seInitialCondition!\n");
} /* end setInitialCondition */


void Incompress_Solver_Smooth_2D_Cartesian_Debug::setInitialCondition_vd(void)
{
    // Incompress_Solver_Smooth_2D_Cartesian::setInitialCondition_vd();
    printf("\nEnter 2D_Debug::setInitialCondition_vd()\n");
    int i,j,l;
    int index;
    COMPONENT comp;
    double coords[MAXD];
        int reflect[MAXD];
        reflect[0] = reflect[1] = YES;

    FT_MakeGridIntfc(front);
    setDomain_vd();
    setComponent();

    m_rho[0] = 1;
    m_rho[1] = 1;
    m_mu[0] = 1;
    m_mu[1] = 1;
    m_comp[0] = iFparams->m_comp1;
    m_comp[1] = iFparams->m_comp2;
    m_smoothing_radius = 0;
    m_sigma = 0;
    mu_min = rho_min = HUGE;

    // for vd
    m_c[0] = 1;
    m_c[1] = 1;
    m_Dcoef[0] = 0;
    m_Dcoef[1] = 0;
    m_rho_old[0] = 1;
    m_rho_old[1] = 1;
    Dcoef_min = c_min = HUGE;
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
    double CFL = 0.1;
    front->max_time = 0.5;
    double dh;
    dh = std::min(front->rect_grid->h[0], front->rect_grid->h[1]);
    double dt = CFL * dh;
    printf("\nBefore adjusting, the dt = %20.16g\n",dt);
    int nStep = int (front->max_time/dt) + 1;
    dt = front->max_time/ nStep;
    printf("\nAfter adjusting, the dt = %20.16g\n",dt);
    front->dt = dt;
    for (i = 0; i < cell_center.size(); i++)
    {
        getRectangleCenter(i, coords);
        cell_center[i].m_state.setZero();
        getExactSolution_vd(coords,front->time,cell_center[i].m_state);
        cell_center[i].m_state.m_q = cell_center[i].m_state.m_P;
    }
    front->dt = 0.0;

    for (l = 0; l < dim; l++)
    {
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index2d(i,j,top_gmax);
            array[index] = cell_center[index].m_state.m_U[l];
        }
        scatMeshArray(reflect);
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index2d(i,j,top_gmax);
            cell_center[index].m_state.m_U[l] = array[index];
        }
    }

    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
        index = d_index2d(i,j,top_gmax);
        array[index] = cell_center[index].m_state.m_P;
    }
    scatMeshArray(reflect);
    for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
    {
        index = d_index2d(i,j,top_gmax);
        cell_center[index].m_state.m_P = array[index];
        cell_center[index].m_state.m_q = array[index];
    }

    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
        index = d_index2d(i,j,top_gmax);
        array[index] = cell_center[index].m_state.m_rho;
    }
    scatMeshArray(reflect);
    for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
    {
        index = d_index2d(i,j,top_gmax);
        cell_center[index].m_state.m_rho = array[index];
        cell_center[index].m_state.m_rho_old = array[index];
    }

    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
        index = d_index2d(i,j,top_gmax);
        array[index] = cell_center[index].m_state.m_c;
    }
    scatMeshArray(reflect);
    for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
    {
        index = d_index2d(i,j,top_gmax);
        cell_center[index].m_state.m_c = array[index];
    }

    computeGradientQ();
    copyMeshStates_vd();
    setAdvectionDt(); //Not used in debug mode
    printf("\nAfter setInitialCondition_vd!\n");
} /* end setInitialCondition_vd */


void Incompress_Solver_Smooth_2D_Cartesian_Debug::solve(double dt)
{

    printf("\nEntering Incompress Solve!! The dt for this time step is %.16g\n",dt);
    //    Incompress_Solver_Smooth_2D_Cartesian::solve(dt);

    //    solve_diffusionOnly(dt);
    //    return;


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
    start_clock("setSmoothedPropertiesOnePhase");
    setSmoProOnePhase();
    stop_clock("setSmoothedPropertiesOnePhase");
    if (debugging("trace"))
	printf("Passed setSmoothedProperties_OnePhase()\n");

    //	// 1) solve for intermediate velocity
    //	start_clock("computeAdvection");
    //	computeAdvection();
    //	stop_clock("computeAdvection");
    //	if (debugging("step_size"))
    //		printf("max_speed after computeAdvection(): %20.14f\n",
    //				max_speed);
    //

    start_clock("compAdvectionTerm");
    compAdvectionTerm_coupled_upgraded();
    stop_clock("compAdvectionTerm");

    start_clock("compDiffWithSmoothProperty");
    compDiffWithSmoothProperty_2nd_coupled();
    //	compDiffWithSmoothProperty();
    stop_clock("compDiffWithSmoothProperty");

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
} /* end solve */


void Incompress_Solver_Smooth_2D_Cartesian_Debug::solve_vd(double dt)
{
    printf("\nEntering Incompress Solve_vd!! The dt for this time step is %.16g\n",dt);

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

    start_clock("solve_vd");
    setDomain_vd();

    setComponent();
    if (debugging("trace"))
        printf("Passed setComponent()\n");
    setGlobalIndex();
    if (debugging("trace"))
        printf("Passed setGlobalIndex()\n");

    // start_clock("setSmoothedPropertiesOnePhase_vd");
    setSmoProOnePhase_vd();
    // stop_clock("setSmoothedPropertiesOnePhase_vd");
    if (debugging("trace"))
        printf("Passed setSmoothedProperties_vd_OnePhase()\n");

    // 1) solve for intermediate U, update density and concentration
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
    //compSGS();        //Subgrid model by Hyunkyun Lim
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
} /* end solve_vd */


void Incompress_Solver_Smooth_2D_Cartesian_Debug::solve_diffusionOnly(double dt)
{
    //    Incompress_Solver_Smooth_2D_Cartesian::solve(dt);

    m_t_old = front->time;
    m_t_int = front->time + front->dt/2;
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


    start_clock("compDiffWithSmoothProperty");
    compDiffWithSmoothProperty_2nd_decoupled();
    //	compDiffWithSmoothProperty();
    stop_clock("compDiffWithSmoothProperty");



    start_clock("copyMeshStates");
    copyMeshStates();
    stop_clock("copyMeshStates");

    setAdvectionDt();
    stop_clock("solve");
} /* end solve_diffusionOnly */


void Incompress_Solver_Smooth_2D_Cartesian_Debug::getVelocity(double *p, double *U)
{
    int i;
    for (i = 0; i < dim; ++i)
	U[i] = 0.0;
}

/**
* This function is written here so that it can be overwritten by a derived
* class.
*/

// for vd
boolean Incompress_Solver_Smooth_2D_Cartesian_Debug::FT_StateStructAtGridCrossing_tmp(
	Front *front,
	int *icoords,
	GRID_DIRECTION dir,
	COMPONENT comp,
	Locstate *state,
	HYPER_SURF **hs,
	double *crx_coords,
	double t)
{
    boolean bRet = FT_StateStructAtGridCrossing(front,icoords,dir,comp,state,hs,crx_coords);

    if(bRet)
    {
	static int debug_index = 0;
	L_STATE exact_state;

        getExactSolution(crx_coords,t,exact_state);
        // for vd
	//getExactSolution_vd(crx_coords,t,exact_state);

	STATE *fstate = (STATE*)(*state);
	fstate->vel[0] = exact_state.m_U[0];
	fstate->vel[1] = exact_state.m_U[1];
        // for vd
        //fstate->dens = exact_state.m_rho;
        //fstate->conc = exact_state.m_c;
    }
    return bRet;
} /* end FT_StateStructAtGridCrossing_tmp */


/**
* see also:
* EBM2D_Liquid_Solution::getExactSolution, using TEST_BROWN.
*/
void Incompress_Solver_Smooth_2D_Cartesian_Debug::computeSourceTerm_Adv(
	double *coords,
	L_STATE &state)
{
    //    Incompress_Solver_Smooth_2D_Cartesian::computeSourceTerm(coords,state);

    double t = front->time;
    double x = coords[0];
    double y = coords[1];

    double mu = m_mu[0];
    double rho = m_rho[0];

    // BROWN
    state.m_U[0] = rho * (0.8e1 * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * 0.3141592654e1 * cos(0.2e1 * 0.3141592654e1 * t * t) * t * (double) (3 * y * y - 2 * y) - 0.2e1 * cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double)  pow((double) (3 * y * y - 2 * y), (double) 2) * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 + 0.2e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) (y * y) * (double) (y - 1) * cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) (6 * y - 2)) - 0.4e1 * cos(0.2e1 * 0.3141592654e1 * t * t) * t * cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * (sin(0.2e1 * 0.3141592654e1 * (double) y) - 0.2e1 * 0.3141592654e1 * (double) y + 0.3141592654e1) + 0.2e1 * mu * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * (-0.2e1 * sin(0.2e1 * 0.3141592654e1 * (double) y) + 0.2e1 * 0.3141592654e1 * (double) y - 0.3141592654e1) - mu * (-0.4e1 * cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * 0.3141592654e1 * (double) (3 * y * y - 2 * y) + 0.6e1 * cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))));
    state.m_U[1] = rho * (-0.16e2 * pow(0.3141592654e1, 0.3e1) * cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * cos(0.2e1 * 0.3141592654e1 * t * t) * t * (double) (y * y) * (double) (y - 1) + 0.4e1 * pow(cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))), 0.2e1) * (double) (3 * y * y - 2 * y) * 0.3141592654e1 * 0.3141592654e1 * (double) (y * y) * (double) (y - 1) + 0.2e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) (y * y) * (double) (y - 1) * (0.4e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) y * (double) (y - 1) + 0.2e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) (y * y))) - 0.2e1 * cos(0.2e1 * 0.3141592654e1 * t * t) * t * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (0.2e1 * cos(0.2e1 * 0.3141592654e1 * (double) y) * 0.3141592654e1 - 0.2e1 * 0.3141592654e1) - mu * cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (-0.4e1 * cos(0.2e1 * 0.3141592654e1 * (double) y) * 0.3141592654e1 + 0.2e1 * 0.3141592654e1) - mu * (-0.8e1 * pow(0.3141592654e1, 0.3e1) * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) (y * y) * (double) (y - 1) + 0.4e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) (y - 1) + 0.8e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) y);

    //	// diffusion
    //	state.m_U[0] = -0.24e2 * rho * sin(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * 0.3141592654e1 * cos(0.2e1 * 0.3141592654e1 * t * t) * t * y * y + 0.16e2 * rho * sin(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * 0.3141592654e1 * cos(0.2e1 * 0.3141592654e1 * t * t) * t * y + 0.12e2 * mu * cos(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * 0.3141592654e1 * y * y - 0.8e1 * mu * cos(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * 0.3141592654e1 * y - 0.6e1 * mu * cos(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t)));
    ////	state.m_U[1] = -0.4e1 * 0.3141592654e1 * (0.4e1 * rho * 0.3141592654e1 * 0.3141592654e1 * cos(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * cos(0.2e1 * 0.3141592654e1 * t * t) * t * pow(y, 0.3e1) - 0.4e1 * rho * 0.3141592654e1 * 0.3141592654e1 * cos(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * cos(0.2e1 * 0.3141592654e1 * t * t) * t * y * y + 0.2e1 * mu * 0.3141592654e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * pow(y, 0.3e1) - 0.2e1 * mu * 0.3141592654e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * y * y - 0.3e1 * mu * sin(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * y + mu * sin(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))));
    //	state.m_U[1] = state.m_U[0];

    state.m_P = HUGE_VAL;
} /* end computeSourceTerm_Adv */


void Incompress_Solver_Smooth_2D_Cartesian_Debug::computeSourceTerm(
	double *coords,
	L_STATE &state)
{
    //    Incompress_Solver_Smooth_2D_Cartesian::computeSourceTerm(coords,state);

    double t = front->time + front->dt/2.0;
    double x = coords[0];
    double y = coords[1];

    double mu = m_mu[0];
    double rho = m_rho[0];

    // BROWN
    state.m_U[0] = rho * (0.8e1 * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * 0.3141592654e1 * cos(0.2e1 * 0.3141592654e1 * t * t) * t * (double) (3 * y * y - 2 * y) - 0.2e1 * cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double)  pow((double) (3 * y * y - 2 * y), (double) 2) * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 + 0.2e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) (y * y) * (double) (y - 1) * cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) (6 * y - 2)) - 0.4e1 * cos(0.2e1 * 0.3141592654e1 * t * t) * t * cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * (sin(0.2e1 * 0.3141592654e1 * (double) y) - 0.2e1 * 0.3141592654e1 * (double) y + 0.3141592654e1) + 0.2e1 * mu * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * (-0.2e1 * sin(0.2e1 * 0.3141592654e1 * (double) y) + 0.2e1 * 0.3141592654e1 * (double) y - 0.3141592654e1) - mu * (-0.4e1 * cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * 0.3141592654e1 * (double) (3 * y * y - 2 * y) + 0.6e1 * cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))));
    state.m_U[1] = rho * (-0.16e2 * pow(0.3141592654e1, 0.3e1) * cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * cos(0.2e1 * 0.3141592654e1 * t * t) * t * (double) (y * y) * (double) (y - 1) + 0.4e1 * pow(cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))), 0.2e1) * (double) (3 * y * y - 2 * y) * 0.3141592654e1 * 0.3141592654e1 * (double) (y * y) * (double) (y - 1) + 0.2e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) (y * y) * (double) (y - 1) * (0.4e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) y * (double) (y - 1) + 0.2e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) (y * y))) - 0.2e1 * cos(0.2e1 * 0.3141592654e1 * t * t) * t * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (0.2e1 * cos(0.2e1 * 0.3141592654e1 * (double) y) * 0.3141592654e1 - 0.2e1 * 0.3141592654e1) - mu * cos(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (-0.4e1 * cos(0.2e1 * 0.3141592654e1 * (double) y) * 0.3141592654e1 + 0.2e1 * 0.3141592654e1) - mu * (-0.8e1 * pow(0.3141592654e1, 0.3e1) * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) (y * y) * (double) (y - 1) + 0.4e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) (y - 1) + 0.8e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (x - 0.1e1 - sin(0.2e1 * 0.3141592654e1 * t * t))) * (double) y);

    //	// diffusion
    //	state.m_U[0] = -0.24e2 * rho * sin(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * 0.3141592654e1 * cos(0.2e1 * 0.3141592654e1 * t * t) * t * y * y + 0.16e2 * rho * sin(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * 0.3141592654e1 * cos(0.2e1 * 0.3141592654e1 * t * t) * t * y + 0.12e2 * mu * cos(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * 0.3141592654e1 * y * y - 0.8e1 * mu * cos(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * 0.3141592654e1 * 0.3141592654e1 * y - 0.6e1 * mu * cos(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t)));
    ////	state.m_U[1] = -0.4e1 * 0.3141592654e1 * (0.4e1 * rho * 0.3141592654e1 * 0.3141592654e1 * cos(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * cos(0.2e1 * 0.3141592654e1 * t * t) * t * pow(y, 0.3e1) - 0.4e1 * rho * 0.3141592654e1 * 0.3141592654e1 * cos(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * cos(0.2e1 * 0.3141592654e1 * t * t) * t * y * y + 0.2e1 * mu * 0.3141592654e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * pow(y, 0.3e1) - 0.2e1 * mu * 0.3141592654e1 * 0.3141592654e1 * sin(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * y * y - 0.3e1 * mu * sin(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))) * y + mu * sin(0.2e1 * 0.3141592654e1 * (-x + sin(0.2e1 * 0.3141592654e1 * t * t))));
    //	state.m_U[1] = state.m_U[0];

    state.m_P = HUGE_VAL;
} /* end computeSourceTerm */


void Incompress_Solver_Smooth_2D_Cartesian_Debug::getExactSolution(
	double *coords,
	double t,
	L_STATE &state)
{
    //printf(" The front->dt is %20.16g\n",front->dt);
    //printf(" The front->time is %20.16g\n",front->time);
    double pi = 3.14159265358979323846264338327950288419716939937510;
    double omega = 1+sin(2*pi*t*t);
    double x = coords[0];
    double y = coords[1];

    double t_int;
    //if (t == 0)
	//t_int = 0.0;
    //else if(t > 0)
	t_int = t - (front->dt)/2.0;
    double mu = m_mu[0];

    // BROWN
    state.m_U[0] = cos(2*pi*(x-omega))*(3*y*y-2*y);
    state.m_U[1] = 2*pi*sin(2*pi*(x-omega))*y*y*(y-1);

    state.m_P = -mu*(cos(2*pi*(-1 + x - sin(2*pi*t_int*t_int)))*
      (-pi + 2*pi*y - 2*sin(2*pi*y))) -
   2*t_int*cos(2*pi*t_int*t_int)*(pi - 2*pi*y + sin(2*pi*y))*
    sin(2*pi*(-1 + x - sin(2*pi*t_int*t_int)));
    // diffusion
    //	state.m_U[0] = cos(2*pi*(x-omega))*(3*y*y-2*y);
    ////	state.m_U[1] = 2*pi*sin(2*pi*(x-omega))*y*y*(y-1);
    //	state.m_U[1] = state.m_U[0];
} /* end getExactSolution */


void Incompress_Solver_Smooth_2D_Cartesian_Debug::getExactSolution_vd(
        double *coords,
        double t,
        L_STATE &state)
{
    //printf(" The front->dt is %20.16g\n",front->dt);
    //printf(" The front->time is %20.16g\n",front->time);
    double pi = 3.14159265358979323846264338327950288419716939937510;
    double omega = 1+sin(2*pi*t*t);
    double x = coords[0];
    double y = coords[1];

    double t_int;
    //if (t == 0)
        //t_int = 0.0;
    //else if(t > 0)
        t_int = t - (front->dt)/2.0;
    double mu = m_mu[0];

    // BROWN
    state.m_U[0] = cos(2*pi*(x-omega))*(3*y*y-2*y);
    state.m_U[1] = 2*pi*sin(2*pi*(x-omega))*y*y*(y-1);

    state.m_P = -mu*(cos(2*pi*(-1 + x - sin(2*pi*t_int*t_int)))*
      (-pi + 2*pi*y - 2*sin(2*pi*y))) -
   2*t_int*cos(2*pi*t_int*t_int)*(pi - 2*pi*y + sin(2*pi*y))*
    sin(2*pi*(-1 + x - sin(2*pi*t_int*t_int)));

    state.m_rho = state.m_rho_old = m_rho[0];
    state.m_c = m_c[0];
} /* end getExactSolution_vd */


void Incompress_Solver_Smooth_2D_Cartesian_Debug::saveStates_Tecplot(
	const char*out_name,
	double t,
	bool bPrintCoords,
	bool bPrintExact,
	bool bPrintError)
{
    bPrintExact = true;
    /*
    char filename[200];
    sprintf(filename,"%s/states_%s.plt",out_name,
	    right_flush(front->step,7));

    FILE *hfile=fopen(filename,"w");

    fprintf(hfile, "VARIABLES = X, Y, U, V, P");
    fprintf(hfile, ", U_exact, V_exact, P_exact");
    fprintf(hfile, ", U_error, V_error, P_error");
    fprintf(hfile, "\n");

    fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n",
	    imax-imin+1, jmax-jmin+1);
*/
    double *coords, cellArea = top_h[0]*top_h[1];

    int max_ij[3][2], index;
    max_ij[0][0] = max_ij[0][1] = -1;
    max_ij[1][0] = max_ij[1][1] = -1;
    max_ij[2][0] = max_ij[2][1] = -1;

    L_STATE error, L1, L2, LInf;

    L_STATE state, state_exact;

    LInf.m_U[0] = L1.m_U[0] = L2.m_U[0] = 0;
    LInf.m_U[1] = L1.m_U[1] = L2.m_U[1] = 0;
    LInf.m_P = L1.m_P = L2.m_P = 0;

    for(int j=jmin; j<=jmax; j++)
	for(int i=imin; i<=imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    coords = cell_center[index].m_coords;
	    state = cell_center[index].m_state;

	    if(bPrintExact)
		getExactSolution(coords, t, state_exact);
	    else
	    {
		state_exact.m_U[0] = state_exact.m_U[1] = state_exact.m_P = 0;
	    }

	    error.m_U[0] = fabs(state.m_U[0] - state_exact.m_U[0]);
	    error.m_U[1] = fabs(state.m_U[1] - state_exact.m_U[1]);
	    error.m_P    = fabs(state.m_P - state_exact.m_P);

	    // LInfinity
	    if(error.m_U[0] > LInf.m_U[0])
	    {
		LInf.m_U[0] = error.m_U[0];
		max_ij[0][0] = i;
		max_ij[0][1] = j;
	    }
	    if(error.m_U[1] > LInf.m_U[1])
	    {
		LInf.m_U[1] = error.m_U[1];
		max_ij[1][0] = i;
		max_ij[1][1] = j;
	    }
	    if(error.m_P > LInf.m_P)
	    {
		LInf.m_P = error.m_P;
		max_ij[2][0] = i;
		max_ij[2][1] = j;
	    }


	    // L1
	    L1.m_U[0] += error.m_U[0]*cellArea;
	    L1.m_U[1] += error.m_U[1]*cellArea;
	    L1.m_P += error.m_P*cellArea;


	    // L2
	    L2.m_U[0] += error.m_U[0]*error.m_U[0]*cellArea;
	    L2.m_U[1] += error.m_U[1]*error.m_U[1]*cellArea;
	    L2.m_P += error.m_P*error.m_P*cellArea;
/*
	    if(bPrintCoords)
		fprintf(hfile, "% 12.8e % 12.8e", coords[0], coords[1]);
	    else
		fprintf(hfile, "%3d %3d", i, j);

	    fprintf(hfile, " % 12.8e % 12.8e % 12.8e",
		    state.m_U[0], state.m_U[1], state.m_P);
	    if(bPrintExact)
		fprintf(hfile, " % 12.8e % 12.8e % 12.8e",
			state_exact.m_U[0], state_exact.m_U[1], state_exact.m_P);
	    if(bPrintError)
		fprintf(hfile, " % 12.8e % 12.8e % 12.8e",
			state_exact.m_U[0] - state.m_U[0],
			state_exact.m_U[1] - state.m_U[1],
			state_exact.m_P - state.m_P);
	    fprintf(hfile, "\n");
	    */
	}
	    pp_global_sum(&L1.m_U[0],1);
	    pp_global_sum(&L1.m_U[1],1);
	    pp_global_sum(&L1.m_P,1);

	    pp_global_sum(&L2.m_U[0],1);
	    pp_global_sum(&L2.m_U[1],1);
	    pp_global_sum(&L2.m_P,1);

	    pp_global_max(&LInf.m_U[0],1);
	    pp_global_max(&LInf.m_U[1],1);
	    pp_global_max(&LInf.m_P,1);

    // LInfinity
    // L1
    // L2
    L2.m_U[0] = sqrt(L2.m_U[0]);
    L2.m_U[1] = sqrt(L2.m_U[1]);
    L2.m_P = sqrt(L2.m_P);

    //fclose(hfile);

    printf("Incompress_Solver_Smooth_2D_Cartesian_Debug::saveStates_Tecplot: \n");
    printf("\t#### L1  ={%12.8f, %12.8f, %12.8f}\n", L1.m_U[0], L1.m_U[1], L1.m_P);
    printf("\t#### L2  ={%12.8f, %12.8f, %12.8f}\n", L2.m_U[0], L2.m_U[1], L2.m_P);
    printf("\t#### LInf={%12.8f, %12.8f, %12.8f}\n", LInf.m_U[0], LInf.m_U[1], LInf.m_P);

} /* saveStates_Tecplot */


void Incompress_Solver_Smooth_2D_Cartesian_Debug::saveStates_Tecplot_vd(
        const char*out_name,
        double t,
        bool bPrintCoords,
        bool bPrintExact,
        bool bPrintError)
{
    bPrintExact = true;
    /*
    char filename[200];
    sprintf(filename,"%s/states_%s.plt",out_name,
            right_flush(front->step,7));

    FILE *hfile=fopen(filename,"w");

    fprintf(hfile, "VARIABLES = X, Y, U, V, P");
    fprintf(hfile, ", U_exact, V_exact, P_exact");
    fprintf(hfile, ", U_error, V_error, P_error");
    fprintf(hfile, "\n");

    fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n",
            imax-imin+1, jmax-jmin+1);
*/
    double *coords, cellArea = top_h[0]*top_h[1];

    int max_ij[5][2], index;
    max_ij[0][0] = max_ij[0][1] = -1;
    max_ij[1][0] = max_ij[1][1] = -1;
    max_ij[2][0] = max_ij[2][1] = -1;
    max_ij[3][0] = max_ij[3][1] = -1;
    max_ij[4][0] = max_ij[4][1] = -1;

    L_STATE error, L1, L2, LInf;

    L_STATE state, state_exact;

    LInf.m_U[0] = L1.m_U[0] = L2.m_U[0] = 0;
    LInf.m_U[1] = L1.m_U[1] = L2.m_U[1] = 0;
    LInf.m_P = L1.m_P = L2.m_P = 0;
    LInf.m_rho = L1.m_rho = L2.m_rho = 0;
    LInf.m_c = L1.m_c = L2.m_c = 0;

    for(int j=jmin; j<=jmax; j++)
        for(int i=imin; i<=imax; i++)
        {
            index = d_index2d(i,j,top_gmax);
            coords = cell_center[index].m_coords;
            state = cell_center[index].m_state;

            if(bPrintExact)
                getExactSolution_vd(coords, t, state_exact);
            else
            {
                state_exact.m_U[0] = state_exact.m_U[1] = state_exact.m_P = 0;
                // for vd
                state_exact.m_rho = state_exact.m_c = 0;
            }

            error.m_U[0] = fabs(state.m_U[0] - state_exact.m_U[0]);
            error.m_U[1] = fabs(state.m_U[1] - state_exact.m_U[1]);
            error.m_P    = fabs(state.m_P - state_exact.m_P);
            // for vd
            error.m_rho  = fabs(state.m_rho - state_exact.m_rho);
            error.m_c    = fabs(state.m_c - state_exact.m_c);

            // LInfinity
            if(error.m_U[0] > LInf.m_U[0])
            {
                LInf.m_U[0] = error.m_U[0];
                max_ij[0][0] = i;
                max_ij[0][1] = j;
            }
            if(error.m_U[1] > LInf.m_U[1])
            {
                LInf.m_U[1] = error.m_U[1];
                max_ij[1][0] = i;
                max_ij[1][1] = j;
            }
            if(error.m_P > LInf.m_P)
            {
                LInf.m_P = error.m_P;
                max_ij[2][0] = i;
                max_ij[2][1] = j;
            }
            // LInfinity for vd
            if(error.m_rho > LInf.m_rho)
            {
                LInf.m_rho = error.m_rho;
                max_ij[3][0] = i;
                max_ij[3][1] = j;
            }
            if(error.m_c > LInf.m_c)
            {
                LInf.m_c = error.m_c;
                max_ij[4][0] = i;
                max_ij[4][1] = j;
            }

            // L1
            L1.m_U[0] += error.m_U[0]*cellArea;
            L1.m_U[1] += error.m_U[1]*cellArea;
            L1.m_P += error.m_P*cellArea;
            // L1 for vd
            L1.m_rho += error.m_rho*cellArea;
            L1.m_c += error.m_c*cellArea;

            // L2
            L2.m_U[0] += error.m_U[0]*error.m_U[0]*cellArea;
            L2.m_U[1] += error.m_U[1]*error.m_U[1]*cellArea;
            L2.m_P += error.m_P*error.m_P*cellArea;
            // L2 for vd
            L2.m_rho += error.m_rho*error.m_rho*cellArea;
            L2.m_c += error.m_c*error.m_c*cellArea;

            //if(bPrintCoords)
            //    fprintf(hfile, "% 12.8e % 12.8e", coords[0], coords[1]);
            //else
                  //fprintf(hfile, "%3d %3d", i, j);
            /*
            printf( "\t#### i=%3d, j= %3d\n", i, j);
            printf("\t#### m_rho  ={%12.8f}\n", state.m_rho);
            */

           // fprintf(hfile, " % 12.8e % 12.8e % 12.8e",
           //         state.m_U[0], state.m_U[1], state.m_P);
           // fprintf(hfile, " % 12.8e", state.m_rho);
            //if(bPrintExact)
            //    fprintf(hfile, " % 12.8e % 12.8e % 12.8e",
            //            state_exact.m_U[0], state_exact.m_U[1], state_exact.m_P);
            //if(bPrintError)
            //   fprintf(hfile, " % 12.8e % 12.8e % 12.8e",
	    //state_exact.m_U[0] - state.m_U[0],
            //            state_exact.m_U[1] - state.m_U[1],
            //            state_exact.m_P - state.m_P);
           // fprintf(hfile, "\n");

        }
            pp_global_sum(&L1.m_U[0],1);
            pp_global_sum(&L1.m_U[1],1);
            pp_global_sum(&L1.m_P,1);
            // for vd
            pp_global_sum(&L1.m_rho,1);
            pp_global_sum(&L1.m_c,1);

            pp_global_sum(&L2.m_U[0],1);
            pp_global_sum(&L2.m_U[1],1);
            pp_global_sum(&L2.m_P,1);
            // for vd
            pp_global_sum(&L2.m_rho,1);
            pp_global_sum(&L2.m_c,1);

            pp_global_max(&LInf.m_U[0],1);
            pp_global_max(&LInf.m_U[1],1);
            pp_global_max(&LInf.m_P,1);
            // for vd
            pp_global_max(&LInf.m_rho,1);
            pp_global_max(&LInf.m_c,1);

    // LInfinity
    // L1
    // L2
    L2.m_U[0] = sqrt(L2.m_U[0]);
    L2.m_U[1] = sqrt(L2.m_U[1]);
    L2.m_P = sqrt(L2.m_P);
    // for vd
    L2.m_rho = sqrt(L2.m_rho);
    L2.m_c = sqrt(L2.m_c);

    //fclose(hfile);

    printf("Incompress_Solver_Smooth_2D_Cartesian_Debug::saveStates_Tecplot_vd: \n");
    printf("\t#### L1  ={%12.8f, %12.8f, %12.8f, %12.8f, %12.8f}\n", L1.m_U[0], L1.m_U[1], L1.m_P, L1.m_rho, L1.m_c);
    printf("\t#### L2  ={%12.8f, %12.8f, %12.8f, %12.8f, %12.8f}\n", L2.m_U[0], L2.m_U[1], L2.m_P, L2.m_rho, L2.m_c);
    printf("\t#### LInf={%12.8f, %12.8f, %12.8f, %12.8f, %12.8f}\n", LInf.m_U[0], LInf.m_U[1], LInf.m_P, LInf.m_rho, LInf.m_c);
    printf("\t#### (i,j) with rho_max={%3d, %3d}\n",max_ij[3][0], max_ij[3][1]);

} /* end saveStates_Tecplot_vd */


void Incompress_Solver_Smooth_2D_Cartesian_Debug::saveDivUPhi_Tecplot(
	const char*out_name,
	double t,
	bool bPrintCoords)
{
    char filename[200];
    sprintf(filename,"%s/DivUPhi_%s.plt",out_name,
	    right_flush(front->step,7));

    FILE *hfile=fopen(filename,"w");

    fprintf(hfile, "VARIABLES = X, Y, DivU, Phi\n");

    fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n",
	    imax-imin+1, jmax-jmin+1);

    double *coords;
    int index;
    L_STATE state;

    double max_div_U = 0, max_phi = 0;

    for(int j=jmin; j<=jmax; j++)
	for(int i=imin; i<=imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    coords = cell_center[index].m_coords;

	    state = cell_center[index].m_state;
	    if(fabs(state.div_U)>max_div_U)
		max_div_U = fabs(state.div_U);
	    if(fabs(state.m_phi)>max_phi)
		max_phi = fabs(state.m_phi);

	    if(bPrintCoords)
		fprintf(hfile, "% 12.8e % 12.8e", coords[0], coords[1]);
	    else
		fprintf(hfile, "%3d %3d", i, j);

	    fprintf(hfile, " % 12.8e % 12.8e",
		    state.div_U, state.m_phi);
	    fprintf(hfile, "\n");
	}

    fclose(hfile);

    printf("Incompress_Solver_Smooth_2D_Cartesian_Debug::saveDivUPhi_Tecplot: max_div_U = %12.8f, max_phi = %12.8f\n",
    	    max_div_U, max_phi);
} /* end saveDivUPhi_Tecplot */


void Incompress_Solver_Smooth_2D_Cartesian_Debug::saveParameters_Tecplot(
	const char*out_name,
	double t,
	bool bPrintCoords)
{
    char filename[200];
    sprintf(filename,"%s/parameters_%s.plt",out_name,
	    right_flush(front->step,7));

    FILE *hfile=fopen(filename,"w");

    fprintf(hfile, "VARIABLES = X, Y, rho, mu\n");

    fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n",
	    imax-imin+1, jmax-jmin+1);

    double *coords;
    int index;
    L_STATE state;

    for(int j=jmin; j<=jmax; j++)
	for(int i=imin; i<=imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    coords = cell_center[index].m_coords;
	    state = cell_center[index].m_state;
	    if(bPrintCoords)
		fprintf(hfile, "% 12.8e % 12.8e", coords[0], coords[1]);
	    else
		fprintf(hfile, "%3d %3d", i, j);

	    fprintf(hfile, " % 12.8e % 12.8e",
		    state.m_mu, state.m_rho);
	    fprintf(hfile, "\n");
	}

    fclose(hfile);
} /* end saveParameters_Tecplot */


void Incompress_Solver_Smooth_2D_Cartesian_Debug::saveParameters_Tecplot_vd(
        const char*out_name,
        double t,
        bool bPrintCoords)
{
    char filename[200];
    sprintf(filename,"%s/parameters_%s.plt",out_name,
            right_flush(front->step,7));

    FILE *hfile=fopen(filename,"w");

    fprintf(hfile, "VARIABLES = X, Y, mu\n");

    fprintf(hfile, "ZONE I=%d, J=%d, F=POINT\n",
            imax-imin+1, jmax-jmin+1);

    double *coords;
    int index;
    L_STATE state;

    for(int j=jmin; j<=jmax; j++)
        for(int i=imin; i<=imax; i++)
        {
            index = d_index2d(i,j,top_gmax);
            coords = cell_center[index].m_coords;
            state = cell_center[index].m_state;
            if(bPrintCoords)
                fprintf(hfile, "% 12.8e % 12.8e", coords[0], coords[1]);
            else
                fprintf(hfile, "%3d %3d", i, j);

            fprintf(hfile, " % 12.8e ", state.m_mu);
            fprintf(hfile, "\n");
        }

    fclose(hfile);
} /* end saveParameters_Tecplot_vd */

