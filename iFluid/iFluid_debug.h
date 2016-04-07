/**
 * iFluid_debug.h
 *
 *  Created on: Jun 23, 2010
 *      Author: shuqiangwang
 *  Merged and Modified on Aug 02, 2010 by Yijie Zhou
 *  
 *
 *  This class is used to check the order of accuracy of the incompressible
 *  code defined in iFluid.h using exact solution.
 *
 * Notable change:
 * 1 FT_StateStructAtGridCrossing_tmp() is defined as a member function which
 *   is easy to be overwritten by derived class.
 * 2 iFluid.cpp needs to be modified.
 */

#ifndef IFLUID_DEBUG_H_
#define IFLUID_DEBUG_H_
#include "iFluid.h"

/*****************************************************************************
 *			Incompress_Solver_Smooth_2D_Cartesian_Debug
 * 2D debug.
 *
 * Using input file in-1C2_10,in-1C2_20, in-1C2_40, ... where 1C for 1 component,
 * 2 for 2D and 10, 20, 40 for mesh size.
 *
 * The parameter for the density, viscosity are hard coded in setInitialCondition().
 * The exact solution and the right hand side are hard coded in computeSourceTerm().
 *
 * Reference:
 * [1] Brown and Cortez and Minion, Accurate Projection Methods for the
 * Incompressible Navier–Stokes Equations, 168, 464–499 (2001).
 *****************************************************************************/
class Incompress_Solver_Smooth_2D_Cartesian_Debug: public Incompress_Solver_Smooth_2D_Cartesian {
public:
	Incompress_Solver_Smooth_2D_Cartesian_Debug(Front &front);
	void setInitialCondition(void);
        // for vd
        void setInitialCondition_vd(void);
        void setInitialCondition_RSRV_vd(LEVEL_FUNC_PACK*) {};
        void setInitialDiffusionVelocity_vd(void) {};

	void solve(double dt);
        // for vd
        void solve_vd(double dt);

	void solve_diffusionOnly(double dt);

	void getVelocity(double *p, double *U); //Rewritten to let the interface to be stable
	boolean FT_StateStructAtGridCrossing_tmp(Front*,int*,GRID_DIRECTION,
				COMPONENT,Locstate*,HYPER_SURF**,double*,double t=0);
protected:

	// time at time step n, n+1/2 and n+1; set in solve()
//	double m_t_old, m_t_int, m_t_new;

	//! solution dependent
	void computeSourceTerm(double *coords,L_STATE &state);
	void computeSourceTerm_Adv(double *coords, L_STATE &state);

	void getExactSolution(double *coords,double t,L_STATE &state);
        // for vd
        void getExactSolution_vd(double *coords,double t,L_STATE &state);
public:
	void saveStates_Tecplot(
			const char*out_name, double t,
			bool bPrintCoords=true,
			bool bPrintExact=true,
			bool bPrintError=true);
        // for vd
        void saveStates_Tecplot_vd(
                        const char*out_name, double t,
                        bool bPrintCoords=true,
                        bool bPrintExact=true,
                        bool bPrintError=true);

	void saveDivUPhi_Tecplot(
		const char*out_name, double t,
		bool bPrintCoords=true);

	void saveParameters_Tecplot(
		const char*out_name, double t,
		bool bPrintCoords=true);
        // for vd
        void saveParameters_Tecplot_vd(
                const char*out_name, double t,
                bool bPrintCoords=true);

	bool m_bStartDebugging;
	L_STATE *m_pStateDebug;
};


/*****************************************************************************
 *			Incompress_Solver_Smooth_3D_Cartesian_Debug
 * 3D debug.
 *
 * Using input file in-1C3_10,in-1C3_20, in-1C3_40, ... where 1C for 1 component,
 * 3 for 3D and 10, 20, 40 for mesh size.
 *
 * The parameter for the density, viscosity are hard coded in setInitialCondition().
 * The exact solution and the right hand side are hard coded in computeSourceTerm().
 *
 * Reference:
 * [1] Ethier and Steinman, Exact Fully 3D Navier-Stokes Solutions for Benchmarking,
 * International Journal for Numerical Methods in Fluids, Vol 19, 369-375 (1994).
 *****************************************************************************/
class Incompress_Solver_Smooth_3D_Cartesian_Debug: public Incompress_Solver_Smooth_3D_Cartesian {
public:
	Incompress_Solver_Smooth_3D_Cartesian_Debug(Front &front);
	void setInitialCondition(void);
        // for vd
        void setInitialCondition_vd(void);
        void setInitialCondition_RSRV_vd(LEVEL_FUNC_PACK*) {};
        void setInitialDiffusionVelocity_vd(void){};

	void solve(double dt);
        // for vd
        void solve_vd(double dt);

	void solve_diffusionOnly(double dt); //Solve the diffusion equation only

	void getVelocity(double *p, double *U);
	boolean FT_StateStructAtGridCrossing_tmp(Front*,int*,GRID_DIRECTION,
			COMPONENT,Locstate*,HYPER_SURF**,double*,double t=0); 
	//Get the boundary states using exact solution
protected:
	// time at time step n, n+1/2 and n+1; set in solve()
	//	double m_t_old, m_t_int, m_t_new;

	//! solution dependent
	void computeSourceTerm(double *coords,L_STATE &state);
	void computeSourceTerm_Adv(double *coords, L_STATE &state);

	void getExactSolution(double *coords,double t,L_STATE &state);
        // for vd
        void getExactSolution_vd(double *coords,double t,L_STATE &state);
        void getExactSolution_MAC_vd(double *coords,double t,L_STATE &state);

public:
	void saveStates_Tecplot(
			const char*out_name, double t,
			bool bPrintCoords=true,
			bool bPrintExact=true,
			bool bPrintError=true);
        // for vd
        void saveStates_Tecplot_vd(
                        const char*out_name, double t,
                        bool bPrintCoords=true,
                        bool bPrintExact=true,
                        bool bPrintError=true);

	void saveDivUPhi_Tecplot(
		const char*out_name, double t,
		bool bPrintCoords=true);

	void saveParameters_Tecplot(
		const char*out_name, double t,
		bool bPrintCoords=true);
        // for vd
        void saveParameters_Tecplot_vd(
                const char*out_name, double t,
                bool bPrintCoords=true);
};

/*****************************************************************************
*			Incompress_Solver_Smooth_3D_Cylindrical_Debug
* 3D testing using exact solution.
*
* Using input file in-1C3_10,in-1C3_20, in-1C3_40, ... where 1C for 1 component,
* 3 for 3D and 10, 20, 40 for mesh size.
*
* The parameter for the density, viscosity are hard coded in setInitialCondition().
* The exact solution and the right hand side are hard coded in computeSourceTerm().
*
* Reference:
* [1] Ethier and Steinman, Exact Fully 3D Navier-Stokes Solutions for Benchmarking,
* International Journal for Numerical Methods in Fluids, Vol 19, 369-375 (1994).
*****************************************************************************/
class Incompress_Solver_Smooth_3D_Cylindrical_Debug: public Incompress_Solver_Smooth_3D_Cylindrical {
public:
    Incompress_Solver_Smooth_3D_Cylindrical_Debug(Front &front);

//    virtual void setComponent(void); //init components;
    void setInitialCondition(void);
    void setInitialDiffusionVelocity_vd(void) {};
    void solve(double dt);
    void solve_diffusionOnly(double dt);

    virtual void getVelocity(double *p, double *U);
    virtual boolean FT_StateStructAtGridCrossing_tmp(Front*,int*,GRID_DIRECTION,
	    COMPONENT,Locstate*,HYPER_SURF**,double*,double t=0);
protected:


    // time at time step n, n+1/2 and n+1; set in solve()
    //	double m_t_old, m_t_int, m_t_new;

    //! solution dependent
    //! the source term does not contain the density.
    void computeSourceTerm(double *coords,L_STATE &state);
    void computeSourceTerm_Adv(double *coords, L_STATE &state);

    void getExactSolution(double *coords,double t,L_STATE &state);
    // for vd
    void getExactSolution_vd(double *coords,double t,L_STATE &state){};

    double getAveragePressure(double t);

public:
    void saveStates_Tecplot(
	    const char*out_name, double t,
	    bool bPrintCoords=true,
	    bool bPrintExact=true,
	    bool bPrintError=true);

    void saveDivUPhi_Tecplot(
	    const char*out_name, double t,
	    bool bPrintCoords=true);

    void saveParameters_Tecplot(
	    const char*out_name, double t,
	    bool bPrintCoords=true);

    void getMaximumSourceTerm(double t);

    enum TEST_CASE
    {
	TEST_EthierSteinman,	// assuming that the density is 1
	TEST_EthierSteinman_Simple,
	TEST_EthierSteinman_Simple2,
	TEST_DIFFUSION,		// assuming that the density is 1
	TEST_ELLIPTIC
    };
    static TEST_CASE m_testingCase;
};

#endif /* IFLUID_DEBUG_H_ */
