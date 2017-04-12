/**********************************************************************
 * 		iFluid.h
 **********************************************************************/

#ifndef _FT_IFLUID_H_
#define _FT_IFLUID_H_

#include <vector>
#include <petscksp.h>
#include <assert.h>

#include "FronTier.h"
#include "solver.h"

#define         SOLID_COMP		1
#define         LIQUID_COMP1		2
#define         LIQUID_COMP2		3

// This is a shortcut to mark FREE SLIP WALL and NO SLIP WALL
// This feature comes with the introduction of REFLECTION BOUNDARY CONDITION
#define SLIP    YES
#define NOSLIP  YES

#define		ifluid_comp(comp)   (((comp) == LIQUID_COMP1 || 	\
		comp == LIQUID_COMP2) ? YES : NO)

enum EBM_COORD
{
    COORD_X = 0,  COORD_Y = 1,  COORD_Z = 2
};

enum intfc_geometry { planar, spherical };

struct _IF_FIELD {
	double *pres;			/* Pressure */
	double **vel;			/* Velocities */
	double *vort;			/* Vorticity in 2D */
	double **vort3d;		/* Vorticity in 3D */
	double *div;
        // for vd
        double *dens;                   /* Density */
        double *dens_old;               /* Old Density */
        double *conc;                   /* Concentration */
        double **grad_pres;             /* Gradient of pressure */
};
typedef struct _IF_FIELD IF_FIELD;

struct _IF_MOVIE_OPTION {
	/* HDF movie options */
	boolean plot_comp;
	boolean plot_pres;
	boolean plot_vort;
	boolean plot_velo;
        boolean plot_dens;
        boolean plot_conc;
	boolean plot_cross_section[MAXD]; /* 3D 0: yz; 1: zx; 2: xy */
	/* VTK movie options */
	boolean plot_vel_vector;	  /* Plot velocity vector field */
};
typedef struct _IF_MOVIE_OPTION IF_MOVIE_OPTION;

enum _NS_SCHEME {
	ERROR_SCHEME		= -1,
        SIMPLE			=  1,
        BELL_COLELLA,
        KIM_MOIN,
        PEROT_BOTELLA
};
typedef enum _NS_SCHEME NS_SCHEME;

typedef struct {
        int dim;
        POINTER level_func_params;
	IF_MOVIE_OPTION *movie_option;
	NS_SCHEME num_scheme;
        double rho1;
        double rho2;
	double mu1;
	double mu2;
	double U1[MAXD];
	double U2[MAXD];
        double bvel[2][MAXD];
	double gravity[MAXD];
	double surf_tension;
	double smoothing_radius;
	COMPONENT m_comp1;
	COMPONENT m_comp2;
	IF_FIELD *field;
        double init_vel[MAXD];
        boolean use_couette_init_vel;
        // for vd
        double c1;
        double c2;
        double Dcoef1;
        double Dcoef2;
        double vol_frac_threshold;
	double width_idl; //width of initial diffusion layer
    int ifluid_type;
} IF_PARAMS;

struct _FLOW_THROUGH_PARAMS {
        POINT *oldp;
        COMPONENT comp;
};
typedef struct _FLOW_THROUGH_PARAMS FLOW_THROUGH_PARAMS;

enum _TIME_FUNC_TYPE {
	CONSTANT		=  1,
	PULSE_FUNC,
	SINE_FUNC
};
typedef enum _TIME_FUNC_TYPE TIME_FUNC_TYPE;

struct _TIME_DEPENDENT_PARAMS {
	TIME_FUNC_TYPE td_type;
	double v_base[MAXD],p_base;
	double v_peak[MAXD],p_peak;
	double v_tail[MAXD],p_tail;
	double v_amp[MAXD],p_amp;
	double omega,phase;
	double T[10];
};
typedef struct _TIME_DEPENDENT_PARAMS TIME_DEPENDENT_PARAMS;

/******************************************************************************
 * 		lcartsn.h
 * A simple incompressible flow solver using the ghost fluid method and the
 * projection method.
 *
 * the main function is
 * 	L_CARTESIAN::solve().
 *
 * References:
 ******************************************************************************/

class SOLVER;
class Incompress_Solver_Basis;

//enum VISITED_TYPE {UNVISITED, VISITED, PARTIAL_VISITED, FULL_VISITED};

//-------------------------------------------------
//		STATES
// NOTE:
//      L_STATE/L_STATE_RECT_EDGE should be put into
// lcartsn.h. However, there are some trouble
// to compile in that way.
//-------------------------------------------------
// states inside a cell

class L_STATE{
public:
	double m_U[MAXD];		// velocity vector
	double m_exactU[MAXD];
	double m_P;			// pressure
	double m_exactP;

	double m_phi;			// (Brown's) phi*dt
	double m_q;
	double m_adv[MAXD];

	double m_mu;		        // smoothed
        double m_mu_old;                // smoothed
	double m_rho;		        // smoothed
	double div_U;			// Velocity divergence
	double grad_q[MAXD];		// Gradient of q
	double f_surf[MAXD];		// Surface force (such as tension)

	L_STATE();
	void setZero(void);

        //for vd
        double m_c;
        double m_c_old;      // concentration at last time step
        double m_rho_old;
        double m_Dcoef;                 // smoothed
        double m_rho_adv;               // for continuity eqn
        double m_c_adv;                 // for concentration eqn
//        double m_Ud[MAXD];              // velocity vector
//        double m_Ud_tmp[MAXD];          // velocity vector
//        double m_Ud_tmp_tmp[MAXD];      // velocity vector
        double m_U_tmp[MAXD];           // velocity vector
        double m_U_velo_var[MAXD];       // velocity vector
//        double m_U_fromBar[MAXD];       // velocity vector
        double div_U_tmp;               // Velocity divergence
//        double div_U_tmp_tmp;           // Velocity divergence
//        double div_U_face[2*MAXD];      // Velocity divergence on cell faces
//        double div_U_edge[18];          // Velocity divergence on cell edges
//        double div_U_vertex[2*MAXD+2];  // Velocity divergence on cell vertices
//        double m_u_hat[2*MAXD];         // for u_hat
//        double m_v_hat[2*MAXD];         // for v_hat
//        double m_w_hat[2*MAXD];         // for w_hat
        double m_u_bar[2*MAXD];         // for u_bar
        double m_v_bar[2*MAXD];         // for v_bar
//        double m_w_bar[2*MAXD];         // for w_bar
        double m_U_center_bar[MAXD];    // for U_center_bar
        double m_U_face_bar[MAXD];      // for U_face_bar
        double m_u_bar_tmp[2*MAXD];     // for u_bar_tmp
        double m_v_bar_tmp[2*MAXD];     // for v_bar_tmp
//        double m_w_bar_tmp[2*MAXD];     // for w_bar_tmp
//        double m_mflux[2*MAXD];         // for mass flux
        double m_rho_bar[2*MAXD];       // for rho_bar
        double m_c_bar[2*MAXD];         // for c_bar
        double m_mu_turbulent[MAXD+1];    // SGS term
        double m_Dcoef_turbulent[MAXD+1]; // SGS term
//        boolean isCutCell;              // is cut cell or not
//        boolean isRegCell;              // is regular cell or not
};
//states on MAC grid

//------------------------------------------------------
//		MESH
//------------------------------------------------------
// note that the following VERTEX2D/RECT_EDGE are different
// from those defined in MESH2D.h

class L_RECTANGLE {
public:
	int m_index;			// rectangle index
	int comp;
	L_STATE m_state;
	double m_coords[MAXD];
	int icoords[MAXD];

	L_RECTANGLE();

	void setCoords(double*,int);
};

class Incompress_Solver_Basis{
public:
       	Incompress_Solver_Basis() {}; // constructor
	virtual ~Incompress_Solver_Basis() {};

	virtual void setInitialCondition(void) = 0;    //Initialization
        virtual void setInitialCondition_vd(LEVEL_FUNC_PACK*) = 0; //for vd, Initialization
        virtual void setInitialCondition_RSRV_vd(LEVEL_FUNC_PACK*) = 0; //for vd, RS_RV case Initialization
        virtual void setInitialCondition_RSSY_vd(LEVEL_FUNC_PACK*) = 0; //for vd, Smeeton-Young case Initialization
        virtual void setInitialDiffusionVelocity_vd(void) = 0; //for vd, Initialization
	virtual void solve(double dt) = 0;             //main step function
        virtual void solve_vd(double dt) = 0;          // for vd, main step function
};

class Incompress_Solver_Smooth_Basis:public Incompress_Solver_Basis{
public:
        //constructor
	Incompress_Solver_Smooth_Basis(Front &front);
	virtual ~Incompress_Solver_Smooth_Basis() {};

	double m_dt;
	double accum_dt;
	double max_speed;
	double max_value; //for debugging
	double max_dt;
	double min_dt;
	double *top_h;
	int dim;
        // for vd
        double max_density;
        double min_density;
        double max_concentration;
        double min_concentration;
        double max_phi;
        double min_phi;

	void initMesh(void);
        // for vd
        void initMesh_vd(void);

	void setAdvectionDt(void);
			//using max speed and hmin to determine max_dt, min_dt
        // for vd
        void setAdvectionDt_vd(void);

	void readFrontInteriorStates(char *state_name,bool,bool);
        // for vd
        void readFrontInteriorStates_vd(char *state_name,bool,bool);

        virtual void printExpandedMesh(char *outname, bool binary) = 0;

	void printFrontInteriorStates(char *state_name,bool);
        void printFrontInteriorStatesRegridRep(char *state_name,bool,bool);
        // for vd
        void printFrontInteriorStates_vd(char *state_name,bool);
        void printFrontInteriorStatesRegridRep_vd(char *state_name,bool,bool);

	void initMovieVariables(void);
	void augmentMovieVariables(void);
	void getVelocity(double *p, double *U);
        void initSampleVelocity(char *in_name);
        // for vd
        void getVelocity_MAC_vd(double *p, double *U);
        void getDensity_vd(double *, double *);
        void getDensityOld_vd(double *, double *);
        void getConcentration_vd(double *, double *);

	//Initialization of States
	double (*getInitialState) (COMPONENT,double*,L_STATE&,int,IF_PARAMS*);
	virtual boolean FT_StateStructAtGridCrossing_tmp(Front*, int*, GRID_DIRECTION, COMPONENT, Locstate*, HYPER_SURF**, double*, double t=0);
        virtual boolean nearest_interface_point_tmp(double*, COMPONENT, INTERFACE*, USE_BOUNDARIES, HYPER_SURF*, double*, double*, HYPER_SURF_ELEMENT**, HYPER_SURF**);

	//User interface
	virtual void printInteriorVelocity(char *outname) = 0;

	virtual void setInitialCondition(void) = 0;
        // for vd
        virtual void setInitialCondition_vd(LEVEL_FUNC_PACK*) = 0;
        virtual void setInitialCondition_RSRV_vd(LEVEL_FUNC_PACK*) = 0;
        virtual void setInitialCondition_RSSY_vd(LEVEL_FUNC_PACK*) = 0;
        virtual void setInitialDiffusionVelocity_vd(void) = 0;

	virtual void solve(double dt) = 0;             // main step function
        // for vd
        virtual void solve_vd(double dt) = 0;          //main step function

	virtual void computeError(void) = 0;

        virtual void computeRTParameters(double, char*, int) {};
	virtual double new_height_at_fraction_vd(double,int,double){};
	virtual void new_accumulate_fractions_in_layer_vd(double,double*){};

protected:
	Front *front;
	// On topological grid
	RECT_GRID *top_grid;
        RECT_GRID *comp_grid; //Computational grid
        PP_GRID *pp_grid;
	double *array;
	double *source;
    //removal tag: HAOZ this is for REFLECTION BOUNDARY CONDITION
    double **vecarray;
        double *source_tmp;
	double *diff_coeff;
        double *diff_coeff_old;
	COMPONENT *top_comp;
	IF_PARAMS *iFparams;
	IF_FIELD  *field;

	int *top_gmax;
	int *lbuf, *ubuf;
	double *top_L, *top_U;
	int **ij_to_I, **I_to_ij;
	int ***ijk_to_I, ***I_to_ijk;

	// Sweeping limites
	int imin, jmin, kmin;
	int imax, jmax, kmax;

	//member data: mesh storage
	std::vector<L_RECTANGLE>   cell_center;

	//member data:
	int    m_comp[2];
	double m_mu[2];
	double m_rho[2];//two components at most
	double m_sigma; //surface tension
	double m_smoothing_radius;//used by smoothing function

	double hmin; //smallest spacing
	double mu_min; //smallest viscocity
	double rho_min;//smallest density
	double m_t;
	double m_t_old, m_t_int, m_t_new;
        // for vd, member data
        double m_c[2];
        double m_rho_old[2];
        double m_Dcoef[2];
        double c_min, c_max;
        double Dcoef_min,Dcoef_max,rho_max,mu_max;
        double z0; // midplane position
        double zmin_intfc,zmax_intfc; //top/bottom of contact interface
        double zmin_vf,zmax_vf; //top/bottom of n% volume fraction contour
        boolean bGhostCell;

	// for parallel partition
	int NLblocks, ilower, iupper;
	int *n_dist;

protected:
	void setComponent(void); //init components;
	void setDomain();
        // for vd
        void setDomain_vd();

	// parallelization related functions
	void scatMeshArray(void);
	void setGlobalIndex(void);
	void setIndexMap(void);

/*  These functions should be rewritten in 2D basis and 3D basis classes */
	virtual double getSmoothingFunction(double r) = 0; //Heaviside function
	virtual double getSmoothingFunctionD(double*, double*) = 0;
		//Heaviside function
	virtual double smoothedDeltaFunction(double*, double*) = 0;
	virtual double smoothedStepFunction(double*, double*, int) = 0;
	virtual void   sampleVelocity() = 0;
	virtual void   setSmoothedProperties(void) = 0;
		//smooth discontinuous properties
        virtual void   setSmoProOnePhase(void) = 0;
        // for vd
        virtual void   setSmoProOnePhase_vd(void) = 0;
        virtual void   setSmoothedProperties_vd(void) = 0;
                //smooth discontinuous properties

/****************  Functions related to solve() *********/

	virtual void copyMeshStates(void) = 0;
        //for vd
        virtual void copyMeshStates_vd(void) = 0;

	virtual void computeAdvection(void) = 0;
	//compute advection step in first order scheme

	virtual void compAdvectionTerm_coupled(void) = 0;
	//get 2nd order advection term for coupled system
	virtual void compAdvectionTerm_coupled_upgraded(void) = 0;
	//Upgraded algorithm for getting 2nd order advection term for coupled system
        // for vd
        virtual void compAdvectionTerm_coupled_vd(int) = 0;
        virtual void compAdvectionTerm_decoupled_vd(int) = 0;

	//The following 2 functions are for decoupled system
	virtual void compAdvectionTerm_decoupled(void) = 0;
	//get 2nd order advection term
	virtual void compAdvectionTerm_decoupled_upgraded(void) = 0;
	//Upgraded algorithm for getting 2nd order advection term


	virtual void compDiffWithSmoothProperty_1st_coupled(void) = 0; //1st order coupled diffusion solver
	virtual void compDiffWithSmoothProperty_1st_decoupled(void) = 0;
	virtual void compDiffWithSmoothProperty_2nd_coupled(void) = 0;
	virtual void compDiffWithSmoothProperty_2nd_decoupled(void) = 0;//2nd order decoupled diffusion solver
        // for vd
        virtual void compDiffWithSmoothProperty_velocity_decoupled_vd(void) = 0;
        virtual void compDiffWithSmoothProperty_velocity_vd(void) = 0; // 2nd order coupled diffusion solver

	virtual void computeProjection(void) = 0;
        // for vd
        virtual void computeProjection_vd(void) = 0;
        virtual void computeProjection_MAC_vd(void) = 0;

	virtual void computePressure(void) = 0;
	virtual void computePressurePmI(void) = 0;
	virtual void computePressurePmII(void) = 0;
	virtual void computePressurePmIII(void) = 0;
	virtual void computeGradientQ(void) = 0;
        // for vd
        virtual void computePressurePmII_vd(void) = 0;
        virtual void computeGradientQ_MAC_vd(void) = 0;

	virtual void computeNewVelocity(void) = 0;
        // for vd
        virtual void computeNewVelocity_vd(void) = 0;
        virtual void computeNewDensity_vd(int) = 0;
        virtual void computeNewConcentration_vd(void) = 0;

	virtual void computeSourceTerm(double *coords, L_STATE &state) = 0;
	virtual void computeSourceTerm(double *coords, double t, L_STATE
		&state) = 0;
        virtual void computeSourceTerm_Adv(double *coords, L_STATE &state) = 0;

        virtual void getExactSolution(double *coords,double t,L_STATE &state) = 0;
        // for vd
        virtual void getExactSolution_vd(double *coords,double t,L_STATE &state) = 0;

	virtual void surfaceTension(double*, HYPER_SURF_ELEMENT*,
		HYPER_SURF*, double*, double) = 0;
        virtual void surfaceTension_Peskin(void) = 0;

	virtual void computeSubgridModel(void) = 0;    // subgrid model

/***********************  Utility functions  *******************/

	void   computeExactSolution(double *coords, L_STATE &state);
	void   getRectangleIndex(int indexRectangle, int &i, int &j);
	void   getRectangleIndex(int indexRectangle, int &i, int &j, int &k);
	int    getRectangleComponent(int index);	// the center component
	void   getRectangleCenter(int index, double *coords);
	double getDistance(double *coords0, double *coords1);
	int    getComponent(int *icoords);
	int    getComponent(double *coords);
	void   save(char *filename);

/************* TMP Functions which are not implemented or used ***********/

	void getNearestInterfacePoint(COMPONENT,double*,double*,double*,
					double*);
};






///////////////Interface for Embedded Boundary Method////////////////////

class Incompress_Solver_EBM:public Incompress_Solver_Basis{
public:
        Incompress_Solver_EBM(Front &front) {};//constructor
	~Incompress_Solver_EBM() {};

	virtual void setInitialCondition() {};
	virtual void solve(double dt) {};
};




/////////////////////////////////////////////////////////////////////////////////

class Incompress_Solver_Smooth_2D_Basis:
public 	Incompress_Solver_Smooth_Basis{


public:
        Incompress_Solver_Smooth_2D_Basis(Front &front):
	Incompress_Solver_Smooth_Basis(front) {};
	virtual ~Incompress_Solver_Smooth_2D_Basis() {};

        virtual void printExpandedMesh(char* outname, bool binary) = 0;
	virtual void printInteriorVelocity(char *outname) = 0;

	virtual void setInitialCondition(void) = 0;
        // for vd
        virtual void setInitialCondition_vd(LEVEL_FUNC_PACK*) = 0;
        virtual void setInitialCondition_RSRV_vd(LEVEL_FUNC_PACK*) = 0;
        virtual void setInitialCondition_RSSY_vd(LEVEL_FUNC_PACK*) = 0;
        virtual void setInitialDiffusionVelocity_vd(void) = 0;

        virtual void solve(double dt) = 0;
        // for vd
        virtual void solve_vd(double dt) = 0;

	void computeError(void) {};
protected:
	double getSmoothingFunction(double r);
	double getSmoothingFunctionD(double*, double*);
	double smoothedDeltaFunction(double*, double*);
	double smoothedStepFunction(double*, double*, int);
	void sampleVelocity();

	void setSmoothedProperties(void);
        void setSmoProOnePhase(void);
        // for vd
        void setSmoProOnePhase_vd(void);
        void setSmoothedProperties_vd(void);

	virtual void copyMeshStates(void) = 0;
        // for vd
        virtual void copyMeshStates_vd(void) = 0;

	virtual void computeAdvection(void) = 0;
	virtual void compAdvectionTerm_coupled(void) = 0;
	virtual void compAdvectionTerm_coupled_upgraded(void) = 0;
	virtual void compAdvectionTerm_decoupled(void) = 0;
	virtual void compAdvectionTerm_decoupled_upgraded(void) = 0;
        // for vd
        virtual void compAdvectionTerm_coupled_vd(int) = 0;
        virtual void compAdvectionTerm_decoupled_vd(int) = 0;

	virtual void compDiffWithSmoothProperty_1st_coupled(void) = 0;
	//1st order coupled diffusion solver
	virtual void compDiffWithSmoothProperty_1st_decoupled(void) = 0;
	virtual void compDiffWithSmoothProperty_2nd_coupled(void) = 0;
	virtual void compDiffWithSmoothProperty_2nd_decoupled(void) = 0;
	//2nd order decoupled diffusion solver
        // for vd
        virtual void compDiffWithSmoothProperty_velocity_decoupled_vd(void) = 0;
        virtual void compDiffWithSmoothProperty_velocity_vd(void) = 0;

	virtual void computeProjection(void) = 0;
        // for vd
        virtual void computeProjection_vd(void) = 0;
        virtual void computeProjection_MAC_vd(void) = 0;

	virtual void computePressure(void) = 0;
	virtual void computePressurePmI(void) = 0;
	virtual void computePressurePmII(void) = 0;
	virtual void computePressurePmIII(void) = 0;
	virtual void computeGradientQ(void) = 0;
        // for vd
        virtual void computePressurePmII_vd(void) = 0;
        virtual void computeGradientQ_MAC_vd(void) = 0;

	virtual void computeNewVelocity(void) = 0;
        // for vd
        virtual void computeNewVelocity_vd(void) = 0;
        virtual void computeNewDensity_vd(int) = 0;
        virtual void computeNewConcentration_vd(void) = 0;

	virtual void computeSourceTerm(double *coords, L_STATE &state) = 0;
	virtual void computeSourceTerm(double *coords, double t, L_STATE
		&state) = 0;
        virtual void computeSourceTerm_Adv(double *coords, L_STATE &state) = 0;

        virtual void getExactSolution(double *coords,double t,L_STATE &state) = 0;
        virtual void getExactSolution_vd(double *coords,double t,L_STATE &state) = 0;

	virtual void surfaceTension(double*, HYPER_SURF_ELEMENT*,
		HYPER_SURF*, double*, double) = 0; // ?????????
        virtual void surfaceTension_Peskin(void) = 0;

	virtual void computeSubgridModel(void) = 0;    // subgrid model
};


class Incompress_Solver_Smooth_3D_Basis:
public 	Incompress_Solver_Smooth_Basis{

public:
        Incompress_Solver_Smooth_3D_Basis(Front &front):
	Incompress_Solver_Smooth_Basis(front) {};
	virtual ~Incompress_Solver_Smooth_3D_Basis() {};

        virtual void printExpandedMesh(char* outname, bool binary) = 0;
	virtual void printInteriorVelocity(char *outname) = 0;

	virtual void setInitialCondition(void) = 0;
        // for vd
        virtual void setInitialCondition_vd(LEVEL_FUNC_PACK*) = 0;
        virtual void setInitialCondition_RSRV_vd(LEVEL_FUNC_PACK*) = 0;
        virtual void setInitialCondition_RSSY_vd(LEVEL_FUNC_PACK*) = 0;
        virtual void setInitialDiffusionVelocity_vd(void) = 0;

	virtual void solve(double dt) = 0;
        // for vd
        virtual void solve_vd(double dt) = 0;

	virtual void computeError(void) = 0;
protected:
	double getSmoothingFunction(double r);
	double getSmoothingFunctionD(double*, double*);
	double smoothedDeltaFunction(double*, double*);
	double smoothedStepFunction(double*, double*, int);
	void sampleVelocity();

	void setSmoothedProperties(void);
        void setSmoProOnePhase(void);
        // for vd
        void setSmoProOnePhase_vd(void);
        void setSmoothedProperties_vd(void);


	virtual void copyMeshStates(void) = 0;
        // for vd
        virtual void copyMeshStates_vd(void) = 0;

	virtual void computeAdvection(void) = 0;
	virtual void compAdvectionTerm_coupled(void) = 0;
	virtual void compAdvectionTerm_coupled_upgraded(void) = 0;
	virtual void compAdvectionTerm_decoupled(void) = 0;
	virtual void compAdvectionTerm_decoupled_upgraded(void) = 0;
        // for vd
        virtual void compAdvectionTerm_coupled_vd(int) = 0;
        virtual void compAdvectionTerm_decoupled_vd(int) = 0;

	virtual void compDiffWithSmoothProperty_1st_coupled(void) = 0;
	//1st order coupled diffusion solver
	virtual void compDiffWithSmoothProperty_1st_decoupled(void) = 0;
	virtual void compDiffWithSmoothProperty_2nd_coupled(void) = 0;
	virtual void compDiffWithSmoothProperty_2nd_decoupled(void) = 0;
	//2nd order decoupled diffusion solver
        // for vd
        virtual void compDiffWithSmoothProperty_velocity_decoupled_vd(void) = 0;
        virtual void compDiffWithSmoothProperty_velocity_vd(void) = 0;

	virtual void computeProjection(void) = 0;
        // for vd
        virtual void computeProjection_vd(void) = 0;
        virtual void computeProjection_MAC_vd(void) = 0;

	virtual void computePressure(void) = 0;
	virtual void computePressurePmI(void) = 0;
	virtual void computePressurePmII(void) = 0;
	virtual void computePressurePmIII(void) = 0;
	virtual void computeGradientQ(void) = 0;
        // for vd
        virtual void computePressurePmII_vd(void) = 0;
        virtual void computeGradientQ_MAC_vd(void) = 0;

	virtual void computeNewVelocity(void) = 0;
        // for vd
        virtual void computeNewVelocity_vd(void) = 0;
        virtual void computeNewDensity_vd(int) = 0;
        virtual void computeNewConcentration_vd(void) = 0;

	virtual void computeSourceTerm(double *coords, L_STATE &state) = 0;
	virtual void computeSourceTerm(double *coords, double t, L_STATE
		&state) = 0;
        virtual void computeSourceTerm_Adv(double *coords, L_STATE &state) = 0;

        virtual void getExactSolution(double *coords,double t,L_STATE &state) = 0;
        // for vd
        virtual void getExactSolution_vd(double *coords,double t,L_STATE &state) = 0;

	virtual void surfaceTension(double*, HYPER_SURF_ELEMENT*,
		HYPER_SURF*, double*, double) = 0;
        virtual void surfaceTension_Peskin(void) = 0;

	virtual void computeSubgridModel(void) = 0;  // subgrid model
};

class Incompress_Solver_Smooth_2D_Cartesian:
public 	Incompress_Solver_Smooth_2D_Basis{

public:
        Incompress_Solver_Smooth_2D_Cartesian(Front &front):
	Incompress_Solver_Smooth_2D_Basis(front) {};
	~Incompress_Solver_Smooth_2D_Cartesian() {};

        void printExpandedMesh(char* outname, bool binary) {};
	void printInteriorVelocity(char* outname) {};
	void setInitialCondition(void);
        // for vd
        void setInitialCondition_vd(LEVEL_FUNC_PACK*);
        void setInitialCondition_RSRV_vd(LEVEL_FUNC_PACK*) {};
        void setInitialCondition_RSSY_vd(LEVEL_FUNC_PACK*) {};
        void setInitialDiffusionVelocity_vd(void) {};

	void solve(double dt);
        // for vd
        void solve_vd(double dt);

protected:
	void copyMeshStates(void);
        // for vd
        void copyMeshStates_vd(void);

	void computeAdvection(void); //first order advection

	void compAdvectionTerm_coupled(void); //not implemented yet????
	void compAdvectionTerm_coupled_upgraded(void); //not implemented yet????
	void compAdvectionTerm_decoupled(void);
	void compAdvectionTerm_decoupled_upgraded(void);
        // for vd
        void compAdvectionTerm_coupled_vd(int);
        void compAdvectionTerm_decoupled_vd(int) {};

	void compDiffWithSmoothProperty_1st_coupled(void);
	void compDiffWithSmoothProperty_1st_decoupled(void);
	void compDiffWithSmoothProperty_2nd_coupled(void); //Not implemented yet????
	void compDiffWithSmoothProperty_2nd_decoupled(void);
        // for vd
        void compDiffWithSmoothProperty_velocity_decoupled_vd(void);
        void compDiffWithSmoothProperty_velocity_vd(void);

	void computeProjection(void);
	void computeProjection_new(void);
        // for vd
        void computeProjection_vd(void);
        void computeProjection_MAC_vd(void){};

	void computePressure(void);
	void computePressurePmI(void);
	void computePressurePmII(void);
	void computePressurePmIII(void);
	void computeGradientQ(void);
        // for vd
        void computePressurePmII_vd(void){};
        void computeGradientQ_MAC_vd(void){};

	void computeNewVelocity(void);
        // for vd
        void computeNewVelocity_vd(void);
        void computeNewDensity_vd(int);
        void computeNewConcentration_vd(void);

	virtual void computeSourceTerm(double *coords, L_STATE &state);
        virtual void computeSourceTerm_Adv(double *coords, L_STATE &state);
	virtual void computeSourceTerm(double *coords, double t, L_STATE &state);

        void getExactSolution(double *coords,double t,L_STATE &state){};
        // for vd
        void getExactSolution_vd(double *coords,double t,L_STATE &state){};

	void surfaceTension(double*, HYPER_SURF_ELEMENT*, HYPER_SURF*, double*, double); //??????????????????
        void surfaceTension_Peskin(void) {};
        //Not implemented

	void computeSubgridModel(void);    // subgrid model

	/***************   Low level computation functions  *************/
	double computeFieldPointDiv(int*, double**);
	void   computeFieldPointGrad(int*, double*, double*);
        void   computeFieldPointGradPhi(int*, double*, double*);
        // for vd
        void   computeFieldPointGradRho(int*, double*, double*);

	double computeFieldPointCurl(int*, double**, double*);
	double getVorticity(int i, int j);

	//---------------------------------------------------------------
	//         utility functions for the advection step
	//---------------------------------------------------------------

	void getAdvectionTerm_decoupled(int *icoords, double convectionTerm[2]);
	//get second-order advection term at each cell
	void getAdvectionTerm_decoupled_upgraded(int *icoords, double convectionTerm[2]);
	//Upgraded version of getAdvectionTerm
	void getAdvectionTerm_coupled(int *icoords, double convectionTerm[2]);
	//get second-order advection term at each cell with cross derivative terms
	void getAdvectionTerm_coupled_upgraded(int *icoords, double convectionTerm[2]);
	//Upgraded version of getAdvectionTerm_coupled
        // for vd
        void getStatesBar_coupled_vd(int*, double);
        void computeMacPhi_vd(int);
        void computeNewUbar_vd(double,int);

	void getFaceVelocity_middleStep(int *icoords,GRID_DIRECTION dir, L_STATE &state_face);
	void getFaceVelocity_middleStep_hat(int *icoords,GRID_DIRECTION dir,L_STATE &state_hat);
	//For the upgraded algorithm
	void getFaceVelocity_middleStep_bar(int *icoords,GRID_DIRECTION dir,L_STATE &state_bar, double transverseD[2], L_STATE state_hat);
	//For the upgraded algorithm

	void getFaceVelocity_middleStep_coupled(int *icoords,GRID_DIRECTION dir, L_STATE &state_face); //for variable mu
	void getFaceVelocity_middleStep_coupled_hat(int *icoords,GRID_DIRECTION dir,L_STATE &state_hat) {};
	//For the upgraded algorithm with variable mu
	void getFaceVelocity_middleStep_coupled_bar(int *icoords,GRID_DIRECTION dir,L_STATE &state_bar, double transverseD[2], L_STATE state_hat);
	//For the upgraded algorithm with variable mu
        // for vd
        void getFaceState_middleStep_hat_vd(int *icoords,GRID_DIRECTION dir,L_STATE &state_hat);
        void getFaceState_middleStep_bar_vd(int *icoords,GRID_DIRECTION dir,L_STATE &state_bar, double transverseD[4], L_STATE state_hat);

	void getDifffusion(int *icoords,double diffusion[2]); //Get the diffusion terms
	void getDiffusion_coupled(int *icoords, double diffusion[2]); //Get the diffusion terms for variable mu
        // for vd
        void getDiffusion_coupled_vd(int *icoords, double diffusion[4]);
        void getDivU_coupled_vd(int *icoords, double diffusion[4], int flag);
        void getDiffusionC_coupled_vd(int *icoords, double diffusion[4]);

	void getDU2(int *icoords,EBM_COORD xyz,double dU2[2]); //Get the U_xx, U_yy, U_zz
	void getLimitedSlope(int *icoords,EBM_COORD xzy,double slope[2]); //mimmod slope limiter
	void getLimitedSlope_Vanleer(int *icoords,EBM_COORD xyz, double slope[2]); //Van Leer slope limiter
        // for vd
        void getLimitedSlope_vd(int *icoords,EBM_COORD xzy,double slope[4]);

	double EBM_minmod(double x, double y); //minmod function

	bool getNeighborOrBoundaryState(int icoords[2],GRID_DIRECTION dir,L_STATE &state,double t); //get the neighbor state or boundary state
        // for vd
        bool getNeighborOrBoundaryState_vd(int icoords[2],GRID_DIRECTION dir,L_STATE &state,double t);

	void getRiemannSolution(EBM_COORD xyz,L_STATE &u_left,L_STATE &state_right,L_STATE &ans); //Compute Riemann solution using left and right state
        // for vd
        void getRiemannSolution_vd(EBM_COORD xyz,L_STATE &u_left,L_STATE &state_right,L_STATE &ans);

	//void   computeVelDivergence(void);
	//void getVelocityGradient(double *p,double *gradU,double *gradV);

    //Smeeton Youngs' Experiment 105
    double computeFieldPointDiv_MAC_vd(int*, double**); // 2D version not implemented yet. TODO && FIXME: In Progress.
    void   computeFieldPointGrad_MAC_vd(int*, double*, double*) {};
    void getDivU_MAC_vd(int *icoords, double *diffusion, int flag, boolean) {};

};

class Incompress_Solver_Smooth_3D_Cartesian:
public 	Incompress_Solver_Smooth_3D_Basis{

public:
        Incompress_Solver_Smooth_3D_Cartesian(Front &front):
	Incompress_Solver_Smooth_3D_Basis(front) {};
	~Incompress_Solver_Smooth_3D_Cartesian() {};
    // enforce Reflection Boundary Condition
    void enforceVecState(double**);
    void checkBoundaryCondition(GRID_DIRECTION,int*,int*,double,COMPONENT);
    bool FT_Reflect(int*,int,int);

        void printExpandedMesh(char* outname,bool binary);
        void printExpandedMesh_big_endian(char* outname);
        void printExpandedMesh_little_endian(char* outname);
        void printExpandedMesh_ascii(char* outname);

	void printInteriorVelocity(char* outname);

	void setInitialCondition(void);
        // for vd
        void setInitialCondition_vd(LEVEL_FUNC_PACK*);
        void setInitialCondition_RSRV_vd(LEVEL_FUNC_PACK*);
        void setInitialCondition_RSSY_vd(LEVEL_FUNC_PACK*);
        void setInitialDiffusionVelocity_vd(void) {};

	void solve(double dt);
        // for vd
        void solve_vd(double dt);

	void computeError(void) {};
protected:
	void copyMeshStates(void);
        // for vd
        void copyMeshStates_vd(void);

	void computeAdvection(void);
	void compAdvectionTerm_coupled(void);
	void compAdvectionTerm_coupled_upgraded(void);
	void compAdvectionTerm_decoupled(void);
	void compAdvectionTerm_decoupled_upgraded(void);
        // for vd
        void compAdvectionTerm_coupled_vd(int);
        void compAdvectionTerm_decoupled_vd(int);
        void compAdvectionTerm_MAC_decoupled_vd(int, boolean);

	void compDiffWithSmoothProperty_1st_coupled(void);
	void compDiffWithSmoothProperty_1st_decoupled(void);
	void compDiffWithSmoothProperty_2nd_coupled(void);
	void compDiffWithSmoothProperty_2nd_decoupled(void);
        // for vd
        void compDiffWithSmoothProperty_velocity_decoupled_vd(void);
        void compDiffWithSmoothProperty_velocity_MAC_decoupled_vd(void);
        void compDiffWithSmoothProperty_velocity_MAC_decoupled_zeroW_vd(void);
        void compDiffWithSmoothProperty_velocity_MAC_decoupled_zeroV_vd(void);
        void compDiffWithSmoothProperty_velocity_MAC_coupled_vd(void);
        void compDiffWithSmoothProperty_velocity_vd(void);

	void computeProjection(void);
        // for vd
        void computeProjection_vd(void);
        void computeProjection_Expand_vd(void);
        void computeProjection_MAC_vd(void);
        void computeProjection_MAC_PPE_vd(void);

	void computePressure(void);
	void computePressurePmI(void);
	void computePressurePmII(void);
	void computePressurePmIII(void);
	void computeGradientQ(void);
        // for vd
        void computePressure_vd(void);
        void computePressure_MAC_vd(void);
        void computePressure_MAC_zeroW_vd(void);
        void computePressure_MAC_zeroV_vd(void);
        void computePressurePmII_vd(void);
        void computeGradientQ_MAC_vd(void);
        void computeGradientQ_MAC_zeroW_vd(void);
        void computeGradientQ_MAC_zeroV_vd(void);

	void computeNewVelocity(void);
        // for vd
        void computeNewVelocity_vd(void);
        void computeNewVelocity_Filter_vd(void);
        void computeNewVelocity_MAC_vd(void);
        void computeNewVelocity_fullMAC_vd(void);
        void computeNewVelocity_fullMAC_zeroW_vd(void);
        void computeNewVelocity_fullMAC_zeroV_vd(void);
        void computeNewDensity_vd(int);
        void computeImmiscibleConcentration_vd(int);
        void updateImmiscibleDensity_vd();
        void computeNewDensityByConcentration_vd();
        void computeNewDensity_MAC_constRho_vd();
        void computeNewConcentration_vd(void);

        // for vd
        void compFilter_Ud_vd(void);
        void compFilter_VertexProjection_vd(void);
        void compFilter_VertexRhoProjection_vd(void);
        void compFilter_EdgeProjection_vd(void);
        void compFilter_FaceProjection_vd(void);
        void compFilter_FaceRhoProjection_vd(void);
        void compFilter_FaceProjection_Asymmetric_vd(void);

	virtual void computeSourceTerm(double *coords, L_STATE &state);
        virtual void computeSourceTerm_Adv(double *coords, L_STATE &state);
	virtual void computeSourceTerm(double *coords, double t, L_STATE &state);

        void getExactSolution(double *coords,double t,L_STATE &state){};
        // for vd
        void getExactSolution_vd(double *coords,double t,L_STATE &state){};

	void surfaceTension(double*, HYPER_SURF_ELEMENT*, HYPER_SURF*, double*, double);
        void surfaceTension_Peskin(void);// transformed Oak Ridge Version in Cylindrical3D
        void compSurfaceTension_Tri(TRI*);

	void computeSubgridModel(void) {};    // subgrid model
        void computeSubgridModel_vd(void);    // subgrid model

	double computeFieldPointDiv(int*, double**);
        double computeFieldPointDiv_Neumann_vd(int*, double**);
        double computeFieldPointDiv_General_vd(int*, double**);
        double computeFieldPointDiv_MAC_vd(int*, double**);
	void   computeFieldPointGrad(int*, double*, double*);
        void   computeFieldPointGradPhi(int*, double*, double*);
        // for vd
        void   computeFieldPointGradRho(int*, double*, double*);
        void   computeFieldPointGradRho_MAC_vd(int*, double*, double*);
        void   computeFieldPointGradRho_GhostDensity(int*, double*, double*);
        void   computeFieldPointGrad_MAC_vd(int*, double*, double*);

        void   computeRTParameters(double, char*, int);
        double new_height_at_fraction_vd(double,int,double);
        void   new_accumulate_fractions_in_layer_vd(double,double*);
	double computeFieldPointCurl(int*, double**, double*);
	double getVorticityX(int i, int j, int k);
	double getVorticityY(int i, int j, int k);
	double getVorticityZ(int i, int j, int k);

	//void   computeVelDivergence(void);
	//void getVelocityGradient(double *p,double *gradU,double *gradV);
	//-------------------------------------------------
	//	utility function for the advection step
	//-------------------------------------------------
	void getAdvectionTerm_decoupled(int *icoords, double convectionTerm[3]);
	void getAdvectionTerm_decoupled_upgraded(int *icoords, double convectionTerm[3]);
	void getAdvectionTerm_coupled(int *icoords, double convectionTerm[3]);
	void getAdvectionTerm_coupled_upgraded(int *icoords, double convectionTerm[3]);
	//not implemented yet
        // for vd
        void getVelocityBar_coupled_vd(int*, double);
        void getVelocityBar_decoupled_vd(int*, double);
        void getCellCenterVelocityBar_MAC_decoupled_vd(int*, double);
        void getCellFaceVelocityBar_MAC_decoupled_vd(int*, double);
        void getCellEdgeVelocityBar_MAC_decoupled_vd(int*, double);
        void getVelocityBar_Normal_coupled_vd(int*, double);
        void getScalarBar_coupled_vd(int*, double);
        void getScalarBar_MAC_vd(int*, double, boolean);
        void getStateBar_coupled_vd(int*, double);
        void updateDensityBar_GhostStates_vd(int*, double);
        void computeMacPhi_vd(int);
        void computeMacPhi_MAC_vd(int);
        void computeNewUbar_vd(double,int);
        void computeNewUbar_Normal_vd(double,int);
        void computeNewUbar_Tangential_vd(double,int);
        void computeNewU_MAC_vd(double);
        void computeNewUFaceBar_MAC_vd(double,int);

	void getFaceVelocity_middleStep(int *icoords,GRID_DIRECTION dir,L_STATE &state_face);
	void getFaceVelocity_middleStep_hat(int *icoords,GRID_DIRECTION dir, L_STATE &state_hat);
	void getFaceVelocity_middleStep_bar(int *icoords,GRID_DIRECTION dir, L_STATE &state_bar, double transverseD[3], L_STATE state_hat);

	void getFaceVelocity_middleStep_coupled(int *icoords,GRID_DIRECTION dir,L_STATE &state_face);
	void getFaceVelocity_middleStep_coupled_hat(int *icoords,GRID_DIRECTION dir, L_STATE &state_hat) {};
	void getFaceVelocity_middleStep_coupled_bar(int *icoords,GRID_DIRECTION dir, L_STATE &state_bar, double transverseD[3], L_STATE state_hat);
        // for vd
        void getFaceState_middleStep_hat_vd(int *icoords,GRID_DIRECTION dir,L_STATE &state_hat);
        void getFaceVelocity_middleStep_hat_vd(int *icoords,GRID_DIRECTION dir,L_STATE &state_hat);
        void getFaceScalar_middleStep_hat_vd(int *icoords,GRID_DIRECTION dir,L_STATE &state_hat);
        void getFaceScalar_MAC_middleStep_hat_vd(int *icoords,EBM_COORD xyz,GRID_DIRECTION dir,L_STATE &state_hat,boolean);
        void getFaceState_middleStep_bar_vd(int *icoords,GRID_DIRECTION dir,L_STATE &state_bar, double transverseD[5], L_STATE state_hat);
        void getFaceVelocity_middleStep_bar_vd(int *icoords,GRID_DIRECTION dir,L_STATE &state_bar, double transverseD[5], L_STATE state_hat);
        void getFaceVelocity_middleStep_bar_decoupled_vd(int *icoords,GRID_DIRECTION dir,L_STATE &state_bar, double transverseD[5], L_STATE state_hat);
        void getCenterVelocity_MAC_middleStep_hat_vd(int *icoords,GRID_DIRECTION dir,L_STATE &state_hat);
        void getFaceVelocity_MAC_middleStep_hat_vd(int *icoords,EBM_COORD xyz,L_STATE &state_hat);
        void getEdgeVelocity_MAC_middleStep_hat_vd(int *icoords,EBM_COORD xyz,GRID_DIRECTION dir,L_STATE &state_hat);
        void getCenterVelocity_MAC_zBdry_middleStep_hat_vd(int *icoords,GRID_DIRECTION dir,L_STATE &state_hat);
        void getVelocity_MAC_middleStep_bar_decoupled_vd(int *icoords,EBM_COORD xyz,GRID_DIRECTION dir,L_STATE &state_bar,L_STATE state_hat);
        void getScalar_MAC_middleStep_bar_decoupled_vd(int *icoords,EBM_COORD xyz,GRID_DIRECTION dir,L_STATE &state_bar,L_STATE state_hat,boolean);
        void getFaceScalar_middleStep_bar_vd(int *icoords,GRID_DIRECTION dir,L_STATE &state_bar, double transverseD[5], L_STATE state_hat);
        void updateFaceDensity_GhostStates_hat_vd(int *icoords,GRID_DIRECTION dir,L_STATE,L_STATE,L_STATE &);

	void getDifffusion(int *icoords,double diffusion[3]);
	void getDiffusion_coupled(int *icoords, double diffusion[3]);
        // for vd
        void getViscousTerm_coupled_vd(int *icoords, double diffusion[5]);
        void getViscousTerm_decoupled_vd(int *icoords, double diffusion[5]);
        void getViscousTerm_MAC_decoupled_vd(int *icoords, EBM_COORD xyz, double diffusion[3]);
        void getViscousTerm_coupled_tmp_vd(int *icoords, double diffusion[5]);
        void getTransverseDTerm_Velocity_MAC_vd(int *icoords, double transverseD[3]);
        void getTransverseDTerm_Scalar_MAC_vd(int *icoords, EBM_COORD xyz, double transverseD[2], boolean);
        void getDivU_coupled_vd(int *icoords, double diffusion[5], int flag);
        void getDivU_MAC_vd(int *icoords, double *diffusion, int flag, boolean);
        void getDivU_GhostStates_coupled_vd(int *icoords, double diffusion[5], int flag);
        void getDiffusionC_coupled_vd(int *icoords, double diffusion[5]);
        void getDiffusionC_MAC_vd(int *icoords, double *diffusion, boolean);

	void getDU2(int *icoords,EBM_COORD xyz,double dU2[3]);
	void getLimitedSlope(int *icoords,EBM_COORD xzy,double slope[3]);
	void getLimitedSlope_Vanleer(int *icoords, EBM_COORD xyz, double slope[3]);
        // for vd
        void getLimitedSlope_vd(int *icoords,EBM_COORD xzy,double slope[6]);
        void getLimitedSlope_Velocity_MAC_vd(int *icoords,EBM_COORD xzy,double slope[3]);
        void getLimitedSlope_Scalar_MAC_vd(int *icoords,EBM_COORD xzy,double slope[2],boolean);

        void getCenterDifference_vd(int *icoords,EBM_COORD xzy,double slope[6]);

	double EBM_minmod(double x, double y);

	bool getNeighborOrBoundaryState(int icoords[3],GRID_DIRECTION dir,L_STATE &state,double t);
        // for vd
        bool getNeighborOrBoundaryState_vd(int icoords[3],GRID_DIRECTION dir,L_STATE &state,double t);
        bool getNeighborOrBoundaryState_tmp_vd(int icoords[3],GRID_DIRECTION dir,L_STATE &state,double t);
        int  getNeighborOrBoundaryScalar_MAC_vd(int icoords[3],GRID_DIRECTION dir,L_STATE &state,double t);
        int  getNeighborOrBoundaryScalar_MAC_GhostCell_vd(int icoords[3],GRID_DIRECTION dir,L_STATE &state,double t);
        void getGhostDensity_vd(int icoords[3],COMPONENT,GRID_DIRECTION dir,L_STATE &state,double t);

	void getRiemannSolution(EBM_COORD xyz,L_STATE &u_left,L_STATE &state_right,L_STATE &ans);
        // for vd
        void getRiemannSolution_vd(EBM_COORD xyz,L_STATE &u_left,L_STATE &state_right,L_STATE &ans);
        void getRiemannSolution_Velocity_vd(EBM_COORD xyz,L_STATE &u_left,L_STATE &state_right,L_STATE &ans, int *icoords, int *ICoords);
        void getRiemannSolution_Velocity_hat_vd(EBM_COORD xyz,L_STATE &u_left,L_STATE &state_right,L_STATE &ans);
        void getRiemannSolution_Scalar_vd(L_STATE &u_left,L_STATE &state_right,L_STATE &ans, int *icoords, GRID_DIRECTION dir);
        void getRiemannSolution_localMAC_vd(L_STATE &u_left,L_STATE &state_right,L_STATE &ans, int *icoords, GRID_DIRECTION dir);
        void getRiemannSolution_MAC_CenterVelocity_vd(EBM_COORD xyz,L_STATE &u_left,L_STATE &u_right,L_STATE &ans, int *icoords);
        void getRiemannSolution_MAC_EdgeVelocity_vd(EBM_COORD xyz1,EBM_COORD xyz2,L_STATE &u_left,L_STATE &u_right,L_STATE &u_LEFT,L_STATE &u_RIGHT,L_STATE &ans);
        void getRiemannSolution_MAC_Scalar_vd(L_STATE &u_left,L_STATE &u_right,L_STATE &ans, int *icoords, GRID_DIRECTION dir);
        //Reflection Boundary Condition
        void ReflectBC(double**);
        void NeumannBC(double**);
};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////





class Incompress_Solver_Smooth_3D_Cylindrical:public Incompress_Solver_Smooth_3D_Basis{
public:
        Incompress_Solver_Smooth_3D_Cylindrical(Front &front):Incompress_Solver_Smooth_3D_Basis(front) {};
	~Incompress_Solver_Smooth_3D_Cylindrical() {};

        void printExpandedMesh(char *out_name,bool binary){};
	void printInteriorVelocity(char *out_name);
	void setInitialCondition(void);
        // for vd
        void setInitialCondition_vd(LEVEL_FUNC_PACK*){};
        void setInitialCondition_RSRV_vd(LEVEL_FUNC_PACK*){};
        void setInitialCondition_RSSY_vd(LEVEL_FUNC_PACK*){};
        void setInitialDiffusionVelocity_vd(void) {};

	void solve(double dt);
        // for vd
        void solve_vd(double dt){};

	void computeError(void);

protected:
	void copyMeshStates(void);
        // for vd
        void copyMeshStates_vd(void){};

	void computeAdvection(void); //First order advection, operator splitting
	void computeAdvection_test(void); //Paper version, first order advection, operator splitting

	void compAdvectionTerm_coupled(void); //not implemented
	void compAdvectionTerm_coupled_upgraded(void) {}; //not implemented
        // for vd
        void compAdvectionTerm_coupled_vd(void){}; //not implemented
        void compAdvectionTerm_decoupled_vd(void){}; //not implemented

	void compAdvectionTerm_decoupled(void);
	void compAdvectionTerm_decoupled_upgraded(void) {}; //not implemented

	void compDiffWithSmoothProperty_1st_coupled(void) {}; //not implemented
	void compDiffWithSmoothProperty_1st_decoupled(void);
	void compDiffWithSmoothProperty_1st_decoupled_test(void); //Paper version
	void compDiffWithSmoothProperty_1st_decoupled_source(void); //half Crank-Nicholson

	void compDiffWithSmoothProperty_2nd_coupled(void); //no implemented
	void compDiffWithSmoothProperty_2nd_decoupled(void);
	//just adjusting coefficients, do not support complext B.C.
        // for vd
        void compDiffWithSmoothProperty_velocity_decoupled_vd(void){}; //no implemented
        void compDiffWithSmoothProperty_velocity_vd(void){}; //no implemented

	//---------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//   utility functions for diffusion and projection
	//----------------------------------------------------------------------------

	void compDiff_CellFace( //Set the coefficient for cell corners
		PETSc *pSolver,
		int I,
		int I_nb[18],
		double U_center[3],
		double U_nb[3][18],
		int flag[6],
		int equation_index,
		int vel_comp,
		int face_index,
		double coeff);


	void compDiff_CellCorner( //Set the coefficient for cell corners
		PETSc *pSolver,
		int I,
		int I_nb[18],
		double U_center[3],
		double U_nb[3][18],
		int flag[6],
		int equation_index,
		int vel_comp,
		int corner_index,
		double coeff);

	double getFaceArea(int *icoords, GRID_DIRECTION dir);
	double getCellVolume(int *icoords);
	void getFaceCenter(int *icoords, GRID_DIRECTION dir, double faceCenter[3]);

	void compDiffWithSmoothProperty_2nd_decoupled_Shuqiang(void);
	void compDiffWithSmoothProperty_cellFace(
		PETSc *pSolver,
		int *icoords,
		int I, int I_nb,
		GRID_DIRECTION dir,
		double dh,
		double faceArea,
		double cellVolume,
		double r,
		double mu,
		double rho,
		L_STATE &U_nb,
		L_STATE &U_nb_new,
		L_STATE &U_center);

	void compDiffWithSmoothProperty_Dirichlet(
		PETSc *pSolver,
		int *icoords,
		int I, int I_nb,
		GRID_DIRECTION dir,
		double dh,
		double faceArea,
		double cellVolume,
		double r,
		double mu,
		double rho,
		L_STATE &U_nb,
		L_STATE &U_nb_new,
		L_STATE &U_center);

	void compDiffWithSmoothProperty_cellInterior(
		PETSc *pSolver,
		int *icoords,
		int I,
		double cellVolume,
		double r,
		double mu,
		double rho,
		L_STATE &U_center);


	void computeProjection_Shuqiang(void);

	//--------------------------------------------------------

	void computeProjection(void);
        // for vd
        void computeProjection_vd(void){};
        void computeProjection_MAC_vd(void){};

	void computePressure(void);
	void computePressurePmI(void);
	void computePressurePmII(void);
	void computePressurePmIII(void);
	void computeGradientQ(void);
        // for vd
        void computePressurePmII_vd(void){};
        void computeGradientQ_MAC_vd(void){};

	void computeNewVelocity(void);
        // for vd
        void computeNewVelocity_vd(void){};
        void computeNewDensity_vd(void){};
        void computeNewConcentration_vd(void){};

	virtual void computeSourceTerm(double *coords, L_STATE &state);
	virtual void computeSourceTerm(double *coords, double t, L_STATE &state);
        void computeSourceTerm_Adv(double *coords, L_STATE &state){};

        void getExactSolution(double *coords,double t,L_STATE &state){};
        // for vd
        void getExactSolution_vd(double *coords,double t,L_STATE &state){};

	void surfaceTension(double*, HYPER_SURF_ELEMENT*, HYPER_SURF*, double*, double);
        void surfaceTension_Peskin(void){};

	void computeSubgridModel(void);    // subgrid model


	double computeFieldPointDiv(int*, double**);
	void   computeFieldPointGrad(int*, double*, double*);
	double computeFieldPointCurl(int*, double**, double*);
	double getVorticityX(int i, int j, int k);
	double getVorticityY(int i, int j, int k);
	double getVorticityZ(int i, int j, int k);

	//void   computeVelDivergence(void);
	//void getVelocityGradient(double *p,double *gradU,double *gradV);
	//-------------------------------------------------
	//	utility function for the advection step
	//-------------------------------------------------
	void getAdvectionTerm_decoupled(int *icoords, double convectionTerm[3]);
	void getAdvectionTerm_decoupled_upgraded(int *icoords, double convectionTerm[3]) {};//not implemented yet

	void getAdvectionTerm_coupled(int *icoords, double convectionTerm[3]);
	void getAdvectionTerm_coupled_upgraded(int *icoords, double convectionTerm[3]) {};
	//not implemented yet

	void getFaceVelocity_middleStep(int *icoords,GRID_DIRECTION dir,L_STATE &state_face);
	void getFaceVelocity_middleStep_hat(int *icoords,GRID_DIRECTION dir, L_STATE &state_hat) {}; //not implemented
	void getFaceVelocity_middleStep_bar(int *icoords,GRID_DIRECTION dir, L_STATE &state_bar, double transverseD[3], L_STATE state_hat) {}; //not implemented

	void getFaceVelocity_middleStep_coupled(int *icoords,GRID_DIRECTION dir,L_STATE &state_face);
	void getFaceVelocity_middleStep_coupled_hat(int *icoords,GRID_DIRECTION dir, L_STATE &state_hat) {}; //not implemented
	void getFaceVelocity_middleStep_coupled_bar(int *icoords,GRID_DIRECTION dir, L_STATE &state_bar, double transverseD[3], L_STATE state_hat) {}; //not implemented


	void getDifffusion(int *icoords,double diffusion[3]);
	void getDiffusion_coupled(int *icoords, double diffusinon[3]);

	void getDU2(int *icoords,EBM_COORD xyz,double dU2[3], double dU[3]);
	void getLimitedSlope(int *icoords,EBM_COORD xzy,double slope[3]);
	void getLimitedSlope_Vanleer(int *icoords, EBM_COORD xyz, double slope[3]);

	double EBM_minmod(double x, double y);

	bool getNeighborOrBoundaryState(int icoords[3],GRID_DIRECTION dir,L_STATE &state,double t);
	void getRiemannSolution(EBM_COORD xyz,L_STATE &u_left,L_STATE &state_right,L_STATE &ans);
};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////




extern double getStatePres(POINTER);
extern double getStateVort(POINTER);
extern double getStateXvel(POINTER);
extern double getStateYvel(POINTER);
extern double getStateZvel(POINTER);
extern double getStateXvort(POINTER);
extern double getStateYvort(POINTER);
extern double getStateZvort(POINTER);
// for vd
extern double getStateDens(POINTER);
extern double getStateDensOld(POINTER);
extern double getStateConc(POINTER);

extern double burger_flux(double,double,double);
extern double linear_flux(double,double,double,double);

extern void fluid_print_front_states(FILE*,Front*,bool);
extern void fluid_read_front_states(FILE*,Front*,bool);
//for vd
extern void fluid_print_front_states_vd(FILE*,Front*,bool);
extern void fluid_read_front_states_vd(FILE*,Front*,bool);

extern void read_iF_movie_options(char*,IF_PARAMS*);
extern void read_iF_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
extern boolean isDirichletPresetBdry(Front*,int*,GRID_DIRECTION,COMPONENT);
extern void convertGridDirectionToDirSide(GRID_DIRECTION, int*, int*);

#endif
