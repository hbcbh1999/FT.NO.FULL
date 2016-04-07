/**********************************************************************
 * 		crystal.h					      *
 **********************************************************************/

#ifndef _FT_CRYSTAL_H_
#define _FT_CRYSTAL_H_

#include <FronTier.h>
#include <vector>
#include <petscksp.h>
#include <assert.h>

#define         CRYSTAL_COMP            1
#define         SOLUTE_COMP             3

typedef struct {
        boolean grow_from_floor;
        boolean grow_from_ceiling;
        boolean grow_from_space;
        int dim;
        int num_floor_seeds;
        int num_ceiling_seeds;
        int num_space_seeds;
        double floor_level;
        double ceiling_level;
        double **floor_center;
        double **ceiling_center;
        double **space_center;
        double seed_radius;
        boolean add_space_seed_pert;  // 2D
        double nu;
        double amp;
        double phase;
        double point;           // 1D crystal boundary
} SEED_PARAMS;

enum _DF_SCHEME {
	UNSPLIT_EXPLICIT = 1,
	UNSPLIT_IMPLICIT,
	CRANK_NICOLSON
};
typedef enum _DF_SCHEME DF_SCHEME;

struct _CRT_FIELD {
        double *solute;
	double **vel;
};
typedef struct _CRT_FIELD CRT_FIELD;

enum _POINT_PROP_SCHEME {
        EXPLICIT_EULER = 1,
        IMPLICIT_EULER,
        MIDDLE_POINT,
	CONSTANT_STATE
};
typedef enum _POINT_PROP_SCHEME POINT_PROP_SCHEME;

struct _CRT_MOVIE_OPTION {
        boolean plot_pres;
        boolean plot_vort;
        boolean plot_velo;
        boolean plot_solute;
        boolean plot_cross_section[MAXD];  /* 3D 0: yz; 1: zx; 2: xy */
};
typedef struct _CRT_MOVIE_OPTION CRT_MOVIE_OPTION;

struct _CRT_PARAMS {
        int dim;
	DF_SCHEME num_scheme;
	POINT_PROP_SCHEME point_prop_scheme;
	CRT_MOVIE_OPTION *movie_option;
	boolean add_curvature;
	CRT_FIELD *field;	// field of solute concentration
	double C_0;		// Ambient concentration
        double C_eq;    /* solute concentration in equilibrium with solid */
        double rho_s;   /* density of the precipitated solid phase */
        double D;       /* diffusion coefficient of the solute concentration */
        double k;       /* local reaction rate coefficient */
	double max_solute,min_solute;
};
typedef struct _CRT_PARAMS CRT_PARAMS;

extern double   crystal_seed_curve(POINTER,double*);

typedef class C_CARTESIAN C_CARTESIAN_EB;

/******************************************************************************
 * 		lcartsn.h
 * A simple incompressible flow solver using the ghost fluid method and the
 * projection method.
 *
 * the main function is 
 * 	C_CARTESIAN::solve().
 *
 * References:
 ******************************************************************************/

class SOLVER;
class C_CARTESIAN;

//-------------------------------------------------
//		STATES
// NOTE:
//      C_STATE/C_STATE_RECT_EDGEshould be put into 
// lcartsn.h. However, there are some trouble
// to compile in that way.
//-------------------------------------------------
// states inside a cell

class C_STATE{
public:
	double C;		// solute concentration
	double D;		// diffusion rate

	void setZero(void);
};
// states on edge

//------------------------------------------------------
//		MESH
//------------------------------------------------------

class C_RECTANGLE {
public:
	int index;			// rectangle index
	int comp;			 
	//C_STATE state;
	double area;
	double coords[MAXD];	
	int icoords[MAXD];

	C_RECTANGLE();

	void setCoords(double*,int);
};


class C_CARTESIAN{
	Front *front;
public:
	C_CARTESIAN(Front &front);

	// member data: RECT_GRID
	int dim;

	// On topological grid
	RECT_GRID *top_grid;
	double *array;		// for scatter states;
	double *top_L,*top_U,*top_h,hmin;
	int *top_gmax;
	COMPONENT *top_comp;
	CRT_PARAMS *cRparams;
	CRT_FIELD *field;

	int *lbuf,*ubuf,*gmax;
	int *i_to_I,*I_to_i;		// Index mapping for 1D
	int **ij_to_I,**I_to_ij;	// Index mapping for 2D
	int ***ijk_to_I,**I_to_ijk;	// Index mapping for 3D

	// Sweeping limites
	int imin,jmin,kmin;
	int imax,jmax,kmax;

	// member data: mesh storage
	std::vector<C_RECTANGLE> 	cell_center;

	double m_t;                     // time
	double m_dt;			// time increment
        double max_dt;

	// constructor
	~C_CARTESIAN();

	// for parallel partition
	int             NLblocks,ilower,iupper;
        int             *n_dist;

	// mesh: full cells mesh
	void initMesh(void);		// setup the cartesian grid
	void setComponent(void);	// init components	
	void setAdvectionDt(void);	// compute time step 	
	void setBoundary(void);		// set up boundary conditions 	
	void printFrontInteriorStates(char*);
	void copySolute(void);

	void computeAdvection(void);

	void computeAdvectionExplicit(void);
	void computeAdvectionImplicit(void);
	void computeAdvectionCN(void);

	// interface functions
	void setDomain();
	void deleteGridIntfc();

	// Extra plot functions
	void oneDimPlot(char*);
	void xgraphOneDimPlot(char*);
#if defined __GD__
	void gdOneDimPlot(char*);
#endif /* defined __GD__ */

	// Extra movie functions
        void initMovieVariables(void);
        void augmentMovieVariables(void);
	void setInitialCondition(SEED_PARAMS);

	void checkStates();
	void sampleStates(SAMPLE);

	// parallelization related functions
	//
	void scatMeshArray();
	void setGlobalIndex();

	// physics calculation
	void readFrontInteriorStates(char*);

	void setIndexMap(void);
		// for compProjWithSmoothProperty(), 
		// should be changed to use setIndexMap() only.

	// -------------------------------------------------------
	// 		incompressible solver functions

	void computeSourceTerm(double *coords, double t, C_STATE &state); 

	void getInitialState(C_RECTANGLE&); 

	// main step function
	void solve(double dt);		

	void getVelocity(double *p, double *U);

	// ----------------------------------------------------------
	// 		utility functions
	// ----------------------------------------------------------

	void getRectangleIndex(int indexRectangle, int &i, int &j);
	void getRectangleIndex(int indexRectangle, int &i, int &j, int &k);
	int getRectangleComponent(int index);	// the center component
	void getRectangleCenter(int index, double *coords);
	void getRectangleCenter(int index0, int index1, double *coords);
	
	double getDistance(double *coords0, double *coords1);
	
			// incompletely implemented
	void getNearestInterfacePoint(double *q,double *p); 
		
	int  getComponent(int *icoords);	
	int  getComponent(double *coords);	
	void save(char *filename);
};

extern double   getStateSolute(POINTER);
extern void 	solute_print_front_states(FILE*,Front*);
extern void 	solute_read_front_states(FILE*,Front*);
extern void     read_crt_movie_options(char*,CRT_PARAMS*);

#endif
