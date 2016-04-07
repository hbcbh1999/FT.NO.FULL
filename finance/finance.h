/**********************************************************************
 * 				finance.h			      *
 **********************************************************************/

#include <FronTier.h>
#include <vector>
#include <petscksp.h>
#include <assert.h>

#define         EXCERCISE_COMP            	2
#define         BLACK_SCHOLES_COMP              3

struct _FIELD {
        double *option_price;                   /* Pressure */
};
typedef struct _FIELD FIELD;

enum _F_TYPE {
	ERROR_TYPE = -1,
	EURO_CALL_OPTION = 1,
	EURO_PUT_OPTION = 2,
	AMRI_CALL_OPTION = 3,
	AMRI_PUT_OPTION = 4
};
typedef enum _F_TYPE F_TYPE;

enum _NUM_SCHEME {
	UNSPLIT_EXPLICIT = 1,
	UNSPLIT_IMPLICIT,
	CRANK_NICOLSON
};
typedef enum _NUM_SCHEME NUM_SCHEME;

struct _PARAMS {
	F_TYPE f_type;
	NUM_SCHEME num_scheme;
	double *rate_of_change;	// Option price on stock value
	double *option_price;	// Option price on stock value
	double *temp_price;     // intermediate price while looking for SF
	double E;		// Excercise price
	double sigma[MAXD];	// Volotility
	double r;		// Interest Rate
	double D;               // Dividend
	double oldroot;
	double a;
	double b;
	double c;
        int idx;                // index for current front
	bool findSF;
	FIELD *field;
};
typedef struct _PARAMS PARAMS;

struct _TIME_DEPENDENT_PARAMS {
};
typedef struct _TIME_DEPENDENT_PARAMS TIME_DEPENDENT_PARAMS;

typedef double STATE;

extern double   price_of_state(POINTER);
extern double	extend_from_put_exc(double,double*,double*,double*,double);
extern double	extend_from_call_exc(double,double*,double*,double*,double);

typedef class CARTESIAN CARTESIAN_EB;

/******************************************************************************
 * 		lcartsn.h
 * A simple incompressible flow solver using the ghost fluid method and the
 * projection method.
 *
 * the main function is 
 * 	CARTESIAN::solve().
 *
 * References:
 ******************************************************************************/

class SOLVER;
class CARTESIAN;

//typedef class CARTESIAN CARTESIAN_EB;
//enum VISITED_TYPE {UNVISITED, VISITED, PARTIAL_VISITED, FULL_VISITED};

//-------------------------------------------------
//		STATES
// NOTE:
//      INC_STATE/INC_STATE_RECT_EDGEshould be put into 
// lcartsn.h. However, there are some trouble
// to compile in that way.
//-------------------------------------------------
// states inside a cell

class INC_STATE{
public:
	double P;		// Option price

	double sigma[MAXD];	// Volotility 
	double r;		// Interest rate

	void setZero(void);
};
// states on edge

//------------------------------------------------------
//		MESH
//------------------------------------------------------

class RECTANGLE {
public:
	int index;			// rectangle index
	int comp;			 
	INC_STATE state;
	double area;
	double coords[MAXD];	
	int icoords[MAXD];

	RECTANGLE();

	void setCoords(double*,int);
};


class CARTESIAN{
	Front *front;
public:
	CARTESIAN(Front &front);

	// member data: RECT_GRID
	int dim;

	// On topological grid
	RECT_GRID *top_grid;
	double *array;		// for scatter states;
	double *top_L,*top_U,*top_h;
	COMPONENT *top_comp;
	PARAMS *eqn_params;
	FIELD *field;

	int *top_gmax;
	int *lbuf,*ubuf;
	int *i_to_I,*I_to_i;		// Index mapping for 1D
	int **ij_to_I,**I_to_ij;	// Index mapping for 2D
	int ***ijk_to_I,**I_to_ijk;	// Index mapping for 3D

	// Sweeping limites
	int imin,jmin,kmin;
	int imax,jmax,kmax;

	enum BC_TYPE { 									// used by ADVECTION
		BC_PERIODIC = 1,
		BC_Extrapolation0 = 2,
		BC_Extrapolation1 = 3,
		BC_InflowOutflow = 4};	
	BC_TYPE m_bc[4];								// down, right, up, left 		

	// member data: mesh storage
	std::vector<RECTANGLE> 	cell_center;

	double m_t;                     // time
	double m_dt;			// time increment

	// constructor
	~CARTESIAN();

	// for parallel partition
	int             NLblocks,ilower,iupper;
        int             *n_dist;

	// mesh: full cells mesh
	void initMesh(void);		// setup the cartesian grid
	void setComponent(void);	// init components	
	void setAdvectionDt(void);	// compute time step 	
	void setBoundary(void);		// set up boundary conditions 	
	void readOptionPrice(char*);
	void printOptionPrice(char*);
	void copyOptionPrice(void);
	void copyBackOptionPrice(void);

	void computeAdvection(void);

	void computeAdvectionExplicit(void);
	void computeAdvectionExplicit1d(void);
	void computeAdvectionExplicit2d(void);
	void computeAdvectionExplicit3d(void);

	void computeAdvectionImplicit(void);
	void computeAdvectionImplicit1d(void);
	void computeAdvectionImplicit2d(void);
	void computeAdvectionImplicit3d(void);

	void computeAdvectionCN(void);
	void computeAdvectionCN1d(void);
	void computeAdvectionCN2d(void);
	void computeAdvectionCN3d(void);

	// Extra plot functions
	void oneDimPlot(char*);
	void xgraphOneDimPlot(char*);
#if defined __GD__
	void gdOneDimPlot(char*);
#endif /* defined __GD__ */

	void checkStates();

	// parallelization related functions
	//
	void scatMeshArray();
	void setGlobalIndex();

	// physics calculation
	void setInitialCondition(void);

	void setDomain();
	void setIndexMap(void);
		// for compProjWithSmoothProperty(), 
		// should be changed to use setIndexMap() only.

	// -------------------------------------------------------
	// 		incompressible solver functions

	void computeSourceTerm(double *coords, double t, INC_STATE &state); 

	void getInitialState(double *coords, INC_STATE &state); 
	void getInitialState(double *coords, STATE *state); 

	// main step function
	void solve(double dt);		

	void getVelocity(double *p, double *U);

	// ----------------------------------------------------------
	// 		utility functions
	// ----------------------------------------------------------

	void getRectangleCenter(int index, double *coords);
	
	int  getComponent(int *icoords);	
	void save(char *filename);
};

