#include <FronTier.h>

enum _AF_PROB_TYPE {
        ERROR_TYPE = -1,
        DOUBLE_VORTEX_VELO,
        FLIP_FLAP,
        LEAF_FALL,
        PARACHUTE,
        REVERSAL_VELO,
	ZERO_VELO
};
typedef enum _AF_PROB_TYPE AF_PROB_TYPE;

enum _INIT_SURF_TYPE {
        ERROR_SURF_TYPE = -1,
        ELLIPTIC,
        SINE_WAVE,
        CIRCULAR_PERT,
	READ_FROM_VTK
};
typedef enum _INIT_SURF_TYPE INIT_SURF_TYPE;

typedef struct {
        int dim;
        INIT_SURF_TYPE init_surf_type;
        POINTER level_func_params;
        IF_MOVIE_OPTION *movie_option;
        NS_SCHEME num_scheme;
        double rho1;
        double rho2;
        double mu1;
        double mu2;
        double U1[MAXD];
        double U2[MAXD];
        double gravity[MAXD];
        double surf_tension;
        double smoothing_radius;
        IF_FIELD *field;
} AF_PARAMS;

/*	rgbody.c functions */

struct _STATE {
	double dens;			/* Density */
        double pres;                    /* Pressure */
        double vel[MAXD];               /* Velocities */
        double vort;                    /* Vorticity in 2D */
        double vort3d[MAXD];            /* Vorticity in 3D */
};
typedef struct _STATE STATE;

typedef struct {
        double i1,i2;
        double cen1[2],cen2[2];
} DOUBLE_VORTEX_PARAMS;

void read_iFparams(char*,IF_PARAMS*);
void read_af_prob_type(char*,AF_PROB_TYPE*);
void read_movie_options(char*,IF_PARAMS*);
void read_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
void restart_set_dirichlet_bdry_function(Front*);
void liquid_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void init_fluid_state_func(Incompress_Solver_Smooth_Basis*,AF_PROB_TYPE);
extern void setInitialIntfc(Front*,LEVEL_FUNC_PACK*,char*,AF_PROB_TYPE);
extern void ifluid_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void initVelocityFunc(Front*,AF_PROB_TYPE);
