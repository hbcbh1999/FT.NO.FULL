/**********************************************************************
 * 		ifluid_basic.h
 **********************************************************************/

#ifndef _FT_IFLUID_BASIC_H_
#define _FT_IFLUID_BASIC_H_

enum _IF_PROB_TYPE {
        ERROR_TYPE = -1,
        TWO_FLUID_BUBBLE = 1,
        TWO_FLUID_RT,
        TWO_FLUID_KH,
        TWO_FLUID_RS_RV,
        TWO_FLUID_RS_SY,
        TWO_FLUID_TC,
        FLUID_SOLID_CIRCLE,
        BUBBLE_SURFACE,
        FLUID_RIGID_BODY,
        ROTOR_ONE_FLUID,
        ROTOR_TWO_FLUID,
	CHANNEL_FLOW,
        FLUID_CRYSTAL
};
typedef enum _IF_PROB_TYPE IF_PROB_TYPE;

struct _STATE {
	double dens;			/* Density */
        double pres;                    /* Pressure */
        double vel[MAXD];               /* Velocities */
        double vort;                    /* Vorticity in 2D */
        double vort3d[MAXD];            /* Vorticity in 3D */
        // for vd
        double conc;                    /* Concentration */
        double dens_old;
};
typedef struct _STATE STATE;

extern void restart_set_dirichlet_bdry_function(Front*);
extern void iF_flowThroughBoundaryState(double*,HYPER_SURF*,Front*,POINTER,
                        POINTER);
extern void iF_timeDependBoundaryState(double*,HYPER_SURF*,Front*,POINTER,
                        POINTER);
extern void ifluid_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
// for vd
extern void ifluid_point_propagate_vd(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);

extern void ifluid_compute_force_and_torque(Front*,CURVE*,double,double*,
                        double*);
extern void setInitialIntfc(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);
extern void init_fluid_state_func(Incompress_Solver_Smooth_Basis*,IF_PROB_TYPE);
extern void read_iFparams(char*,IF_PARAMS*);
extern void read_iF_prob_type(char*,IF_PROB_TYPE*);
extern void recordBdryEnergyFlux(Front*,char*);
#endif
