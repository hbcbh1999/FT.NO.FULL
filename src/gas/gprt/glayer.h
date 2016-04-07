/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 
whose solutions have discontinuities.  


Copyright (C) 1999 by The University at Stony Brook. 
 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

******************************************************************************/


/*
*				glayer.h
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Declarations and definitions for layer statistics and extrema
*/

#if !defined(_GLAYER_H)
#define _GLAYER_H

#include <gdecs/gdecs.h>

#define	DEBUG_STRING	"glayer"

enum {
	MIN_POS_ID   = FIRST_PHYSICS_MESSAGE_ID+1,
	MIN_COMP_ID,
	MIN_STATE_ID,
	MAX_POS_ID,
	MAX_COMP_ID,
	MAX_STATE_ID,

	G_IMIN_ID     = FIRST_PHYSICS_MESSAGE_ID+10,
	G_IMAX_ID,
	G_HMIN_ID,
	G_HMAX_ID,
	LAYER_SUM_ID,
	INTFC_SUM_ID,
	LAYER_FRAC_ID,

/* This ID should be last, because it uses subsequent ID's 
   in the case of parallel processors: */
	L_MIN_MAX_ID = LAYER_FRAC_ID
};

/* mapping from (i,j) entry of a symmetric matrix (j >= i) 
   to an array of dimension MAXD*(MAXD-1)/2 */
#define sym_index(i,j) (i*MAXD-i*(i+1)/2+j)

enum intfc_geometry { planar, spherical };

typedef struct {
    double 		d;		/* density */
    double 		v[MAXD];	/* velocity */
    double 		p;		/* pressure */
    double               k;              /* kinetic energy */
    double 		e;		/* internal energy */
    double		frac;		/* layer fraction */
    int			count;		/* number of pts summed within layer */

  /* Following, including comment lines, added 24/10/2003 - egeorge.*/
  /* Foll. added 21 Aug. 2003 by egeorge to compute an effective and local 
     Atwood no. Currently computed in gas/gprt/gintext.c's 
     accumulate_state_in_layer() only for planar 3d case.
  */
    double               min_d;
    double               max_d;
    double               global_min_d;
    double               global_max_d;
  /* Following added 17 Nov. 2003 to help compute atwood5. It is the density 
     averaged over the largest 50% of density values in the heavy phase and 
     the smallest 50% of density values in the light phase. - egeorge
  */
    double               partial_d;
  /* Following added 6 Apr. 2004 to help compute atwood6. It is the median
     density in a given fluid phase (heavy or light).
  */
    double               median_d;
} Big_State;

typedef struct {
    double		dd;             /* density*density */
    double		de;             /* internal energy density */
    double 		dv[MAXD];       /* momentum */
    double		dkv[MAXD];      /* kinetic energy flux */
    double		dev[MAXD];      /* internal energy flux */
    double 		pv[MAXD];       /* pressure work */
    double 		vv[MAXD*(MAXD+1)/2];   /* velocity covariances  */
                                               /*   XX, XY, XZ, YY, YZ, ZZ */
    double 		dvv[MAXD*(MAXD+1)/2];  /* momentum flux tensor  */
                                               /*   XX, XY, XZ, YY, YZ, ZZ */
} Turbulence_moments;


typedef struct {
    OUTPUT_DATA         odata;
    FILE		*out_min, *out_max, *out_amp;
    FILE		*out_min1, *out_max1, *out_amp1;
    FILE		*out_min5, *out_max5, *out_amp5;
    boolean		do_01, do_05;
    char                *min_fname, *max_fname, *amp_fname;
    char		*min_dname, *max_dname, *amp_dname;
    enum intfc_geometry geom;
    double 		origin[3];
    double		rfactor;
    double		pos[2][3];       /* 0 <--> min; 1 <--> max; */
    int			index_at_min, index_at_max;
} Intfc_extrema;

#define Intfc_extrema_data(odata)	(Intfc_extrema*)odata

typedef struct {
	OUTPUT_DATA    	odata;
	FILE		**torque_files,**omega_files;
	FILE		**force_files,**com_files;	

	char		**torque_fname,**omega_fname;
	char		**force_fname,**com_fname;

	char		**torque_dname,**omega_dname;
	char		**force_dname,**com_dname;
} RGB_DATA;

#define	RGB_data(odata)	(RGB_DATA*)odata

typedef struct {	 

    /* quantities that are summed at interfaces within each layer: */

    double nor[MAXD];	 /* gradient of characteristic function, */
    double vdotn;         /* interface velocity dotted with normal, */
    double pre[MAXD];     /* pressure times normal */
    double pvdotn;        /* pressure times velocity dotted with normal */

    double nor_fac;       /* For 2D, number of points of intersection of 
			  * contact surface with x axis at a particular z.
			  * For 3D, total length of lines of intersection of
			  * contact surface with (x,y) plane at particular z */

    int count;           /* Number of occurrences of intersection of contact
			  * surface with x axis or (x,y) plane at a particular
			  * z.  Redundant for 2D, but needed as a flag for 3D
			  * to avoid division by 0 during normalization. */
} Intfc_data;


typedef struct {
    OUTPUT_DATA	odata;
    char	column_header[Gets_BUF_SIZE];

    enum intfc_geometry geom;           /* layer geometry */
    double 		origin[MAXD];   /* needed for spherical geometry */
    
    double 	h_min;         /* range of height or radius over which to */
    double 	h_max;         /* compute intfc stats */
    			       
    int		n_layers;      /* number of layers ranging from
                                * h_min to h_max inclusive */

    double 	dh;            /* (h_max - h_min)/(n_layers - 1) */
    			       
    Intfc_data *data;    /* array containing the data for each layer */
	                                   
} Intfc_stats;


typedef struct {
    OUTPUT_DATA         odata;
    FILE		**bst_out;  	/* array of output files of */
    char		**bst_fname;    /* layer-averaged Big_States, one */
    char		**bst_dname;    /* for each distinct material, with */
                                        /* associated file and dir names */

    boolean		include_turb;
    FILE		**turb_out;     /* array of output files of */
    char                **turb_fname;   /* layer-averaged Turbulent_moments, */
    char		**turb_dname;   /* one for each material, with */
                                        /* associated file and dir names */

    enum intfc_geometry geom;           /* layer geometry */
    int			n_params;	/* number of distinct materials */
    double 		origin[MAXD];   /* needed for spherical geometry */
    double		rfactor;        /* sub-grid refinement factor */

    char	bst_col_header[Gets_BUF_SIZE];
    char	turb_col_header[Gets_BUF_SIZE];
    
    double 	h_min;         /* range of height or radius over which to */
    double 	h_max;         /* compute layer stats */
    			       
    int		n_layers;      /* number of layers ranging from
                                * h_min to h_max inclusive */

    double 	dh;            /* (h_max - h_min)/(n_layers - 1) */
    			       
    Big_State   *bst_data;    /* arrays containing the data for each layer */
    Turbulence_moments* turb_data;
	                                   
} Layer_stats;


#if defined(TWOD) || defined(THREED)

IMPORT	void accumulate_state_in_layer(Wave*,Front*,double,Big_State*,
				        Turbulence_moments*,double,const double*,
					boolean,boolean,boolean,boolean);
IMPORT	void accumulate_state_totals(const Big_State*,const Turbulence_moments*,
				     Big_State*,Turbulence_moments*);
IMPORT	void normalize_state_totals(Big_State*,Turbulence_moments*,int);

#endif  /* defined(TWOD) || defined(THREED) */

#endif /* !defined(_GLAYER_H) */
