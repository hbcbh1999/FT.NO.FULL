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
*				dprt.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains External Declarations for the Driver Code.
*/

#if !defined(_DPRT_H)
#define _DPRT_H


#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

enum _OUTPUT_TYPE {Float=1,Int=2,Pointer=3,Long=4,Unsigned=5,ULong=6};
typedef	enum _OUTPUT_TYPE	OUTPUT_TYPE;

typedef struct {
	union {
	    int          ival;
	    double        fval;
	    POINTER      pval;
	    long int     lval;
	    unsigned int uval;
	    uint64_t     ulval;
	} uval;
	OUTPUT_TYPE utype;
} OUTPUT_VALUE;


enum {
	CONSTANT_FIT = 01,
	LINEAR_FIT   = 02,
	SMOOTH	     = 04,
	SINGULAR     = 010
};

struct _INPUT_SOLN {
	const char	*name;
	int		fit,smoothness;

	double		*states1d;
	double		**states2d;
	double		***states3d;
	RECT_GRID	grid;

	INTERFACE	*intfc;			/* only needed for SINGULAR */
	void		(*set_intfc_states)(double*,double*,int,POINT*,
					    HYPER_SURF_ELEMENT*,HYPER_SURF*);
};
typedef struct _INPUT_SOLN INPUT_SOLN;


struct _OUTPUT_SOLN {
	const char	*name;
	int		fit,smoothness;

	OUTPUT_VALUE	*(*solution)(struct _OUTPUT_SOLN*,double*,int*);
	RECT_GRID	*grid;		/* only needed if intfc == NULL */

				     /* only needed for LINEAR_FIT | SINGULAR */
	void		(*intfc_solution)(struct _OUTPUT_SOLN*,POINT*,
					  HYPER_SURF_ELEMENT*,HYPER_SURF*,
					  OUTPUT_VALUE*,OUTPUT_VALUE*);
	INTERFACE	*intfc;		/* only needed for SINGULAR */

	int		var;		/* passed to solution functions */
	POINTER		extra;
	boolean            repart_at_end;  /* flag for re-partition output */
	int             *icrds_new;     /* position of new nodes */
};
typedef struct _OUTPUT_SOLN OUTPUT_SOLN;

		/* Types of Printing and Pausing Modes */

enum _PrtMode {
		PRINTING_OFF    =  0,
		EXACT_TIME	=  1,	/* Prints on Exact Time Intervals */
		MESH_TIME	=  2,	/* Prints on Mesh Time Intervals */
		CONSTANT_TIME	=  3,	/* Prints in Constant Time Intervals, */
					/* does not change dt */
		WALL_TIME       =  4,
		SET_PRT_MODE	=  5
};
typedef	enum _PrtMode	PrtMode;

#define	mesh_time_output(mode)		((mode) == MESH_TIME)
#define	exact_time_output(mode)		((mode) == EXACT_TIME)
#define	constant_time_output(mode)	((mode) == CONSTANT_TIME)
#define	real_time_output(mode)						\
		(exact_time_output(mode) || constant_time_output(mode))

struct	_PRINT_OPTIONS {
	/* Real time controls */
	double	_print_time_interval; /* Time interval for printing */
	double	_print_start_time;    /* Time to start printing */

	PrtMode	_prt_mode;	      /* Real vs. mesh time printing */

	/*Output type*/
	boolean		_print_in_binary;

	/* Mesh time controls */
	int	_print_step_interval; /* Interval for printing */
	int	_print_start_step;    /* Time step to start printing */

	/* Wall time controls */
	double	_print_wall_time_interval; /* Wall time interval for printing */
	double	_last_wall_time_dump;      /* Wall time of last wall time dump*/
	int	_wall_time_dump_number;	   /* Parity of wall time dump */

	/*File name specifications*/
	char		_print_directory[1024];
	char		_print_filename[1024];
};
typedef	struct _PRINT_OPTIONS PRINT_OPTIONS;

#define	Prt_mode(Pto)		    (Pto)._prt_mode
#define	prt_mode(pto)		    (pto)->_prt_mode
#define	Print_time_interval(Pto)    (Pto)._print_time_interval
#define	print_time_interval(pto)    (pto)->_print_time_interval
#define	Print_start_time(Pto)       (Pto)._print_start_time
#define	print_start_time(pto)       (pto)->_print_start_time
#define	Print_step_interval(Pto)    (Pto)._print_step_interval
#define	print_step_interval(pto)    (pto)->_print_step_interval
#define	Print_start_step(Pto)       (Pto)._print_start_step
#define	print_start_step(pto)       (pto)->_print_start_step
#define Print_directory(Pto)	    (Pto)._print_directory
#define print_directory(pto)	    (pto)->_print_directory
#define Print_filename(Pto)	    (Pto)._print_filename
#define print_filename(pto)	    (pto)->_print_filename
#define Print_in_binary(Pto)	    (Pto)._print_in_binary
#define print_in_binary(pto)	    (pto)->_print_in_binary
#define	Print_wall_time_interval(Pto)    (Pto)._print_wall_time_interval
#define	print_wall_time_interval(pto)    (pto)->_print_wall_time_interval
#define	Last_wall_time_dump(Pto)    (Pto)._last_wall_time_dump
#define	last_wall_time_dump(pto)    (pto)->_last_wall_time_dump
#define	Wall_time_dump_number(Pto)    (Pto)._wall_time_dump_number
#define	wall_time_dump_number(pto)    (pto)->_wall_time_dump_number

struct	_PRINT_CONTROL {
	PRINT_OPTIONS	_print_options;

	/* Real Time Controls */
	double	_next_print_time;     /* Time of next printout */
	double	_print_time_tolerance;/* TOLERANCE*/

	/* Mesh Time Controls */
	int	_next_print_step;     /* Time step of next printout */
};
typedef	struct _PRINT_CONTROL PRINT_CONTROL;

#define	Print_options(Ptc)		(Ptc)._print_options
#define	print_options(ptc)		(ptc)->_print_options
#define	Next_print_time(Ptc)		(Ptc)._next_print_time
#define	next_print_time(ptc)		(ptc)->_next_print_time
#define	Print_time_tolerance(Ptc)	(Ptc)._print_time_tolerance
#define	print_time_tolerance(ptc)	(ptc)->_print_time_tolerance
#define	Next_print_step(Ptc)		(Ptc)._next_print_step
#define	next_print_step(ptc)		(ptc)->_next_print_step

#define	PrtCtrl(odo)		(odo)->_PrtCtrl

struct _OUTPUT_DATA {
	struct _Printplot		*_prt;
	FILE				*_logfile;
	PRINT_CONTROL			_PrtCtrl;
};
typedef struct _OUTPUT_DATA OUTPUT_DATA;

	/* OUTPUT_DATA access macros */
#define	Output_printplot(data) 	     ((data)->_prt)
#define Output_file(data)	     ((data)->_logfile)
#define Output_PrtCtrl(data)	     ((data)->_PrtCtrl)
#define Print_control(data)	     (&Output_PrtCtrl(data))

#define	PrtOpts(data)	     	     Print_options(Output_PrtCtrl(data))

#define Output_dir(data)	     Print_directory(PrtOpts(data))
#define Output_filename(data)	     Print_filename(PrtOpts(data))
#define	Output_in_binary(data)	     Print_in_binary(PrtOpts(data))
#define Output_mode(data)	     Prt_mode(PrtOpts(data))
#define Output_time_freq(data)	     Print_time_interval(PrtOpts(data))
#define Output_start_time(data)	     Print_start_time(PrtOpts(data))
#define Output_step_freq(data)	     Print_step_interval(PrtOpts(data))
#define Output_start_step(data)	     Print_start_step(PrtOpts(data))

#define Output_next_print_time(data) Next_print_time(Output_PrtCtrl(data))
#define	Output_time_tolerance(data)  Print_time_tolerance(Output_PrtCtrl(data))
#define Output_next_print_step(data) Next_print_step(Output_PrtCtrl(data))

typedef struct {

	OUTPUT_DATA	data;

	/* The following are set in the physics */
	int		nfloats;		/* Num Conserved Quantities */
	double		(*stat_var)(int,int*,Wave*,Front*,double);
	double		(*stat_flux)(int,int,int*,Wave*);
	double		(*stat_inhom_source)(int,int*,Wave*,Front*,double);
	double		(*stat_point_flux)(int,int,Wave*);

	/* These may be set outside the physics */
	int		num_point_sources;
	double		*initial_present,*present,*incremented;
	double		*inhom_source_present,**source_incremented;

	char		*col_dname; /* For optional columnar output of data */
	char	        *col_filename;
	FILE		*col_file;
} Grid_Stats_data;
#define GS_data(data)	((Grid_Stats_data *)data)

struct _PROSTAR_PRINT_OPTIONS {
	double		_prostar_L0[3];	/* plotting window is defined    */
	double		_prostar_U0[3];	/* by the rectangle with corners */
	double		_prostar_len[3];	/* U0 - L0	 	 */
	int		_prostar_pixels[3];	/* pixels in direction   */
	int		_prostar_vertices[3];	/* pixels in direction   */
	int		_prostar_num_vars;
};
typedef struct _PROSTAR_PRINT_OPTIONS PROSTAR_PRINT_OPTIONS;

#define	PROSTAR_L0(Prostar_po)	             (Prostar_po)._prostar_L0
#define	prostar_L0(prostar_po)	             (prostar_po)->_prostar_L0
#define	PROSTAR_U0(Prostar_po)	             (Prostar_po)._prostar_U0
#define	prostar_U0(prostar_po)	             (prostar_po)->_prostar_U0
#define	PROSTAR_len(Prostar_po)	             (Prostar_po)._prostar_len
#define	prostar_len(prostar_po)	             (prostar_po)->_prostar_len
#define	PROSTAR_pixels(Prostar_po)           (Prostar_po)._prostar_pixels
#define	prostar_pixels(prostar_po)           (prostar_po)->_prostar_pixels
#define	PROSTAR_vertices(Prostar_po)           (Prostar_po)._prostar_vertices
#define	prostar_vertices(prostar_po)           (prostar_po)->_prostar_vertices
#define	PROSTAR_num_vars(Prostar_po)         (Prostar_po)._prostar_num_vars
#define	prostar_num_vars(prostar_po)	     (prostar_po)->_prostar_num_vars

struct _PROSTAR_FRAME_OPTS {
	boolean	_prostar_dflt_scale; /* Use default scaling?          */
	char	_prostar_selector[256];
	char	_prostar_plot_name[1024];
	double	_prostar_scale_min;  /* for non-default scaling       */
	double	_prostar_scale_max;	 /* for non-default scaling       */
	double	(*_prostar_plot_function)(double*,Front*,POINTER,
				      COMPONENT,Locstate);
};
typedef struct _PROSTAR_FRAME_OPTS PROSTAR_FRAME_OPTS;

#define	PROSTAR_dflt_scale(Fopts)		(Fopts)._prostar_dflt_scale
#define	prostar_dflt_scale(fopts)		(fopts)->_prostar_dflt_scale
#define	PROSTAR_selector(Fopts)		        (Fopts)._prostar_selector
#define	prostar_selector(fopts)		        (fopts)->_prostar_selector
#define	PROSTAR_plot_name(Fopts)		(Fopts)._prostar_plot_name
#define	prostar_plot_name(fopts)		(fopts)->_prostar_plot_name
#define	PROSTAR_scale_min(Fopts)		(Fopts)._prostar_scale_min
#define	prostar_scale_min(fopts)		(fopts)->_prostar_scale_min
#define	PROSTAR_scale_max(Fopts)		(Fopts)._prostar_scale_max
#define	prostar_scale_max(fopts)		(fopts)->_prostar_scale_max
#define	PROSTAR_plot_function(Fopts)	        (Fopts)._prostar_plot_function
#define	prostar_plot_function(fopts)	        (fopts)->_prostar_plot_function

struct _PROSTAR_plot_data {

	PROSTAR_PRINT_OPTIONS	_PROSTAR_print_opts;

	int		dim;
	double		step[3];	/* step = (U0 - L0)/pixels	 */
	double		*scale[4];
        double		*vert_scale[4]; /*vertices scaling */
        char            *v_file_name;
        char            *c_file_name;
        char            *p_file_name;
	char            *base_file_name;
        boolean            vrt_cel_created; /* flag to see if .vrt and .cel 
					    files have been created */
	struct _PROSTAR_frame_data {

		PROSTAR_FRAME_OPTS	_PROSTAR_frame_opts;

		boolean            first;
		boolean	        append;
	}	*frame_data;

};
typedef struct _PROSTAR_plot_data PROSTAR_plot_data;
#if defined(__cplusplus)
typedef struct PROSTAR_plot_data::_PROSTAR_frame_data PROSTAR_frame_data;
#else /* defined(__cplusplus) */
typedef struct _PROSTAR_frame_data PROSTAR_frame_data;
#endif /* defined(__cplusplus) */

#define	PROSTAR_print_opts(prostar_pdata)  (prostar_pdata)->_PROSTAR_print_opts

#define	PROSTAR_frame_opts(Fdata)		(Fdata)._PROSTAR_frame_opts
#define	prostar_frame_opts(fdata)		(fdata)->_PROSTAR_frame_opts

#define PROSTAR_frame_dflt_scale(Fdata)      	PROSTAR_dflt_scale(PROSTAR_frame_opts(Fdata))
#define prostar_frame_dflt_scale(fdata) 	PROSTAR_dflt_scale(prostar_frame_opts(fdata))
#define PROSTAR_frame_scale_min(Fdata)	        PROSTAR_scale_min(PROSTAR_frame_opts(Fdata))
#define prostar_frame_scale_min(fdata)	        PROSTAR_scale_min(prostar_frame_opts(fdata))
#define PROSTAR_frame_scale_max(Fdata)	        PROSTAR_scale_max(PROSTAR_frame_opts(Fdata))
#define prostar_frame_scale_max(fdata)         	PROSTAR_scale_max(prostar_frame_opts(fdata))
#define PROSTAR_frame_plot_name(Fdata)      	PROSTAR_plot_name(PROSTAR_frame_opts(Fdata))
#define prostar_frame_plot_name(fdata)	        PROSTAR_plot_name(prostar_frame_opts(fdata))
#define PROSTAR_frame_plot_function(Fdata)	PROSTAR_plot_function(PROSTAR_frame_opts(Fdata))
#define prostar_frame_plot_function(fdata)	PROSTAR_plot_function(prostar_frame_opts(fdata))


struct _ScalarPlotItem {
	const char *name;
	double	   (*plot_fn)(double*,Front*,POINTER,COMPONENT,Locstate);
	struct _ScalarPlotItem	*next;
};
typedef	struct _ScalarPlotItem ScalarPlotItem;

struct _Plot_choice {
	const char           *prompt;
	const char	     *selector;
	int	             (*_SelectedScalarItems)(const char*,
			          		   struct  _Plot_choice*,
				          	   ScalarPlotItem**);
	struct	_Plot_choice *next;
};
typedef	struct _Plot_choice Plot_choice;

#define	SelectedScalarItems(c,pc,spi)	(*(pc)->_SelectedScalarItems)(c,pc,spi)

		/* rect state printing variable list */

struct _PRINTING_LIST {
	const char            *name;
	int	              var;
	struct _PRINTING_LIST *next;
};
typedef struct _PRINTING_LIST PRINTING_LIST;


	/* Currently supported output formats */
enum{
    PRT_FRONTS,
    RECT_STATES,
    TRI_STATES,
    HDF_STATES,
    SDS_STATES,
    /* needed for VTK */
    VTK_STATES,
    /* needed for VTK*/
    PROSTAR_STATES,
    GD_MOVIE,
    NUM_OUTPUT_FORMATS
};

    typedef double (VTK_PLOT_FILTER)(double);


typedef struct _VTK_FRAME_OPTS {
	char	_vtk_selector[256];
	char	_vtk_plot_name[1024];
	double	(*_vtk_plot_function)(double*,Front*,void*,
				      COMPONENT,Locstate);
	VTK_PLOT_FILTER *_vtk_plot_filter;
} VTK_FRAME_OPTS;

typedef struct _VTK_PRINT_OPTIONS {
        double           _vtk_L0[3];     /* plotting window is defined    */
        double           _vtk_U0[3];     /* by the rectangle with corners */
        double           _vtk_len[3];    /* U0 - L0                       */
        double           _vtk_V[3];      /* L0 + V*t, U0 + V*t            */
        int             _vtk_pixels[3]; /* pixels in direction           */
        int             _vtk_num_vars;
	char		vars_string[1024];
	char		base_string[1024];
	char            base_name[1024];
        int             numvects;
	int 		subsample_factor; 
        int             vecselect[3][20]; 
	boolean		_print_in_binary; /* print vtk files in binary if set to true */
} VTK_PRINT_OPTIONS;

typedef struct _VTK_frame_data {
        VTK_FRAME_OPTS  _VTK_frame_opts;
        } VTK_frame_data;

typedef struct _VTK_plot_data {

        VTK_PRINT_OPTIONS       _VTK_print_opts;
        int             dim;
        COMPONENT       *comp;
        char		sel_plots[1024];
	int             num_values;
        int             num_raster_data;
        VTK_frame_data       *frame_data;

} VTK_plot_data;

#if defined(__GD__)
typedef struct _GD_FRAME_OPTS {
	char	_gd_selector[256];
	char	_gd_plot_name[1024];
	char    _gd_var_name[1024];
	double	(*_gd_plot_function)(double*,Front*,void*,
				      COMPONENT,Locstate);
} GD_FRAME_OPTS;

typedef struct _GD_PRINT_OPTIONS {
        int             _gd_num_vars;
	char		vars_string[1024];
	char		base_string[1024];
	char            base_name[1024];
} GD_PRINT_OPTIONS;

typedef struct _GD_plot_data {
        GD_PRINT_OPTIONS       _GD_print_opts;
        int             dim;
        COMPONENT       *comp;
        char		sel_plots[1024];
	int             num_values;
        int             num_raster_data;
        GD_FRAME_OPTS   *frame_opts;
} GD_plot_data;
#endif /* defined(__GD__) */


struct _Printplot {
	
	/* Output control structure for main output formats */
	OUTPUT_DATA *main_output_data[NUM_OUTPUT_FORMATS];
	OUTPUT_DATA *store_main_output_data[NUM_OUTPUT_FORMATS];

	/* Output control structure for additional printing/plotting */
	OUTPUT_DATA **user_output_data;
	void	    (**user_output_funcs)(Grid*,Wave*,Front*,
					  struct _Printplot*,
					  OUTPUT_DATA*,boolean);

	/* Computes statistics */
	void	    (*grid_statistics)(Grid_Stats_data*,Grid*,
				       Wave*,Front*,int);
	Grid_Stats_data *gs_data;

	/* Prompts for initialization or changes in the Prinplot data */
	void (*_init_printplot)(INIT_DATA*,Grid*,Front*,struct _Printplot*);
	void (*_init_statistics)(Front*,Grid*,Wave*,
				 struct _Printplot*,INIT_DATA*);

	/* Prints the initial data */
	const char *title;
	void (*print_initial_data)(FILE*,CHART*,struct _Printplot*);

	/* Initializes Before each Printout */
	void (*initialize_for_printout)(Grid*,Wave*,Front*,struct _Printplot*);

	/* Prints the data */
	void (*printout)(CHART*,struct _Printplot*,boolean,int);

	/* Prints state data */
	void (*_print_states[3])(FILE*,struct _Printplot*);

	/* Printout ellip solution variables */
	void (*plot_ellip_vars)(Grid*);

	FILE		*file;		/* Optional output file */
	char		*outfile;	/* Optional name of output file */
	int		compress;	/* If YES compress outfile */

	/* Info for restart input, see dinout.c */
	int		n_restart_vars;
	INPUT_SOLN	**restart_soln;

	/* Info for RECT_STATES printout, see dinout.c */
	int		n_rect_state_vars;
	OUTPUT_SOLN 	**output_soln;

	/* Extreme solution values */
	void (*_print_extreme_values)(FILE*,CHART*,struct _Printplot*);

	/* Info for tri plots */
	size_t		n_tri_vars;
	const char	**tri_plot_name;
	double		(**tri_plot_function)(double*,Front*,POINTER,
					     COMPONENT,Locstate);
        /* Info for totals
	double           *total_Upper
	double           *total_Lower
	int             *total_pix; */

#if defined(USE_HDF)
	/* Info for HDF raster plots */
	int		n_HDF_vars;
	HDF_plot_data	*HDF_data;
	/* Info for SDS plots */
	int		n_SDS_vars;
	HDF_plot_data	*SDS_data;
	/* needed for VTK */
#endif /* defined(USE_HDF) */

	/* Info for VTK plots */
	int             n_VTK_vars;
	VTK_plot_data   *vtk_data;
	/* needed for VTK */

#if defined(__GD__)
	/* Info for GD plots */
	int             n_GD_vars;
	GD_plot_data   *GD_data;
	/* needed for GD */
#endif /* defined(__GD__) */

        /* Info for PROSTAR plots */
	int		n_PROSTAR_vars;
        PROSTAR_plot_data *PROSTAR_data;

	PRINT_OPTIONS	*_wall_time_dump_options;

	struct _INIT_DATA	*init;
        boolean   surface_area;
};
typedef struct _Printplot Printplot;

#define	init_printplot(init,grid,fr,prt)				  \
	(*(prt)->_init_printplot)(init,grid,fr,prt)
#define	init_statistics(fr,grid,wave,prt,init)				   \
	(*(prt)->_init_statistics)(fr,grid,wave,prt,init)
#define	print_states(file,prt,dim)	(*(prt)->_print_states[dim-1])(file,prt)
#define print_extreme_values(file,chart,prt)				\
	if (prt->_print_extreme_values != NULL) 			\
	    (*prt->_print_extreme_values)(file,chart,prt)

#define	wall_time_dump_options(prt)	(prt)->_wall_time_dump_options

#define output_format_on(format,prt)					   \
	((prt)->main_output_data[format] != NULL)

#define is_ts_for_output_format(format,grid,prt)			   \
	(output_format_on(format,prt)					   \
	 		 	&&					   \
	 (is_print_time((grid)->time,					   \
			Print_control((prt)->main_output_data[format])) || \
	  is_print_step((grid)->step,					   \
			Print_control((prt)->main_output_data[format]))))


	/* Structures for printing cross sectional data */
struct _CROSS_SECTION {
	double                 **pcoords;
	int                   num_pts;
	const char            *message;
	struct _CROSS_SECTION *prev;
	struct _CROSS_SECTION *next;
};

typedef struct _CROSS_SECTION CROSS_SECTION;

typedef struct {
	OUTPUT_DATA	odata;
	CROSS_SECTION	*cross_sections;/* Cross sections for linear plots */
} Cross_Sections_data;

#define Cr_Data(data)	((Cross_Sections_data *)data)

#endif /* !defined(_DPRT_H) */



#if defined(__cplusplus)
/*typedef struct VTK_plot_data::_VTK_frame_data VTK_frame_data; */
#else /* defined(__cplusplus) */
/*typedef struct _VTK_frame_data VTK_frame_data; */
#endif /* defined(__cplusplus) */

/********************** end blah *****************************/
/*#define VTK_frame_options(init)         init->_VTK_frame_options */

#define	VTK_plot_function(Fopts)	(Fopts)._vtk_plot_function
#define	vtk_plot_function(fopts)	(fopts)->_vtk_plot_function
#define	VTK_plot_filter(Fopts)		(Fopts)._vtk_plot_filter
#define	vtk_plot_filter(fopts)		(fopts)->_vtk_plot_filter
#define	VTK_plot_name(Fopts)		(Fopts)._vtk_plot_name
#define	vtk_plot_name(fopts)		(fopts)->_vtk_plot_name
#define VTK_frame_plot_name(Fdata)	VTK_plot_name(VTK_frame_opts(Fdata))
#define vtk_frame_plot_name(fdata)	VTK_plot_name(vtk_frame_opts(fdata))
#define	VTK_print_opts(hdf_pdata)	(hdf_pdata)->_VTK_print_opts
#define VTK_numvects(Hdf_po)        (Hdf_po).numvects
#define vtk_numvects(hdf_po)        (hdf_po)->numvects
#define VTK_vecselect(Hdf_po)        (Hdf_po).vecselect
#define vtk_vecselect(hdf_po)        (hdf_po)->vecselect
#define	VTK_selector(Fopts)		(Fopts)._vtk_selector
#define	vtk_selector(fopts)		(fopts)->_vtk_selector
#define	VTK_frame_opts(Fdata)		(Fdata)._VTK_frame_opts
#define	vtk_frame_opts(fdata)		(fdata)->_VTK_frame_opts
#define	VTK_num_vars(Vtk_po)         (Vtk_po)._vtk_num_vars
#define	vtk_num_vars(vtk_po)	     (vtk_po)->_vtk_num_vars


#define VTK_vars(Hdf_po)        (Hdf_po).vars_string
#define vtk_vars(hdf_po)        (hdf_po)->vars_string
#define VTK_base(Hdf_po)        (Hdf_po).base_string
#define vtk_base(hdf_po)        (hdf_po)->base_string
/******************************************************/
#define	VTK_L0(Hdf_po)	             (Hdf_po)._vtk_L0
#define	vtk_L0(hdf_po)	             (hdf_po)->_vtk_L0
#define	VTK_U0(Hdf_po)	             (Hdf_po)._vtk_U0
#define	vtk_U0(hdf_po)	             (hdf_po)->_vtk_U0
#define	VTK_len(Hdf_po)	             (Hdf_po)._vtk_len
#define	vtk_len(hdf_po)	             (hdf_po)->_vtk_len
#define	VTK_V(Hdf_po)	             (Hdf_po)._vtk_V
#define	vtk_V(hdf_po)	             (hdf_po)->_vtk_V
#define	VTK_pixels(Hdf_po)           (Hdf_po)._vtk_pixels
#define	vtk_pixels(hdf_po)           (hdf_po)->_vtk_pixels
#define VTK_frame_plot_function(Fdata)	VTK_plot_function(VTK_frame_opts(Fdata))
#define VTK_frame_sel(Fdata)	VTK_selector(VTK_frame_opts(Fdata))
#define vtk_frame_plot_function(fdata)	VTK_plot_function(vtk_frame_opts(fdata))
/*******************************************************/

#if defined(__GD__)
#define	GD_selector(Fopts)	     (Fopts)._gd_selector
#define	GD_print_opts(gd_pdata)      (gd_pdata)->_GD_print_opts
#define	gd_plot_function(fopts)	     (fopts)->_gd_plot_function
#define	gd_plot_name(fopts)	     (fopts)->_gd_plot_name
#define	gd_var_name(fopts)	     (fopts)->_gd_var_name
#define	gd_num_vars(gd_po)	     (gd_po)->_gd_num_vars
#define GD_frame_plot_function(Fdata) (Fdata)._gd_plot_function
#endif /* defined(__GD__) */

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif
