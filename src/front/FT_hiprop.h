/*! 
 * \file FT_hiprop.h
 *
 * \author Chenzhe Diao
 * \date 2012.11.10
 *
 *
 */

//#include <stdio.h>
//#include "FronTier.h"
//
#include <front/fdecs.h>
//#include "hiprop.h"
#include "hiprop_grid.h"
#include "util.h"


#ifndef __FT_HIPROP_H__
#define __FT_HIPROP_H__

#ifndef EXTERN_C

#ifdef __cplusplus
  #define EXTERN_C extern "C"
#else
  #define EXTERN_C extern
#endif

#endif

//EXTERN_C void print_pinfo(hiPropMesh *mesh);

EXTERN_C int tri_hit_box(TRI* tri, double* L, double* U, double* h);

EXTERN_C void ImportMeshToHiprop(INTERFACE* intfc, hiPropMesh* mesh, PP_GRID* ppgrid, int debug_step);


EXTERN_C boolean read_hiprop_surface(INTERFACE   *intfc, COMPONENT   neg_comp, COMPONENT   pos_comp, 
		hiPropMesh  *mesh, SURFACE **ps);

EXTERN_C void ExportMeshToFronTier(
		hiPropMesh* mesh, 
		Front *front,
		COMPONENT neg_comp,
		COMPONENT pos_comp);

EXTERN_C void print_tri_normal(INTERFACE* intfc);
EXTERN_C void print_intfc_tmp_vtk(
	Front *front,
        char *out_name,
        boolean print_in_binary);
EXTERN_C void print_hpmesh(
	hiPropMesh* mesh,
	char* out_name,	// output dir name
	const char* special_name,
	int step);
#endif
