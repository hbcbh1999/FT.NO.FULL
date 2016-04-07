/*!
 * \file FT_hiprop.c
 * \brief Implementation of functions as interface between FronTier and hiprop library
 *
 * \author Chenzhe Diao
 * \date 2012.10.25
 *
 */

#include <front/FT_hiprop.h>
/*
 * Check if the bounding box of the triangle has intersection with
 * the box bounded by L-h and U+h.
 * h is the relaxation here for safty.
 */
int tri_hit_box(TRI* tri, double* L, double* U, double* h)
{
    int i;
    double eps = 1e-8;

    double x1 = Point_of_tri(tri)[0]->_coords[0];
    double y1 = Point_of_tri(tri)[0]->_coords[1];
    double z1 = Point_of_tri(tri)[0]->_coords[2];

    double x2 = Point_of_tri(tri)[1]->_coords[0];
    double y2 = Point_of_tri(tri)[1]->_coords[1];
    double z2 = Point_of_tri(tri)[1]->_coords[2];

    double x3 = Point_of_tri(tri)[2]->_coords[0];
    double y3 = Point_of_tri(tri)[2]->_coords[1];
    double z3 = Point_of_tri(tri)[2]->_coords[2];

    double tri_lx = min3(x1, x2, x3);
    double tri_ux = max3(x1, x2, x3);
    double tri_ly = min3(y1, y2, y3);
    double tri_uy = max3(y1, y2, y3);
    double tri_lz = min3(z1, z2, z3);
    double tri_uz = max3(z1, z2, z3);

    double boxLx = L[0] - (0.5+eps)*h[0];
    double boxUx = U[0] + (0.5+eps)*h[0];
    double boxLy = L[1] - (0.5+eps)*h[1];
    double boxUy = U[1] + (0.5+eps)*h[1];
    double boxLz = L[2] - (0.5+eps)*h[2];
    double boxUz = U[2] + (0.5+eps)*h[2];

    /* find the intersection box, if it is a valid box, 
     * then the bounding box of the triangle
     * and the bounding box of L, U have intersection, 
     * otherwise, they don't.
     * This method requires less comparison */
    double inter_box[6];
    inter_box[0] = max(tri_lx, boxLx);	// lower x
    inter_box[1] = min(tri_ux, boxUx);	// upper x
    inter_box[2] = max(tri_ly, boxLy);	// lower y
    inter_box[3] = min(tri_uy, boxUy);	// upper y
    inter_box[4] = max(tri_lz, boxLz);	// lower z
    inter_box[5] = min(tri_uz, boxUz);	// upper z

    if((inter_box[0] <= inter_box[1]) 
	    && (inter_box[2] <= inter_box[3]) 
	    && (inter_box[4] <= inter_box[5]))
	return 1;
    else
	return 0;

}


void ImportMeshToHiprop(INTERFACE* intfc, hiPropMesh* mesh, PP_GRID* ppgrid, int debug_step)
{
    int i,j,k;

    RECT_GRID grid = ppgrid->Zoom_grid;

    /* 1. fill mesh->tris and mesh->ps by bounding box */
    POINT *p;
    SURFACE **s;
    TRI *tri;
    double* U = ppgrid->Zoom_grid.U;
    double* L = ppgrid->Zoom_grid.L;
    double* h = ppgrid->Zoom_grid.h;

    int num_tris = 0;
    int num_pts = 0;
    

    /* 1.1 count the number of tris and label the points and triangles */
    for (s = intfc->surfaces; s && *s; ++s)
    {
        if (Boundary(*s))
            continue;
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	    if(tri_hit_box(tri,L,U,h))
	    {
		tri->local_index = 1;	// indicate that the triangle should be saved
		for (k = 0; k < 3; k++)
		{
		    p = Point_of_tri(tri)[k];
		    p->local_index = -1;
		}
		num_tris++;
	    }
	    else
		tri->local_index = 0;	// indicate that the triangle should not be saved
    }
    mesh->tris = emxCreate_int32_T(num_tris, 3);

    /* 1.2 count the number of points to allocate the memory */
    for (s = intfc->surfaces; s && *s; ++s)
    {
        if (Boundary(*s))
            continue;
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	    if(tri->local_index)
	    {
		for (k = 0; k < 3; k++)
		{
		    p = Point_of_tri(tri)[k];
		    if(p->local_index == -1)
		    {
			num_pts++;
			p->local_index = -2;
		    }
		}
	    }
    }

    mesh->ps = emxCreate_real_T(num_pts, 3);
    emxArray_real_T* pt = mesh->ps;

    /* 1.3 fill in mesh->tris and mesh->ps */
    int pt_index=1;
    int tri_index=1;
    for (s = intfc->surfaces; s && *s; ++s)
    {
        if (Boundary(*s))
            continue;
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{
	    if(tri->local_index)
	    {
		for (k = 1; k <= 3; k++)
		{
		    p = Point_of_tri(tri)[k-1];
		    if (p->local_index == -2)
		    {
			pt->data[I2dm(pt_index,1,pt->size)] = p->_coords[0];
			pt->data[I2dm(pt_index,2,pt->size)] = p->_coords[1];
			pt->data[I2dm(pt_index,3,pt->size)] = p->_coords[2];
			p->local_index = pt_index;
			pt_index++;
		    }

		    mesh->tris->data[I2dm(tri_index, k, mesh->tris->size)] = p->local_index;
		}
		tri_index++;
	    }
	}
    }

    double domain[3];
    domain[0] = grid.GU[0]-grid.GL[0];
    domain[1] = grid.GU[1]-grid.GL[1];
    domain[2] = grid.GU[2]-grid.GL[2];

    boolean_T bdry[3];

    bdry[0] = 1;
    bdry[1] = 1;
    bdry[2] = 0;

    hpInitDomainBoundaryInfo(mesh, domain, bdry);
    hpGetNbProcListAuto(mesh);

    /* 3. build pinfo for triangles and points, done */
    hpInitPInfo(mesh);
    hpBuildPInfoWithOverlappingTris(mesh);
    
//    hpPrint_pinfo(mesh);

    /* 4. cut the mesh to be no overlapping triangles, done */
    hpCleanMeshByPinfo(mesh);
    /*
    printf("After CleanMesh\n");
    printf("Number of points after cleanmesh: %d\n", mesh->ps->size[0]);
    printf("Number of tris after cleanmseh: %d\n", mesh->tris->size[0]);
    fflush(stdout);
    pp_gsync();
    */


//    hpPrint_pinfo(mesh);

    hpBuildOppositeHalfEdge(mesh);
    hpBuildIncidentHalfEdge(mesh);
    hpBuildNRingGhost(mesh,1);
    printf("Leave hpBuildNRingGhost1\n");
    fflush(stdout);
    pp_gsync();


    hpBuildOppositeHalfEdge(mesh);
    printf("Leave hpBuildOppositeHalfEdge2\n");
    hpBuildIncidentHalfEdge(mesh);
    printf("Leave hpBuildIncidentHalfEdge2\n");
    hpBuildNRingGhost(mesh,2); // TODO: check if this means 1+2 ring or max(1,2) ring?
    printf("Leave hpBuildNRingGhost2\n");
    fflush(stdout);
    pp_gsync();


    /* 5. build PUpdate info, done */
    hpBuildPUpdateInfo(mesh);
    hpMeshSmoothing(mesh, 2, "CMF");

    /* 6. Reconstruct the buffers by bounding box  */
//    printf("Before Bounding Box\n");
//    fflush(stdout);
//

    double bd_box[6];
    bd_box[0] = ppgrid->Zoom_grid.L[0] - ppgrid->Zoom_grid.h[0];
    bd_box[1] = ppgrid->Zoom_grid.U[0] + ppgrid->Zoom_grid.h[0];
    bd_box[2] = ppgrid->Zoom_grid.L[1] - ppgrid->Zoom_grid.h[1];
    bd_box[3] = ppgrid->Zoom_grid.U[1] + ppgrid->Zoom_grid.h[1];
    bd_box[4] = ppgrid->Zoom_grid.L[2] - ppgrid->Zoom_grid.h[2];
    bd_box[5] = ppgrid->Zoom_grid.U[2] + ppgrid->Zoom_grid.h[2];

    hpGetAllTrisFromBoundingBox(mesh,bd_box);
    printf("Passed ImportMeshToHiprop()\n");
}


/* read_surface, modified from read_vtk_surface */
boolean read_hiprop_surface(
	INTERFACE   *intfc,
	COMPONENT   neg_comp,
	COMPONENT   pos_comp,
	hiPropMesh  *mesh,
	SURFACE	    **ps)
{
    	printf("Getting into read_hiprop_surface()\n");

	if (mesh->tris->size[0] == 0)
	{
	    printf("No triangles in the surface.\n");
	    return YES;
	}

	SURFACE *surf;
	double coords[MAXD];
	double *vertex;
	int *index;
	int i,j,k,l,m,i1,i2,i3,num_pts,num_point_tris,num_ptris;
	int max_num_pts,num_tris;
	POINT **points,*p;
	TRI **tris,**ptris;
	INTERFACE *sav_intfc;
	double L[MAXD],U[MAXD];
	double dist,max_side,min_side,ave_side;
	int N,num_edges;

	sav_intfc = current_interface();
	set_current_interface(intfc);
	surf = make_surface(neg_comp,pos_comp,NULL,NULL);


	max_side = -HUGE;
	min_side =  HUGE;
	ave_side = 0.0;
	for (i = 0; i < 3; ++i)
	{
	    L[i] =  HUGE;
	    U[i] = -HUGE;
	}

	num_pts = mesh->ps->size[0];
	printf("In the HiPropMesh, num_pts = %d\n",num_pts);
	uni_array(&vertex,3*num_pts,DOUBLE);

	N = num_pts;
	num_pts = 0;
	for (k = 0; k < N; ++k)
	{
	    for (i = 1; i<=3; i++)
		coords[i-1] = mesh->ps->data[I2dm(k+1, i, mesh->ps->size)];

	    for (i = 0; i < 3; ++i)
	    {
		vertex[3*num_pts+i] = coords[i];
		if (L[i] > coords[i]) L[i] = coords[i];
		if (U[i] < coords[i]) U[i] = coords[i];
	    }
	    num_pts++;
	}
	printf("In the FT interface num_pts = %d\n", num_pts);
	uni_array(&points,num_pts,sizeof(POINT*));
	for (i = 0; i < num_pts; ++i)
	{
	    points[i] = Point(vertex+3*i);
	    points[i]->num_tris = 0;
	}

	num_tris = mesh->tris->size[0];
	uni_array(&tris,num_tris,sizeof(TRI*));
	num_point_tris = 0;
	for (i = 0; i < num_tris; ++i)
	{

	    i1 = mesh->tris->data[I2dm(i+1, 1, mesh->tris->size)] - 1;
	    i2 = mesh->tris->data[I2dm(i+1, 2, mesh->tris->size)] - 1;
	    i3 = mesh->tris->data[I2dm(i+1, 3, mesh->tris->size)] - 1;
	    tris[i] = make_tri(points[i1],points[i2],points[i3],
	    			NULL,NULL,NULL,NO);

	    tris[i]->surf = surf;
	    (points[i1]->num_tris)++;
	    (points[i2]->num_tris)++;
            (points[i3]->num_tris)++;
            num_point_tris += 3;
	    for (j = 0; j < 3; ++j)
	    {
	    	dist = distance_between_positions(
			Coords(Point_of_tri(tris[i])[j]),
			Coords(Point_of_tri(tris[i])[(j+1)%3]),3);
	    	ave_side += dist;
		if (max_side < dist) max_side = dist;
		if (min_side > dist) min_side = dist;
	    }
	}
	printf("Passed reading the triangles\n");
	intfc->point_tri_store = (TRI**)store(num_point_tris*sizeof(TRI*));
        ptris = intfc->point_tri_store;
	for (i = 0; i < num_pts; ++i)
        {
	    points[i]->tris = ptris;
            ptris += points[i]->num_tris;
            points[i]->num_tris = 0;
        }
	for (i = 0; i < num_tris; ++i)
        {
	    if (i != 0) 
	    {
	    	tris[i]->prev = tris[i-1];
	    	tris[i-1]->next = tris[i];
	    }
	    for (j = 0; j < 3; ++j)
            {
                p = Point_of_tri(tris[i])[j];
                p->tris[(p->num_tris)++] = tris[i];	// find the triangles connected to each point, save in intfc->point_tri_store
            }
        }
        for (i = 0; i < num_pts; ++i)
        {
            ptris = points[i]->tris;
            num_ptris = points[i]->num_tris;
            for (j = 0; j < num_ptris; ++j)	// set the neighbour triangles for each triangle
            for (k = 0; k < j; ++k)
            {
                TRI *tri1 = ptris[j];
                TRI *tri2 = ptris[k];
                for (m = 0; m < 3; ++m)
                for (l = 0; l < 3; ++l)
                {
                    if ((Point_of_tri(tri1)[m] == Point_of_tri(tri2)[(l+1)%3]) &&                         
				(Point_of_tri(tri1)[(m+1)%3] == Point_of_tri(tri2)[l]))
                    {
                        Tri_on_side(tri1,m) = tri2;
                        Tri_on_side(tri2,l) = tri1;
                    }
                }
            }
        }


	for (i = 0; i<num_tris; i++)
	{
	    for(j = 0; j<3; j++)
	    {
		if(Tri_on_side(tris[i],j)==NULL)
		    set_side_bdry(Boundary_tri(tris[i]),j,YES);
		else
		    set_side_bdry(Boundary_tri(tris[i]),j,NO);
	    }
	}

	wave_type(surf) = FIRST_PHYSICS_WAVE_TYPE;	

	ave_side /= num_pts;
	surf->num_tri = num_tris;
	first_tri(surf) = tris[0];
	last_tri(surf) = tris[num_tris-1];
	last_tri(surf)->next = tail_of_tri_list(surf);
	first_tri(surf)->prev = head_of_tri_list(surf);
	reset_intfc_num_points(surf->interface);


	*ps = surf;
	set_current_interface(sav_intfc);
	return YES;
}	/* end read_hiprop_surface */


void ExportMeshToFronTier(
	hiPropMesh* mesh, 
	Front *front,
	COMPONENT neg_comp,
	COMPONENT pos_comp)

{
    INTERFACE *intfc = front->interf;
    SURFACE *surf = NULL;
    const double eps = 10.0*MACH_EPS;
    SURFACE **s;
    POINT **p;
    NODE **node;
    CURVE **curve;
    INTERFACE *sav_intfc;


    printf("Getting into ExportMeshToFronTier\n");

    /*
     * if original intfc doesn't have triangles, intfc->surfaces==NULL, so we cannot get two comp info here,
     * we need to provide component info at the input.
     * For CONTACTOR, neg_comp=2, pos_comp=3
     */

    /*
    if(intfc->surfaces)
    {
	printf("intfc->surfaces = %p\n", intfc->surfaces);
	if (*(intfc->surfaces))
	{
	    printf("*(intfc->surfaces) = %p\n", *(intfc->surfaces));
	    if((*(intfc->surfaces))->hs)
    		printf("(*(intfc->surfaces))->hs = %p\n",(*(intfc->surfaces))->hs );
	    else
		printf("(*(intfc->surfaces))->hs = NULL\n");
	}
	else
	    printf("*(intfc->surfaces) = NULL\n");
    }
    else
	printf("intfc->surfaces = NULL\n");


    COMPONENT neg_comp = (*(intfc->surfaces))->hs->neg_comp;
    COMPONENT pos_comp = (*(intfc->surfaces))->hs->pos_comp;
    printf("neg_comp = %d, pos_comp = %d\n", neg_comp, pos_comp);
    */

    /*
    for (p = intfc->points; p && *p; p++)
	delete_point(*p);
    intfc->points = NULL;
    for (node = intfc->nodes; node && *node; node++)
	delete_node(*node);
    intfc->nodes = NULL;
    for (curve = intfc->curves; curve && *curve; curve++)
	delete_curve(*curve);
    intfc->curves = NULL;

    for (s = intfc->surfaces; s && *s; ++s)
	delete_surface(*s);
    intfc->surfaces = NULL; */

    // delete physical surface and its dependencies.
    // Need to keep bdry surfaces
    for (s = intfc->surfaces; s && *s; ++s) {
        if (Boundary(*s))
            continue;
        if (!f_delete_surface(*s))
	    clean_up(ERROR);
    }


    if (!read_hiprop_surface(intfc,neg_comp,pos_comp,
		mesh,&surf))
    {
	screen("read_hiprop_surface() failed!\n");
	clean_up(ERROR);
    }

    printf("Passed read_hiprop_surface\n");

//    print_tri_normal(intfc);


//    char intfc_after_export[256];
//    sprintf(intfc_after_export,"output10/after_export");
//    print_intfc_tmp_vtk(front, intfc_after_export,NO);

//    set_topological_grid(intfc,computational_grid(intfc));
//    printf("Passed set_topological_grid\n");


//    print_intfc_number(intfc);
    sav_intfc = current_interface();
//    print_intfc_number(sav_intfc);

    set_current_interface(intfc);
    install_subdomain_bdry_curves(intfc);
    reset_intfc_num_points(intfc);
    printf("Passed install_subdomain_bdry_curves\n");

    reset_normal_on_intfc(intfc);
    init_intfc_curvature3d(front, intfc);

    set_current_interface(sav_intfc);

    if (!consistent_interface(front->interf))
	printf("inconsistent interface!\n");

    interface_reconstructed(front->interf) = NO;

    printf("Passed ExportMeshToFronTier\n");
}

void print_tri_normal(INTERFACE* intfc)
{
    SURFACE** s;
    TRI* tri;
    int i,k;
    POINT *p[3];
    double v1[3], v2[3], v[3];

    printf("\n\nThe normal of the triangles:\n\n");
    i = 1;
    for (s = intfc->surfaces; s && *s; ++s)
    {
	for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
	{
	    for (k = 0; k < 3; k++)
		p[k] = Point_of_tri(tri)[k];

	    for(k = 0; k<3; k++)
	    {
		v1[k] = p[1]->_coords[k] - p[0]->_coords[k];
		v2[k] = p[2]->_coords[k] - p[0]->_coords[k];
	    }

	    v[0] = v1[1]*v2[2] - v1[2]*v2[1];
	    v[1] = v1[2]*v2[0] - v1[0]*v2[2];
	    v[2] = v1[0]*v2[1] - v1[1]*v2[0];

	    printf("For triangle %d, normal = {%f, %f, %f}\n",i, v[0], v[1], v[2]);
	    i++;
	}
    }
}

void print_intfc_tmp_vtk(			// Chenzhe
	Front *front,
        char *out_name,
        boolean print_in_binary)
{
	char dirname[256];
        int step = front->step; 
	int dim = front->rect_grid->dim;

	if (dim == 1) return;

	/* Create vtk directories */

        //sprintf(dirname,"%s/vtk.ts%s",out_name,right_flush(step,7));
	//if (pp_numnodes() > 1)
        //    sprintf(dirname,"%s-nd%s",dirname,right_flush(pp_mynode(),4));

        if (pp_numnodes() > 1)
            sprintf(dirname,"%s/P-%s",out_name,right_flush(pp_mynode(),4));
	else
	    sprintf(dirname,"%s",out_name);

	if (!create_directory(dirname,YES))
	{
	    screen("Cannot create directory %s\n",dirname);
	    clean_up(ERROR);
	}
        vtk_interface_plot(dirname,front->interf,print_in_binary,
		           front->time,front->step,front->coordinate);
	/*
	if (front->vtk_movie_var != NULL)
	    vtk_plot_vector_field(dirname,front);
	    */
}

void print_hpmesh(
	hiPropMesh* mesh,
	char* out_name,	// output dir name
	const char* special_name,
	int step)
{
    char dirname[128];
    char mesh_name[256];
    char fname[256];

//    sprintf(dirname,"%s/hpmesh",out_name);
//    if (!create_directory(dirname,NO))
//    {
//	screen("Cannot create directory %s\n",dirname);
//	clean_up(ERROR);
//    }

    if(special_name!=NULL)
	sprintf(mesh_name,"hpmesh-%s-t%s",special_name,right_flush(step,7));
    else
	sprintf(mesh_name,"hpmesh-t%s",right_flush(step,7));
//    printf("%s\n", mesh_name);

    sprintf(fname,"%s/%s-p%s.vtk",out_name,mesh_name,right_flush(pp_mynode(),4));
//    printf("%s\n", fname);
//    fflush(stdout);
//    sprintf(fname,"%s/%s",dirname,mesh_name);
//    printf("created the file name %s\n", fname);
    hpWriteUnstrMeshWithPInfo(fname, mesh);
//    hpWritePolyMeshVtk3d(fname, mesh);
}
