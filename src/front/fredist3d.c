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
*				fredist3d.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This file contains the high-level routine for the re-distribution 
*	of triangles on given surfaces according to given criteria,
*
*			redistribute_surface()
*
*	and several elementary routines for splitting and joining
*	triangles.  Consider two adjacent triangles with a common side.
*	The two elementary operations, which are regarded as inverses of
*	each other, are to either split the common side by insertion of
*	a new point (and thus to split the triangles) or to shrink the
*	common side, and reduce it to a single point (and thus to reduce the
*	two triangles to a pair of edges).  A further elementary operation is
*	to flip the diagonal of a pair of adjacent triangles, ie of a
*	elementary quadralateral.
*
*	Also in this file are the three D additions to insert and delete_
*	point_in_bond.	With these three D additions, these two D routines
*	applied to curves of a three D interface, will modify the triangles
*	of the surfaces bounded by the curves and bonds, so that a valid
*	interface is produced.	These subroutines must not be called 
*	independently.
*
*	Each has the functionality as stated above.
*/


#define DEBUG_STRING	"redist3d"
#include <front/fdecs.h>		/* includes int.h, table.h */

typedef  struct {
        TRI             *tri1;
        int             side1;
        TRI             *tri2;
        int             side2;
        double           dist;
        int             connect;
}       TRI_PAIR;

typedef  struct {
        TRI             **tris;
        int             *sides;
        POINT           **pts;
        int             n_sides;
        double           *tnor;
        ORIENTATION     orient;
        int             max_side;
}       NULL_LOOP;

enum _SPQ_FLAG {
	SHORTEST = 0,
	LONGEST	 = 1,
	BD_DIST	 = 2
};
typedef enum _SPQ_FLAG SPQ_FLAG;

/* enum _TRI_STATUS {
	SMALL = 1,
	LARGE,
	BAD_ANGLE,
	GOOD_ANGLE
};
typedef enum _TRI_STATUS TRI_STATUS; */

	/* LOCAL Function Declarations */
LOCAL	POINTER_Q  *alloc_and_add_to_queue(TRI*,SURFACE*,POINTER_Q*,int,RECT_GRID*);
LOCAL	POINTER_Q  *dequeue(TRI*,POINTER_Q*);
//LOCAL	TRI_STATUS tri_status(TRI*,RECT_GRID*);
LOCAL	boolean	delete_min_side_of_tri(TRI*,int,SURFACE*,POINTER_Q**,Front*);
LOCAL	boolean	flip_max_side_of_tri(TRI*,int,Front*,POINTER_Q**);
LOCAL	boolean	is_tri_in_queue(TRI*,POINTER_Q*);
LOCAL	boolean	redistribute_surface(SURFACE*,Front*);
LOCAL	boolean	redistribute_surface_temp(SURFACE*,Front*);
LOCAL   boolean redistribute_surface_init(SURFACE*,Front*);
LOCAL 	boolean 	is_critical_side(TRI*,int);
LOCAL	double	angle_weight_of_tri(double*,double*,double*,RECT_GRID*);
LOCAL	double	scaled_tri_area(double*,double*,double*,RECT_GRID*);
LOCAL	int	find_scaled_extrem_edge(TRI*,RECT_GRID*,SPQ_FLAG);
LOCAL	void	exchange_queues(POINTER_Q*,POINTER_Q*);
LOCAL	void	redistribute_curve3d(CURVE*,Front*);
LOCAL	void	sort_pointer_queue(POINTER_Q*,INTERFACE*,SPQ_FLAG);
LOCAL 	void 	print_tri_surf_queue(POINTER_Q*);  
LOCAL   void  	set_tri_status_order(int);
LOCAL   void    tecplot_tri_queue(const char*, FILE*, POINTER_Q*);
LOCAL	boolean	check_and_rm_tetrahedron(TRI*,int,POINTER_Q**,SURFACE*);
LOCAL	void	tri_queue_test(const char*, POINTER_Q*);
LOCAL	void	change_buffer_for_intfc(INTERFACE*);
LOCAL	boolean	small_area_tri_on_intfc(INTERFACE*);


#if defined(DEBUG_STRING)
#define	DEBUG_FRONT(mesg,fr)	 debug_front(DEBUG_STRING,mesg,fr);
#else /* defined(DEBUG_STRING) */
#define	DEBUG_FRONT(mesg,fr)
#endif /* defined(DEBUG_STRING) */

/*
*			redistribute3d():
*
*	Returns one of the following values
*
*	BAD_REDISTRIBUTION  	
*		if Surface_redistribute() fails.
*	      	if intersections() fails after Surface_redistribute() succeeds. 
*		if restart_init is set and the restart interface is tangled.
*						
*	UNABLE_TO_UNTANGLE	
*		if redistributed intfc is tangled and scalar_unravel fails.
*		if Surface_redistribute() fails after scalar_unravel succeeds.
*		if intersections() fails after Surface_redistribute() succeeds.
*		if new tangle found after old tangle was successfully removed.
*
*	GOOD_REDISTRIBUTION	
*		otherwise.
*
*	If the do_redist flag is NO, the redistribution step
*	is omitted.
*/

#define MAX_SMOOTH_PARA  50000



double	smooth_min_tri_area = 0.0;

LOCAL	void	set_tol_for_smooth(Front *fr)
{
	smooth_min_tri_area = 1.0e-2*sqrt(Min_tri_sqr_area(fr,GENERAL_WAVE));
	//DEBUG_TMP printf("smooth_min_tri_area = %24.16e\n", smooth_min_tri_area);
}

EXPORT	boolean	compute_smooth_para(
	SMOOTH_PARA	*smooth_para,
	POINT		*p,
	TRI		*tri,
	SURFACE		*s,
	SMOOTH_TOL	*stol)
{
	TRI	*tri2, *tri1, **ptris;
	int	n, j, k;
	POINT	*p1, *p2, *p3;
	double	hv1[3], hv2[3], ang, tri_area;
	double	dist, lenk, max_cos, avep[3];
	
	if(Boundary_point(p))
	    return NO;

	for(k=0; k<3; k++)
	    avep[k] = 0.0;
	
	lenk = 0.0;
	max_cos = -1.0;

	n = set_tri_list_around_point(p, tri, &ptris, s->interface);

	for(j=0; j<n; j++)
	{
	    tri1 = ptris[j];
	    p1 = Point_of_tri(tri1)[Next_m3(Vertex_of_point(tri1,p))];
	    p2 = Point_of_tri(tri1)[Prev_m3(Vertex_of_point(tri1,p))];
		        
	    for(k=0; k<3; k++)
		avep[k] += Coords(p1)[k];
			
	    //compute cone dist
	    lenk += distance_between_positions(Coords(p1), Coords(p2), 3);

	    //compute the smallest angle between two triangles in an edge.
	    tri2 = ptris[(j-1+n)%n];
	    p3 = Point_of_tri(tri2)[Next_m3(Vertex_of_point(tri2,p))];
	    triangle_height_vec(hv1, Coords(p), Coords(p2), Coords(p1));
	    
	    //if a very small tri is found, just skip it, 
	    //delete_min_side_of_tri will deal with this case.
	    tri_area = 0.5*distance_between_positions(Coords(p),Coords(p1),
	    			3)*Mag3d(hv1);
	    if(tri_area < smooth_min_tri_area)
		return NO;

	    triangle_height_vec(hv2, Coords(p), Coords(p3), Coords(p1));
	    ang = Dot3d(hv1,hv2)/(Mag3d(hv1)*Mag3d(hv2));
	    if(ang > max_cos)
		max_cos = ang;
	}
		    
	for(k=0; k<3; k++)
	    avep[k] /= n;
		    
	dist = distance_between_positions(avep, Coords(p), 3);
		    
	// dist/lenk, max_cos to test bad point
	if(dist/lenk > stol->cone_ratio || max_cos > stol->max_cos)
	{
	    smooth_para->pt = p;
	    smooth_para->tri = tri;
	    ft_assign(smooth_para->avep, avep, 3*FLOAT);
	    smooth_para->cor = dist/lenk > stol->cone_ratio ? dist/lenk : -1.0;
	    smooth_para->cos = max_cos > stol->max_cos ? max_cos : -1.0;
	    
	    return YES;
	}

	return NO;
}	/* end compute_smooth_para */

EXPORT	boolean	compute_average_point(
	SMOOTH_PARA	*smooth_para,
	POINT		*p,
	TRI		*tri,
	SURFACE		*s,
	SMOOTH_TOL	*stol)
{
	TRI	*tri1, **ptris;
	int	j, k, n;
	POINT	*p1;
	double	avep[3];
	
	if(Boundary_point(p))
	    return NO;

	for(k=0; k<3; k++)
	    avep[k] = 0.0;
	
	n = set_tri_list_around_point(p, tri, &ptris, s->interface);

	for(j=0; j<n; j++)
	{
	    tri1 = ptris[j];
	    p1 = Point_of_tri(tri1)[Next_m3(Vertex_of_point(tri1,p))];
		        
	    for(k=0; k<3; k++)
		avep[k] += Coords(p1)[k];
	}
		    
	for(k=0; k<3; k++)
	    avep[k] /= n;
		    
	smooth_para->pt = p;
	smooth_para->tri = tri;
	ft_assign(smooth_para->avep, avep, 3*FLOAT);
	smooth_para->cor = -1.0;
	smooth_para->cos = -1.0;
	    
	return YES;
}

EXPORT  void    detect_and_move_points(
	SURFACE		*s)
{
	TRI		*tri;
	POINT		*p;
	int		i, num;
	SMOOTH_PARA	smooth_que[MAX_SMOOTH_PARA];
	SMOOTH_TOL	stol;

	if (!(first_tri(s)))
	    return;

	//smooth paramaters.
	stol.cone_ratio = 0.24;  // h/r = 1.5
	//stol.cone_ratio = 0.2;  // h/r = 1.5
	stol.max_cos = 0.939;      //20 deg
	//stol.max_cos = 0.7;      //20 deg
	stol.alpha = sqrt(0.25);  //0.65 prev

	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		Index_of_point(Point_of_tri(tri)[i]) = -1;
	    }
	}

	num = 0;
	//Compute the the parameters in each points
	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		if(Boundary_point(p) || Index_of_point(p) != -1)
		    continue;

		Index_of_point(p) = 1;
		if(!compute_smooth_para(&smooth_que[num], p,tri,s,&stol))
		    continue;
		
		num++;
		if(num >= MAX_SMOOTH_PARA)
		{
	 	    printf("ERROR detect_and_move_points, array is too large.\n");
		    clean_up(ERROR);
		}
	    }
	}

	//DEBUG_TMP printf("#check am  %d\n", num);

	if(num > 0)
	    s->interface->modified = YES;
	
	//Apply Laplacian smooth
	for(i=0; i<num; i++)
	    compute_point_smooth(&smooth_que[i], &stol, s->interface);

}	/* end tecplot_surface_states */

#define	MAX_CURVE_PTS	2000

//must be called after sep_common_pt_for_open_bdry is called

EXPORT  void	smooth_curve(
	CURVE	*c)
{
        POINT   *p, *prevp, *nextp, *pt[MAX_CURVE_PTS];
        BOND    *b;
        int     i, k, num;
        double  nor[4], smo_fac = sqrt(0.6);
        double  newpt[MAX_CURVE_PTS][3];

        num = 0;
        for(b=c->first; b != NULL; b=b->next)
        {
            p = b->start;
            if(!point_outside_open_bdry(&k, nor, p, c->interface))
                continue;

            nextp = b->end;

            prevp = b->prev != NULL ? b->prev->start : NULL;
            if(prevp == NULL)
                if(is_closed_curve(c))
                    prevp = c->last->start;
                else
                    continue;

            for(i=0; i<3; i++)
                newpt[num][i] = Coords(p)[i]*smo_fac +
                    (Coords(prevp)[i] + Coords(nextp)[i])*(1.0-smo_fac)*0.5;
            pt[num] = p;

            num++;
            if(num >= MAX_CURVE_PTS)
            {
                printf("ERROR smooth_curves, too many points.\n");
                clean_up(ERROR);
            }
        }

        for(i=0; i<num; i++)
            ft_assign(Coords(pt[i]), newpt[i], 3*FLOAT);
}

#define	MAX_SUBDOMAIN_TRIS	3000

int	collect_bdry_tris(
		double		*m_coord,
		TRI		**tris,
		int		*sides,
		double		plane,
		int		dir,
		int		nor,
		SURFACE		*s)
{
	TRI	*tri;
	POINT	*p0, *p1;
	double	pc0, pc1, pm;
	int	i, nt;

	pm = nor == 0 ? HUGE_VAL : -HUGE_VAL;
	nt = 0;
	for(tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for(i=0; i<3; i++)
	    {
		if(is_side_bdry(tri, i))	
		    continue;
		if(Tri_on_side(tri,i) != NULL)
		    continue;
		
		p0 = Point_of_tri(tri)[i];
		p1 = Point_of_tri(tri)[Next_m3(i)];
		pc0 = Coords(p0)[dir];
		pc1 = Coords(p1)[dir];
		pm = nor==0 ? min3(pc0,pc1,pm) : max3(pc0,pc1,pm);

		//the null side should be outside the plane.
		if(nor==0 && (pc0>plane || pc1>plane))
		    continue;
		if(nor==1 && (pc0<plane || pc1<plane))
		    continue;

		tris[nt] = tri;
		sides[nt] = i;
		nt++;

		if(nt >= MAX_SUBDOMAIN_TRIS)
		{
		    printf("ERROR collect_bdry_tris "
			   "too many tris.\n");
		    clean_up(ERROR);
		}
	    }
	}
	
	*m_coord = pm;
	return nt;
}

//assume no curves are attached to s, need to be upgraded.
void	invert_surface(
	SURFACE		*s)
{
	TRI		*t;
	POINT		*p;
	TRI_NEIGHBOR	tn;
	int		i;
	Locstate	st;

	for(t = first_tri(s); !at_end_of_tri_list(t,s); t = t->next)
	{
	    for(i=0; i<3; i++)
	    {
		p = Point_of_tri(t)[i];
		sorted(p) = NO;
	    }
	}

	for(t = first_tri(s); !at_end_of_tri_list(t,s); t = t->next)
	{
	    p = Point_of_tri(t)[1];
	    Point_of_tri(t)[1] = Point_of_tri(t)[2];
	    Point_of_tri(t)[2] = p;

	    tn = Tri_neighbor(t)[0];
	    Tri_neighbor(t)[0] = Tri_neighbor(t)[2];
	    Tri_neighbor(t)[2] = tn;

	    set_normal_of_tri(t);
	/*    
	    for(i=0; i<3; i++)
	    {
		p = Point_of_tri(t)[i];
		if(!sorted(p))
		{
		    st = left_state(p);
		    left_state(p) = right_state(p);
		    right_state(p) = st;
		    sorted(p) = YES;
		}
	    }
	 */
	}
}

LOCAL  boolean  check_adjecant_constrain(
        POINT   *p0,
        POINT   *p1,
        POINT   *p2,
        POINT   *p3,
        TRI     **new_tris,
        int     num_new_tris)
{
        int     i, j;
        TRI     *tri;
        POINT   *pt, *ptn;
        boolean    status_in, status_out;

        status_in = NO;
        status_out = NO;
        for(i=0; i<num_new_tris; i++)
        {
            tri = new_tris[i];
            for(j=0; j<3; j++)
            {
                if(Tri_on_side(tri, j) != NULL)
                    continue;

                pt = Point_of_tri(tri)[j];
                ptn = Point_of_tri(tri)[(j+1)%3];

                //(1) adjecant to a null side tri is valid.
                if( ((p0 == pt) && (p2 == ptn)) || ((p3 == pt) && (p1 == ptn)) )
                    return YES;

                //(2) only the out_side or the in_side have common points with new_tris, 
                //this constraint can avoid tangled surface.
                if(p0 == pt || p0 == ptn || p1 == pt || p1 == ptn)
                    status_out = YES;
                if(p2 == pt || p2 == ptn || p3 == pt || p3 == ptn)
                    status_in = YES;
            }
        }

        if(status_in && status_out)
            return NO;
        else
            return YES;
}


LOCAL   boolean    make_tri_pairs_with_constraint(
        TRI_PAIR        *tri_pair,
        TRI             **new_tris,
        int             *n_new)
{
        TRI     *new_tri1,*new_tri2;
        TRI     *in_tri, *out_tri;
        int     in_side, out_side, num_new, i;
        SURFACE *surf;
        POINT   *p0, *p1, *p2, *p3;

        out_tri = tri_pair->tri1;
        out_side = tri_pair->side1;
        in_tri = tri_pair->tri2;
        in_side = tri_pair->side2;

        //in_tri or out_tri is already connectted, return.
        if(Tri_on_side(out_tri, out_side) != NULL)
            return NO;
        if(Tri_on_side(in_tri, in_side) != NULL)
            return NO;

        num_new = *n_new;
        surf = out_tri->surf;

        p0 = Point_of_tri(out_tri)[Next_m3(out_side)];
        p1 = Point_of_tri(out_tri)[out_side];
        p2 = Point_of_tri(in_tri)[in_side];
        p3 = Point_of_tri(in_tri)[Next_m3(in_side)];

        //logical constraint is voilated, return.
        if(!check_adjecant_constrain(p0,p1,p2,p3, new_tris, num_new))
            return NO;

        if (tri_pair->connect == 0)
        {
            //distance from p1 to p2 is less than that from p0 to p3.
            new_tri1 = make_tri(p0,p1,p2,NULL,NULL,NULL,NO);
            new_tris[num_new++] = new_tri1;
            insert_tri_at_tail_of_list(new_tri1,surf);
            link_neighbor_null_side_tris(new_tri1,out_tri);

            new_tri2 = make_tri(p3,p2,p1,NULL,NULL,NULL,NO);
            new_tris[num_new++] = new_tri2;
            insert_tri_at_tail_of_list(new_tri2,surf);
            link_neighbor_null_side_tris(new_tri2,in_tri);

            link_neighbor_null_side_tris(new_tri1,new_tri2);
        }
        else
        {
            //distance from p0 to p3 is less than that from p1 to p2.
            new_tri1 = make_tri(p0,p1,p3,NULL,NULL,NULL,NO);
            new_tris[num_new++] = new_tri1;
            insert_tri_at_tail_of_list(new_tri1,surf);
            link_neighbor_null_side_tris(new_tri1,out_tri);

            new_tri2 = make_tri(p3,p2,p0,NULL,NULL,NULL,NO);
            new_tris[num_new++] = new_tri2;
            insert_tri_at_tail_of_list(new_tri2,surf);
            link_neighbor_null_side_tris(new_tri2,in_tri);

            link_neighbor_null_side_tris(new_tri1,new_tri2);
        }

        for(i=0; i<num_new-2; i++)
        {
            link_neighbor_null_side_tris(new_tri1,new_tris[i]);
            link_neighbor_null_side_tris(new_tri2,new_tris[i]);
        }

        *n_new = num_new;
        return YES;
}

void	connect_tris(TRI**,int*,TRI**,int*,int,INTERFACE*);

void    connect_tris(
        TRI             **otris,
        int             *osides,
        TRI             **rtris,
        int             *rsides,
        int             nt,
        INTERFACE       *intfc)
{
        int             i, side, n_new, nt1;
        TRI             *tri, **ptris, *new_tris[5];
        TRI_PAIR        tri_pair;
        POINT           *p;
        int             side_map[3] = {2, 1, 0};

        for(i=0; i<nt; i++)
        {
            //printf("i=%d\n", i);
            tri = otris[i];
            side = osides[i];

            tri_pair.tri1 = tri;
            tri_pair.side1 = side;
            tri_pair.tri2 = rtris[i];
            //after invert_surface, sides index are changed.
            tri_pair.side2 = side_map[rsides[i]];
            tri_pair.connect = 0;

            n_new = 0;
            make_tri_pairs_with_constraint(&tri_pair, new_tris, &n_new);
            //printf("n_new = %d\n", n_new);

            p = Point_of_tri(tri)[side];
            nt1 = set_tri_list_around_point(p, tri, &ptris, intfc);
            link_neighbor_null_side_tris(ptris[nt1-1], new_tris[1]);

            p = Point_of_tri(rtris[i])[side_map[rsides[i]]];
            nt1 = set_tri_list_around_point(p, rtris[i], &ptris, intfc);
            link_neighbor_null_side_tris(ptris[nt1-1], new_tris[0]);
        }
}

#define MAX_TRIS_IN_POINT  100

boolean	connect_reflect_surface(
	int		dir,
	int		nor,
	double		pm,
	double		plane,
	SURFACE		*os,
	INTERFACE	*intfc)
{
	INTERFACE	*sav_intfc;
	SURFACE		*rs;
	TRI		*t;
	POINT		*p;
	TRI		*otris[MAX_SUBDOMAIN_TRIS], *rtris[MAX_SUBDOMAIN_TRIS];
	int		osides[MAX_SUBDOMAIN_TRIS], rsides[MAX_SUBDOMAIN_TRIS];
	int		i, nt;
	double		*h, pm1;
	boolean		sav_copy;
	
	DEBUG_ENTER(connect_reflect_surface)

	sav_intfc = current_interface();
	set_current_interface(intfc);
	
	sav_copy = copy_intfc_states();
	set_copy_intfc_states(YES);
	rs = copy_surface(os, NULL, NULL, YES);
	set_copy_intfc_states(sav_copy);
	
	nt = collect_bdry_tris(&pm1, otris, osides, plane, dir, nor, os);
	nt = collect_bdry_tris(&pm1, rtris, rsides, plane, dir, nor, rs);

	printf("#reflect bf  %d %24.15e\n", nt, pm);

	//reflect rs wrt pm, do not reflect states
	for(t = first_tri(rs); !at_end_of_tri_list(t,rs); t = t->next)
	{
	    for(i=0; i<3; i++)
	    {
		p = Point_of_tri(t)[i];
		sorted(p) = NO;
	    }
	}
	for(t = first_tri(rs); !at_end_of_tri_list(t,rs); t = t->next)
	{
	    for(i=0; i<3; i++)
	    {
		p = Point_of_tri(t)[i];
		if(!sorted(p))
		{
		    Coords(p)[dir] = 2*pm - Coords(p)[dir];
		    sorted(p) = YES;
		}
	    }
	}
	invert_surface(rs);

	//adjoin the reflected surface.
	last_tri(os)->next = first_tri(rs);
	first_tri(rs)->prev = last_tri(os);
	link_tri_list_to_surface(first_tri(os),last_tri(rs),os);
	(void) delete_surface(rs);

	connect_tris(otris, osides, rtris, rsides, nt, intfc);

	set_current_interface(sav_intfc);
	DEBUG_LEAVE(connect_reflect_surface)
	
	return YES;
}

double	intfc_bdry_coord(
	int		dir,
	int		nor,
	INTERFACE	*intfc)
{
	HYPER_SURF              *hs;
        HYPER_SURF_ELEMENT      *hse;
        POINT                   *p;
	double			pm;

	pm = nor == 0 ? HUGE_VAL : -HUGE_VAL;
	next_point(intfc,NULL,NULL,NULL);
        while(next_point(intfc,&p,&hse,&hs))
	{
	    pm = nor == 0 ? min(pm,Coords(p)[dir]) : max(pm,Coords(p)[dir]);
	}

	return pm;
}

void    sep_common_pt_for_open_bdry(INTERFACE*);
void sep_common_pt_for_open_bdry(
        INTERFACE       *intfc)
{
        SURFACE                 **s;
        POINT                   *p;
        TRI                     *tri, *tris[MAX_SUBDOMAIN_TRIS];
        int                     i, k, nt;
        double                   nor[4], tmp_fac;

        if(!debugging("sep_for_open"))
            return;

        nt = 0;
        for(s=intfc->surfaces; s && *s; s++)
        {
            for(tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
            {
                for(i=0; i<3; i++)
                {
                    p = Point_of_tri(tri)[i];
                    if(point_outside_open_bdry(&k, nor, p, intfc))
                    {
                        tris[nt++] = tri;
                        if(nt >= MAX_SUBDOMAIN_TRIS)
                        {
                            printf("ERROR sep_common_pt_for_open_bdry ",
                                   "too many tris.\n");
                            clean_up(ERROR);
                        }
                        break;
                    }
                }
            }
        }

        tmp_fac = get_comm_pt_fac();
        set_comm_pt_fac(0.999);
        if(nt > 0)
            sep_common_point_from_loop(tris, nt, NULL, NULL, intfc);
        set_comm_pt_fac(tmp_fac);
}

//IMPORT void sep_common_pt_for_open_bdry(INTERFACE*);

void	reflect_extend_intfc(
	Front	*fr)
{
	INTERFACE	*intfc;
	SURFACE		**s;
	RECT_GRID	*gr;
	double		*h, pm, plane;
	boolean		status;
	int		dir, nor;
	double		tol = 0.01, ref_shift = 0.1;
	
	DEBUG_ENTER(reflect_extend_intfc)
	
	intfc = fr->interf;
	gr = computational_grid(intfc);
	h = gr->h;
	
	dir = 2;
	nor = 0;
	
	//for nor==0 case ONLY, pm: the global min coords in direction dir,
	//plane: the global min reflect plane.
	pm = intfc_bdry_coord(dir,nor,intfc);
	pp_global_min(&pm, 1L);
	printf("#pm %24.15e\n", pm);

	plane = nor == 0 ? gr->VL[dir]+tol*h[dir] : gr->VU[dir]-tol*h[dir];
	pp_global_min(&plane, 1L);
	printf("#plane %24.15e\n", plane);

	//no curve formed.
	if((nor==0 && pm > plane) || (nor==1 && pm < plane))
	{
	    DEBUG_LEAVE(reflect_extend_intfc)
	    return;
	}
	
	strip_subdomain_bdry_curves(intfc);
	
	add_to_debug("sep_for_open");
	sep_common_pt_for_open_bdry(intfc);
	remove_from_debug("sep_for_open");

	printf("#ref enter %24.15e %24.15e\n", pm, plane);
	//the reflect plane.
	pm = nor == 0 ? pm - ref_shift*h[dir] : pm + ref_shift*h[dir];
	
	//null sides outside plane should be connectted with the reflected 
	//surface.
	if(rect_boundary_type(intfc, dir, nor) == OPEN_BOUNDARY)
	{
	    for(s=intfc->surfaces; s && *s; s++)
	    {
		connect_reflect_surface(dir, nor, pm, plane, *s, intfc);
	    }
	}
	
	null_sides_are_consistent();
	check_print_intfc("After connect ref intfc", "ref_af", 'g', 
			intfc, fr->step, -1, YES);

	fflush(NULL);	
	status = scatter_front(fr);
	
	check_print_intfc("After scat af connect ref intfc", "ref_scat_af", 
			'g', intfc, 1, -1, NO);

	DEBUG_LEAVE(reflect_extend_intfc)
}

EXPORT	double	min_null_pair_angle(
	double	*p0,
	double	*p1,
	double	*p2,
	double	*p3)
{
        double	hv1[3], hv2[3], ang;
	
	triangle_height_vec(hv1, p1, p0, p2);
	triangle_height_vec(hv2, p1, p3, p2);
	ang = Dot3d(hv1,hv2)/(Mag3d(hv1)*Mag3d(hv2));
	return	ang;
}

EXPORT	void	set_increased_buffer_grid(
	RECT_GRID       *rgr,
	const RECT_GRID *gr,
	int		inc,
	INTERFACE	*intfc)
{
	int		dlbuf[MAXD], dubuf[MAXD];
	int		i;

	if (rgr == NULL || gr == NULL)
	    return;
	rgr->dim = gr->dim;
	for (i = 0; i < gr->dim; ++i)
	{
	    dlbuf[i] = gr->lbuf[i];
	    rgr->GL[i] = gr->GL[i];
	    if(dlbuf[i] != 0 && 
	       rect_boundary_type(intfc,i,0) == SUBDOMAIN_BOUNDARY)
	    {
		dlbuf[i] += inc;
		rgr->GL[i] = gr->GL[i] - inc*gr->h[i];
	    }

	    dubuf[i] = gr->ubuf[i];
	    rgr->GU[i] = gr->GU[i];
	    if(dubuf[i] != 0 &&
	       rect_boundary_type(intfc,i,1) == SUBDOMAIN_BOUNDARY)
	    {
		dubuf[i] += inc;
		rgr->GU[i] = gr->GU[i] + inc*gr->h[i];
	    }
	    
	    rgr->gmax[i] = gr->gmax[i];
	    rgr->L[i] = gr->L[i];
	    rgr->U[i] = gr->U[i];
	}
	set_rect_grid(rgr->L,rgr->U,rgr->GL,rgr->GU,dlbuf,dubuf,
		      rgr->gmax,gr->dim,&gr->Remap,rgr);
}

//ref:  clip_intfc_at_grid_bdry1 for 4
EXPORT	void	set_increased_open_bdry_grid(
	RECT_GRID       *rgr,
	const RECT_GRID *gr,
	int		inc,
	INTERFACE	*intfc)
{
	int		dlbuf[MAXD], dubuf[MAXD];
	int		i;

	if (rgr == NULL || gr == NULL)
	    return;
	rgr->dim = gr->dim;
        for (i = 0; i < gr->dim; ++i)
        {
            dlbuf[i] = gr->lbuf[i];
            rgr->GL[i] = gr->GL[i];
            if(dlbuf[i] != 0 &&
               rect_boundary_type(intfc,i,0) == OPEN_BOUNDARY)
            {
                dlbuf[i] += inc;
                rgr->GL[i] = gr->GL[i] - inc*gr->h[i];
            }

            dubuf[i] = gr->ubuf[i];
            rgr->GU[i] = gr->GU[i];
            if(dubuf[i] != 0 &&
               rect_boundary_type(intfc,i,1) == OPEN_BOUNDARY)
            {
                dubuf[i] += inc;
                rgr->GU[i] = gr->GU[i] + inc*gr->h[i];
            }

            rgr->gmax[i] = gr->gmax[i];
            rgr->L[i] = gr->L[i];
            rgr->U[i] = gr->U[i];
        }
	set_rect_grid(rgr->L,rgr->U,rgr->GL,rgr->GU,dlbuf,dubuf,
		      rgr->gmax,gr->dim,&gr->Remap,rgr);
}


LOCAL	void	change_buffer_for_intfc(
	INTERFACE	*intfc)
{
	RECT_GRID	Dual_grid, c_gr, *comp_grid;

	c_gr = Computational_grid(intfc);
	comp_grid = computational_grid(intfc);
	set_increased_buffer_grid(comp_grid, &c_gr, 3, intfc);
        set_dual_grid(&Dual_grid, comp_grid);
	set_expanded_grid(&Dual_grid, &topological_grid(intfc));

	//DEBUG_TMP printf("#change rect grid.\n");
	if (debugging("change_buffer"))
	    print_rectangular_grid(&topological_grid(intfc));
}

LOCAL	int	repeat_count = 0;

EXPORT	void	set_repeat(int	cnt)
{
	repeat_count = cnt;
}
EXPORT	int	recon_repeat()
{
	return repeat_count;
}

int smooth_redistribute(Front*);
int smooth_redistribute_temp(Front*);
int smooth_redistribute_init(Front*);

LOCAL   boolean use_bd_dist = NO;
void    set_use_bd_dist(boolean);
void    set_use_bd_dist(boolean fg)
{
        use_bd_dist = fg;
        if (fg == YES)
            printf("Set use_bd_dist to YES\n");
}

LOCAL   boolean use_bd_dist_smooth = NO;
void    set_use_bd_dist_smooth(boolean);
void    set_use_bd_dist_smooth(boolean fg)
{
        use_bd_dist_smooth = fg;
        if (fg == YES)
            printf("Set use_bd_dist_smooth to YES\n");
}

EXPORT int smooth_redistribute_init(
	Front		*fr)
{
	INTERFACE	*intfc;
	boolean		status;
	boolean		force_redistribute = YES;
	int		k, i;

	DEBUG_ENTER(smooth_redistribute_init)

	intfc = fr->interf;
	/*
	if(fr->step > 0)
	{
	    SURFACE		**s;
	    CURVE		**c;

	    set_tol_for_smooth(fr);

	    for(i=0; i<1; i++)
	    {
		for (s = intfc->surfaces; s && *s; ++s)
		{
	    	    if(wave_type(*s) != FIRST_SCALAR_PHYSICS_WAVE_TYPE)
			continue;
	    
		    printf("#detect_and_move_points\n");
		    detect_and_move_points(*s);
		}
	
		//printf("\nAfter moving the points\n");
		start_clock("scatter_front");
		status = scatter_front(fr);
		stop_clock("scatter_front");


		if(!status)
		{
		    printf("ERROR in redistribute3d, scatter_front fails "
		    	   "after detect_and_move_points.\n");
                    return NO;
		    //clean_up(ERROR);
		}
	
	    }
	}
	*/
	{
	    int		num_iter = 0;
            int         max_num_iter = 10;
            fr->is_bad_tri = YES;	    
            Interface_redistributed(fr) = NO;

            while(fr->is_bad_tri == YES)
	    {
                num_iter++;
		//printf("\nNum of iter in repairing: %d\n", num_iter);

                if(num_iter > max_num_iter)
                    break;
                if(num_iter > 1)
                    set_use_bd_dist(NO);
                fr->is_bad_tri = NO;
                force_redistribute = YES;
	       
		status = (Interface_redistributed(fr) == YES) ?
		        FUNCTION_SUCCEEDED :
		        surface_redistribute_init(fr,&force_redistribute);

                if(fr->is_bad_tri == YES)
		    Interface_redistributed(fr) = NO;
		
		//There are two cases, 
		//    1. redistribute_surface fails for one proc, intfc is ok,
		//    2. the scatter_front in surface_redistribute fails, 
		//    intfc is changed and is not consistent.
	    	if (status == FUNCTION_FAILED)
		{
		    printf("WARNING in redistribute3d, Surface_redistribute "
		    	   "fails. %d  %d\n", fr->step, i);
		    DEBUG_LEAVE(smooth_redistribute)
		    return BAD_REDISTRIBUTION;
		}

	    	check_print_intfc("After surface redist", "sur_redist_af", 'g', 
	        	      fr->interf, fr->step, -1, YES);
	
		//printf("#redist finish %d\n\n", i);
	
		if(debugging("pt_surface"))
		    clean_up(0);
	    }

	    if(debugging("pt_surface"))
		clean_up(0);
	    
            set_use_bd_dist(NO);
	    intfc = fr->interf;
	}
	
	DEBUG_LEAVE(smooth_redistribute_init)
	return FUNCTION_SUCCEEDED;
}


EXPORT int smooth_redistribute_temp(
	Front		*fr)
{
	INTERFACE	*intfc;
	boolean		status;
	boolean		force_redistribute = YES;
	int		k, i;

	DEBUG_ENTER(smooth_redistribute_temp)

	intfc = fr->interf;
	/*
	if(fr->step > 0)
	{
	    SURFACE		**s;
	    CURVE		**c;

	    set_tol_for_smooth(fr);

	    for(i=0; i<1; i++)
	    {
		for (s = intfc->surfaces; s && *s; ++s)
		{
	    	    if(wave_type(*s) != FIRST_SCALAR_PHYSICS_WAVE_TYPE)
			continue;
	    
		    printf("#detect_and_move_points\n");
		    detect_and_move_points(*s);
		}
	
		//printf("\nAfter moving the points\n");
		start_clock("scatter_front");
		status = scatter_front(fr);
		stop_clock("scatter_front");


		if(!status)
		{
		    printf("ERROR in redistribute3d, scatter_front fails "
		    	   "after detect_and_move_points.\n");
                    return NO;
		    //clean_up(ERROR);
		}
	
	    }
	}
	*/
	{
	    int		num_iter = 0;
            int         max_num_iter;
            fr->is_bad_tri = YES;	    
            Interface_redistributed(fr) = NO;

            if(fr->step == 0)
                max_num_iter = 20;
            else
            {
                if(use_bd_dist_smooth == YES)
                    max_num_iter = 1;
                else
                    max_num_iter = 10;
            }

            while(fr->is_bad_tri == YES)
	    {
                num_iter++;
		//printf("\nNum of iter in repairing: %d\n", num_iter);

                if(num_iter > max_num_iter)
                    break;
                if(num_iter > 1)
                    set_use_bd_dist(NO);
                fr->is_bad_tri = NO;
                force_redistribute = YES;
	       
		status = (Interface_redistributed(fr) == YES) ?
		        FUNCTION_SUCCEEDED :
		        surface_redistribute_temp(fr,&force_redistribute); /*surface_redistribute_temp called*/

                if(fr->is_bad_tri == YES)
		    Interface_redistributed(fr) = NO;
		
		//There are two cases, 
		//    1. redistribute_surface fails for one proc, intfc is ok,
		//    2. the scatter_front in surface_redistribute fails, 
		//    intfc is changed and is not consistent.
	    	if (status == FUNCTION_FAILED)
		{
		    printf("WARNING in redistribute3d, surface_redistribute_temp "
		    	   "fails. %d  %d\n", fr->step, i);
		    DEBUG_LEAVE(smooth_redistribute)
		    return BAD_REDISTRIBUTION;
		}

	    	check_print_intfc("After surface redist", "sur_redist_af", 'g', 
	        	      fr->interf, fr->step, -1, YES);
	
		//printf("#redist finish %d\n\n", i);
	
		if(debugging("pt_surface"))
		    clean_up(0);
	    }

	    if(debugging("pt_surface"))
		clean_up(0);
	    
            set_use_bd_dist(NO);
	    intfc = fr->interf;
	}
	
	DEBUG_LEAVE(smooth_redistribute_temp)
	return FUNCTION_SUCCEEDED;
}

EXPORT int smooth_redistribute(
	Front		*fr)
{
	INTERFACE	*intfc;
	boolean		status;
	boolean		force_redistribute = YES;
	int		k, i;

	DEBUG_ENTER(smooth_redistribute)

	intfc = fr->interf;
	/*
	if(fr->step > 0)
	{
	    SURFACE		**s;
	    CURVE		**c;

	    set_tol_for_smooth(fr);

	    for(i=0; i<1; i++)
	    {
		for (s = intfc->surfaces; s && *s; ++s)
		{
	    	    if(wave_type(*s) != FIRST_SCALAR_PHYSICS_WAVE_TYPE)
			continue;
	    
		    printf("#detect_and_move_points\n");
		    detect_and_move_points(*s);
		}
	
		//printf("\nAfter moving the points\n");
		start_clock("scatter_front");
		status = scatter_front(fr);
		stop_clock("scatter_front");


		if(!status)
		{
		    printf("ERROR in redistribute3d, scatter_front fails "
		    	   "after detect_and_move_points.\n");
                    return NO;
		    //clean_up(ERROR);
		}
	
	    }
	}
        */
	
	{
	    int		num_iter = 0;
            int         max_num_iter;
            fr->is_bad_tri = YES;	    
            Interface_redistributed(fr) = NO;

            if(fr->step == 0)
                max_num_iter = 20;
            else
                max_num_iter = 10;

            while(fr->is_bad_tri == YES)
	    {
                num_iter++;
		//printf("\nNum of iter in repairing: %d\n", num_iter);

                if(num_iter > max_num_iter)
                    break;
                if(num_iter > 1)
                    set_use_bd_dist(NO);
                fr->is_bad_tri = NO;
                force_redistribute = YES;
	       
		status = (Interface_redistributed(fr) == YES) ?
		        FUNCTION_SUCCEEDED :
		        Surface_redistribute(fr,&force_redistribute); /*surface_redistribute called*/

                if(fr->is_bad_tri == YES)
		    Interface_redistributed(fr) = NO;
		
		//There are two cases, 
		//    1. redistribute_surface fails for one proc, intfc is ok,
		//    2. the scatter_front in surface_redistribute fails, 
		//    intfc is changed and is not consistent.
	    	if (status == FUNCTION_FAILED)
		{
		    printf("WARNING in redistribute3d, Surface_redistribute "
		    	   "fails. %d  %d\n", fr->step, i);
		    DEBUG_LEAVE(smooth_redistribute)
		    return BAD_REDISTRIBUTION;
		}

	    	check_print_intfc("After surface redist", "sur_redist_af", 'g', 
	        	      fr->interf, fr->step, -1, YES);
	
		//printf("#redist finish %d\n\n", i);
	
		if(debugging("pt_surface"))
		    clean_up(0);
	    }

	    if(debugging("pt_surface"))
		clean_up(0);
	    
            set_use_bd_dist(NO);
	    intfc = fr->interf;
	}
	
	DEBUG_LEAVE(smooth_redistribute)
	return FUNCTION_SUCCEEDED;
}


LOCAL	boolean  small_area_tri_on_intfc(
	INTERFACE	*intfc)
{
	const double	*nor;
	SURFACE		**s;
	TRI		*t;
	RECT_GRID	*gr = &topological_grid(intfc);
	double		tol = min3(gr->h[0],gr->h[1],gr->h[2]);
	double		sqr_area;

	printf("In small_area_tri_on_intfc(), use_bd_dist_smooth = %d\n", use_bd_dist_smooth);
        printf("In small_area_tri_on_intfc(), use_bd_dist = %d\n", use_bd_dist);


	//tol = tol*tol*1.0e-10;
        //tol = 2.0*smooth_min_tri_area;
        tol = tol*tol*1.0e-8;

	for (s = intfc->surfaces; s && *s; ++s)
	    for (t = first_tri(*s); !at_end_of_tri_list(t,*s); t = t->next)
	    {
		nor = Tri_normal(t);
		sqr_area = 0.25*Dot3d(nor,nor);
		//if ( (use_bd_dist_smooth && Mag3d(nor) < tol) ||
		//     (!use_bd_dist_smooth && sqr_area <= min_sqr_area) )
                if (Mag3d(nor) < tol)
		{
		    //printf("WARNING small_area_tri_on_intfc,  "
		    //	   "bad normal appears.\n");
		    //print_tri(t, intfc);
		    //print_general_vector("nor=", nor, 3, "\n");
		    return NO;
		}
	    }

	//DEBUG_TMP printf("#small_area_tri_on_intfc, tri area is valid.\n");
	return YES;
}

boolean	set_use_rect_tris(boolean);
//LOCAL   boolean    use_rect_tris = YES;
//boolean    set_use_rect_tris(boolean);
//boolean    set_use_rect_tris(boolean fg)
//{
//        use_rect_tris = fg;
//}

LOCAL   boolean    no_redist = NO;
void    set_no_redist(boolean);
void    set_no_redist(boolean flag)
{
        no_redist = flag;
        if (flag == YES)
            printf("set_no_redist to YES\n");
}

LOCAL   boolean    no_redist_small = NO;
void    set_no_redist_small(boolean);
void    set_no_redist_small(boolean flag)
{
        no_redist_small = flag;
        if (flag == YES)
            printf("set_no_redist_small to YES\n");
}

EXPORT int redistribute3d(
	Front		*fr,
	boolean		do_redist,
	boolean		restart_init)
{
	CROSS		*cross;
	INTERFACE	*intfc;
	boolean		status;
	boolean		force_redistribute = YES;
	int		k, i;
	static int	cnt = 0;

        SURFACE **s;
        TRI *tri;
        double min = HUGE;
        double max = 0.0;
        int p_in_domain = 0;
        POINT *p;
        double min_diameter_tri, max_diameter_tri;
        double prec[3][3];
        double edge[3][3];
        int m;

	DEBUG_ENTER(redistribute3d)
	start_clock("redist_and_untangle");

	intfc = fr->interf;
	
	if(fr->step == -1)
	{
	    char  sn[30];
	    sprintf(sn, "pt_surface_bf%d", fr->step);
	    tecplot_interface_in_ball(sn, fr->interf);
	}
        set_use_bd_dist_smooth(NO);
        if(no_redist == YES)
            do_redist = NO;

        if(no_redist_small == YES)
        {
            if(no_redist == YES)
                printf("No redistribution\n");
            else
            {
                status = smooth_redistribute_init(fr);
                status = pp_min_status(status);
                if(!status)
                {
                    DEBUG_LEAVE(redistribute3d)
                    set_no_redist(YES);
                    return MODIFY_TIME_STEP;
                }
            }
        }
        else
        {
            /*
            start_clock("Mesh_Smoothing");
	    if ((do_redist && fr->step > 0) || fr->step == 0)
	    {
	        SURFACE		**s;
	        CURVE		**c;


	        set_tol_for_smooth(fr);

	        cnt++;

	        for(i=0; i<1; i++)
	        {

		    for (s = intfc->surfaces; s && *s; ++s)
		    {
	    	        if(wave_type(*s) != FIRST_SCALAR_PHYSICS_WAVE_TYPE)
			    continue;
	    
		        detect_and_move_points(*s);
		    }
		    start_clock("scatter_front");
		    status = scatter_front(fr);
		    stop_clock("scatter_front");

		    if(!status)
		    {
		        printf("ERROR in redistribute3d, scatter_front fails "
		    	       "after detect_and_move_points.\n");
		        stop_clock("Mesh_Smoothing");
		        stop_clock("redist_and_untangle");
                        return NO;
		        //clean_up(ERROR);
		    }
	        }
	    }
	    stop_clock("Mesh_Smoothing");
            */
	    if(fr->step == 0)
	        null_sides_are_consistent();


	    start_clock("Mesh_Reparing");

	    if ((do_redist && fr->step > 0) || fr->step == 0)
	    {
	        int		num_iter = 0;
                fr->is_bad_tri = YES;
                Interface_redistributed(fr) = NO;

                while(fr->is_bad_tri == YES)
                {
                    num_iter++;
                    if(num_iter > 10)
                        break;
                    if(num_iter > 1)
                        set_use_bd_dist(NO);
                    fr->is_bad_tri = NO;
	            force_redistribute = YES;

		    status = (Interface_redistributed(fr) == YES) ?
		            FUNCTION_SUCCEEDED :
		            Surface_redistribute(fr,&force_redistribute);/* surface_redistribute called */

                    if(fr->is_bad_tri == YES)
                        Interface_redistributed(fr) = NO;
		
		    //There are two cases, 
		    //    1. redistribute_surface fails for one proc, intfc is ok,
		    //    2. the scatter_front in surface_redistribute fails, 
		    //    intfc is changed and is not consistent.
	    	    if (status == FUNCTION_FAILED)
		    {
                        set_no_redist_small(YES);
		        printf("WARNING in redistribute3d, Surface_redistribute "
		    	       "fails. %d  %d\n", fr->step, i);
		        stop_clock("Mesh_Reparing");
		        stop_clock("redist_and_untangle");
		        DEBUG_LEAVE(redistribute3d)
		        return BAD_REDISTRIBUTION;
		    }
	
		    if(debugging("pt_surface"))
		        clean_up(0);
	        }

	        set_use_bd_dist(NO);
	        intfc = fr->interf;

                {
                    status = smooth_redistribute_temp(fr);
                    set_use_bd_dist(NO);
                    status = pp_min_status(status);
                    if(!status)
                    {
                        DEBUG_LEAVE(redistribute3d)
                        set_no_redist_small(YES);
                        return MODIFY_TIME_STEP;
                    }
                }
	    }
	    stop_clock("Mesh_Reparing");
        }

	if (Tracking_algorithm(fr) == LOCALLY_GRID_BASED_TRACKING)
	{
	    RECT_GRID   ngr, sav_c_gr, sav_t_gr, *sav_fr_gr;

	    if(fr->step == 0)
	        null_sides_are_consistent();

	    //set_use_rect_tris(NO);

	    if(recon_repeat() >= 1)
	    {
		printf("\ncecon_repeat() >= 1\n");
                RECT_GRID	sav_c_gr, sav_t_gr, *sav_fr_gr;

	        sav_c_gr = Computational_grid(intfc);
	        sav_t_gr = topological_grid(intfc);
	    
	        //use a large grid, scatter_front will construct a large 
	        //buffer for intfc.
	        change_buffer_for_intfc(intfc);
	    
	        //after changing the grid, do not comm component because component3d 
	        //will be called and make_tri_lists is needed.
	        //set_default_comp(YES);
	        start_clock("scatter_front");
	        status = scatter_front(fr);
	        stop_clock("scatter_front");
	        //set_default_comp(NO);
	
	        if(fr->step == 0)
	            null_sides_are_consistent();
	        //DEBUG_TMP check_print_intfc("Before surface recon in repeat", 
	    	//DEBUG_TMP "repair_bf", 'g', fr->interf, fr->step, -1, YES);
	
	        //add_to_debug("db_un");
	        //add_to_debug("box_intfc");
	
	        sav_fr_gr = fr->rect_grid;
	        fr->rect_grid = computational_grid(intfc);
   
	        /* repair_intfc_at_crossings3d */


	        start_clock("repair_front_at_grid_crossing");
	        status = repair_front_at_grid_crossing(fr); /* repair_intfc_at_crossings3d */
	        stop_clock("repair_front_at_grid_crossing");
	
	        //recover the previous grid.
	        fr->rect_grid = sav_fr_gr;
	        Computational_grid(intfc) = sav_c_gr;
	        topological_grid(intfc) = sav_t_gr;
	
	        Computational_grid(fr->interf) = sav_c_gr;
	        topological_grid(fr->interf) = sav_t_gr;

	        set_use_rect_tris(YES);
	    }
	    else
	    {
		if(debugging("ref_scat"))
		{
		    set_increased_open_bdry_grid(&ngr,
				    fr->rect_grid,4,fr->interf);
		    sav_fr_gr = fr->rect_grid;
		    fr->rect_grid = &ngr;
		}

		/* repair_intfc_at_crossings3d */
		start_clock("repair_front_at_grid_crossing");
		status = repair_front_at_grid_crossing(fr); /* repair_intfc_at_crossings3d */
		stop_clock("repair_front_at_grid_crossing");


		if(debugging("ref_scat"))
		    fr->rect_grid = sav_fr_gr;
	    }

	    if(fr->step == 0)
	        null_sides_are_consistent();
	    //DEBUG_TMP check_print_intfc("After surface recon", "repair_af", 'g', 
	        	//DEBUG_TMP fr->interf, fr->step, -1, NO);

            status = pp_min_status(status);
	    if (!status)
            {
                fprintf(stdout, "after untangling\n");
                //set_increased_recon_box(YES);
                return UNABLE_TO_UNTANGLE;
            }

	    //DEBUG_TMP if(fr->step == 458)
	    //DEBUG_TMP {
	        //add_to_debug("step_out");
    	        //add_to_debug("out_surf");
		//add_to_debug("err_stop");
	    //DEBUG_TMP }

            //if ((do_redist && fr->step > 0) || fr->step == 0)
            {
            //check for bad triangles.
                //if(use_bd_dist_smooth == YES)
                //    set_use_bd_dist(YES);
                status = smooth_redistribute_init(fr);
                //set_use_bd_dist(NO);
                status = pp_min_status(status);
                if(!status)
                {
                    DEBUG_LEAVE(redistribute3d)
                    set_no_redist_small(YES); 
                    //set_increased_recon_box(YES);
                    //set_use_bd_dist_smooth(YES);
                    return MODIFY_TIME_STEP;
                }
            }

	    //check for small triangles.
            set_tol_for_smooth(fr);
	    status = small_area_tri_on_intfc(fr->interf);
            //pp_gsync();
	    if(!pp_min_status(status))
            {
		printf("After lgb, small triangle appears.\n");
                //set_use_bd_dist(YES);
                set_use_bd_dist_smooth(YES);
		status = smooth_redistribute_temp(fr);
                set_use_bd_dist_smooth(NO);
                set_use_bd_dist(NO);
                status = pp_min_status(status);
		if(!status)
		{
		    DEBUG_LEAVE(redistribute3d)
                    set_no_redist(YES);
                    //set_increased_recon_box(YES);
                    //set_use_bd_dist_smooth(YES);
		    return MODIFY_TIME_STEP;
                    //return status;
		}
	    }
            set_no_redist(NO);
            set_no_redist_small(NO);
            //set_increased_recon_box(NO);
            //set_use_bd_dist_smooth(NO);
	    
	    start_clock("scatter_front");
	    
	    add_to_debug("sep_for_open");
	    remove_from_debug("ref_scat");
	    
	    status = scatter_front(fr);
	    
	    remove_from_debug("sep_for_open");
	    stop_clock("scatter_front");
	    
	    
	    if (!status)
	    {
		printf("WARNING, scatter_front fails after "
		       "locally reconstruction in redistribute3d.\n");
		
		set_use_rect_tris(NO);

		stop_clock("redist_and_untangle");
		DEBUG_LEAVE(redistribute3d)
		return INCONSISTENT_RECONSTRUCTION;
	    }
	    
	    if(fr->step == 0)
	        null_sides_are_consistent();

	    stop_clock("redist_and_untangle");

	    DEBUG_LEAVE(redistribute3d)
	    return GOOD_REDISTRIBUTION;
	}

	//for LGB belowing is never reached.

	if (debugging("gvrdst3d"))
	{
	    char s[120];

	    (void) sprintf(s,"before-intersect-ts%d",fr->step);
	    gview_plot_interface(s,fr->interf);
	}
	DEBUG_FRONT("before intersections check",fr)

		/* Check for Intersections in Front */

	  /* intersections does one of  the following:
	     1) returns NO if topology construction fails;
	     2) crashes if add_to_pointers() fails;
	     3) returns YES otherwise.
	  */

	if (pp_min_status(intersections(intfc,&cross,YES)) == FUNCTION_FAILED)
	{
	    stop_clock("redistribute");
	    (void) printf("WARNING in redistribute3d(), "
	                  "intersections() failed\n");
	    print_interface(intfc);
	    DEBUG_LEAVE(redistribute3d)
	    return BAD_REDISTRIBUTION;
	}

	if (debugging("gvrdst3d"))
	{
	    char s[120];

	    (void) sprintf(s,"before-untangle-ts%d",fr->step);
	    gview_plot_interface(s,fr->interf);
	}

	if (interface_is_tangled(cross))
	{
	    static const int Max_nattemps = 3;
	    int              nattemps;
	    (void) print_number_of_tangles("",intfc,cross);
	    start_clock("untangle");
	    if (restart_init) 
	    {
		stop_clock("untangle");
		stop_clock("redistribute");
		(void) printf("WARNING in redistribute(), "
		              "Restart interface tangled, cannot continue\n");
		DEBUG_LEAVE(redistribute3d)
		return BAD_REDISTRIBUTION;
	    }

	    nattemps = 0;
	    while (cross) 
	    {
		++nattemps;
		if (!scalar_unravel_3d(fr,&cross))
		{
		    stop_clock("untangle");
		    stop_clock("redistribute");
		    (void) printf("WARNING in redistribute3d(), "
		                  "scalar_unravel_3d() failed\n");
		    DEBUG_LEAVE(redistribute3d)
		    return UNABLE_TO_UNTANGLE;
		}
		force_redistribute = YES;
		if (!Surface_redistribute(fr,&force_redistribute))
		{
		    stop_clock("untangle");
		    stop_clock("redistribute");
		    (void) printf("WARNING in redistribute3d(), after "
		                  "untangling Surface_redistribute() failed\n");
		    DEBUG_LEAVE(redistribute3d)
		    return UNABLE_TO_UNTANGLE;
		}
		intfc = fr->interf;
		if (!pp_min_status(intersections(intfc,&cross,YES)))
		{
		    (void) printf("WARNING in redistribute3d(), "
		                  "After untangle, intersections() failed\n");
		    DEBUG_LEAVE(redistribute3d)
		    return UNABLE_TO_UNTANGLE;
		}
		if (interface_is_tangled(cross))
		{
		    if (nattemps>=Max_nattemps)
		    {
		        (void) printf("WARNING in redistribute3d(), "
				      "After untangle, intersections() finds "
				      "new tangle, too many attemps\n");
		        DEBUG_LEAVE(redistribute3d)
		        return UNABLE_TO_UNTANGLE;
		    }
		    else
		    {
		        (void) printf("WARNING in redistribute3d(), "
				      "After untangle, intersections() finds "
				      "new tangle, will attempt to untangle\n");
		    }
		}
		
	    }
	    stop_clock("untangle");
	}

	if (debugging("gvrdst3d"))
	{
	    char s[120];

	    (void) sprintf(s,"after-untangle-ts%d",fr->step);
	    gview_plot_interface(s,fr->interf);
	}

	stop_clock("redist_and_untangle");
	DEBUG_LEAVE(redistribute3d)
	return GOOD_REDISTRIBUTION;
}		/*end redistribute3d*/

void	set_use_bd_dist(boolean);
void    set_use_bd_dist_smooth(boolean);
boolean	is_use_bd_dist();

/*
* 		surface_redistribute()
*
* 	The choice of name surface_redistribute() clashes in an obvious 
* 	way with the name redistribute_surface().  
*/

EXPORT boolean surface_redistribute_init(
	Front		*fr,
	boolean		*force_redist)
{
	SURFACE		**s;
	CURVE		**c;
	boolean		status, min_status, delete_flag;
	boolean		redist_non_vec_cur, redist_vec_cur;
	INTERFACE	*intfc = fr->interf;
	RECT_GRID	ngr, sav_c_gr, sav_t_gr, *sav_fr_gr;

	DEBUG_ENTER(surface_redistribute_init)

		// Check on redistribution conditions

	if (*force_redist)
	{
	    redist_non_vec_cur = redist_vec_cur = YES;
	    *force_redist = NO;
	}
	else if (Redistribution_count(fr) < 0)
	{
	    redist_non_vec_cur = redist_vec_cur = NO;
	}
	else
	{
	    redist_vec_cur = redist_needed(fr,VECTOR_WAVE);
	    redist_non_vec_cur = redist_needed(fr,GENERAL_WAVE);
	}
	++Redistribution_count(fr);

		// Redistribute vector surfaces
	
	set_size_of_intfc_state(size_of_state(fr->interf));
	set_copy_intfc_states(YES);
	set_current_interface(fr->interf);


	status = YES;
	if (redist_vec_cur == YES)
	{
	    for (s = fr->interf->surfaces; s && *s ; ++s) 
	    {
	    	if ((!omit_redistribution(*s)) &&
		    (wave_type(*s) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE))
	    	    if (!redistribute_surface_init(*s,fr))
			status = NO;
	    }
	}
	
	start_clock("Add_and_Delete_Points");
	if(is_use_bd_dist())
	{

	    printf("#extend intfc in redist.\n");
	    sav_c_gr = Computational_grid(intfc);
	    sav_t_gr = topological_grid(intfc);
	    
	    //use a large grid, scatter_front will construct a large 
	    //buffer for intfc.
	    change_buffer_for_intfc(intfc);
	   
	    //after changing the grid, do not comm component because component3d 
	    //will be called and make_tri_lists is needed.
	    //set_default_comp(YES);
	    start_clock("scatter_front");
	    status = scatter_front(fr);
	    stop_clock("scatter_front");
	    //set_default_comp(NO);
	    
	    sav_fr_gr = fr->rect_grid;
	    fr->rect_grid = computational_grid(intfc);
   
	}
	
	if (redist_non_vec_cur == YES)
	{
	    //printf("#redist curve3d bf\n");
	    //for (c = fr->interf->curves; c && *c; ++c)
	    //{
		//if (hsbdry_type(*c) != SUBDOMAIN_HSBDRY)
		//redistribute_curve3d(*c,fr);
	    //}

	    //DEBUG_TMP printf("#redist surface bf\n");
	    for (s = fr->interf->surfaces; s && *s ; ++s) 
	    {
	    	if ((!omit_redistribution(*s)) &&
	    	    (wave_type(*s) >= FIRST_PHYSICS_WAVE_TYPE) &&
	    	    (wave_type(*s) < FIRST_VECTOR_PHYSICS_WAVE_TYPE))
		{
		    start_clock("redistribute_surface_init");
		    if (!redistribute_surface_init(*s,fr))
			status = NO;
		    stop_clock("redistribute_surface_init");
		}
	    }

	    delete_flag = YES;
	    while(delete_flag)
	    {
		delete_flag = NO;
		for (s = fr->interf->surfaces; s && *s; ++s)
		    if ((*s)->num_tri == 0)
		    {
			printf("#delete af redis.\n");
			delete_surface(*s);
			delete_flag = YES;
			break;
		    }
	    }
	}
	stop_clock("Add_and_Delete_Points");

        boolean status_is_bad_tri;
        status_is_bad_tri = pp_max_status(fr->is_bad_tri);
	
        if(status_is_bad_tri == YES)
            fr->is_bad_tri = YES;

	//insert_point_in_tri_side fails or delete_min_side_of_tri fails.
	if(status == NO)
	{
	    printf("WARNING after #surface_redistribute_init, "
	    	   "status is NO in step %d\n", fr->step);
	}
	min_status = pp_min_status(status);
	if (min_status == NO)
	{
	    (void) printf("WARNING in surface_redistribute_init(), "
		          "redistribute_surface_init(), failed\n");
	    DEBUG_LEAVE(surface_redistribute)
	    return min_status;
	}

	if(is_use_bd_dist())
	{
	    //recover the previous grid.
	    fr->rect_grid = sav_fr_gr;
	    Computational_grid(intfc) = sav_c_gr;
	    topological_grid(intfc) = sav_t_gr;
	
	    Computational_grid(fr->interf) = sav_c_gr;
	    topological_grid(fr->interf) = sav_t_gr;
	}
	else
	{
	    //DEBUG_TMP check_print_intfc("After redistribute_surface", "redist_sf_af", 'g', 
		//DEBUG_TMP	fr->interf, fr->step, -1, NO);
	    //DEBUG_TMP make_interface_topology_lists(fr->interf);
	    //DEBUG_TMP check_valid_intfc("after redistribute_surface", fr->interf);
	}

	start_clock("scatter_front");
	status = scatter_front(fr);
	stop_clock("scatter_front");
	
	if (!status)
	{
	    set_use_bd_dist(YES);
	    (void) printf("WARNING in surface_redistribute_init(), "
	    	          "scatter_front() failed\n");
	    DEBUG_LEAVE(surface_redistribute_init)
	    return status;
	}

	set_current_interface(fr->interf);

	Interface_redistributed(fr) =
	    ((redist_vec_cur==YES) || (redist_non_vec_cur==YES)) ? YES : NO;

	DEBUG_LEAVE(surface_redistribute_init)
	return status;
}		/*end surface_redistribute_init*/


EXPORT boolean surface_redistribute_temp(
	Front		*fr,
	boolean		*force_redist)
{
	SURFACE		**s;
	CURVE		**c;
	boolean		status, min_status, delete_flag;
	boolean		redist_non_vec_cur, redist_vec_cur;
	INTERFACE	*intfc = fr->interf;
	RECT_GRID	ngr, sav_c_gr, sav_t_gr, *sav_fr_gr;

	DEBUG_ENTER(surface_redistribute_temp)

		// Check on redistribution conditions

	if (*force_redist)
	{
	    redist_non_vec_cur = redist_vec_cur = YES;
	    *force_redist = NO;
	}
	else if (Redistribution_count(fr) < 0)
	{
	    redist_non_vec_cur = redist_vec_cur = NO;
	}
	else
	{
	    redist_vec_cur = redist_needed(fr,VECTOR_WAVE);
	    redist_non_vec_cur = redist_needed(fr,GENERAL_WAVE);
	}
	++Redistribution_count(fr);

		// Redistribute vector surfaces
	
	set_size_of_intfc_state(size_of_state(fr->interf));
	set_copy_intfc_states(YES);
	set_current_interface(fr->interf);


	status = YES;
	if (redist_vec_cur == YES)
	{
	    for (s = fr->interf->surfaces; s && *s ; ++s) 
	    {
	    	if ((!omit_redistribution(*s)) &&
		    (wave_type(*s) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE))
                {

                    printf("redist_vec_cur == YES\n");

                    start_clock("redistribute_surface_temp");
	    	    if (!redistribute_surface_temp(*s,fr))
			status = NO;
                    stop_clock("redistribute_surface_temp");
                }
	    }
	}
	
	start_clock("Add_and_Delete_Points");
	if(is_use_bd_dist())
	{

	    printf("#extend intfc in redist.\n");
	    sav_c_gr = Computational_grid(intfc);
	    sav_t_gr = topological_grid(intfc);
	    
	    //use a large grid, scatter_front will construct a large 
	    //buffer for intfc.
	    change_buffer_for_intfc(intfc);
	   
	    //after changing the grid, do not comm component because component3d 
	    //will be called and make_tri_lists is needed.
	    //set_default_comp(YES);
	    start_clock("scatter_front");
	    status = scatter_front(fr);
	    stop_clock("scatter_front");
	    //set_default_comp(NO);
	    
	    sav_fr_gr = fr->rect_grid;
	    fr->rect_grid = computational_grid(intfc);
   
	}
	
	if (redist_non_vec_cur == YES)
	{
	    //printf("#redist curve3d bf\n");
	    //for (c = fr->interf->curves; c && *c; ++c)
	    //{
		//if (hsbdry_type(*c) != SUBDOMAIN_HSBDRY)
		//redistribute_curve3d(*c,fr);
	    //}

	    //DEBUG_TMP printf("#redist surface bf\n");
	    for (s = fr->interf->surfaces; s && *s ; ++s) 
	    {
	    	if ((!omit_redistribution(*s)) &&
	    	    (wave_type(*s) >= FIRST_PHYSICS_WAVE_TYPE) &&
	    	    (wave_type(*s) < FIRST_VECTOR_PHYSICS_WAVE_TYPE))
		{
                    printf("redist_non_vec_cur == YES\n");

		    start_clock("redistribute_surface_temp");
		    if (!redistribute_surface_temp(*s,fr))
			status = NO;
		    stop_clock("redistribute_surface_temp");
		}
	    }

	    delete_flag = YES;
	    while(delete_flag)
	    {
		delete_flag = NO;
		for (s = fr->interf->surfaces; s && *s; ++s)
		    if ((*s)->num_tri == 0)
		    {
			printf("#delete af redis.\n");
			delete_surface(*s);
			delete_flag = YES;
			break;
		    }
	    }
	}
	stop_clock("Add_and_Delete_Points");

        boolean status_is_bad_tri;
        status_is_bad_tri = pp_max_status(fr->is_bad_tri);
	
        if(status_is_bad_tri == YES)
            fr->is_bad_tri = YES;

	//insert_point_in_tri_side fails or delete_min_side_of_tri fails.
	if(status == NO)
	{
	    printf("WARNING after #surface_redistribute_temp, "
	    	   "status is NO in step %d\n", fr->step);
	}
	min_status = pp_min_status(status);
	if (min_status == NO)
	{
	    (void) printf("WARNING in surface_redistribute()_temp, "
		          "redistribute_surface(), failed\n");
	    DEBUG_LEAVE(surface_redistribute)
	    return min_status;
	}

	if(is_use_bd_dist())
	{
	    //recover the previous grid.
	    fr->rect_grid = sav_fr_gr;
	    Computational_grid(intfc) = sav_c_gr;
	    topological_grid(intfc) = sav_t_gr;
	
	    Computational_grid(fr->interf) = sav_c_gr;
	    topological_grid(fr->interf) = sav_t_gr;
	}
	else
	{
	    //DEBUG_TMP check_print_intfc("After redistribute_surface", "redist_sf_af", 'g', 
		//DEBUG_TMP	fr->interf, fr->step, -1, NO);
	    //DEBUG_TMP make_interface_topology_lists(fr->interf);
	    //DEBUG_TMP check_valid_intfc("after redistribute_surface", fr->interf);
	}

	start_clock("scatter_front");
	status = scatter_front(fr);
	stop_clock("scatter_front");
	
	if (!status)
	{
	    set_use_bd_dist(YES);
	    (void) printf("WARNING in surface_redistribute_temp(), "
	    	          "scatter_front() failed\n");
	    DEBUG_LEAVE(surface_redistribute)
	    return status;
	}

	set_current_interface(fr->interf);

	Interface_redistributed(fr) =
	    ((redist_vec_cur==YES) || (redist_non_vec_cur==YES)) ? YES : NO;

	DEBUG_LEAVE(surface_redistribute_temp)
	return status;
}		/*end surface_redistribute_temp*/

EXPORT boolean surface_redistribute(
	Front		*fr,
	boolean		*force_redist)
{
	SURFACE		**s;
	CURVE		**c;
	boolean		status, min_status, delete_flag;
	boolean		redist_non_vec_cur, redist_vec_cur;
	INTERFACE	*intfc = fr->interf;
	RECT_GRID	ngr, sav_c_gr, sav_t_gr, *sav_fr_gr;

	DEBUG_ENTER(surface_redistribute)

		// Check on redistribution conditions

	if (*force_redist)
	{
	    redist_non_vec_cur = redist_vec_cur = YES;
	    *force_redist = NO;
	}
	else if (Redistribution_count(fr) < 0)
	{
	    redist_non_vec_cur = redist_vec_cur = NO;
	}
	else
	{
	    redist_vec_cur = redist_needed(fr,VECTOR_WAVE);
	    redist_non_vec_cur = redist_needed(fr,GENERAL_WAVE);
	}
	++Redistribution_count(fr);

		// Redistribute vector surfaces
	
	set_size_of_intfc_state(size_of_state(fr->interf));
	set_copy_intfc_states(YES);
	set_current_interface(fr->interf);


	status = YES;
	if (redist_vec_cur == YES)
	{
	    for (s = fr->interf->surfaces; s && *s ; ++s) 
	    {
	    	if ((!omit_redistribution(*s)) &&
		    (wave_type(*s) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE))
	    	    if (!redistribute_surface(*s,fr))
			status = NO;
	    }
	}
	
	start_clock("Add_and_Delete_Points");
	if(is_use_bd_dist())
	{

	    printf("#extend intfc in redist.\n");
	    sav_c_gr = Computational_grid(intfc);
	    sav_t_gr = topological_grid(intfc);
	    
	    //use a large grid, scatter_front will construct a large 
	    //buffer for intfc.
	    change_buffer_for_intfc(intfc);
	   
	    //after changing the grid, do not comm component because component3d 
	    //will be called and make_tri_lists is needed.
	    //set_default_comp(YES);
	    start_clock("scatter_front");
	    status = scatter_front(fr);
	    stop_clock("scatter_front");
	    //set_default_comp(NO);
	    
	    sav_fr_gr = fr->rect_grid;
	    fr->rect_grid = computational_grid(intfc);
   
	}
	
	if (redist_non_vec_cur == YES)
	{
	    //printf("#redist curve3d bf\n");
	    //for (c = fr->interf->curves; c && *c; ++c)
	    //{
		//if (hsbdry_type(*c) != SUBDOMAIN_HSBDRY)
		//redistribute_curve3d(*c,fr);
	    //}

	    //DEBUG_TMP printf("#redist surface bf\n");
	    for (s = fr->interf->surfaces; s && *s ; ++s) 
	    {
	    	if ((!omit_redistribution(*s)) &&
	    	    (wave_type(*s) >= FIRST_PHYSICS_WAVE_TYPE) &&
	    	    (wave_type(*s) < FIRST_VECTOR_PHYSICS_WAVE_TYPE))
		{
		    start_clock("redistribute_surface");
		    if (!redistribute_surface(*s,fr))
			status = NO;
		    stop_clock("redistribute_surface");
		}
	    }

	    delete_flag = YES;
	    while(delete_flag)
	    {
		delete_flag = NO;
		for (s = fr->interf->surfaces; s && *s; ++s)
		    if ((*s)->num_tri == 0)
		    {
			printf("#delete af redis.\n");
			delete_surface(*s);
			delete_flag = YES;
			break;
		    }
	    }
	}
	stop_clock("Add_and_Delete_Points");

        boolean status_is_bad_tri;
        status_is_bad_tri = pp_max_status(fr->is_bad_tri);
	
        if(status_is_bad_tri == YES)
            fr->is_bad_tri = YES;

	//insert_point_in_tri_side fails or delete_min_side_of_tri fails.
	if(status == NO)
	{
	    printf("WARNING after #surface_redistribute, "
	    	   "status is NO in step %d\n", fr->step);
	}
	min_status = pp_min_status(status);
	if (min_status == NO)
	{
	    (void) printf("WARNING in surface_redistribute(), "
		          "redistribute_surface(), failed\n");
	    DEBUG_LEAVE(surface_redistribute)
	    return min_status;
	}

	if(is_use_bd_dist())
	{
	    //recover the previous grid.
	    fr->rect_grid = sav_fr_gr;
	    Computational_grid(intfc) = sav_c_gr;
	    topological_grid(intfc) = sav_t_gr;
	
	    Computational_grid(fr->interf) = sav_c_gr;
	    topological_grid(fr->interf) = sav_t_gr;
	}
	else
	{
	    //DEBUG_TMP check_print_intfc("After redistribute_surface", "redist_sf_af", 'g', 
		//DEBUG_TMP	fr->interf, fr->step, -1, NO);
	    //DEBUG_TMP make_interface_topology_lists(fr->interf);
	    //DEBUG_TMP check_valid_intfc("after redistribute_surface", fr->interf);
	}

	start_clock("scatter_front");
	status = scatter_front(fr);
	stop_clock("scatter_front");
	
	if (!status)
	{
	    set_use_bd_dist(YES);
	    (void) printf("WARNING in surface_redistribute(), "
	    	          "scatter_front() failed\n");
	    DEBUG_LEAVE(surface_redistribute)
	    return status;
	}

	set_current_interface(fr->interf);

	Interface_redistributed(fr) =
	    ((redist_vec_cur==YES) || (redist_non_vec_cur==YES)) ? YES : NO;

	DEBUG_LEAVE(surface_redistribute)
	return status;
}		/*end surface_redistribute*/

LOCAL	const int	Num_pqs_in_block = 1000;

struct _TRI_SURF {
	TRI	*tri;
	SURFACE *surf;
	CURVE   *c01, *c12, *c20;
	double   sqr_norm, dist, bd_dist;
	int     side;
};
typedef struct _TRI_SURF	TRI_SURF; 

#define Tri_surf(p)			((TRI_SURF *) (p)->pointer)
#define PQ_for_tri(tri)			((POINTER_Q *) Tri_workspace(tri))
#define tri_surface_from_queue(tri)	(Tri_surf(PQ_for_tri(tri)))
#define Bond_of_q(pq)		        ((BOND *)(pq)->pointer)
#define Tri_of_q(pq)		        (Tri_surf(pq)->tri)


LOCAL	double	max_sqr_area;	/* maximum triangle area	*/
LOCAL   double	min_sqr_area;	/* minimum triangle area	*/
LOCAL	double	max_sqr_length; /* maximun triangle side length */
LOCAL	double	aspect_tol2;	/* square of aspect ratio tolerance */
LOCAL   double  min_sqr_norm;   //min sqr area of a tri

double   Min_sqr_norm(RECT_GRID *gr)
{
	double  *h = gr->h;
        double	hmin = min3(h[0], h[1], h[2]);

	//DEBUG_TMP printf("min_sqr_norm = %24.16e\n", hmin*hmin*1.0e-8);

	return  hmin*hmin*1.0e-8;
	//1e-4 area of min block
}

double   reset_min_sqr_norm(double  new_min)
{
	if(new_min < min_sqr_norm)
	    min_sqr_norm = new_min;
}

LOCAL POINTER_Q *dequeue(
	TRI       *tri,
	POINTER_Q *pq)
{
	if (PQ_for_tri(tri))
	{
	    if (head_of_pointer_queue(PQ_for_tri(tri))
				    == head_of_pointer_queue(pq))
	    {
	        pq = delete_from_pointer_queue(PQ_for_tri(tri));
	        Tri_workspace(tri) = NULL;
	    }
	}
	return pq;
}		/*end dequeue*/

LOCAL POINTER_Q *alloc_and_add_to_queue(
	TRI       *t,
	SURFACE   *s,
	POINTER_Q *pq,
	int       nside,
	RECT_GRID *gr)
{
	TRI_SURF    *ts;
	const double *nor;
	int	    i, j;
	double	    cen[3];
	double	    coef[3]={1.414, 0.618, 1.323};

	pq = add_to_pointer_queue(NULL,pq);
	Tri_workspace(t) = (POINTER) pq;
	ts = tri_surface_from_queue(t);
	if (ts == NULL)
	{
	    screen("ERROR in alloc_and_add_to_queue(), "
	           "tri_surface_from_queue() returns NULL\n");
	    clean_up(ERROR);
	}
	ts->tri = t;
	nor = Tri_normal(t);
	ts->sqr_norm = Dot3d(nor,nor);

	centroid_of_tri(cen, t);
	// a distant for further check, should be shift invariant
	ts->dist = 0.0;
	ts->bd_dist = HUGE_VAL;
	for(i=0; i<3; i++)
	{
	    for(j=0; j<3; j++)
	        ts->dist += coef[j]*Coords(Point_of_tri(t)[i])[j];
	    
	    ts->bd_dist = min3(ts->bd_dist, 
			    fabs(cen[i] - gr->L[i]), fabs(cen[i] - gr->U[i]));
	}


	ts->surf = s;
	ts->c01 = ts->c12 = ts->c20 = NULL;
	ts->side = nside;
	return pq;
}		/*end alloc_and_add_to_queue*/

LOCAL  void  tecplot_tri_queue(
        const char      *msg,
        FILE            *file,
        POINTER_Q       *p_q)
{
        POINTER_Q       *q;
        TRI_SURF        *t_surf;
        POINT           *p;
        TRI             *tri;
        int             k, i, cnt = 0;

        double  pt[3] = { 0.9733961893900565,     -0.301512040795133,     -15.64756169180028  };
        double  pt1[3] = { -1.026603810609944,     -0.301512040795133,     -15.64756169180028 };


        if (p_q == NULL)
        {
            (void) printf("tecplot_tri_queue NULL POINTER_Q %s\n", msg);
            return;
        }

        q = head_of_pointer_queue(p_q);
        while (q != tail_of_pointer_queue(p_q))
        {
            t_surf = Tri_surf(q);
            tri = t_surf->tri;
            q = q->next;
            cnt++;
        }
        (void) fprintf(file, "ZONE T=\"%s\" N=%d E=%d\nF=FEPOINT, ET=TRIANGLE\n",
                             msg, 3*cnt, cnt);

        k = 0;
        q = head_of_pointer_queue(p_q);
        while (q != tail_of_pointer_queue(p_q))
        {
            t_surf = Tri_surf(q);
            tri = t_surf->tri;
            for (i = 0; i < 3; i++)
            {
                p = Point_of_tri(tri)[i];
                fprintf(file,"%-9g %-9g %-9g\n",Coords(p)[0],
                                 Coords(p)[1],Coords(p)[2]);

                if(distance_between_positions(Coords(Point_of_tri(tri)[i]), pt, 3)<0.05  ||
                   distance_between_positions(Coords(Point_of_tri(tri)[i]), pt1, 3)<0.05 )
                {
                        //printf("#de tri sort  k = %d\n", k);
                        //printf("sqr_norm = %24.16e, dist = %24.16e\n", t_surf->sqr_norm, t_surf->dist);
                        //print_tri(tri, t_surf->surf->interface);
                }
            }
            k++;
            q = q->next;
        }

        for(i=0; i<cnt; i++)
        {
            fprintf(file, "%d %d %d\n", 3*i+1, 3*i+2, 3*i+3);
        }
}

LOCAL  void  tri_queue_test(
	const char	*msg,
	POINTER_Q 	*p_q)
{
	POINTER_Q 	*q;
	TRI_SURF 	*t_surf;
	POINT		*p;
	TRI 		*tri;
	int		k, i;
	double		tst_pt[3] = {  0.0,     0.8985,                 1.47066};
	double		tst_pt1[3] = { 1.0,     0.8985,                 1.47066};
	double		tol = 1.0/40.0;
 
	if (p_q == NULL) 
	{
	    (void) printf("tri_queue_test NULL POINTER_Q %s\n", msg);
	    return;
        }
	
	k = 0;
	q = head_of_pointer_queue(p_q);
	while (q != tail_of_pointer_queue(p_q))
        {
	    t_surf = Tri_surf(q);
	    tri = t_surf->tri;

	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		
		if(distance_between_positions(Coords(Point_of_tri(tri)[i]), tst_pt, 3)<tol  || 
		   distance_between_positions(Coords(Point_of_tri(tri)[i]), tst_pt1, 3)<tol )
		{
		    printf("#de tri sort  k = %d\n", k);
		    printf("sqr_norm = %24.15e, dist = %24.15e\n", t_surf->sqr_norm, t_surf->dist);
		    print_tri(tri, t_surf->surf->interface);
		    
		    break;
		}
	    }

	    k++;
	    q = q->next;
	}
}

boolean	is_use_bd_dist();
boolean	is_use_bd_dist()
{
	return use_bd_dist;
}

LOCAL   boolean    use_modi_max_sqr_length = NO;
EXPORT  void    modi_max_sqr_length(boolean flag)
{
        use_modi_max_sqr_length = flag;
}

LOCAL  boolean redistribute_surface_init(
	SURFACE		*s,
	Front		*fr)
{
	INTERFACE *intfc;
	POINT	  *midp;
	POINTER_Q *insert_queue, *delete_queue;
	RECT_GRID *gr;
	TRI	  *tri, *oppt;
	boolean      status;
	int       nside, nf, nl, ns, nt, i;
	int	  dim, wc;
	FILE	  *db_file;
	double	  coords[3];

	DEBUG_ENTER(redistribute_surface_init)

        printf("In redistribute_surface_init(), use_bd_dist_smooth = %d\n", use_bd_dist_smooth);

	set_pointer_queue_opts(PQ_BLOCK_SIZE,Num_pqs_in_block,PQ_ALLOC_TYPE,
		               "vmalloc",PQ_ALLOC_SIZE_FOR_POINTERS,
			       sizeof(TRI_SURF),0);

	intfc = s->interface;
	dim = intfc->dim;
	gr = fr->rect_grid;
	
	//set the tolerance for tri_status
	wc = wave_type(s)< FIRST_VECTOR_PHYSICS_WAVE_TYPE ?
				GENERAL_WAVE : VECTOR_WAVE;
	max_sqr_area = Max_tri_sqr_area(fr,wc);
	min_sqr_area = Min_tri_sqr_area(fr,wc);
	max_sqr_length = Max_scaled_tri_side_sqr_length(fr);
	aspect_tol2 = sqr(Aspect_ratio_tolerance(fr,wc));
	min_sqr_norm = Min_sqr_norm(gr);
	    
	status = YES;
	insert_queue = delete_queue = NULL;
	nt = nf = nl = ns = 0;
	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{

	    tri->local_index = nt;

	    double  *p0 = Coords(Point_of_tri(tri)[0]);
	    double  *p1 = Coords(Point_of_tri(tri)[1]);
	    double  *p2 = Coords(Point_of_tri(tri)[2]);	

	    nt++;
	    switch (tri_status(tri,gr))
	    {
	    case BAD_ANGLE:
		//printf("\nBAD ANGLE tri index = %d\n", nt-1);	
		//printf(" Point 1: (%10.8g, %10.8g, %10.8g)\n", p0[0], p0[1], p0[2]);
		//printf(" Point 2: (%10.8g, %10.8g, %10.8g)\n", p1[0], p1[1], p1[2]);
		//printf(" Point 3: (%10.8g, %10.8g, %10.8g)\n", p2[0], p2[1], p2[2]);
		insert_queue = alloc_and_add_to_queue(tri,s,insert_queue,nside,gr);
                fr->is_bad_tri = YES;
		++nf;
		break;
	    case LARGE:
		//printf("\nLARGE TRIANGLE tri index = %d\n", nt-1);	
		//printf(" Point 1: (%10.8g, %10.8g, %10.8g)\n", p0[0], p0[1], p0[2]);
		//printf(" Point 2: (%10.8g, %10.8g, %10.8g)\n", p1[0], p1[1], p1[2]);
		//printf(" Point 3: (%10.8g, %10.8g, %10.8g)\n", p2[0], p2[1], p2[2]);
		++nl;
		insert_queue = alloc_and_add_to_queue(tri,s,insert_queue,nside,gr);
                fr->is_bad_tri = YES;
		break;
                /*
	    case SMALL:
		//printf("\nSMALL TRIANGLE tri index = %d\n", nt-1);	
		//printf(" Point 1: (%10.8g, %10.8g, %10.8g)\n", p0[0], p0[1], p0[2]);
		//printf(" Point 2: (%10.8g, %10.8g, %10.8g)\n", p1[0], p1[1], p1[2]);
		//printf(" Point 3: (%10.8g, %10.8g, %10.8g)\n", p2[0], p2[1], p2[2]);
		++ns;
		delete_queue = alloc_and_add_to_queue(tri,s,delete_queue,nside,gr);
                fr->is_bad_tri = YES;
		break;
                */
	    case GOOD_ANGLE:
	    default:
		Tri_workspace(tri) = NULL;
		break;
	    }
	}
        //printf("\nNumber of large tris in redistribution : %d %d %d %d\n", nt,nf,nl,ns);

	if(insert_queue==NULL && delete_queue==NULL)
	{
	    DEBUG_ENTER(redistribute_surface_init)
	    return YES;
	}

	if(debugging("prt_que"))
	{
	    char  ch[60];

	    sprintf(ch, "queue_tri%d_%d.plt", fr->step, pp_mynode());
	    db_file = fopen(ch, "w");
	    (void) fprintf(db_file,"TITLE = \"tecplot surface\"\n"
		   	"VARIABLES = \"x\", \"y\", \"z\"\n");
	    printf("#queue_tri open \n");

	    tecplot_tri_queue("insert", db_file, insert_queue);
	    tecplot_tri_queue("delete", db_file, delete_queue);
	    fclose(db_file);
	}

	if(use_bd_dist)
	    sort_pointer_queue(insert_queue,intfc,BD_DIST);
	else
	    sort_pointer_queue(insert_queue,intfc,LONGEST);
	if(debugging("pt_surface"))
	{
	    tecplot_surface_in_ball("in_bf", s);

	}

	while (insert_queue)
	{
	    insert_queue = head_of_pointer_queue(insert_queue);
	    tri = Tri_of_q(insert_queue);
	    
	    nside = find_scaled_extrem_edge(tri,gr,LONGEST);
	    insert_queue = dequeue(tri,insert_queue);
	    
	    //if (is_side_bdry(tri,nside))
		//continue;
	    
	    oppt = Tri_on_side(tri,nside);	    
	    if(is_tri_in_queue(oppt,insert_queue))
		insert_queue = dequeue(oppt,insert_queue);
	    if(is_tri_in_queue(oppt,delete_queue))
		delete_queue = dequeue(oppt,delete_queue);

	    //if(skip_bdry_tri(oppt) || skip_bdry_tri(tri))
		//continue;
            
	    // find and make tri side mid point
	    for(i = 0; i < dim; ++i)
		coords[i] = 0.5*(Coords(Point_of_tri(tri)[nside])[i] +
		                 Coords(Point_of_tri(tri)[Next_m3(nside)])[i]);
	    midp = Point(coords);
            //printf("\nNumber of large tris in redistribution count : %e coords %e %e %e\n", count,Coords(Point_of_tri(tri)[nside])[0],Coords(Point_of_tri(tri)[nside])[1],Coords(Point_of_tri(tri)[nside])[2]);

	    if (!insert_point_in_tri_side(midp,nside,tri,s))
	    {
		printf("WARNING redistribute_surface, "
		       "insert_point_in_tri_side fails.\n");
		status = NO;
	    }
	}
        /*
	nt = nf = nl = ns = 0;
	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{

	    tri->local_index = nt;

	    double  *p0 = Coords(Point_of_tri(tri)[0]);
	    double  *p1 = Coords(Point_of_tri(tri)[1]);
	    double  *p2 = Coords(Point_of_tri(tri)[2]);	

	    nt++;
	    switch (tri_status(tri,gr))
	    {
	    case LARGE:
		//printf("\nLARGE TRIANGLE tri index = %d\n", nt-1);	
		//printf(" Point 1: (%10.8g, %10.8g, %10.8g)\n", p0[0], p0[1], p0[2]);
		//printf(" Point 2: (%10.8g, %10.8g, %10.8g)\n", p1[0], p1[1], p1[2]);
		//printf(" Point 3: (%10.8g, %10.8g, %10.8g)\n", p2[0], p2[1], p2[2]);
		++nl;
		insert_queue = alloc_and_add_to_queue(tri,s,insert_queue,nside,gr);
                fr->is_bad_tri = YES;
		break;
	    case GOOD_ANGLE:
	    default:
		Tri_workspace(tri) = NULL;
		break;
	    }
	}
	
        printf("\nafter Number of large tris in redistribution : %d %d %d %d\n", nt,nf,nl,ns);
        */
	    
	if(use_bd_dist)
	    sort_pointer_queue(delete_queue,intfc,BD_DIST);
	else
            //sort_pointer_queue(delete_queue,intfc,LONGEST);
	    sort_pointer_queue(delete_queue,intfc,SHORTEST);
	if(debugging("pt_surface"))
	{
	    tecplot_surface_in_ball("de_bf", s);
	    tri_queue_test("delete_bf_test",  delete_queue);
	}

	while (delete_queue)
	{
	    delete_queue = head_of_pointer_queue(delete_queue);
	    tri = Tri_of_q(delete_queue);

	    nside = find_scaled_extrem_edge(tri,gr,SHORTEST);
		
	    delete_queue = dequeue(tri, delete_queue);

	    if(!delete_min_side_of_tri(tri,nside,s,&delete_queue,fr))
	    {
		printf("WARNING, redistribute_surface_init, "
		       "delete_min_side_of_tri fails.\n");
		status = NO;
	    }

	}


	DEBUG_LEAVE(redistribute_surface_init)
	return status;
}		/*end redistribute_surface_init*/

LOCAL void fprint_tri_temp(
        FILE      *file,
        TRI       *tri,
        INTERFACE *intfc)
{
        POINT              *p[3];
        int                i;

        if (tri == NULL)
        {
            (void) fprintf(file,"NULL Triangle\n");
            return;
        }
        (void) fprintf(file,"Tri = %llu\n",tri_number(tri,intfc));
        for (i = 0; i < 3; ++i)
        {
            p[i] = Point_of_tri(tri)[i];
            (void) fprintf(file,"p%d: ",i);
            if (p[i] != NULL)
            {
                (void) fprintf(file,"( %"FFMT" %"FFMT" %"FFMT" ) ",
                           Coords(p[i])[0],Coords(p[i])[1],Coords(p[i])[2]);
                (void) fprintf(file,"%llu Boundary %d Boundary_point %d\n",
                               point_number(p[i]),Boundary(p[i]),
                               Boundary_point(p[i]));
            }
            else
                (void) fprintf(file,"NULL ");
        }
}

LOCAL isSameTri(TRI *s, TRI *t)
{
	if (s == NULL && t == NULL)	return YES;
	if (s == NULL || t == NULL)	return NO;
	int i;
	for (i=0; i<3; ++i)
	{
	    double diff = Point_of_tri(s)[i] - Point_of_tri(t)[i];
	    if (fabs(diff) > 1e-12)	return NO;
	}
	return YES;
}

LOCAL  boolean redistribute_surface_temp(
	SURFACE		*s,
	Front		*fr)
{
	INTERFACE *intfc;
	POINT	  *midp;
	POINTER_Q *insert_queue, *delete_queue;
	RECT_GRID *gr;
	TRI	  *tri, *oppt;
	boolean      status;
	int       nside, nf, nl, ns, nt, i;
	int	  dim, wc;
	FILE	  *db_file;
	double	  coords[3];

	DEBUG_ENTER(redistribute_surface_temp)

        printf("In redistribute_surface_temp(), use_bd_dist_smooth = %d\n", use_bd_dist_smooth);

	set_pointer_queue_opts(PQ_BLOCK_SIZE,Num_pqs_in_block,PQ_ALLOC_TYPE,
		               "vmalloc",PQ_ALLOC_SIZE_FOR_POINTERS,
			       sizeof(TRI_SURF),0);

	intfc = s->interface;
	dim = intfc->dim;
	gr = fr->rect_grid;


        if(debugging("prt_small_tri_que") && (pp_mynode() == 950 || pp_mynode() == 951) )
        {
            char  ch[256];
            sprintf(ch, "init_tri_redistribute_surface_temp%d_%d.txt", fr->step, pp_mynode());
            db_file = fopen(ch, "w");
            for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
                fprint_tri_temp(db_file, tri, intfc);
            fclose(db_file);
        }



	//set the tolerance for tri_status
	wc = wave_type(s)< FIRST_VECTOR_PHYSICS_WAVE_TYPE ?
				GENERAL_WAVE : VECTOR_WAVE;
	max_sqr_area = Max_tri_sqr_area(fr,wc);
	min_sqr_area = Min_tri_sqr_area(fr,wc);
	max_sqr_length = Max_scaled_tri_side_sqr_length(fr);
	aspect_tol2 = sqr(Aspect_ratio_tolerance(fr,wc));
	min_sqr_norm = Min_sqr_norm(gr);
	    
	status = YES;
	insert_queue = delete_queue = NULL;
	nt = nf = nl = ns = 0;
	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{

	    tri->local_index = nt;

	    double  *p0 = Coords(Point_of_tri(tri)[0]);
	    double  *p1 = Coords(Point_of_tri(tri)[1]);
	    double  *p2 = Coords(Point_of_tri(tri)[2]);	

	    nt++;
	    switch (tri_status(tri,gr))
	    {
		
		/*
	    case BAD_ANGLE:
		printf("\nBAD ANGLE tri index = %d\n", nt-1);	
		printf(" Point 1: (%10.8g, %10.8g, %10.8g)\n", p0[0], p0[1], p0[2]);
		printf(" Point 2: (%10.8g, %10.8g, %10.8g)\n", p1[0], p1[1], p1[2]);
		printf(" Point 3: (%10.8g, %10.8g, %10.8g)\n", p2[0], p2[1], p2[2]);
		delete_queue = alloc_and_add_to_queue(tri,s,delete_queue,nside,gr);
                fr->is_bad_tri = YES;
		++nf;
		break;
	    case LARGE:
		printf("\nLARGE TRIANGLE tri index = %d\n", nt-1);	
		printf(" Point 1: (%10.8g, %10.8g, %10.8g)\n", p0[0], p0[1], p0[2]);
		printf(" Point 2: (%10.8g, %10.8g, %10.8g)\n", p1[0], p1[1], p1[2]);
		printf(" Point 3: (%10.8g, %10.8g, %10.8g)\n", p2[0], p2[1], p2[2]);
		++nl;
		insert_queue = alloc_and_add_to_queue(tri,s,insert_queue,nside,gr);
                fr->is_bad_tri = YES;
		break;
	
		*/
	    case SMALL:
		/*
		printf("\nSMALL TRIANGLE tri index = %d\n", nt-1);	
		printf(" Point 1: (%10.8g, %10.8g, %10.8g)\n", p0[0], p0[1], p0[2]);
		printf(" Point 2: (%10.8g, %10.8g, %10.8g)\n", p1[0], p1[1], p1[2]);
		printf(" Point 3: (%10.8g, %10.8g, %10.8g)\n", p2[0], p2[1], p2[2]);
		*/
		++ns;
		delete_queue = alloc_and_add_to_queue(tri,s,delete_queue,nside,gr);
                fr->is_bad_tri = YES;


        	if(debugging("prt_small_tri_que") && (pp_mynode() == 950 || pp_mynode() == 951) ) {
            	    char  ch[256];
            	    sprintf(ch, "delete_queue_small_tri%d_%d.txt", fr->step, pp_mynode());
            	    db_file = fopen(ch, "w");
                    fprint_tri_temp(db_file, tri, intfc);
            	    fclose(db_file);
        	}


		break;
	    case GOOD_ANGLE:
	    default:
		Tri_workspace(tri) = NULL;
		break;
	    }
	}
        printf("\nNumber of small tris in redistribution_temp: %d %d %d %d\n", nt,nf,nl,ns);

	if(insert_queue==NULL && delete_queue==NULL)
	{
	    DEBUG_ENTER(redistribute_surface_temp)
	    return YES;
	}

	if(debugging("prt_que"))
	{
	    char  ch[60];

	    sprintf(ch, "queue_tri%d_%d.plt", fr->step, pp_mynode());
	    db_file = fopen(ch, "w");
	    (void) fprintf(db_file,"TITLE = \"tecplot surface\"\n"
		   	"VARIABLES = \"x\", \"y\", \"z\"\n");
	    printf("#queue_tri open \n");

	    tecplot_tri_queue("insert", db_file, insert_queue);
	    tecplot_tri_queue("delete", db_file, delete_queue);
	    fclose(db_file);
	}


	if(use_bd_dist)
	    sort_pointer_queue(insert_queue,intfc,BD_DIST);
	else
	    sort_pointer_queue(insert_queue,intfc,LONGEST);
	if(debugging("pt_surface"))
	{
	    tecplot_surface_in_ball("in_bf", s);
	}

	while (insert_queue)
	{
	    insert_queue = head_of_pointer_queue(insert_queue);
	    tri = Tri_of_q(insert_queue);
	    
	    nside = find_scaled_extrem_edge(tri,gr,LONGEST);
	    insert_queue = dequeue(tri,insert_queue);
	    
	    //if (is_side_bdry(tri,nside))
		//continue;
	    
	    oppt = Tri_on_side(tri,nside);	    
	    if(is_tri_in_queue(oppt,insert_queue))
		insert_queue = dequeue(oppt,insert_queue);
	    if(is_tri_in_queue(oppt,delete_queue))
		delete_queue = dequeue(oppt,delete_queue);

	    //if(skip_bdry_tri(oppt) || skip_bdry_tri(tri))
		//continue;
            
	    // find and make tri side mid point
	    for(i = 0; i < dim; ++i)
		coords[i] = 0.5*(Coords(Point_of_tri(tri)[nside])[i] +
		                 Coords(Point_of_tri(tri)[Next_m3(nside)])[i]);
	    midp = Point(coords);

	    if (!insert_point_in_tri_side(midp,nside,tri,s))
	    {
		printf("WARNING redistribute_surface_temp, "
		       "insert_point_in_tri_side fails.\n");
		status = NO;
	    }
	}
	    
	if(use_bd_dist)
	    sort_pointer_queue(delete_queue,intfc,BD_DIST);
	else
	    sort_pointer_queue(delete_queue,intfc,SHORTEST);
	if(debugging("pt_surface"))
	{
	    tecplot_surface_in_ball("de_bf", s);
	    tri_queue_test("deletei_bf_test",  delete_queue);
	}

	while (delete_queue)
	{
	    delete_queue = head_of_pointer_queue(delete_queue);
	    tri = Tri_of_q(delete_queue);

	    nside = find_scaled_extrem_edge(tri,gr,SHORTEST);
		
	    delete_queue = dequeue(tri, delete_queue);

	    if(!delete_min_side_of_tri(tri,nside,s,&delete_queue,fr))
	    {
		printf("WARNING, redistribute_surface_temp, "
		       "delete_min_side_of_tri fails.\n");
		status = NO;
	    }

	}




        if(debugging("prt_small_tri_que") && (pp_mynode() == 950 || pp_mynode() == 951) )
	{
            char  ch[256];
            sprintf(ch, "after_delete_small_tri%d_%d.txt", fr->step, pp_mynode());
            db_file = fopen(ch, "w");

            double  pt_p0951_0[3] = { 24.0125502876005541,      1.92332475341791631,       10.096417590337051 };
            for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	    {
                for (i = 0; i < 3; i++)
                {
                    if(distance_between_positions(Coords(Point_of_tri(tri)[i]), pt_p0951_0, 3)<1e-2)
                    {
                        fprint_tri_temp(db_file, tri, intfc);
                        break;
                    }
                }
	    }
            fclose(db_file);
	}





	DEBUG_LEAVE(redistribute_surface_temp)
	return status;
}		/*end redistribute_surface_temp*/

LOCAL  boolean redistribute_surface(
	SURFACE		*s,
	Front		*fr)
{
	INTERFACE *intfc;
	POINT	  *midp;
	POINTER_Q *insert_queue, *delete_queue;
	RECT_GRID *gr;
	TRI	  *tri, *oppt;
	boolean      status;
	int       nside, nf, nl, ns, nt, i;
	int	  dim, wc;
	FILE	  *db_file;
	double	  coords[3];

	DEBUG_ENTER(redistribute_surface)

        printf("In redistribute_surface(), use_bd_dist_smooth = %d\n", use_bd_dist_smooth);

	set_pointer_queue_opts(PQ_BLOCK_SIZE,Num_pqs_in_block,PQ_ALLOC_TYPE,
		               "vmalloc",PQ_ALLOC_SIZE_FOR_POINTERS,
			       sizeof(TRI_SURF),0);

	intfc = s->interface;
	dim = intfc->dim;
	gr = fr->rect_grid;
	
	//set the tolerance for tri_status
	wc = wave_type(s)< FIRST_VECTOR_PHYSICS_WAVE_TYPE ?
				GENERAL_WAVE : VECTOR_WAVE;
	max_sqr_area = Max_tri_sqr_area(fr,wc);
	min_sqr_area = Min_tri_sqr_area(fr,wc);
	max_sqr_length = Max_scaled_tri_side_sqr_length(fr);
	aspect_tol2 = sqr(Aspect_ratio_tolerance(fr,wc));
	min_sqr_norm = Min_sqr_norm(gr);
	    
	status = YES;
	insert_queue = delete_queue = NULL;
	nt = nf = nl = ns = 0;
	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{

	    tri->local_index = nt;
	    nt++;
	    switch (tri_status(tri,gr))
	    {
	    case BAD_ANGLE:
		// delete_queue = alloc_and_add_to_queue(tri,s,delete_queue,nside,gr);
                insert_queue = alloc_and_add_to_queue(tri,s,insert_queue,nside,gr);
                fr->is_bad_tri = YES;
		++nf;
		break;
	    case LARGE:
		++nl;
		insert_queue = alloc_and_add_to_queue(tri,s,insert_queue,nside,gr);
                fr->is_bad_tri = YES;
		break;
            /*
	    case SMALL:
		++ns;
		delete_queue = alloc_and_add_to_queue(tri,s,delete_queue,nside,gr);
                fr->is_bad_tri = YES;
		break;
            */
	    case GOOD_ANGLE:
	    default:
		Tri_workspace(tri) = NULL;
		break;
	    }
	}
        printf("\nNumber of bad tris in redistribution for proc #%d ts=%d: %d %d %d %d\n", pp_mynode(),fr->step,nt,nf,nl,ns);

	if(insert_queue==NULL && delete_queue==NULL)
	{
	    DEBUG_ENTER(redistribute_surface)
	    return YES;
	}

	if(debugging("prt_que"))
	{
	    char  ch[60];

	    sprintf(ch, "queue_tri%d_%d.plt", fr->step, pp_mynode());
	    db_file = fopen(ch, "w");
	    (void) fprintf(db_file,"TITLE = \"tecplot surface\"\n"
		   	"VARIABLES = \"x\", \"y\", \"z\"\n");
	    printf("#queue_tri open \n");

	    tecplot_tri_queue("insert", db_file, insert_queue);
	    tecplot_tri_queue("delete", db_file, delete_queue);
	    fclose(db_file);
	}

	if(use_bd_dist)
	    sort_pointer_queue(insert_queue,intfc,BD_DIST);
	else
	    sort_pointer_queue(insert_queue,intfc,LONGEST);
	if(debugging("pt_surface"))
	{
	    tecplot_surface_in_ball("in_bf", s);
	}

	while (insert_queue)
	{
	    insert_queue = head_of_pointer_queue(insert_queue);
	    tri = Tri_of_q(insert_queue);
	    
	    nside = find_scaled_extrem_edge(tri,gr,LONGEST);
	    insert_queue = dequeue(tri,insert_queue);
	    
	    //if (is_side_bdry(tri,nside))
		//continue;
	    
	    oppt = Tri_on_side(tri,nside);	    
	    if(is_tri_in_queue(oppt,insert_queue))
		insert_queue = dequeue(oppt,insert_queue);
	    if(is_tri_in_queue(oppt,delete_queue))
		delete_queue = dequeue(oppt,delete_queue);

	    //if(skip_bdry_tri(oppt) || skip_bdry_tri(tri))
		//continue;
            
	    // find and make tri side mid point
	    for(i = 0; i < dim; ++i)
		coords[i] = 0.5*(Coords(Point_of_tri(tri)[nside])[i] +
		                 Coords(Point_of_tri(tri)[Next_m3(nside)])[i]);
	    midp = Point(coords);

	    if (!insert_point_in_tri_side(midp,nside,tri,s))
	    {
		printf("WARNING redistribute_surface, "
		       "insert_point_in_tri_side fails.\n");
		status = NO;
	    }
	}
	    
	if(use_bd_dist)
	    sort_pointer_queue(delete_queue,intfc,BD_DIST);
	else
            //sort_pointer_queue(delete_queue,intfc,LONGEST);
	    sort_pointer_queue(delete_queue,intfc,SHORTEST);
	if(debugging("pt_surface"))
	{
	    tecplot_surface_in_ball("de_bf", s);
	    tri_queue_test("delete_bf_test",  delete_queue);
	}

	while (delete_queue)
	{
	    delete_queue = head_of_pointer_queue(delete_queue);
	    tri = Tri_of_q(delete_queue);

	    nside = find_scaled_extrem_edge(tri,gr,SHORTEST);
		
	    delete_queue = dequeue(tri, delete_queue);

	    if(!delete_min_side_of_tri(tri,nside,s,&delete_queue,fr))
	    {
		printf("WARNING, redistribute_surface, "
		       "delete_min_side_of_tri fails.\n");
		status = NO;
	    }

	}

	
	DEBUG_LEAVE(redistribute_surface)
	return status;
}		/*end redistribute_surface*/

LOCAL  int find_scaled_extrem_edge(
	TRI		*tri,
	RECT_GRID	*grid,
	SPQ_FLAG	to_find)
{
	const double* const *s;
	double	h0 = grid->h[0], h1 = grid->h[1], h2 = grid->h[2];
	double	s00, s01, s02;
	double  s10, s11, s12;
	double	s20, s21, s22;
	double	len0, len1, len2;

	double h_min = min3(grid->h[0], grid->h[1], grid->h[2]);
	h0 = h_min;
	h1 = h_min;
	h2 = h_min;

	s = side_vector(tri);
	s00 = s[0][0]/h0; s01 = s[0][1]/h1; s02 = s[0][2]/h2;
	s10 = s[1][0]/h0; s11 = s[1][1]/h1; s12 = s[1][2]/h2;
	s20 = s[2][0]/h0; s21 = s[2][1]/h1; s22 = s[2][2]/h2;
	len0 = QDot3d(s0,s0); len1 = QDot3d(s1,s1); len2 = QDot3d(s2,s2);

	switch (to_find)
	{
	case LONGEST:
	    return (len0<len1) ? ((len1<len2) ? 2:1) : ((len0<len2) ? 2:0);
	case SHORTEST:
	    return (len0>len1) ? ((len1>len2) ? 2:1) : ((len0>len2) ? 2:0);
	default:
	    return -1;
	}
}		/*end find_scaled_extrem_edge*/


/*
*				tri_status():
*
*	Determines whether redistribution of a triangle is needed by
*	comparing the triangles normalized area with two tolerances.
*	In addition, this routine also checks the squared edge lengths and
*	determines their aspect ratios.
*	This routine has the following return values:
*
*	if (norm_area > max_sqr_area) return LARGE;
*	if (norm_area < min_sqr_area) return SMALL;
*	if (aspect ratio < aspect_tol) return BAD_ANGLE;
*	return GOOD_ANGLE;
*
*	The aspect ratio of a triangle is defined as A/(l0^2+l1^2+l2^2)
*	where A is the area of the triangle in the grid scaled metric
*	and the li are the lengths of the sides of the triangle in
*	grid scaled metric.
*/

//LOCAL  TRI_STATUS tri_status(
EXPORT	TRI_STATUS tri_status(
	TRI *tri,
	RECT_GRID *gr)
{
	double	*p0 = Coords(Point_of_tri(tri)[0]);
	double  *p1 = Coords(Point_of_tri(tri)[1]);
	double	*p2 = Coords(Point_of_tri(tri)[2]);
	double	s00, s01, s02;
	double	s10, s11, s12;
	double	s20, s21, s22;
	double	N0, N1, N2;
	double	h0 = gr->h[0], h1 = gr->h[1], h2 = gr->h[2];
	double	sqr_area;
	double	a2,len[3],tri_area;
	const double *nor = Tri_normal(tri);
        double  h_min = min3(gr->h[0], gr->h[1], gr->h[2]);
        h0 = h_min;
        h1 = h_min;
        h2 = h_min;

	sqr_area = 0.25*Dot3d(nor,nor); 


        if (use_bd_dist_smooth == YES)
        {
            if (Mag3d(nor) < h_min*h_min*1.0e-8)
                return SMALL;
        }
        else
        {
	    if (sqr_area <= min_sqr_area)
	        return SMALL;
        }


	s00 = (p1[0]-p0[0])/h0; s01 = (p1[1]-p0[1])/h1; s02 = (p1[2]-p0[2])/h2;
	s10 = (p2[0]-p1[0])/h0; s11 = (p2[1]-p1[1])/h1; s12 = (p2[2]-p1[2])/h2;
	s20 = (p0[0]-p2[0])/h0; s21 = (p0[1]-p2[1])/h1; s22 = (p0[2]-p2[2])/h2;
	QCross3d(s0,s2,N);

	tri_area = 0.5*sqrt(QDot3d(N,N));

			/* Check aspect ratio	*/
	len[0] = QDot3d(s0,s0);
	len[1] = QDot3d(s1,s1);
	len[2] = QDot3d(s2,s2);

	a2 = len[0]+len[1]+len[2];
	a2 = tri_area/sqr(a2);


			/* Check max/min angle	*/
        double	max_cos, min_cos;
	double	max_angle, min_angle;
	double	l[3], cos_angle[3];

	l[2] = sqrt( sqr(p1[0]-p0[0]) + sqr(p1[1]-p0[1]) + sqr(p1[2]-p0[2]) );
        l[0] = sqrt( sqr(p2[0]-p1[0]) + sqr(p2[1]-p1[1]) + sqr(p2[2]-p1[2]) );
        l[1] = sqrt( sqr(p0[0]-p2[0]) + sqr(p0[1]-p2[1]) + sqr(p0[2]-p2[2]) );

	cos_angle[0] = (sqr(l[1]) + sqr(l[2]) - sqr(l[0])) / (2*l[1]*l[2]);
        cos_angle[1] = (sqr(l[0]) + sqr(l[2]) - sqr(l[1])) / (2*l[0]*l[2]);
        cos_angle[2] = (sqr(l[0]) + sqr(l[1]) - sqr(l[2])) / (2*l[0]*l[1]);

	max_cos = max3(cos_angle[0], cos_angle[1], cos_angle[2]);
	max_cos = min(1, max_cos);
        min_cos = min3(cos_angle[0], cos_angle[1], cos_angle[2]);
	min_cos = max(-1, min_cos);

        max_angle = acos(min_cos) / 3.1415926535897931 * 180.0;
	min_angle = acos(max_cos) / 3.1415926535897931 * 180.0;


        if (sqr_area >= max_sqr_area)
            return LARGE;

	if ((len[0] > max_sqr_length) || (len[1] > max_sqr_length) ||
	    (len[2] > max_sqr_length))
	    return LARGE;


        if (max_angle > 165 || min_angle < 15) {
	    // printf("Detection: max_angle = %lf, min_angle = %lf\n", max_angle, min_angle);
            return BAD_ANGLE;
	}

        if (a2 < aspect_tol2)
            return BAD_ANGLE;

	return GOOD_ANGLE;
}		/*end tri_status*/

LOCAL  int status_order = 0;
LOCAL   void  set_tri_status_order(int order)
{
	status_order = order;
}

boolean    null_side_loop(TRI*,int,ORIENTATION,TRI***,int**,POINT***,
                                int*,double**);
void	seal_null_loop(TRI**,int*,POINT**,int);
LOCAL	boolean	is_critical_tri(TRI *, int);


//COND: there should be no boundary point in the bound of tris
EXPORT	int	remove_tris_and_seal(
	TRI		**new_tris,
	TRI		**tris,
	int		nt,
	SURFACE		*s,
	POINTER_Q	**pq,
	INTERFACE	*intfc)
{
	TRI	*out_tris[500], **new_out_tris;
	int	i, num_out_tris, num_new_tris;

	DEBUG_ENTER(remove_tris_and_seal);
	
	num_out_tris = bound_tris_set(out_tris, tris, nt);

        if (num_out_tris >= 500) {
            printf("Stack overflow: num_out_tris >= 500!\n");
            clean_up(ERROR);
        }

	//tris are already dequeue in delete_min_side_of_tri
	for(i=0; i<nt; i++)
	    remove_tri_from_surface(tris[i], s, NO);
	    
	if(num_out_tris == 0)
	{
	    DEBUG_LEAVE(remove_tris_and_seal);
	    return 0;
	}

	sep_common_point_from_loop(out_tris, num_out_tris, NULL, NULL, intfc);
	
	//new null side tris can be added into out_tris
	num_out_tris = sep_common_edge_from_tris(&new_out_tris, 
				out_tris, num_out_tris, intfc);

	//since a smooth_null_loop is applied above, the positions of 
	//3 vertics of a tri is changed, all the bound tris should be 
	//removed from the que.
	
	for(i=0; i<num_out_tris; i++)
	    *pq = dequeue(new_out_tris[i], *pq);

	num_new_tris = 0;
	nt = seal_all_loops_wo_constraint(new_tris, &num_new_tris, 
			new_out_tris, num_out_tris, 0, NO);


        if (num_new_tris >= 500) {
            printf("Stack overflow: num_new_tris >= 500!\n");
            clean_up(ERROR);
        }


	nt = merge_tris_set(new_tris, num_new_tris, new_out_tris, nt);

	DEBUG_LEAVE(remove_tris_and_seal);
	return nt;
}


LOCAL boolean delete_min_side_of_tri(
	TRI	  *tri,
	int	  side,
	SURFACE	  *s,
	POINTER_Q **pq,
	Front     *fr)
{
	INTERFACE	*intfc;
	TRI	*nbtri, *t, *nbt, **tmp_tris;
	TRI	*new_tris[500], *in_tris[200], *tris[2][100];
	POINT	*p[4], *pt, *pmid, *plist[2][100];
	int	i, j, k, nt, np[2], nside, ntris[2];
	boolean	rm_flag;
	static	int	cnt = 0;
	FILE	*file;
	char	fname[100];

	DEBUG_ENTER(delete_min_side_of_tri)

	//printf("#tri side %d %d\n", tri, side);

	intfc = fr->interf;
	p[0] = Point_of_tri(tri)[side];
	p[1] = Point_of_tri(tri)[Next_m3(side)];

        if(debugging("prt_small_tri_que") && (pp_mynode() == 950 || pp_mynode() == 951) ) {
	    printf("tri #%llu side %d, proc %d    ts %d\n", tri_number(tri,intfc),side,pp_mynode(),fr->step);
	}

	if(Boundary_point(p[0]) || Boundary_point(p[1]))
	{
	    DEBUG_LEAVE(delete_min_side_of_tri)
	    return YES;
	}

	nbtri = Tri_on_side(tri,side);
	for(nside=0; nside<3; nside++)
	    if (Tri_on_side(nbtri,nside) == tri)
		break;

        if(debugging("prt_small_tri_que") && (pp_mynode() == 950 || pp_mynode() == 951) ) {
           printf("nbtri #%llu nside %d, proc %d    ts %d\n", tri_number(nbtri,intfc),nside,pp_mynode(),fr->step);
        }

	p[2] = Point_of_tri(tri)[Prev_m3(side)];
	p[3] = Point_of_tri(nbtri)[Prev_m3(nside)];


        if(debugging("prt_small_tri_que") && (pp_mynode() == 950 || pp_mynode() == 951) ) { 
	    for (i = 0; i < 4; ++i)
            {
                (void) printf("p%d: ",i);
                if (p[i] != NULL)
                {
                    (void) printf("( %"FFMT" %"FFMT" %"FFMT" ) ",
                               Coords(p[i])[0],Coords(p[i])[1],Coords(p[i])[2]);
                    (void) printf("%llu\n", point_number(p[i]));
                }
                else
                    (void) printf("NULL\n");
            }
	}


	for(k=0; k<2; k++)
	{
	    ntris[k] = set_tri_list_around_point(p[k],tri,&tmp_tris,intfc);

	    if (ntris[k] >= 100) {
		printf("Stack overflow: ntris[%d] >= 100!\n", k);
		clean_up(ERROR);
	    }

	    for(i=0; i<ntris[k]; i++)
	    {
		tris[k][i] = tmp_tris[i];
	        *pq = dequeue(tris[k][i],*pq);
	    }

	    np[k] = 0;
	    //finding bounding points except the 4 common points.
	    for(i=0; i<ntris[k]; i++)
	    {
		t = tris[k][i];
		j = Vertex_of_point(t, p[k]);
		pt = Point_of_tri(t)[Prev_m3(j)];
		
		for(j=0; j<4; j++)
		    if(pt == p[j])
			break;
		if(j < 4)
		    continue;

		plist[k][np[k]] = pt;
		np[k]++;

        	if (np[k] >= 100) {
                    printf("Stack overflow: np[%d] >= 100!\n", k);
                    clean_up(ERROR);
                }

	    }
	}

	//skip the bdry case.
	if(Boundary_point(p[2]) || Boundary_point(p[3]))
	{
	    DEBUG_LEAVE(delete_min_side_of_tri)
	    return YES;
	}
	for(i=0; i<np[0]; i++)
	    if(Boundary_point(plist[0][i]))
	    {
		DEBUG_LEAVE(delete_min_side_of_tri)
		return YES;
	    }
	for(i=0; i<np[1]; i++)
	    if(Boundary_point(plist[1][i]))
	    {
		DEBUG_LEAVE(delete_min_side_of_tri)
		return YES;
	    }

	//check if there are duplicate points in the bounding tris.
	rm_flag = NO;
	if(np[0] > 0 && np[1] > 0)
	{
	    //the general case, test duplicate points
	    for(i=0; i<np[0]; i++)
		for(j=0; j<np[1]; j++)
		    if(plist[0][i] == plist[1][j])
			rm_flag = YES;
	}
	else if(np[0] == 0 && np[1] == 0)
	{
	    //the tetrahedron case
	    rm_flag = YES;
	}

	//printf("#shape np %3d  %3d %3d\n", rm_flag, np[0], np[1]);
	if(rm_flag)
	{
	    //make sure, after removing, no boundary point appears in the 
	    //null loop
	    //DEBUG_TMP printf("#shape np %3d  %3d %3d\n", rm_flag, np[0], np[1]);

	    nt = 0;
	    for(k=0; k<2; k++)
		nt = merge_tris_set(in_tris, nt, tris[k], ntris[k]);

            if (nt >= 200) {
                printf("Stack overflow: nt >= 200!\n", k);
                clean_up(ERROR);
            }

	    //DEBUG_TMP printf("#delete_min,  duplicate point, nt = %d\n", nt);	

	    if(debugging("delete_dup"))
	    {
		sprintf(fname,"dup_min%d_%d.plt",pp_mynode(),cnt);
		cnt++;
		printf("debug file %s\n", fname);
		
		file = fopen(fname,"w");
		tecplot_show_tris("in_tris", in_tris, nt, file);
		fclose(file);
	    }

	    nt = remove_tris_and_seal(new_tris, in_tris, nt, s, pq, intfc);

	    if(debugging("delete_dup"))
	    {
		file = fopen(fname,"a");
		tecplot_show_tris("new_tris", new_tris, nt, file);
		fclose(file);
	    }

	    //DEBUG_TMP printf("#rm_seal finish\n");

	    DEBUG_LEAVE(delete_min_side_of_tri)
	    return   nt==-1 ? NO : YES;
	}

	//collapse two tris.
	pmid = average_points(YES,p[0],Hyper_surf_element(tri),Hyper_surf(s),
				  p[1],Hyper_surf_element(tri),Hyper_surf(s));
	
	//change the point for the surrounding tris.
	for (i = 0; i < 2; ++i)
	{
	    for (j = 0; j < ntris[i]; ++j)
	    {
		t = tris[i][j];
		k = Vertex_of_point(t,p[i]);
		Point_of_tri(t)[k] = pmid;
		if ((t != tri) && (t != nbtri))
		    set_normal_of_tri(t);
	    }
	}

	//change tri neighbor for tri.
	nbt = Tri_on_side(tri,Next_m3(side));
	t = Tri_on_side(tri,Prev_m3(side));
	for(i=0; i<3; i++)
	{
	    if (Tri_on_side(t,i) == tri)
		Tri_on_side(t,i) = nbt;
	    if (Tri_on_side(nbt,i) == tri)
		Tri_on_side(nbt,i) = t;
	}

	//change tri neighbor for nbtri.
	nbt = Tri_on_side(nbtri,Next_m3(nside));
	t = Tri_on_side(nbtri,Prev_m3(nside));
	for (i = 0; i < 3; ++i)
	{
	    if (Tri_on_side(t,i) == nbtri)
		Tri_on_side(t,i) = nbt;
	    if (Tri_on_side(nbt,i) == nbtri)
		Tri_on_side(nbt,i) = t;
	}

	remove_tri_from_surface(tri,s,YES);
	remove_tri_from_surface(nbtri,s,YES);

	DEBUG_LEAVE(delete_min_side_of_tri)
	return YES;
}


/*
*			angle_weight_of_tri():
*
*	Returns 2*(1 + cos(theta)) where theta is the maximum internal
*	vertex of the triangle when mapped into the scaled metric space
*	x[i] -> x[i]/h[i].
*	
*	Note the formulas are based on the law of cosines which states
*	that
*	
*	sqr(si) + sqr(sj) - sqr(sk) = 2*si*sj*cos(theta_k)
*
*	where si, sj, and sk are the lengths of the three sides of the
*	triangle and theta_k is the interior angle of the triangle at the
*	vertex opposite to side k.
*/

LOCAL double angle_weight_of_tri(
	double     *p0,
	double     *p1,
	double     *p2,
	RECT_GRID *grid)
{
	double s00, s01, s02;
	double s10, s11, s12;
	double s20, s21, s22;
	double l0,l1,l2,angle;
	double t0,t1,t2;
	double h0 = grid->h[0], h1 = grid->h[1], h2 = grid->h[2];

	double h_min = min3(grid->h[0], grid->h[1], grid->h[2]);
	h0 = h_min;
	h1 = h_min;
	h2 = h_min;

	s00 = (p1[0]-p0[0])/h0; s01 = (p1[1]-p0[1])/h1; s02 = (p1[2]-p0[2])/h2;
	s10 = (p2[0]-p1[0])/h0; s11 = (p2[1]-p1[1])/h1; s12 = (p2[2]-p1[2])/h2;
	s20 = (p0[0]-p2[0])/h0; s21 = (p0[1]-p2[1])/h1; s22 = (p0[1]-p2[1])/h2;
	l0 = QDot3d(s0,s0); l1 = QDot3d(s1,s1); l2 = QDot3d(s2,s2);

	if (l0 == 0 || l1 == 0 || l2 == 0)
	    return 0;
	t0 = (l0+l1-l2)/sqrt(l0*l1);
	t1 = (l1+l2-l0)/sqrt(l1*l2);
	t2 = (l0+l2-l1)/sqrt(l0*l2);
	angle = 2.0 + (t0 < t1) ? min(t0,t2) : min(t1,t2);

	return angle;
}		/*end angle_weight_of_tri*/


LOCAL boolean flip_max_side_of_tri(
	 TRI		*tri,
	 int		side,
	 Front          *fr,
	 POINTER_Q	**pq)
{
	BOND_TRI     *bt;
	TRI          *otri;
	TRI_NEIGHBOR Nbr[3], ONbr[3];
	POINT        *pt1, *pt2;
	RECT_GRID    *gr = fr->rect_grid;
	double        *p, *op, *p1, *p2;
	double        f1,f2,l1,l2;
	int          i, nside, pside;
	int          oside, noside, poside;
	boolean         bdry[3], obdry[3];

	DEBUG_ENTER(flip_max_side_of_tri)

	if (is_side_bdry(tri,side))
	{
	    DEBUG_LEAVE(flip_max_side_of_tri)
	    return NO;
	}

	otri = Tri_on_side(tri,side);
	nside = Next_m3(side);
	pside = Prev_m3(side);
	p = Coords(Point_of_tri(tri)[side]);
	op = Coords(Point_of_tri(tri)[nside]);
	pt1 = Point_of_tri(tri)[pside];
	p1 = Coords(pt1);

	for (oside = 0; oside < 3; ++oside)
	    if (Tri_on_side(otri,oside) == tri)
		break;

	noside = Next_m3(oside);
	poside = Prev_m3(oside);
	pt2 = Point_of_tri(otri)[poside];
	p2 = Coords(pt2);

	if (two_points_share_side(pt1,tri,pt2,fr->interf) == 1)
	{
	    screen("WARNING: Cannot do flipping, side already existed\n");
	    *pq = dequeue(tri,*pq);
	    *pq = dequeue(Tri_on_side(tri,side),*pq);
	    DEBUG_LEAVE(flip_max_side_of_tri)
	    return NO;
	}

	l1 = angle_weight_of_tri(p,op,p1,gr);
	l2 = angle_weight_of_tri(p,p2,op,gr);

	if (DEBUG)
	{
	    (void) printf("flipping common side of tri %llu and otri %llu\n",
			  tri_number(tri,fr->interf),
			  tri_number(otri,fr->interf));
	    (void) printf("side = %d, oside = %d\n",side,oside);
	    (void) printf("tri - ");
	    print_tri(tri,fr->interf);
	    (void) printf("otri - ");
	    print_tri(otri,fr->interf);
	    (void) printf("angle weight of tri = %g\n",l1);
	    (void) printf("angle weight of otri = %g\n",l2);
	}

	if (l2 < l1)/* we may have bad otri*/
	{
	    if (find_scaled_extrem_edge(otri,gr,LONGEST) == oside)
	    	f1 = l2;
	    else
	    {
		*pq = dequeue(tri,*pq);
		*pq = dequeue(Tri_on_side(tri,side),*pq);
		DEBUG_LEAVE(flip_max_side_of_tri)
		return NO;
	    }
	}
	else
	    f1 = l1;

	l1 = angle_weight_of_tri(p,p2,p1,gr);
	l2 = angle_weight_of_tri(op,p1,p2,gr);

	if (DEBUG)
	{
	    (void) printf("angle weight of flipped tri = %g\n",l1);
	    (void) printf("angle weight of flipped otri = %g\n",l2);
	}

	f2 = min(l1,l2);
	if (f2 <= f1)
	{
	    /* The flipped tris with an max angle bigger than the
	     * original max angle. In this case we need to compare
	     * the area of the tris.  If the min_area of the
	     * original two will be smaller than that of the flipped
	     * two, do flipping.  The area business is good for
	     * these tris in delete_queue since in some case we
	     * couldn't delete_side. We do flipping instead.
	     */

	    double atri, otri, antri, aontri, amin, namin;

	    atri = scaled_tri_area(p,op,p1,gr);
	    otri = scaled_tri_area(p,p2,op,gr);
	    amin = min(atri,otri);
	    antri = scaled_tri_area(p,p2,p1,gr);
	    aontri = scaled_tri_area(op,p1,p2,gr);
	    namin = min(antri,aontri);
	    if (amin >= namin)
	    {
		*pq = dequeue(tri,*pq);
		*pq = dequeue(Tri_on_side(tri,side),*pq);
		DEBUG_LEAVE(flip_max_side_of_tri)
	        return NO;
	    }
	}
	
	if(skip_bdry_tri(tri) || skip_bdry_tri(otri))
	{
	    *pq = dequeue(tri,*pq);
	    *pq = dequeue(otri,*pq);
	    DEBUG_LEAVE(flip_max_side_of_tri)
	    return NO;
	}

	for (i = 0; i < 3; ++i)
	{
	    bdry[i] = is_side_bdry(tri,i) ? YES : NO;
	    Nbr[i] = Tri_neighbor(tri)[i];
	    obdry[i] = is_side_bdry(otri,i) ? YES : NO;
	    ONbr[i] = Tri_neighbor(otri)[i];
	}

	Point_of_tri(tri)[nside] = pt2;
	Tri_on_side(tri,nside) = otri;
	set_side_bdry(Boundary_tri(tri),nside,0);
	if (obdry[noside])
	{
	    bt = ONbr[noside].btri;
	    (void) link_tri_to_bond(bt,tri,bt->surface,bt->bond,bt->curve);
	}
	else
	{
	    set_side_bdry(Boundary_tri(tri),side,0);
	    Tri_on_side(tri,side) = ONbr[noside].tri;
	    for (i = 0; i < 3; ++i)
	    {
		if (Tri_on_side(ONbr[noside].tri,i) == otri)
		{
		    Tri_on_side(ONbr[noside].tri,i) = tri;
		    break;
		}
	    }
	}

	Point_of_tri(otri)[noside] = pt1;
	Tri_on_side(otri,noside) = tri;
	set_side_bdry(Boundary_tri(otri),noside,0);
	if (bdry[nside])
	{
	    bt = Nbr[nside].btri;
	    (void) link_tri_to_bond(bt,otri,bt->surface,bt->bond,bt->curve);
	}
	else
	{
	    set_side_bdry(Boundary_tri(otri),oside,0);
	    Tri_on_side(otri,oside) = Nbr[nside].tri;
	    for (i = 0; i < 3; ++i)
	    {
		if (Tri_on_side(Nbr[nside].tri,i) == tri)
		{
		    Tri_on_side(Nbr[nside].tri,i) = otri;
		    break;
		}
	    }
	}
	set_normal_of_tri(tri);
	set_normal_of_tri(otri);
	*pq = dequeue(tri,*pq);
	*pq = dequeue(otri,*pq);
	DEBUG_LEAVE(flip_max_side_of_tri)
	return YES;
}		/*end flip_max_side_of_tri*/


LOCAL boolean is_tri_in_queue(
	TRI		*tri,
	POINTER_Q		*pq)
{
	POINTER_Q		*tri_q;

	for (tri_q = head_of_pointer_queue(pq); tri_q; tri_q = tri_q->next)
	{
	    if (Tri_of_q(tri_q) == tri)
		return YES;
	}
	return NO;
}		/*end is_tri_in_queue*/

LOCAL boolean is_critical_side(
	TRI	*tri,
	int	side)
{
	TRI *nbtri = Tri_on_side(tri,side);
	TRI *nbp,*nbn;
	int i,nside;

	for (nside = 0; nside < 3; ++nside)
	    if (Tri_on_side(nbtri,nside) == tri) break;

	if (!is_side_bdry(tri,Prev_m3(side)) &&
	    !is_side_bdry(tri,Next_m3(side)))
	{
	    nbp = Tri_on_side(tri,Prev_m3(side));
	    nbn = Tri_on_side(tri,Next_m3(side));
	    for (i = 0; i < 3; ++i)
		if (Tri_on_side(nbp,i) == nbn) return YES;
	}

	if (!is_side_bdry(nbtri,Prev_m3(nside)) &&
	    !is_side_bdry(nbtri,Next_m3(nside)))
	{
	    nbp = Tri_on_side(nbtri,Prev_m3(nside));
	    nbn = Tri_on_side(nbtri,Next_m3(nside));
	    for (i = 0; i < 3; ++i)
		if (Tri_on_side(nbp,i) == nbn) return YES;
	}

	
	return NO;
}	/* end is_critical_side */

LOCAL boolean is_critical_tri(
	TRI	*tri,
	int	side)
{
	TRI *nbp,*nbn;
	int i,nside;

	if (!is_side_bdry(tri,Prev_m3(side)) &&
	    !is_side_bdry(tri,Next_m3(side)))
	{
	    nbp = Tri_on_side(tri,Prev_m3(side));
	    nbn = Tri_on_side(tri,Next_m3(side));
	    for (i = 0; i < 3; ++i)
		if (Tri_on_side(nbp,i) == nbn) return YES;
	}

	return NO;
}



LOCAL boolean check_and_rm_tetrahedron(
	TRI		*tri,
	int		side,
	POINTER_Q	**pq,
	SURFACE		*s)
{
	TRI	*nbp, *nbn, *nbtri;
	POINT	*p;
	int	nside;

	nbtri = Tri_on_side(tri,side);
	
	if (is_side_bdry(tri,Prev_m3(side)) ||
	    is_side_bdry(tri,Next_m3(side)))
	    return NO;
	    
	nbp = Tri_on_side(tri,Prev_m3(side));
	nbn = Tri_on_side(tri,Next_m3(side));
	p = Point_of_tri(tri)[Prev_m3(side)];

	nside = Vertex_of_point(nbp, p); 
	if(is_side_bdry(nbp,nside) || Tri_on_side(nbp,nside) != nbn)
	    return NO;
	
	nside = Next_m3(nside);
	if(is_side_bdry(nbp,nside) || Tri_on_side(nbp,nside) != nbtri)
	    return NO;
	
	for (nside = 0; nside < 3; ++nside)
	    if (Tri_on_side(nbtri,nside) == tri) break;
	
	nside = Prev_m3(nside);
	if(is_side_bdry(nbtri,nside) || Tri_on_side(nbtri,nside) != nbn)
	    return NO;

	printf("#check_and_rm_tetrahedron, removing tetrahedron.\n");
	
	if(pq == NULL)
	    return YES;
	
	*pq = dequeue(tri,*pq);
	*pq = dequeue(nbtri,*pq);
	*pq = dequeue(nbp,*pq);
	*pq = dequeue(nbn,*pq);

	remove_tri_from_surface(tri, s, YES);
	remove_tri_from_surface(nbtri, s, YES);
	remove_tri_from_surface(nbp, s, YES);
	remove_tri_from_surface(nbn, s, YES);
	
	return YES;
}

LOCAL	void	exchange_queues(
	POINTER_Q *pq1,
	POINTER_Q *pq2)
{
	TRI_SURF *ts1, *ts2, T;

	ts1 = Tri_surf(pq1);
	ts2 = Tri_surf(pq2);
	T = *ts1;
	*ts1 = *ts2;
	*ts2 = T;
	Tri_workspace(ts1->tri) = (POINTER) pq1;
	Tri_workspace(ts2->tri) = (POINTER) pq2;
}		/*end exchange_queues*/

boolean	compare_pointers(
	POINTER_Q	*pq1,
	POINTER_Q	*pq2)
{
	double	ave_norm, norm1, norm2;
	double	tol = 1.0e-8;

	norm1 = Tri_surf(pq1)->sqr_norm;
	norm2 = Tri_surf(pq2)->sqr_norm;
	ave_norm = (norm1 + norm2)*0.5;

	//two tris have very similar area, compare the postions.
	if( fabs(norm1 - norm2) < ave_norm*tol ) 
	    return   Tri_surf(pq1)->dist > Tri_surf(pq2)->dist;
	else
	    return   norm1 > norm2;
}

boolean	compare_bd_dist(
	POINTER_Q	*pq1,
	POINTER_Q	*pq2)
{
	double	ave_dist, dist1, dist2;
	double	tol = 1.0e-8;

	dist1 = Tri_surf(pq1)->bd_dist;
	dist2 = Tri_surf(pq2)->bd_dist;
	ave_dist = (dist1 + dist2)*0.5;

	//two tris have very similar area, compare the postions.
	if( fabs(dist1 - dist2) < ave_dist*tol ) 
	    return   Tri_surf(pq1)->dist > Tri_surf(pq2)->dist;
	else
	    return   dist1 > dist2;
}

/*ARGSUSED*/
LOCAL void sort_pointer_queue(
	POINTER_Q	*pq,
	INTERFACE	*intfc,
	SPQ_FLAG	flag)
{
	POINTER_Q	*pq1,*pq2;
	
	if (pq == NULL)
	    return;

	pq1 = head_of_pointer_queue(pq);
	while (pq1 != tail_of_pointer_queue(pq))
	{
	    pq2 = pq1->next;
	    if (flag == SHORTEST)
	    {
		if (compare_pointers(pq1, pq2))
	    	    exchange_queues(pq1,pq2);
	        
		while (pq2 != tail_of_pointer_queue(pq))
	        {
	    	    pq2 = pq2->next;
		    if (compare_pointers(pq1, pq2))
		    {
	    	        exchange_queues(pq1,pq2);
	    	    }
	        }
	    }
	    else if (flag == LONGEST)
	    {
		if (!compare_pointers(pq1, pq2))
	    	    exchange_queues(pq1,pq2);
	        
		while (pq2 != tail_of_pointer_queue(pq))
	        {
	    	    pq2 = pq2->next;
		    if (!compare_pointers(pq1, pq2))
		    {
	    	        exchange_queues(pq1,pq2);
	    	    }
	        }
	    }
	    else if (flag == BD_DIST)
	    {
		if (compare_bd_dist(pq1, pq2))
	    	    exchange_queues(pq1,pq2);
	        
		while (pq2 != tail_of_pointer_queue(pq))
	        {
	    	    pq2 = pq2->next;
		    if (compare_bd_dist(pq1, pq2))
		    {
	    	        exchange_queues(pq1,pq2);
	    	    }
	        }
	    }

	    pq1 = pq1->next;
	}
}		/*end sort_pointer_queue*/


LOCAL void redistribute_curve3d(
	CURVE		*c,
	Front		*fr)
{

	BOND		*b;
	static boolean	first = YES;
	static double	max_b_length, min_b_length;


	if (first == YES)
	{
	    first = NO;
	    max_b_length = Max_bond_len(fr,GENERAL_WAVE);
	    min_b_length = Min_bond_len(fr,GENERAL_WAVE);
	}
       
	if (debugging("b_length"))
	{
	    (void) printf("max_b_length = %g, min_b_length = %g\n",
	                  max_b_length,min_b_length);
	    detail_of_curve(c);
	}
	
	for (b = c->first; b != NULL; b = b->next)
	{
	    if (b->length >= max_b_length)
	    {
	        POINT *pm = Point(NULL);
		int   i;
	        for (i = 0; i < 3; ++i)
		    Coords(pm)[i] = 0.5*(Coords(b->start)[i]+Coords(b->end)[i]);
	        (void) insert_point_in_bond(pm,b,c);
		b = b->next;
	    }
	    else if (hsbdry_type(c) < FIRST_PHYSICS_HSBDRY_TYPE)
	    {
	        if (b->length <= min_b_length)
		{
		    BOND *bp, *bn;
		    if ((bp = b->prev))
		    {
			if (!delete_start_of_bond(b,c))
		        {
			    screen("ERROR in redistribute_curve3d(), "
			           "delete_start_of_bond() failed\n");
			    clean_up(ERROR);
		        }
			b = bp;
		    }
		    else if ((bn = b->next))
		    {
			if (!delete_end_of_bond(b,c))
		        {
			    screen("ERROR in redistribute_curve3d(), "
			           "delete_end_of_bond() failed\n");
			    clean_up(ERROR);
		        }
			b = bn;
		    }
		}
	    }
	}
	if (debugging("redist_curve"))
	{
	    detail_of_curve(c);
	    summarize_interface("redist_curve3d","exit",c->interface,
				XY_PLANE,"redist_curve3d","exit");
	}
}		/*end redistribute_curve3d */

LOCAL double scaled_tri_area(
	double     *p0,
	double     *p1,
	double     *p2,
	RECT_GRID *gr)
{
	double   	s00, s01, s02;
	double          s20, s21, s22;
	double          N0, N1, N2;
	double		h0 = gr->h[0], h1 = gr->h[1], h2 = gr->h[2];
	double		sqr_area;

	double h_min = min3(gr->h[0], gr->h[1], gr->h[2]);
	h0 = h_min;
	h1 = h_min;
	h2 = h_min;

	s00 = (p1[0]-p0[0])/h0; s01 = (p1[1]-p0[1])/h1; s02 = (p1[2]-p0[2])/h2;
	s20 = (p0[0]-p2[0])/h0; s21 = (p0[1]-p2[1])/h1; s22 = (p0[2]-p2[2])/h2;
	QCross3d(s0,s2,N);

	sqr_area = 0.25*QDot3d(N,N);

	return	sqr_area;
}		/*end scaled_tri_area*/

LOCAL 	void	print_tri_surf_queue(
	POINTER_Q *p_q)
{
  	POINTER_Q 	*q;
	TRI_SURF 	*t_surf;
	TRI 		*tri;
	int		cnt = 0;
  
	if (p_q == NULL) 
	{
	    (void) printf("NULL POINTER_Q\n");
	    return;
        }

	q = head_of_pointer_queue(p_q);
	while (q != tail_of_pointer_queue(p_q))
        {
	    t_surf = Tri_surf(q);
	    tri = t_surf->tri;
	    (void) printf("%3d ( %g %g %g ) ( %g %g %g ) "
			  "( %g %g %g )  %p  side = %d\n", ++cnt,
			  Coords(Point_of_tri(tri)[0])[0],
			  Coords(Point_of_tri(tri)[0])[1],
			  Coords(Point_of_tri(tri)[0])[2],
			  Coords(Point_of_tri(tri)[1])[0],
			  Coords(Point_of_tri(tri)[1])[1],
			  Coords(Point_of_tri(tri)[1])[2],
			  Coords(Point_of_tri(tri)[2])[0],
			  Coords(Point_of_tri(tri)[2])[1],
			  Coords(Point_of_tri(tri)[2])[2],tri,t_surf->side);
	    q = q->next;
	}
}		/*end print_tri_surf_queue*/

EXPORT	boolean	point_outside_open_bdry(
	int		*k,
	double		*nor,
	POINT		*p,
	INTERFACE	*intfc)
{
	int		i;
	RECT_GRID	*gr = computational_grid(intfc);
	double   	tol = 1.0e-4;

	zero_scalar(nor, 3*FLOAT);
	for(i=0; i<3; i++)
	{
	    if(rect_boundary_type(intfc,i,0) == OPEN_BOUNDARY && 
	       Coords(p)[i] < gr->VL[i] + gr->h[i])
	    {
		*k = i;
		nor[i] = -1.0;
		nor[3] = (gr->VL[i] - Coords(p)[i])/gr->h[i];
		return YES;
	    }
	    if(rect_boundary_type(intfc,i,1) == OPEN_BOUNDARY && 
	       Coords(p)[i] > gr->VU[i] - gr->h[i])
	    {
		*k = i;
		nor[i] = 1.0;
		nor[3] = (Coords(p)[i] - gr->VU[i])/gr->h[i];
		return YES;
	    }
	}
	return NO;
}	/* end point_outside_open_bdry */


//height vector in edge p p2
EXPORT	void	triangle_height_vec(
	double		*hv,
	double		*p,
	double   	*p1,
	double		*p2)
{
	double	v[3], nor[3], len;
	int	k;

	difference(p1, p, v, 3);
	difference(p2, p, nor, 3);
	len = Dot3d(v,nor)/Dot3d(nor,nor);
	for(k=0; k<3; k++)
	    hv[k] = v[k] - len*nor[k];
}


EXPORT  void    tecplot_interface_in_ball(
	const char	*bname,
	INTERFACE	*intfc)
{
	SURFACE	**s;

	for (s = intfc->surfaces; s && *s; ++s)
	{
	    if(wave_type(*s) < FIRST_PHYSICS_WAVE_TYPE)
	        continue;
	    printf("#show surface in ball\n");
	    tecplot_surface_in_ball(bname,*s);
	    break;
	}

}	/* end tecplot_interface */

