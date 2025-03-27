//
// Created by Michele on 13/06/24.
//

#ifndef GITBOOLEANS_INTERSECT_CUSTOM_H
#define GITBOOLEANS_INTERSECT_CUSTOM_H

#include <rationals.h>
#include <numerics.h>
#include <type_traits>
#include <code/triangles_intersections/orientTemplated.h>

typedef enum
{
    STRICTLY_OUTSIDE = 0,  // strictly outside the input simplex
    STRICTLY_INSIDE  = 1,  // strictly inside  the input simplex
    ON_VERT0         = 2,  // used for segs, tris and tets
    ON_VERT1         = 3,  // used for segs, tris and tets
    ON_VERT2         = 4,  // used for tris and tets
    ON_VERT3         = 5,  // used for tets
    ON_EDGE0         = 6,  // used for tris and tets
    ON_EDGE1         = 7,  // used for tris and tets
    ON_EDGE2         = 8,  // used for tris and tets
    ON_EDGE3         = 9,  // used for tets
    ON_EDGE4         = 10, // used for tets
    ON_EDGE5         = 11, // used for tets
    ON_FACE0         = 12, // used for tets
    ON_FACE1         = 13, // used for tets
    ON_FACE2         = 14, // used for tets
    ON_FACE3         = 15, // used for tets
}
        PointInSimplex;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// intersection types
typedef enum
{
    DO_NOT_INTERSECT   = 0, // simplices do not intersect
    SIMPLICIAL_COMPLEX = 1, // simplices form a valid simplicial complex (i.e. they are coincident or share a sub-simplex)
    INTERSECT          = 2, // simplices intersect in a non conforming way
    OVERLAP            = 3, // for corner cases: simplices intersect and partially overlap
}                           // (e.g. colinear segments or coplanar triangles)
SimplexIntersection;



bool points_are_colinear_2d(const bigrational * p0, const bigrational * p1, const bigrational * p2);

bool points_are_colinear_3d(const bigrational * p0, const bigrational * p1, const bigrational * p2);

bool triangle_is_degenerate_3d(const bigrational * t0, const bigrational * t1, const bigrational * t2);

bool vec_equals_3d(const bigrational * v0, const bigrational * v1);

bool segment_is_degenerate_3d(const bigrational * s0, const bigrational * s1);

PointInSimplex point_in_segment_3d(const bigrational * p, const bigrational * s0, const bigrational * s1);

bool vec_equals_2d(const bigrational * v0, const bigrational * v1);

PointInSimplex point_in_triangle_2d(const bigrational * p, const bigrational * t0, const bigrational * t1, const bigrational * t2);

PointInSimplex point_in_triangle_3d(const bigrational * p, const bigrational * t0, const bigrational * t1, const bigrational * t2);

bool points_are_coplanar_3d(const bigrational * p0, const bigrational * p1, const bigrational * p2, const bigrational * p3);

SimplexIntersection segment_segment_intersect_2d(const bigrational * s00, const bigrational * s01, const bigrational * s10, const bigrational * s11);

SimplexIntersection segment_segment_intersect_3d(const bigrational * s00, const bigrational * s01, const bigrational * s10, const bigrational * s11);

SimplexIntersection segment_triangle_intersect_3d(const bigrational * s0, const bigrational * s1, const bigrational * t0, const bigrational * t1, const bigrational * t2);

void cross(const bigrational* va, const bigrational* vb, bigrational* res);

void triangle_normal(const bigrational* pa, const bigrational* pb, const bigrational* pc, bigrational* n);

bigrational dot(const bigrational * pa, const bigrational * pb);

void plane_line_intersection(const bigrational* p0, const bigrational* p1, const bigrational* p2, const bigrational* l0, const bigrational* l1, bigrational* res);


#endif //GITBOOLEANS_INTERSECT_CUSTOM_H
