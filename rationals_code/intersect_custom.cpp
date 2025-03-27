//
// Created by Michele on 13/06/24.
//

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "intersect_custom.h"
#include <stack>
#include "cinolib/predicates.h"

bigrational zero_rat = bigrational(0,0,0);

// returns:
// DO_NOT_INTERSECT     if s and t are fully disjoint
// SIMPLICIAL_COMPLEX   if s is an edge of t, or s is degenerate and coincides with a vertex of t
// INTERSECT            if s and t intersect and do not form a valid simplex

bool points_are_colinear_2d(const bigrational * p0,
                                   const bigrational * p1,
                                   const bigrational * p2)
{
    //return (cinolib::orient2d(p0,p1,p2).sgn()==0);
    return orient2dT(p0,p1,p2) == 0;

}


bool points_are_colinear_3d(const bigrational * p0,
                                   const bigrational * p1,
                                   const bigrational * p2)
{
    bigrational p0_dropX[2] = { p0[1], p0[2] };
    bigrational p0_dropY[2] = { p0[0], p0[2] };
    bigrational p0_dropZ[2] = { p0[0], p0[1] };
    bigrational p1_dropX[2] = { p1[1], p1[2] };
    bigrational p1_dropY[2] = { p1[0], p1[2] };
    bigrational p1_dropZ[2] = { p1[0], p1[1] };
    bigrational p2_dropX[2] = { p2[1], p2[2] };
    bigrational p2_dropY[2] = { p2[0], p2[2] };
    bigrational p2_dropZ[2] = { p2[0], p2[1] };

    // check if all the 2d orthogonal projections of p are colinear
    if(points_are_colinear_2d(p0_dropX, p1_dropX, p2_dropX) &&
       points_are_colinear_2d(p0_dropY, p1_dropY, p2_dropY) &&
       points_are_colinear_2d(p0_dropZ, p1_dropZ, p2_dropZ))
    {
        return true;
    }
    return false;
}

bool triangle_is_degenerate_3d(const bigrational * t0,
                                      const bigrational * t1,
                                      const bigrational * t2)
{
    return points_are_colinear_3d(t0, t1, t2);
}

bool vec_equals_3d(const bigrational * v0,
                          const bigrational * v1)
{
    return ((v0[0] == v1[0]) &&
            (v0[1] == v1[1]) &&
            (v0[2] == v1[2]));
}

bool segment_is_degenerate_3d(const bigrational * s0,
                                     const bigrational * s1)
{
    return vec_equals_3d(s0, s1);
}

PointInSimplex point_in_segment_3d(const bigrational * p,
                                          const bigrational * s0,
                                          const bigrational * s1)
{
    if(vec_equals_3d(p, s0)) return ON_VERT0;
    if(vec_equals_3d(p, s1)) return ON_VERT1;

    if(!points_are_colinear_3d(s0, s1, p)) return STRICTLY_OUTSIDE;

    if((p[0] > std::min(s0[0], s1[0]) && p[0] < std::max(s0[0], s1[0])) ||
       (p[1] > std::min(s0[1], s1[1]) && p[1] < std::max(s0[1], s1[1])) ||
       (p[2] > std::min(s0[2], s1[2]) && p[2] < std::max(s0[2], s1[2])))
    {
        return STRICTLY_INSIDE;
    }

    return STRICTLY_OUTSIDE;
}

bool vec_equals_2d(const bigrational * v0,
                          const bigrational * v1)
{
    return ((v0[0] == v1[0]) &&
            (v0[1] == v1[1]));
}

PointInSimplex point_in_triangle_2d(const bigrational * p,
                                           const bigrational * t0,
                                           const bigrational * t1,
                                           const bigrational * t2)
{
    if(vec_equals_2d(p, t0)) return ON_VERT0;
    if(vec_equals_2d(p, t1)) return ON_VERT1;
    if(vec_equals_2d(p, t2)) return ON_VERT2;

    //bigrational e0p_area = cinolib::orient2d(&t0[0], &t1[0], &p[0]);
    //bigrational e1p_area = cinolib::orient2d(&t1[0], &t2[0], &p[0]);
    //bigrational e2p_area = cinolib::orient2d(&t2[0], &t0[0], &p[0]);

    int e0p_area = orient2dT(&t0[0], &t1[0], &p[0]);
    int e1p_area = orient2dT(&t1[0], &t2[0], &p[0]);
    int e2p_area = orient2dT(&t2[0], &t0[0], &p[0]);

    //bool hit = (((e0p_area > zero_rat || e0p_area.sgn() == 0) && (e1p_area > zero_rat || e1p_area.sgn() == 0) && (e2p_area > zero_rat || e2p_area.sgn() == 0)) ||
    //            ((e0p_area < zero_rat || e0p_area.sgn() == 0) && (e1p_area > zero_rat || e1p_area.sgn() == 0) && (e2p_area > zero_rat || e2p_area.sgn() == 0)));

    bool hit = (((e0p_area >= 0) && (e1p_area >= 0) && (e2p_area >= 0)) ||
                ((e0p_area <= 0) && (e1p_area <= 0) && (e2p_area <= 0)));

    if(hit)
    {
        //if(e0p_area.sgn() == 0) return ON_EDGE0;
        //if(e1p_area.sgn() == 0) return ON_EDGE1;
        //if(e2p_area.sgn() == 0) return ON_EDGE2;

        if(e0p_area == 0) return ON_EDGE0;
        if(e1p_area == 0) return ON_EDGE1;
        if(e2p_area == 0) return ON_EDGE2;

        return STRICTLY_INSIDE;
    }

    return STRICTLY_OUTSIDE;
}

PointInSimplex point_in_triangle_3d(const bigrational * p,
                                           const bigrational * t0,
                                           const bigrational * t1,
                                           const bigrational * t2)
{
    // test for point in vert
    if(vec_equals_3d(p, t0)) return ON_VERT0;
    if(vec_equals_3d(p, t1)) return ON_VERT1;
    if(vec_equals_3d(p, t2)) return ON_VERT2;

    // test for point in edge in 3D
    if(point_in_segment_3d(p, t0, t1) == STRICTLY_INSIDE) return ON_EDGE0;
    if(point_in_segment_3d(p, t1, t2) == STRICTLY_INSIDE) return ON_EDGE1;
    if(point_in_segment_3d(p, t2, t0) == STRICTLY_INSIDE) return ON_EDGE2;

    // test for the interior: project t on XYZ and, if the check is never false in
    // any of the projections, then p must be inside t

    bigrational p_dropX[2]  = { p[1],  p[2]};
    bigrational t0_dropX[2] = {t0[1], t0[2]}, t1_dropX[2] = {t1[1], t1[2]}, t2_dropX[2] = {t2[1], t2[2]};

    if(point_in_triangle_2d(p_dropX, t0_dropX, t1_dropX, t2_dropX) == STRICTLY_OUTSIDE) return STRICTLY_OUTSIDE;

    bigrational p_dropY[2]  = { p[0],  p[2]};
    bigrational t0_dropY[2] = {t0[0], t0[2]}, t1_dropY[2] = {t1[0], t1[2]}, t2_dropY[2] = {t2[0], t2[2]};

    if(point_in_triangle_2d(p_dropY, t0_dropY, t1_dropY, t2_dropY) == STRICTLY_OUTSIDE) return STRICTLY_OUTSIDE;

    bigrational p_dropZ[2]  = { p[0],  p[1]};
    bigrational t0_dropZ[2] = {t0[0], t0[1]}, t1_dropZ[2] = {t1[0], t1[1]}, t2_dropZ[2] = {t2[0], t2[1]};

    if(point_in_triangle_2d(p_dropZ, t0_dropZ, t1_dropZ, t2_dropZ) == STRICTLY_OUTSIDE) return STRICTLY_OUTSIDE;

    return STRICTLY_INSIDE;
}

bool points_are_coplanar_3d(const bigrational * p0,
                                   const bigrational * p1,
                                   const bigrational * p2,
                                   const bigrational * p3)
{
    //return (cinolib::orient3d(&p0[0],&p1[0],&p2[0],&p3[0]).sgn()==0);
    return (orient3dT(&p0[0],&p1[0],&p2[0],&p3[0]) == 0);
}

SimplexIntersection segment_segment_intersect_2d(const bigrational * s00,
                                                        const bigrational * s01,
                                                        const bigrational * s10,
                                                        const bigrational * s11)
{
    // https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
    //bigrational det_s00 = cinolib::orient2d(s10, s11, s00);
    //bigrational det_s01 = cinolib::orient2d(s10, s11, s01);
    //bigrational det_s10 = cinolib::orient2d(s00, s01, s10);
    //bigrational det_s11 = cinolib::orient2d(s00, s01, s11);

    int det_s00 = orient2dT(s10, s11, s00);
    int det_s01 = orient2dT(s10, s11, s01);
    int det_s10 = orient2dT(s00, s01, s10);
    int det_s11 = orient2dT(s00, s01, s11);

    // Shewchuk's orient predicates return a rough approximation of the determinant.
    // I am converting values to { -1, 0, 1 } for a simpler check of intersection cases
    //int s00_wrt_s1 = (det_s00 > zero_rat) ? 1 : ((det_s00 < zero_rat) ? -1 : 0);
    //int s01_wrt_s1 = (det_s01 > zero_rat) ? 1 : ((det_s01 < zero_rat) ? -1 : 0);
    //int s10_wrt_s0 = (det_s10 > zero_rat) ? 1 : ((det_s10 < zero_rat) ? -1 : 0);
    //int s11_wrt_s0 = (det_s11 > zero_rat) ? 1 : ((det_s11 < zero_rat) ? -1 : 0);

    int s00_wrt_s1 = (det_s00 > 0) ? 1 : ((det_s00 < 0) ? -1 : 0);
    int s01_wrt_s1 = (det_s01 > 0) ? 1 : ((det_s01 < 0) ? -1 : 0);
    int s10_wrt_s0 = (det_s10 > 0) ? 1 : ((det_s10 < 0) ? -1 : 0);
    int s11_wrt_s0 = (det_s11 > 0) ? 1 : ((det_s11 < 0) ? -1 : 0);

    // segments intersect at a single point
    if(s00_wrt_s1 != s01_wrt_s1 && s10_wrt_s0 != s11_wrt_s0)
    {
        // edges share an endpoint
        if(vec_equals_2d(s00, s10) || vec_equals_2d(s00, s11) || vec_equals_2d(s01, s10) || vec_equals_2d(s01, s11))
            return SIMPLICIAL_COMPLEX;

        // at least one segment endpoint is involved in the intersection
        return INTERSECT;
    }

    // degenerate case: colinear segments
    if(s00_wrt_s1 == 0 && s01_wrt_s1 == 0 && s10_wrt_s0 == 0 && s11_wrt_s0 == 0)
    {
        // coincident segments
        if((vec_equals_2d(s00, s10) && vec_equals_2d(s01, s11)) || (vec_equals_2d(s00, s11) && vec_equals_2d(s01, s10)))
            return SIMPLICIAL_COMPLEX;

        bigrational Xmin_s1 = std::min(s10[0], s11[0]);
        bigrational Xmax_s1 = std::max(s10[0], s11[0]);
        bigrational Ymin_s1 = std::min(s10[1], s11[1]);
        bigrational Ymax_s1 = std::max(s10[1], s11[1]);
        bigrational Xmin_s0 = std::min(s00[0], s01[0]);
        bigrational Xmax_s0 = std::max(s00[0], s01[0]);
        bigrational Ymin_s0 = std::min(s00[1], s01[1]);
        bigrational Ymax_s0 = std::max(s00[1], s01[1]);

        if(// test s0 endpoints against s1 range
                (s00[0] > Xmin_s1 && s00[0] < Xmax_s1) ||
                (s00[1] > Ymin_s1 && s00[1] < Ymax_s1) ||
                (s01[0] > Xmin_s1 && s01[0] < Xmax_s1) ||
                (s01[1] > Ymin_s1 && s01[1] < Ymax_s1) ||
                // test s1 endpoints against s0 range
                (s10[0] > Xmin_s0 && s10[0] < Xmax_s0) ||
                (s10[1] > Ymin_s0 && s10[1] < Ymax_s0) ||
                (s11[0] > Xmin_s0 && s11[0] < Xmax_s0) ||
                (s11[1] > Ymin_s0 && s11[1] < Ymax_s0))
        {
            return INTERSECT;
        }
    }
    return DO_NOT_INTERSECT;
}

SimplexIntersection segment_segment_intersect_3d(const bigrational * s00,
                                                        const bigrational * s01,
                                                        const bigrational * s10,
                                                        const bigrational * s11)
{
    assert(!segment_is_degenerate_3d(s00, s01) && !segment_is_degenerate_3d(s10, s11));

    if(!points_are_coplanar_3d(s00, s01, s10, s11)) return DO_NOT_INTERSECT;

    // check for coincident segments
    bool s00_is_shared = (vec_equals_3d(s00, s10) || vec_equals_3d(s00, s11));
    bool s01_is_shared = (vec_equals_3d(s01, s10) || vec_equals_3d(s01, s11));
    bool s10_is_shared = (vec_equals_3d(s10, s00) || vec_equals_3d(s10, s01));
    bool s11_is_shared = (vec_equals_3d(s11, s00) || vec_equals_3d(s11, s01));

    // s0 and s1 are coincident or one edge is degenerate and coincides with one vertex of the other
    if(s00_is_shared && s01_is_shared && s10_is_shared && s11_is_shared) return SIMPLICIAL_COMPLEX;

    // check 2D projections of the segments

    bigrational s00_dropX[2] = {s00[1], s00[2]}, s01_dropX[2] = {s01[1], s01[2]};
    bigrational s10_dropX[2] = {s10[1], s10[2]}, s11_dropX[2] = {s11[1], s11[2]};
    int x_res = segment_segment_intersect_2d(s00_dropX, s01_dropX, s10_dropX, s11_dropX);
    if (x_res == DO_NOT_INTERSECT) return DO_NOT_INTERSECT;

    bigrational s00_dropY[2] = {s00[0], s00[2]}, s01_dropY[2] = {s01[0], s01[2]};
    bigrational s10_dropY[2] = {s10[0], s10[2]}, s11_dropY[2] = {s11[0], s11[2]};
    int y_res = segment_segment_intersect_2d(s00_dropY, s01_dropY, s10_dropY, s11_dropY);
    if (y_res == DO_NOT_INTERSECT) return DO_NOT_INTERSECT;

    bigrational s00_dropZ[2] = {s00[0], s00[1]}, s01_dropZ[2] = {s01[0], s01[1]};
    bigrational s10_dropZ[2] = {s10[0], s10[1]}, s11_dropZ[2] = {s11[0], s11[1]};
    int z_res = segment_segment_intersect_2d(s00_dropZ, s01_dropZ, s10_dropZ, s11_dropZ);
    if (z_res == DO_NOT_INTERSECT) return DO_NOT_INTERSECT;

    // segments are deemed overlapping if they are so in at least two projections our of three
    // (overlapping axis aligned segments will look like a valid simplcial complex in one projection)
    if((x_res == OVERLAP && y_res == OVERLAP) ||
       (x_res == OVERLAP && z_res == OVERLAP) ||
       (y_res == OVERLAP && z_res == OVERLAP))
    {
        return OVERLAP;
    }

    return INTERSECT;

}



SimplexIntersection segment_triangle_intersect_3d(const bigrational * s0,
                                                         const bigrational * s1,
                                                         const bigrational * t0,
                                                         const bigrational * t1,
                                                         const bigrational * t2) {
    assert(!segment_is_degenerate_3d(s0, s1));
    assert(!triangle_is_degenerate_3d(t0, t1, t2));

    if((vec_equals_3d(s0, t0) || vec_equals_3d(s0, t1) || vec_equals_3d(s0, t2)) &&
       (vec_equals_3d(s1, t0) || vec_equals_3d(s1, t1) || vec_equals_3d(s1, t2)))
    {
        return SIMPLICIAL_COMPLEX;
    }
    std::array<bigrational,3 > s0_tmp = {bigrational(3074457343958179937,153576,1),
                          bigrational(3572095111,235216,1),
                          bigrational(9320677326,108497,1)
    };

    std::array<bigrational,3 > s1_tmp = {bigrational(116913819,89206,-1),
                         bigrational(2305843008117596711,17776,1),
                         bigrational()
    };

    std::array<bigrational,3 > s2_tmp = {bigrational(5504579610,1,-1),
                         bigrational(3074457344089985827,40257,-1),
                         bigrational()
    };

    std::array<bigrational,3 > s3_tmp = {bigrational(18446744070889277111,918766,-1),
                         bigrational(5870338011,958118,-1),
                         bigrational(18446744071096704124,682923,1)
    };
    std::cout << "test" << std::endl;

    //bigrational vol_test = cinolib::orient3d(&s0_tmp[0],&s1_tmp[0],&s2_tmp[0],&s3_tmp[0]);
    int vol_test = orient3dT(&s0_tmp[0],&s1_tmp[0],&s2_tmp[0],&s3_tmp[0]);
    std::cout << "Values before first orient3d: " << std::endl;
    std::cout << "s0: " << s0[0] << " " << s0[1] << " " << s0[2] << std::endl;
    std::cout << "t0: " << t0[0] << " " << t0[1] << " " << t0[2] << std::endl;
    std::cout << "t1: " << t1[0] << " " << t1[1] << " " << t1[2] << std::endl;
    std::cout << "t2: " << t2[0] << " " << t2[1] << " " << t2[2] << std::endl;
    //bigrational vol_s0_t = cinolib::orient3d(s0, t0, t1, t2);
    int vol_s0_t = orient3dT(s0, t0, t1, t2);

    std::cout << "Values before second orient3d: " << std::endl;
    std::cout << "s1: " << s1[0] << " " << s1[1] << " " << s1[2] << std::endl;
    std::cout << "t0: " << t0[0] << " " << t0[1] << " " << t0[2] << std::endl;
    std::cout << "t1: " << t1[0] << " " << t1[1] << " " << t1[2] << std::endl;
    std::cout << "t2: " << t2[0] << " " << t2[1] << " " << t2[2] << std::endl;
    //bigrational vol_s1_t = cinolib::orient3d(s1, t0, t1, t2);
    int vol_s1_t = orient3dT(s1, t0, t1, t2);

    //if(vol_s0_t > zero_rat && vol_s1_t > zero_rat) return DO_NOT_INTERSECT; // s is above t
    //if(vol_s0_t < zero_rat && vol_s1_t < zero_rat) return DO_NOT_INTERSECT; // s is below t

    if(vol_s0_t > 0 && vol_s1_t > 0) return DO_NOT_INTERSECT; // s is above t
    if(vol_s0_t < 0 && vol_s1_t < 0) return DO_NOT_INTERSECT; // s is below t

    //if(vol_s0_t.sgn() == 0 && vol_s1_t.sgn() == 0)                        // s and t are coplanar
    if(vol_s0_t == 0 && vol_s1_t == 0)                        // s and t are coplanar
    {
        // same code as the 2D version, I just copied it here....

        if(point_in_triangle_3d(s0, t0, t1, t2) || point_in_triangle_3d(s1, t0, t1, t2)) return INTERSECT;

        uint simpl_complex = 0;

        switch(segment_segment_intersect_3d(s0, s1, t0, t1))
        {
            case SIMPLICIAL_COMPLEX : ++simpl_complex; break;
            case INTERSECT          : return INTERSECT;
            case OVERLAP            : return INTERSECT;
            case DO_NOT_INTERSECT   : break;
        }

        switch(segment_segment_intersect_3d(s0, s1, t1, t2))
        {
            case SIMPLICIAL_COMPLEX : ++simpl_complex; break;
            case INTERSECT          : return INTERSECT;
            case OVERLAP            : return INTERSECT;
            case DO_NOT_INTERSECT   : break;
        }

        switch(segment_segment_intersect_3d(s0, s1, t2, t0))
        {
            case SIMPLICIAL_COMPLEX : ++simpl_complex; break;
            case INTERSECT          : return INTERSECT;
            case OVERLAP            : return INTERSECT;
            case DO_NOT_INTERSECT   : break;
        }

        // if it is a simplicial complex from any view, then it really is...
        if(simpl_complex == 3) return SIMPLICIAL_COMPLEX;
        return DO_NOT_INTERSECT;
    }

    // s intersects t (borders included), if the signs of the three tetrahedra
    // obtained combining s with the three edges of t are all equal

    // if one point is coplanar and coincides with a triangle vertex, then they form a valid complex
    if(vec_equals_3d(s0, t0) || vec_equals_3d(s0, t1) || vec_equals_3d(s0, t2) ||
       vec_equals_3d(s1, t0) || vec_equals_3d(s1, t1) || vec_equals_3d(s1, t2))
    {
        return SIMPLICIAL_COMPLEX;
    }

    // bigrational vol_s_t01 = cinolib::orient3d(s0, s1, t0, t1);
    // bigrational vol_s_t12 = cinolib::orient3d(s0, s1, t1, t2);
    // bigrational vol_s_t20 = cinolib::orient3d(s0, s1, t2, t0);

    int vol_s_t01 = orient3dT(s0, s1, t0, t1);
    int vol_s_t12 = orient3dT(s0, s1, t1, t2);
    int vol_s_t20 = orient3dT(s0, s1, t2, t0);

    //if((vol_s_t01 > zero_rat && vol_s_t12 < zero_rat) || (vol_s_t01 < zero_rat && vol_s_t12 > zero_rat)) return DO_NOT_INTERSECT;
    //if((vol_s_t12 > zero_rat && vol_s_t20 < zero_rat) || (vol_s_t12 < zero_rat && vol_s_t20 > zero_rat)) return DO_NOT_INTERSECT;
    //if((vol_s_t20 > zero_rat && vol_s_t01 < zero_rat) || (vol_s_t20 < zero_rat && vol_s_t01 > zero_rat)) return DO_NOT_INTERSECT;

    if((vol_s_t01 > 0 && vol_s_t12 < 0) || (vol_s_t01 < 0 && vol_s_t12 > 0)) return DO_NOT_INTERSECT;
    if((vol_s_t12 > 0 && vol_s_t20 < 0) || (vol_s_t12 < 0 && vol_s_t20 > 0)) return DO_NOT_INTERSECT;
    if((vol_s_t20 > 0 && vol_s_t01 < 0) || (vol_s_t20 < 0 && vol_s_t01 > 0)) return DO_NOT_INTERSECT;

    return INTERSECT;
}

void cross(const bigrational* va,
           const bigrational* vb,
           bigrational* res)
{
    res[0] = va[1] * vb[2] - va[2] * vb[1];
    res[1] = va[2] * vb[0] - va[0] * vb[2];
    res[2] = va[0] * vb[1] - va[1] * vb[0];
}


void triangle_normal(const bigrational* pa,
                     const bigrational* pb,
                     const bigrational* pc,
                     bigrational* n) // n is the normal of triangle abc
{
    //print the values above
    std::cout << "pb [0] : " << pb[0] << std::endl;
    std::cout << "pb[1] : " << pb[1] << std::endl;
    std::cout << "pb[2] : " << pb[2] << std::endl;
    bigrational v0[3] = { pb[0]-pa[0], pb[1]-pa[1], pb[2]-pa[2] };


    bigrational v1[3] = { pc[0]-pa[0], pc[1]-pa[1], pc[2]-pa[2] };
    cross(v0,v1,n);
}

bigrational dot(const bigrational * pa,
                             const bigrational * pb)
{
    std::cout << "pa[0] : " << pa[0] << std::endl;
    std::cout << "pa[1] : " << pa[1] << std::endl;
    std::cout << "pa[2] : " << pa[2] << std::endl;
    std::cout << "pb[0] : " << pb[0] << std::endl;
    std::cout << "pb[1] : " << pb[1] << std::endl;
    std::cout << "pb[2] : " << pb[2] << std::endl;

    const bigrational first = pa[0] * pb[0];
    const bigrational second = pa[1] * pb[1];
    const bigrational third = pa[2] * pb[2];
    const bigrational result = first + second + third;
    //nfgMemoryPool;

    return result;
    //return pa[0] * pb[0] + pa[1] * pb[1] + pa[2] * pb[2];
}


void plane_line_intersection(const bigrational* p0,
                             const bigrational* p1,
                             const bigrational* p2,
                             const bigrational* l0,
                             const bigrational* l1,
                             bigrational* res)
{
    // https://en.wikipedia.org/wiki/Line–plane_intersection

    bigrational n[3];
    //print p0
    std::cout << "Dentro plane line intersection " << std::endl;
    std::cout << "p0 : " << p0[0] << " " << p0[1] << " " << p0[2] << std::endl;
    std::cout << "p1 : " << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
    std::cout << "p2 : " << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;

    triangle_normal(&p0[0],&p1[0],&p2[0],&n[0]);

    bigrational l[3];
    l[0] = l1[0] - l0[0];
    l[1] = l1[1] - l0[1];
    l[2] = l1[2] - l0[2];

    bigrational pl[3];
    pl[0] = p0[0] - l0[0];
    pl[1] = p0[1] - l0[1];
    pl[2] = p0[2] - l0[2];

    bigrational d = dot(&pl[0],&n[0])/dot(&l[0],&n[0]);

    res[0] = l0[0] + l[0] * d;
    res[1] = l0[1] + l[1] * d;
    res[2] = l0[2] + l[2] * d;

    std::cout << "Print before assert" << std::endl;
    std::cout << "p0 : " << p0[0] << " " << p0[1] << " " << p0[2] << std::endl;
    std::cout << "p1 : " << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
    std::cout << "p2 : " << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
    std::cout << "res: " << res[0] << " " << res[1] << " " << res[2] << std::endl;
    //assert(cinolib::orient3d(p0,p1,p2,res).sgn() == 0);
    assert(orient3dT(p0,p1,p2,res) == 0);
}




