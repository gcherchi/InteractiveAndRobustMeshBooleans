#include <array>
#include <vector>
#include <algorithm>
#include <cmath>
#include "tri_tri_intersect.h"

//--------------------------------------------------------------
// Funzioni ausiliarie per la 2D e 3D
//--------------------------------------------------------------

template <typename T>
inline bool points_are_coincident_2d(const T* p0, const T* p1) {
    using ElementType = typename std::remove_extent<T>::type;
    if constexpr (std::is_same<ElementType, genericPoint>::value)
        return genericPoint::coincident(*p0, *p1);
    else
        return (p0[0] == p1[0] && p0[1] == p1[1]);
}

template <typename T>
inline bool points_are_coincident_3d(const T* p0, const T* p1) {
    using ElementType = typename std::remove_extent<T>::type;
    if constexpr (std::is_same<ElementType, genericPoint>::value)
        return genericPoint::coincident(*p0, *p1);
    else
        return (p0[0] == p1[0] && p0[1] == p1[1] && p0[2] == p1[2]);
}

template <typename T>
inline bool points_are_colinear_3d(const T* p0, const T* p1, const T* p2) {
    // Verifica se le proiezioni sui tre piani sono collineari
    T proj[2];
    proj[0] = p0[1]; proj[1] = p0[2];
    T proj1[2] = { p1[1], p1[2] };
    T proj2[2] = { p2[1], p2[2] };
    if (orient2dT(proj, proj1, proj2) != 0) return false;

    proj[0] = p0[0]; proj[1] = p0[2];
    proj1[0] = p1[0]; proj1[1] = p1[2];
    proj2[0] = p2[0]; proj2[1] = p2[2];
    if (orient2dT(proj, proj1, proj2) != 0) return false;

    proj[0] = p0[0]; proj[1] = p0[1];
    proj1[0] = p1[0]; proj1[1] = p1[1];
    proj2[0] = p2[0]; proj2[1] = p2[1];
    return (orient2dT(proj, proj1, proj2) == 0);
}

template <typename T>
inline bool point_in_inner_segment(const T* p, const T* v0, const T* v1) {
    if (points_are_coincident_3d(p, v0) || points_are_coincident_3d(p, v1))
        return false;
    using ElementType = typename std::remove_extent<T>::type;
    if constexpr (std::is_same<ElementType, genericPoint>::value)
        return genericPoint::pointInInnerSegment(*p, *v0, *v1);
    else {
        if (!points_are_colinear_3d(v0, v1, p)) return false;
        // Calcola min e max per ogni coordinata una volta sola
        bool between0 = (p[0] > (v0[0] < v1[0] ? v0[0] : v1[0]) &&
                         p[0] < (v0[0] > v1[0] ? v0[0] : v1[0]));
        bool between1 = (p[1] > (v0[1] < v1[1] ? v0[1] : v1[1]) &&
                         p[1] < (v0[1] > v1[1] ? v0[1] : v1[1]));
        bool between2 = (p[2] > (v0[2] < v1[2] ? v0[2] : v1[2]) &&
                         p[2] < (v0[2] > v1[2] ? v0[2] : v1[2]));
        return (between0 || between1 || between2);
    }
}

template <typename T>
inline bool point_in_inner_triangle_2d(const T* p, const T* t0, const T* t1, const T* t2) {
    double a = orient2dT(t0, t1, p);
    double b = orient2dT(t1, t2, p);
    double c = orient2dT(t2, t0, p);
    return ((a >= 0 && b >= 0 && c >= 0) || (a <= 0 && b <= 0 && c <= 0));
}

template <typename T>
inline bool point_in_inner_triangle_3d(const T* p, const T* v0, const T* v1, const T* v2) {
    auto orient = orient3dT(p, v0, v1, v2);
    if (orient != 0) return false;
    using ElementType = typename std::remove_extent<T>::type;
    if constexpr (std::is_same<ElementType, genericPoint>::value)
        return genericPoint::pointInInnerTriangle(*p, *v0, *v1, *v2);
    else {
        if (points_are_coincident_3d(p, v0) || points_are_coincident_3d(p, v1) || points_are_coincident_3d(p, v2))
            return false;
        if (point_in_inner_segment(p, v0, v1) ||
            point_in_inner_segment(p, v1, v2) ||
            point_in_inner_segment(p, v2, v0))
            return false;
        // Verifica in 3 proiezioni: elimina X, Y, Z a turno
        T proj[2], tproj[2];
        // Proiezione su piano YZ (drop X)
        proj[0] = p[1]; proj[1] = p[2];
        tproj[0] = v0[1]; tproj[1] = v0[2];
        T tproj1[2] = { v1[1], v1[2] };
        T tproj2[2] = { v2[1], v2[2] };
        if (!point_in_inner_triangle_2d(proj, tproj, tproj1, tproj2)) return false;
        // Proiezione su piano XZ (drop Y)
        proj[0] = p[0]; proj[1] = p[2];
        tproj[0] = v0[0]; tproj[1] = v0[2];
        tproj1[0] = v1[0]; tproj1[1] = v1[2];
        tproj2[0] = v2[0]; tproj2[1] = v2[2];
        if (!point_in_inner_triangle_2d(proj, tproj, tproj1, tproj2)) return false;
        // Proiezione su piano XY (drop Z)
        proj[0] = p[0]; proj[1] = p[1];
        tproj[0] = v0[0]; tproj[1] = v0[1];
        tproj1[0] = v1[0]; tproj1[1] = v1[1];
        tproj2[0] = v2[0]; tproj2[1] = v2[1];
        if (!point_in_inner_triangle_2d(proj, tproj, tproj1, tproj2)) return false;
        return true;
    }
}

template <typename T>
bool segment_cross_inner_segment_2d (const T * s00, const T * s01, const T * s10, const T * s11){

    // https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
    double det_s00 = orient2dT(s10, s11, s00);
    double det_s01 = orient2dT(s10, s11, s01);
    double det_s10 = orient2dT(s00, s01, s10);
    double det_s11 = orient2dT(s00, s01, s11);

    // Shewchuk's orient predicates return a rough approximation of the determinant.
    // I am converting values to { -1, 0, 1 } for a simpler check of intersection cases
    int s00_wrt_s1 = (det_s00 > 0) ? 1 : ((det_s00 < 0) ? -1 : 0);
    int s01_wrt_s1 = (det_s01 > 0) ? 1 : ((det_s01 < 0) ? -1 : 0);
    int s10_wrt_s0 = (det_s10 > 0) ? 1 : ((det_s10 < 0) ? -1 : 0);
    int s11_wrt_s0 = (det_s11 > 0) ? 1 : ((det_s11 < 0) ? -1 : 0);

    if(s00_wrt_s1 != s01_wrt_s1 && s10_wrt_s0 != s11_wrt_s0)
    {
        // edges share an endpoint
        if(points_are_coincident_2d(s00, s10) || points_are_coincident_2d(s00, s11) || points_are_coincident_2d(s01, s10) || points_are_coincident_2d(s01, s11))
            return true;

        // at least one segment endpoint is involved in the intersection
        return true;
    }

    // degenerate case: colinear segments
    if(s00_wrt_s1 == 0 && s01_wrt_s1 == 0 && s10_wrt_s0 == 0 && s11_wrt_s0 == 0)
    {
        // coincident segments
        if(points_are_coincident_2d(s00, s10) && points_are_coincident_2d(s01, s11) || points_are_coincident_2d(s00, s11) && points_are_coincident_2d(s01, s10))
        {
            return true;
        }


        T Xmin_s1 = std::min(s10[0], s11[0]);
        T Xmax_s1 = std::max(s10[0], s11[0]);
        T Ymin_s1 = std::min(s10[1], s11[1]);
        T Ymax_s1 = std::max(s10[1], s11[1]);
        T Xmin_s0 = std::min(s00[0], s01[0]);
        T Xmax_s0 = std::max(s00[0], s01[0]);
        T Ymin_s0 = std::min(s00[1], s01[1]);
        T Ymax_s0 = std::max(s00[1], s01[1]);

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
            return true;
        }
    }
    return false;
}

template <typename T>
bool segment_cross_inner_segment_3d (const T * s00, const T * s01, const T * s10, const T * s11){
    //check if the segments intersect
    using ElementType = typename std::remove_extent<T>::type;
    if constexpr (std::is_same<ElementType, genericPoint >::value){
        if(orient3dT(s00,s01,s10,s11) != 0) return false;
        return genericPoint::innerSegmentsCross(*s00, *s01, *s10, *s11);
    } else {

        if(orient3dT(s00, s01, s10, s11)!= 0)
            return false;

        if (points_are_coincident_3d(s00, s10) || points_are_coincident_3d(s00, s11) || points_are_coincident_3d(s01, s10) || points_are_coincident_3d(s01, s11))
            return false;

        if (point_in_inner_segment(s00, s10, s11) || point_in_inner_segment(s01, s10, s11) ||
            point_in_inner_segment(s10, s00, s01) || point_in_inner_segment(s11, s00, s01))
        {
            return false;
        }

        // check for coincident segments
        bool s00_is_shared = (points_are_coincident_3d(s00, s10) || points_are_coincident_3d(s00, s11));
        bool s01_is_shared = (points_are_coincident_3d(s01, s10) || points_are_coincident_3d(s01, s11));
        bool s10_is_shared = (points_are_coincident_3d(s10, s00) || points_are_coincident_3d(s10, s01));
        bool s11_is_shared = (points_are_coincident_3d(s11, s00) || points_are_coincident_3d(s11, s01));

        // s0 and s1 are coincident or one edge is degenerate and coincides with one vertex of the other
        if(s00_is_shared && s01_is_shared && s10_is_shared && s11_is_shared)
            return false;

        // check 2D projections of the segments

        T s00_dropX[2] = {s00[1], s00[2]}, s01_dropX[2] = {s01[1], s01[2]};
        T s10_dropX[2] = {s10[1], s10[2]}, s11_dropX[2] = {s11[1], s11[2]};

        int x_res = segment_cross_inner_segment_2d(s00_dropX, s01_dropX, s10_dropX, s11_dropX);
        if (x_res == false)
            return false;

        T s00_dropY[2] = {s00[0], s00[2]}, s01_dropY[2] = {s01[0], s01[2]};
        T s10_dropY[2] = {s10[0], s10[2]}, s11_dropY[2] = {s11[0], s11[2]};
        int y_res = segment_cross_inner_segment_2d(s00_dropY, s01_dropY, s10_dropY, s11_dropY);
        if (y_res == false)
            return false;

        T s00_dropZ[2] = {s00[0], s00[1]}, s01_dropZ[2] = {s01[0], s01[1]};
        T s10_dropZ[2] = {s10[0], s10[1]}, s11_dropZ[2] = {s11[0], s11[1]};
        int z_res = segment_cross_inner_segment_2d(s00_dropZ, s01_dropZ, s10_dropZ, s11_dropZ);
        if (z_res == false)
            return false;


        return true;
    }
}

template <typename T>
inline bool segment_cross_inner_triangle(const T* v0, const T* v1, const T* t1, const T* t2, const T* t3) {
    if (orient3dT(v0, t1, t2, t3) == 0 || orient3dT(v1, t1, t2, t3) == 0)
        return false;
    if (points_are_coincident_3d(v0, t1) || points_are_coincident_3d(v0, t2) || points_are_coincident_3d(v0, t3) ||
        points_are_coincident_3d(v1, t1) || points_are_coincident_3d(v1, t2) || points_are_coincident_3d(v1, t3))
        return false;
    if ( point_in_inner_segment(v0, t1, t2) || point_in_inner_segment(v0, t2, t3) || point_in_inner_segment(v0, t3, t1) ||
         point_in_inner_segment(v1, t1, t2) || point_in_inner_segment(v1, t2, t3) || point_in_inner_segment(v1, t3, t1) ||
         point_in_inner_segment(t1, v0, v1) || point_in_inner_segment(t2, v0, v1) || point_in_inner_segment(t3, v0, v1) )
        return false;
    if ( point_in_inner_triangle_3d(v0, t1, t2, t3) || point_in_inner_triangle_3d(v1, t1, t2, t3) )
        return false;
    if ( segment_cross_inner_segment_3d(v0, v1, t1, t2) ||
         segment_cross_inner_segment_3d(v0, v1, t2, t3) ||
         segment_cross_inner_segment_3d(v0, v1, t3, t1) )
        return false;
    auto vol0 = orient3dT(v0, t1, t2, t3);
    auto vol1 = orient3dT(v1, t1, t2, t3);
    if (vol0 == 0 || vol1 == 0 || (vol0 > 0 && vol1 > 0) || (vol0 < 0 && vol1 < 0))
        return false;
    double vol01 = orient3dT(v0, v1, t1, t3);
    double vol12 = orient3dT(v0, v1, t2, t1);
    double vol20 = orient3dT(v0, v1, t3, t2);
    if ((vol01 > 0 && vol12 < 0) || (vol01 < 0 && vol12 > 0) ||
        (vol12 > 0 && vol20 < 0) || (vol12 < 0 && vol20 > 0) ||
        (vol20 > 0 && vol01 < 0) || (vol20 < 0 && vol01 > 0))
        return false;
    return true;
}

//--------------------------------------------------------------
// Funzione principale: intersezione tra triangoli
//--------------------------------------------------------------

template <typename T>
intersectionResult tri_intersect_tri(const T* t00, const T* t01, const T* t02,
                                       const T* t10, const T* t11, const T* t12) {
    intersectionResult intersection;
    std::array<const T*, 3> t0 = { t00, t01, t02 };
    std::array<const T*, 3> t1 = { t10, t11, t12 };

    // Calcola gli orientamenti dei vertici di un triangolo rispetto al piano dell'altro
    int o_t0_t10 = orient3dT(t00, t01, t02, t10);
    int o_t0_t11 = orient3dT(t00, t01, t02, t11);
    int o_t0_t12 = orient3dT(t00, t01, t02, t12);
    int o_t1_t00 = orient3dT(t10, t11, t12, t00);
    int o_t1_t01 = orient3dT(t10, t11, t12, t01);
    int o_t1_t02 = orient3dT(t10, t11, t12, t02);

    // Early rejection: se un triangolo è interamente da un lato del piano dell'altro
    if ((o_t0_t10 > 0 && o_t0_t11 > 0 && o_t0_t12 > 0) ||
        (o_t0_t10 < 0 && o_t0_t11 < 0 && o_t0_t12 < 0) ||
        (o_t1_t00 > 0 && o_t1_t01 > 0 && o_t1_t02 > 0) ||
        (o_t1_t00 < 0 && o_t1_t01 < 0 && o_t1_t02 < 0))
        return intersection;

    bool coplanar = (o_t0_t10 == 0 && o_t0_t11 == 0 && o_t0_t12 == 0 &&
                     o_t1_t00 == 0 && o_t1_t01 == 0 && o_t1_t02 == 0);

    bool coincidentPoints[3][3] = { {false} };
    bool t1VerticesInT0Edges[3][3] = { {false} };
    bool t0VerticesInT1Edges[3][3] = { {false} };
    bool t1VerticesInT0[3] = { false, false, false };
    bool t0VerticesInT1[3] = { false, false, false };
    bool t0EdgesCrossT1Edges[3][3] = { {false} };
    bool t0EdgesCrossT1InnerTriangle[3] = { false, false, false };
    bool t1EdgesCrossT0InnerTriangle[3] = { false, false, false };

    if (coplanar) {
        intersection.coplanar = true;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                coincidentPoints[i][j]     = points_are_coincident_3d(t0[i], t1[j]);
                t0VerticesInT1Edges[i][j]   = point_in_inner_segment(t0[i], t1[j], t1[(j+1)%3]);
                t1VerticesInT0Edges[i][j]   = point_in_inner_segment(t1[i], t0[j], t0[(j+1)%3]);
                t0EdgesCrossT1Edges[i][j]   = segment_cross_inner_segment_3d(t0[i], t0[(i+1)%3], t1[j], t1[(j+1)%3]);
            }
            t1VerticesInT0[i] = point_in_inner_triangle_3d(t1[i], t0[0], t0[1], t0[2]);
            t0VerticesInT1[i] = point_in_inner_triangle_3d(t0[i], t1[0], t1[1], t1[2]);
        }

        std::array<std::vector<int>, 6> intersectionPointEdges;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (coincidentPoints[i][j])
                    intersection.addIntersectionPoint(VA_ON_VB, i, j);
                if (t1VerticesInT0Edges[i][j])
                    intersection.addIntersectionPoint(VB_ON_EA, i, j);
                if (t0VerticesInT1Edges[i][j])
                    intersection.addIntersectionPoint(VA_ON_EB, i, j);
                if (t0EdgesCrossT1Edges[i][j])
                    intersection.addIntersectionPoint(EA_CROSS_EB, i, j);
            }
            if (t1VerticesInT0[i])
                intersection.addIntersectionPoint(VB_ON_TA, i, -1);
            if (t0VerticesInT1[i])
                intersection.addIntersectionPoint(VA_ON_TB, i, -1);
        }
        //check the edge where the intersection point is located
        for(int point = 0; point < intersection.getIntersectionPointsSize(); point ++) {
            int type = intersection.intersectionsPoints[point*3];
            int i = intersection.intersectionsPoints[point*3 + 1];
            int j = intersection.intersectionsPoints[point*3 + 2];

            switch (type) {
                case VA_ON_VB:
                    intersectionPointEdges[i].push_back(point);
                intersectionPointEdges[(i+2)%3].push_back(point);
                intersectionPointEdges[j+3].push_back(point);
                intersectionPointEdges[(j+2)%3+3].push_back(point);
                break;
                case VB_ON_EA:
                    intersectionPointEdges[i+3].push_back(point);
                intersectionPointEdges[(i+2)%3+3].push_back(point);
                intersectionPointEdges[j].push_back(point);
                break;
                case VA_ON_EB:
                    intersectionPointEdges[i].push_back(point);
                intersectionPointEdges[(i+2)%3].push_back(point);
                intersectionPointEdges[j+3].push_back(point);
                break;
                case EA_CROSS_EB:
                    intersectionPointEdges[i].push_back(point);
                intersectionPointEdges[j+3].push_back(point);
                break;
                case VB_ON_TA:
                    intersectionPointEdges[i+3].push_back(point);
                intersectionPointEdges[(i+2)%3+3].push_back(point);
                break;
                case VA_ON_TB:
                    intersectionPointEdges[i].push_back(point);
                intersectionPointEdges[(i+2)%3].push_back(point);
                break;
                default:
                    break;
            }
        }
        // Se in qualche bordo ci sono due punti, si aggiunge l'arco di intersezione
        for (int e = 0; e < 6; e++) {
            if (intersectionPointEdges[e].size() >= 2)
                intersection.addIntersectionEdge(intersectionPointEdges[e].front(), intersectionPointEdges[e].back());
        }
        return intersection;
    } else {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                coincidentPoints[i][j]     = points_are_coincident_3d(t0[i], t1[j]);
                t0VerticesInT1Edges[i][j]   = point_in_inner_segment(t0[i], t1[j], t1[(j+1)%3]);
                t1VerticesInT0Edges[i][j]   = point_in_inner_segment(t1[i], t0[j], t0[(j+1)%3]);
                t0EdgesCrossT1Edges[i][j]   = segment_cross_inner_segment_3d(t0[i], t0[(i+1)%3], t1[j], t1[(j+1)%3]);
            }
            t1VerticesInT0[i] = point_in_inner_triangle_3d(t1[i], t0[0], t0[1], t0[2]);
            t0VerticesInT1[i] = point_in_inner_triangle_3d(t0[i], t1[0], t1[1], t1[2]);
            t0EdgesCrossT1InnerTriangle[i] = segment_cross_inner_triangle(t0[i], t0[(i+1)%3], t1[0], t1[1], t1[2]);
            t1EdgesCrossT0InnerTriangle[i] = segment_cross_inner_triangle(t1[i], t1[(i+1)%3], t0[0], t0[1], t0[2]);
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (coincidentPoints[i][j])
                    intersection.addIntersectionPoint(VA_ON_VB, i, j);
                if (t1VerticesInT0Edges[i][j])
                    intersection.addIntersectionPoint(VB_ON_EA, i, j);
                if (t0VerticesInT1Edges[i][j])
                    intersection.addIntersectionPoint(VA_ON_EB, i, j);
                if (t0EdgesCrossT1Edges[i][j])
                    intersection.addIntersectionPoint(EA_CROSS_EB, i, j);
                if (intersection.getIntersectionPointsSize() == 2) {
                    intersection.addIntersectionEdge(0, 1);
                    return intersection;
                }
            }
            if (t1VerticesInT0[i])
                intersection.addIntersectionPoint(VB_ON_TA, i, -1);
            if (t0VerticesInT1[i])
                intersection.addIntersectionPoint(VA_ON_TB, i, -1);
            if (t0EdgesCrossT1InnerTriangle[i])
                intersection.addIntersectionPoint(EA_CROSS_TB, i, -1);
            if (t1EdgesCrossT0InnerTriangle[i])
                intersection.addIntersectionPoint(EB_CROSS_TA, i, -1);
            if (intersection.getIntersectionPointsSize() == 2) {
                intersection.addIntersectionEdge(0, 1);
                return intersection;
            }
        }
        return intersection;
    }
}


//--------------------------------------------------------------
// Conversioni: tri_to_double e to_genericPoint
//--------------------------------------------------------------

std::vector<double> tri_to_double(const CGAL::Gmpq* t0, const CGAL::Gmpq* t1, const CGAL::Gmpq* t2) {
    std::vector<double> res;
    res.reserve(9);
    res.push_back(CGAL::to_double(t0[0]));
    res.push_back(CGAL::to_double(t0[1]));
    res.push_back(CGAL::to_double(t0[2]));
    res.push_back(CGAL::to_double(t1[0]));
    res.push_back(CGAL::to_double(t1[1]));
    res.push_back(CGAL::to_double(t1[2]));
    res.push_back(CGAL::to_double(t2[0]));
    res.push_back(CGAL::to_double(t2[1]));
    res.push_back(CGAL::to_double(t2[2]));
    return res;
}

std::vector<double> tri_to_double(const genericPoint &t0, const genericPoint &t1, const genericPoint &t2) {
    std::vector<double> res;
    res.reserve(9);
    double x, y, z;
    t0.getApproxXYZCoordinates(x, y, z);
    res.push_back(x); res.push_back(y); res.push_back(z);
    t1.getApproxXYZCoordinates(x, y, z);
    res.push_back(x); res.push_back(y); res.push_back(z);
    t2.getApproxXYZCoordinates(x, y, z);
    res.push_back(x); res.push_back(y); res.push_back(z);
    return res;
}

template <typename T>
std::vector<genericPoint*> to_genericPoint(T* t00, T* t01, T* t02, T* t10, T* t11, T* t12,
                                             intersectionResult &intersections) {
    std::vector<genericPoint*> points;
    std::vector<genericPoint*> t0_vec, t1_vec;
    using ElementType = typename std::remove_extent<T>::type;
    if constexpr (std::is_same<ElementType, CGAL::Gmpq>::value) {
        auto tmp = tri_to_double(t00, t01, t02);
        for (int i = 0; i < 3; i++) {
            t0_vec.push_back(new explicitPoint3D(tmp[i*3], tmp[i*3+1], tmp[i*3+2]));
            t1_vec.push_back(new explicitPoint3D(tmp[i*3+9], tmp[i*3+10], tmp[i*3+11]));
        }
    } else if constexpr (std::is_same<ElementType, double>::value) {
        std::vector<const double*> tmp0 = { t00, t01, t02 };
        std::vector<const double*> tmp1 = { t10, t11, t12 };
        for (int i = 0; i < 3; i++) {
            t0_vec.push_back(new explicitPoint3D(tmp0[i][0], tmp0[i][1], tmp0[i][2]));
            t1_vec.push_back(new explicitPoint3D(tmp1[i][0], tmp1[i][1], tmp1[i][2]));
        }
    } else if constexpr (std::is_same<ElementType, genericPoint>::value) {
        t0_vec.push_back(t00);
        t0_vec.push_back(t01);
        t0_vec.push_back(t02);
        t1_vec.push_back(t10);
        t1_vec.push_back(t11);
        t1_vec.push_back(t12);
    }
    for (int i = 0; i < intersections.getIntersectionPointsSize(); i++) {
        int type = intersections.intersectionsPoints[i*3];
        int idx0 = intersections.intersectionsPoints[i*3 + 1];
        int idx1 = intersections.intersectionsPoints[i*3 + 2];
        switch (type) {
            case VA_ON_VB:
                points.push_back(t0_vec[idx0]);
                break;
            case VB_ON_EA:
                points.push_back(t1_vec[idx0]);
                break;
            case VA_ON_EB:
                points.push_back(t0_vec[idx0]);
                break;
            case EA_CROSS_EB:
                points.push_back(new implicitPoint3D_LPI(
                    t0_vec[0]->toExplicit3D(), t0_vec[1]->toExplicit3D(), t0_vec[2]->toExplicit3D(),
                    t1_vec[idx1]->toExplicit3D(), t1_vec[(idx1+1)%3]->toExplicit3D()
                ));
                break;
            case VB_ON_TA:
                points.push_back(t1_vec[idx0]);
                break;
            case VA_ON_TB:
                points.push_back(t0_vec[idx0]);
                break;
            case EA_CROSS_TB:
                points.push_back(new implicitPoint3D_LPI(
                    t0_vec[idx0]->toExplicit3D(), t0_vec[(idx0+1)%3]->toExplicit3D(),
                    t1_vec[0]->toExplicit3D(), t1_vec[1]->toExplicit3D(), t1_vec[2]->toExplicit3D()
                ));
                break;
            case EB_CROSS_TA:
                points.push_back(new implicitPoint3D_LPI(
                    t1_vec[idx0]->toExplicit3D(), t1_vec[(idx0+1)%3]->toExplicit3D(),
                    t0_vec[0]->toExplicit3D(), t0_vec[1]->toExplicit3D(), t0_vec[2]->toExplicit3D()
                ));
                break;
            default:
                ip_error("Error in the intersection type");
                break;
        }
    }
    return points;
}
