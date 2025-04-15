#ifndef TRIANGLE_INTERSECTIONS_TRI_TRI_INTERSECT_H
#define TRIANGLE_INTERSECTIONS_TRI_TRI_INTERSECT_H

#include "orientTemplated.h"

typedef enum {
   VA_ON_VB = 0,
   VA_ON_EB = 1,
   VB_ON_EA = 2,
   VA_ON_TB = 3,
   VB_ON_TA = 4,
   EA_CROSS_EB = 5,
   EA_CROSS_TB = 6,
   EB_CROSS_TA = 7,
} intersectionType;

class intersectionResult {
public:
    std::vector<int> intersectionsPoints = {};
    std::vector<int> intersectionsEdges = {};

    bool intersect = false;
    bool coplanar = false;

    void addIntersectionPoint(intersectionType type, int index0, int index1){
        intersectionsPoints.push_back(type);
        intersectionsPoints.push_back(index0);
        intersectionsPoints.push_back(index1);

        intersect = true;
    };

    void addIntersectionEdge(int index0, int index1){
        //check if edge already exists
        for (int i = 0; i < intersectionsEdges.size(); i+=2){
            if ((intersectionsEdges[i] == index0 && intersectionsEdges[i+1] == index1) ||
                (intersectionsEdges[i] == index1 && intersectionsEdges[i+1] == index0)){
                return;
            }
        }

        intersectionsEdges.push_back(index0);
        intersectionsEdges.push_back(index1);
    };

    void printIntersectionPoints(){
        for (int i = 0; i < intersectionsPoints.size(); i+=3){
            std::cout << "Intersection type: " << parseIntersectionType((intersectionType)intersectionsPoints[i]) << ": " << intersectionsPoints[i+1] << " and " << intersectionsPoints[i+2] << std::endl;
        }
    };

    void printIntersectionEdges(){
        for (int i = 0; i < intersectionsEdges.size(); i+=2){
            std::cout << "Intersection edge: " << intersectionsEdges[i] << " and " << intersectionsEdges[i+1] << std::endl;
        }
    };

    string parseIntersectionType(intersectionType t){
        switch (t) {
            case VA_ON_VB: return "VA_ON_VB";
            case VA_ON_EB: return "VA_ON_EB";
            case VB_ON_EA: return "VB_ON_EA";
            case VA_ON_TB: return "VA_ON_TB";
            case VB_ON_TA: return "VB_ON_TA";
            case EA_CROSS_EB: return "EA_CROSS_EB";
            case EA_CROSS_TB: return "EA_CROSS_TB";
            case EB_CROSS_TA: return "EB_CROSS_TA";
        }
        return "Unknown type";
    }

    const int getIntersectionPointsSize(){
        return intersectionsPoints.size()/3;
    }

    const int getIntersectionEdgesSize(){
        return intersectionsEdges.size()/2;
    }


};


template <typename T>
bool points_are_coincident_3d(const T * p0, const T * p1);

template <typename T>
bool points_are_coincident_2d(const T * p0, const T * p1);

template <typename T>
bool points_are_colinear_3d(const T * p0, const T * p1, const T * p2);

template <typename T>
bool point_in_inner_segment(const T * p, const T * a, const T * b);

template <typename T>
bool point_in_inner_triangle_2d(const T * p,const T * t0,const T * t1,const T * t2);

template <typename T>
bool point_in_inner_triangle_3d(const T * p, const T * v0, const T * v1, const T * v2);

template <typename T>
bool segment_cross_inner_segment_2d (const T * a0, const T * a1, const T * b0, const T * b1);

template <typename T>
bool segment_cross_inner_segment_3d (const T * a0, const T * a1, const T * b0, const T * b1);

template <typename T>
bool segment_cross_inner_triangle(const T * v1, const T * v2, const T * t1, const T * t2, const T * t3);

std::pair<int, int> getOppositeEdge(int vertexIndex);

template <typename T>
intersectionResult tri_intersect_tri(const T * t00, const T * t01, const T * t02, const T * t10, const T * t11, const T * t12);

inline
std::vector<double> tri_to_double(const CGAL::Gmpq * t0, const CGAL::Gmpq * t1, const CGAL::Gmpq * t2);

inline
std::vector<double> tri_to_double(const genericPoint & t0, const genericPoint & t1, const genericPoint & t2);

template <typename T>
std::vector<genericPoint *> to_genericPoint( T* t00,  T* t01,  T *t02,  T * t10,  T * t11,  T * t12,  intersectionResult &intersections ) ;

#include "tri_tri_intersect.cpp"

#endif //TRIANGLE_INTERSECTIONS_TRI_TRI_INTERSECT_H
