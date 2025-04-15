
#ifndef TRIANGLE_INTERSECTIONS_PRINTINTERSECTIONS_H
#define TRIANGLE_INTERSECTIONS_PRINTINTERSECTIONS_H
#undef Split
#include <cinolib/drawable_triangle_soup.h>
#include <cinolib/drawable_segment_soup.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/geometry/segment_utils.h>
#include <cinolib/meshes/meshes.h>
#include "tri_tri_intersect.h"
typedef struct{
    std::vector<cinolib::Marker> *markers;
    cinolib::GLcanvas *gui;
    cinolib::DrawableTriangleSoup *triangles;
    cinolib::DrawableSegmentSoup *segments;

    bool viewMarkers = true;
    bool viewTriangles = true;
    bool viewSegments = true;

} printableIntersections;

void printTrianglesIntersections( const double * t00, const double * t01, const double * t02, const double * t10, const double * t11, const double * t12,  intersectionResult &result, printableIntersections &printableIntersections);
void printTrianglesIntersections(const std::vector<double> &t0, const std::vector<double> &t1, intersectionResult &intersections, printableIntersections &printableIntersections);

void drawIntersections(cinolib::DrawableTrimesh<> &m ,cinolib::DrawableSegmentSoup &s ,cinolib::GLcanvas &gui, uint t0, uint t1, intersectionResult &intersection );


#endif //TRIANGLE_INTERSECTIONS_PRINTINTERSECTIONS_H
