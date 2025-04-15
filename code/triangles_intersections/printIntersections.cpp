#include "printIntersections.h"


/**
 * Print the triangles, the intersection points and the intersection edges on the GUI
 * @param t00 
 * @param t01
 * @param t02
 * @param t10
 * @param t11
 * @param t12
 * @param intersections - the intersections result struct
 * @param printableIntersections - the struct containing the drawable objects and the GUI
 */
void printTrianglesIntersections(const double * t00, const double * t01, const double * t02, const double * t10, const double * t11, const double * t12, intersectionResult &intersections, printableIntersections &printableIntersections){
    std::array<const double*, 3> t0 = {t00, t01, t02};
    std::array<const double*, 3> t1 = {t10, t11, t12};

    int n = intersections.getIntersectionPointsSize();

    //calculate the intersection points
    std::vector<std::vector<double>> explicitIntersectionPoints;

    for(int i = 0; i < n; i ++){
        switch (intersections.intersectionsPoints[i * 3]){
            case VA_ON_VB:{
                int index = intersections.intersectionsPoints[i * 3 + 1];
                std::vector<double> tmp = {t0[index][0], t0[index][1], t0[index][2]};
                explicitIntersectionPoints.push_back(tmp);
                break;}
            case VA_ON_EB:{
                int index = intersections.intersectionsPoints[i * 3 + 1];
                std::vector<double> tmp = {t0[index][0], t0[index][1], t0[index][2]};
                explicitIntersectionPoints.push_back(tmp);
                break;}
            case VB_ON_EA:{
                int index = intersections.intersectionsPoints[i * 3 + 1];
                std::vector<double> tmp = {t1[index][0], t1[index][1], t1[index][2]};
                explicitIntersectionPoints.push_back(tmp);
                break;
            }
            case VA_ON_TB:{
                int index = intersections.intersectionsPoints[i * 3 + 1];
                std::vector<double> tmp = {t0[index][0], t0[index][1], t0[index][2]};
                explicitIntersectionPoints.push_back(tmp);
                break;}
            case VB_ON_TA:{
                int index = intersections.intersectionsPoints[i * 3 + 1];
                std::vector<double> tmp = {t1[index][0], t1[index][1], t1[index][2]};
                explicitIntersectionPoints.push_back(tmp);
                break;}
            case EA_CROSS_EB:{
                int p0 = intersections.intersectionsPoints[i * 3 + 1];
                int p1 = intersections.intersectionsPoints[i * 3 + 2];

                cinolib::vec3d v00 = cinolib::vec3d(t0[p0][0], t0[p0][1], t0[p0][2]);
                cinolib::vec3d v01 = cinolib::vec3d(t0[(p0+1)%3][0], t0[(p0+1)%3][1], t0[(p0+1)%3][2]);
                cinolib::vec3d v10 = cinolib::vec3d(t1[p1][0], t1[p1][1], t1[p1][2]);
                cinolib::vec3d v11 = cinolib::vec3d(t1[(p1+1)%3][0], t1[(p1+1)%3][1], t1[(p1+1)%3][2]);
                cinolib::vec3d intersectionPoint = cinolib::segment_intersection(v00, v01, v10, v11);

                std::vector<double> tmp = {intersectionPoint[0], intersectionPoint[1], intersectionPoint[2]};
                explicitIntersectionPoints.push_back(tmp);
                break;
            }

            case EA_CROSS_TB:{
                int e_index = intersections.intersectionsPoints[i * 3 + 1];
                double v0[3] = {t0[e_index][0], t0[e_index][1], t0[e_index][2]};
                double v1[3] = {t0[(e_index+1)%3][0], t0[(e_index+1)%3][1], t0[(e_index+1)%3][2]};


                // plane vectors
                std::array<double, 3> vec1 = {t1[1][0] - t1[0][0], t1[1][1] - t1[0][1], t1[1][2] - t1[0][2]};
                std::array<double, 3> vec2 = {t1[2][0] - t1[0][0], t1[2][1] - t1[0][1], t1[2][2] - t1[0][2]};

                // plane normal
                std::array<double, 3> normal = {
                        vec1[1] * vec2[2] - vec1[2] * vec2[1],
                        vec1[2] * vec2[0] - vec1[0] * vec2[2],
                        vec1[0] * vec2[1] - vec1[1] * vec2[0]
                };

                // plane d
                double d = -(normal[0] * t1[0][0] + normal[1] * t1[0][1] + normal[2] * t1[0][2]);

                // segment direction
                std::array<double, 3> dir = { v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};

                // calculate t param
                double num = -(normal[0] * v0[0] + normal[1] * v0[1] + normal[2] * v0[2] + d);
                double den = normal[0] * dir[0] + normal[1] * dir[1] + normal[2] * dir[2];
                double t = num / den;

                // calculate intersection point
                std::vector<double> p = {
                        v0[0] + t * dir[0],
                        v0[1] + t * dir[1],
                        v0[2] + t * dir[2]
                };
                explicitIntersectionPoints.push_back(p);
                break;}
            case EB_CROSS_TA: {
                int e_index = intersections.intersectionsPoints[i * 3 + 1];
                double v0[3] = {t1[e_index][0], t1[e_index][1], t1[e_index][2]};
                double v1[3] = {t1[(e_index+1)%3][0], t1[(e_index+1)%3][1], t1[(e_index+1)%3][2]};


                // plane vectors
                std::array<double, 3> tmp1 = {t0[1][0] - t0[0][0], t0[1][1] - t0[0][1], t0[1][2] - t0[0][2]};
                std::array<double, 3> tmp2 = {t0[2][0] - t0[0][0], t0[2][1] - t0[0][1], t0[2][2] - t0[0][2]};

                // plane normal
                std::array<double, 3> normal = {
                        tmp1[1] * tmp2[2] - tmp1[2] * tmp2[1],
                        tmp1[2] * tmp2[0] - tmp1[0] * tmp2[2],
                        tmp1[0] * tmp2[1] - tmp1[1] * tmp2[0]
                };

                // d plane
                double d = -(normal[0] * t0[0][0] + normal[1] * t0[0][1] + normal[2] * t0[0][2]);

                // segment direction
                std::array<double, 3> dir = { v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};

                // calculate t param
                double num = -(normal[0] * v0[0] + normal[1] * v0[1] + normal[2] * v0[2] + d);
                double den = normal[0] * dir[0] + normal[1] * dir[1] + normal[2] * dir[2];
                double t = num / den;

                //intersection point
                std::vector<double> p = {
                        v0[0] + t * dir[0],
                        v0[1] + t * dir[1],
                        v0[2] + t * dir[2]
                };
                explicitIntersectionPoints.push_back(p);
                break;
            }
        }

    }

    //print the triangles on the GUI
    for(int i = 0; i < 3 ; i++){
        cinolib::Marker m0,m1;
        m0.pos_3d = cinolib::vec3d(t0[i][0], t0[i][1], t0[i][2]);
        m0.text = "t0" + std::to_string(i);
        printableIntersections.markers->push_back(m0);

        m1.pos_3d = cinolib::vec3d(t1[i][0], t1[i][1], t1[i][2]);
        m1.text = "t1" + std::to_string(i);
        printableIntersections.markers->push_back(m1);
    }
    printableIntersections.gui->push(printableIntersections.triangles);


    //print the intersection points
    n = intersections.getIntersectionPointsSize();
    for(int i = 0; i < n; i++){
        cinolib::Marker m;
        m.pos_3d = cinolib::vec3d(explicitIntersectionPoints[i][0], explicitIntersectionPoints[i][1], explicitIntersectionPoints[i][2]);
        m.text = "p" + std::to_string(i);
        m.color = cinolib::Color::RED();
        printableIntersections.markers->push_back(m);
    }
    printableIntersections.gui->markers = *printableIntersections.markers;

    //print the intersection edges
    n = (int)intersections.intersectionsEdges.size()/2;
    for(int e =0; e < n; e++){
        int v0_index = intersections.intersectionsEdges[e * 2];
        int v1_index = intersections.intersectionsEdges[e * 2 + 1];

        cinolib::vec3d v0 = cinolib::vec3d(explicitIntersectionPoints[v0_index][0], explicitIntersectionPoints[v0_index][1], explicitIntersectionPoints[v0_index][2]);
        cinolib::vec3d v1 = cinolib::vec3d(explicitIntersectionPoints[v1_index][0], explicitIntersectionPoints[v1_index][1], explicitIntersectionPoints[v1_index][2]);

        printableIntersections.segments->push_seg(v0, v1);
    }
    printableIntersections.gui->push(printableIntersections.segments);

}

void printTrianglesIntersections(const std::vector<double> &t0, const std::vector<double> &t1, intersectionResult &intersections, printableIntersections &printableIntersections){
    printTrianglesIntersections(t0.data(), t0.data() + 3, t0.data() + 6, t1.data(), t1.data() + 3, t1.data() + 6, intersections, printableIntersections);
}



void drawIntersections(cinolib::DrawableTrimesh<> &m ,cinolib::DrawableSegmentSoup &s ,cinolib::GLcanvas &gui, uint t0, uint t1, intersectionResult &intersection ){
    using namespace cinolib;
    std::vector<vec3d> t0_verts = m.poly_verts(t0);
    std::vector<vec3d> t1_verts = m.poly_verts(t1);

    int n = intersection.getIntersectionPointsSize();

    //calculate the intersection points
    std::vector<std::vector<double>> explicitIntersectionPoints;

    for(int i = 0; i < n; i ++){
        switch (intersection.intersectionsPoints[i * 3]){
            case VA_ON_VB:{
                int index = intersection.intersectionsPoints[i * 3 + 1];
                std::vector<double> tmp = {t0_verts[index][0], t0_verts[index][1], t0_verts[index][2]};
                explicitIntersectionPoints.push_back(tmp);
                break;}
            case VA_ON_EB:{
                int index = intersection.intersectionsPoints[i * 3 + 1];
                std::vector<double> tmp = {t0_verts[index][0], t0_verts[index][1], t0_verts[index][2]};
                explicitIntersectionPoints.push_back(tmp);
                break;}
            case VB_ON_EA:{
                int index = intersection.intersectionsPoints[i * 3 + 1];
                std::vector<double> tmp = {t1_verts[index][0], t1_verts[index][1], t1_verts[index][2]};
                explicitIntersectionPoints.push_back(tmp);
                break;
            }
            case VA_ON_TB:{
                int index = intersection.intersectionsPoints[i * 3 + 1];
                std::vector<double> tmp = {t0_verts[index][0], t0_verts[index][1], t0_verts[index][2]};
                explicitIntersectionPoints.push_back(tmp);
                break;}
            case VB_ON_TA:{
                int index = intersection.intersectionsPoints[i * 3 + 1];
                std::vector<double> tmp = {t1_verts[index][0], t1_verts[index][1], t1_verts[index][2]};
                explicitIntersectionPoints.push_back(tmp);
                break;}
            case EA_CROSS_EB:{
                int p0 = intersection.intersectionsPoints[i * 3 + 1];
                int p1 = intersection.intersectionsPoints[i * 3 + 2];

                cinolib::vec3d v00 = cinolib::vec3d(t0_verts[p0][0], t0_verts[p0][1], t0_verts[p0][2]);
                cinolib::vec3d v01 = cinolib::vec3d(t0_verts[(p0+1)%3][0], t0_verts[(p0+1)%3][1], t0_verts[(p0+1)%3][2]);
                cinolib::vec3d v10 = cinolib::vec3d(t1_verts[p1][0], t1_verts[p1][1], t1_verts[p1][2]);
                cinolib::vec3d v11 = cinolib::vec3d(t1_verts[(p1+1)%3][0], t1_verts[(p1+1)%3][1], t1_verts[(p1+1)%3][2]);
                cinolib::vec3d intersectionPoint = cinolib::segment_intersection(v00, v01, v10, v11);

                std::vector<double> tmp = {intersectionPoint[0], intersectionPoint[1], intersectionPoint[2]};
                explicitIntersectionPoints.push_back(tmp);
                break;
            }

            case EA_CROSS_TB:{
                int e_index = intersection.intersectionsPoints[i * 3 + 1];
                double v0[3] = {t0_verts[e_index][0], t0_verts[e_index][1], t0_verts[e_index][2]};
                double v1[3] = {t0_verts[(e_index+1)%3][0], t0_verts[(e_index+1)%3][1], t0_verts[(e_index+1)%3][2]};


                // plane vectors
                std::array<double, 3> vec1 = {t1_verts[1][0] - t1_verts[0][0], t1_verts[1][1] - t1_verts[0][1], t1_verts[1][2] - t1_verts[0][2]};
                std::array<double, 3> vec2 = {t1_verts[2][0] - t1_verts[0][0], t1_verts[2][1] - t1_verts[0][1], t1_verts[2][2] - t1_verts[0][2]};

                // plane normal
                std::array<double, 3> normal = {
                        vec1[1] * vec2[2] - vec1[2] * vec2[1],
                        vec1[2] * vec2[0] - vec1[0] * vec2[2],
                        vec1[0] * vec2[1] - vec1[1] * vec2[0]
                };

                // plane d
                double d = -(normal[0] * t1_verts[0][0] + normal[1] * t1_verts[0][1] + normal[2] * t1_verts[0][2]);

                // segment direction
                std::array<double, 3> dir = { v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};

                // calculate t param
                double num = -(normal[0] * v0[0] + normal[1] * v0[1] + normal[2] * v0[2] + d);
                double den = normal[0] * dir[0] + normal[1] * dir[1] + normal[2] * dir[2];
                double t = num / den;

                // calculate intersection point
                std::vector<double> p = {
                        v0[0] + t * dir[0],
                        v0[1] + t * dir[1],
                        v0[2] + t * dir[2]
                };
                explicitIntersectionPoints.push_back(p);
                break;}
            case EB_CROSS_TA: {
                int e_index = intersection.intersectionsPoints[i * 3 + 1];
                double v0[3] = {t1_verts[e_index][0], t1_verts[e_index][1], t1_verts[e_index][2]};
                double v1[3] = {t1_verts[(e_index+1)%3][0], t1_verts[(e_index+1)%3][1], t1_verts[(e_index+1)%3][2]};


                // plane vectors
                std::array<double, 3> tmp1 = {t0_verts[1][0] - t0_verts[0][0], t0_verts[1][1] - t0_verts[0][1], t0_verts[1][2] - t0_verts[0][2]};
                std::array<double, 3> tmp2 = {t0_verts[2][0] - t0_verts[0][0], t0_verts[2][1] - t0_verts[0][1], t0_verts[2][2] - t0_verts[0][2]};

                // plane normal
                std::array<double, 3> normal = {
                        tmp1[1] * tmp2[2] - tmp1[2] * tmp2[1],
                        tmp1[2] * tmp2[0] - tmp1[0] * tmp2[2],
                        tmp1[0] * tmp2[1] - tmp1[1] * tmp2[0]
                };

                // d plane
                double d = -(normal[0] * t0_verts[0][0] + normal[1] * t0_verts[0][1] + normal[2] * t0_verts[0][2]);

                // segment direction
                std::array<double, 3> dir = { v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};

                // calculate t param
                double num = -(normal[0] * v0[0] + normal[1] * v0[1] + normal[2] * v0[2] + d);
                double den = normal[0] * dir[0] + normal[1] * dir[1] + normal[2] * dir[2];
                double t = num / den;

                //intersection point
                std::vector<double> p = {
                        v0[0] + t * dir[0],
                        v0[1] + t * dir[1],
                        v0[2] + t * dir[2]
                };
                explicitIntersectionPoints.push_back(p);
                break;
            }
        }



    }
    n = intersection.getIntersectionPointsSize();
    for(int i = 0; i < n; i++){
        cinolib::Marker marker;
        marker.pos_3d = cinolib::vec3d(explicitIntersectionPoints[i][0], explicitIntersectionPoints[i][1], explicitIntersectionPoints[i][2]);;
        marker.color = cinolib::Color::RED();
        gui.markers.push_back(marker);
    }
    //print the intersection edges
    n = (int)intersection.intersectionsEdges.size()/2;
    for(int e =0; e < n; e++){
        int v0_index = intersection.intersectionsEdges[e * 2];
        int v1_index = intersection.intersectionsEdges[e * 2 + 1];

        cinolib::vec3d v0 = cinolib::vec3d(explicitIntersectionPoints[v0_index][0], explicitIntersectionPoints[v0_index][1], explicitIntersectionPoints[v0_index][2]);
        cinolib::vec3d v1 = cinolib::vec3d(explicitIntersectionPoints[v1_index][0], explicitIntersectionPoints[v1_index][1], explicitIntersectionPoints[v1_index][2]);

        s.push_seg(v0, v1);
    }




}
