#include "../../tri_tri_intersect.h"
#include "../../printIntersections.h"
#undef Split
#include <cinolib/meshes/meshes.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/profiler.h>

using namespace cinolib;

typedef struct{
    intersectionResult intersction;
    std::vector< genericPoint *> points;
    uint tri0;
    uint tri1;
}meshIntersection;

bool hasEdgeChain(std::vector<pair<genericPoint &,genericPoint &>> edges){
    int n = edges.size();
    genericPoint *current = &edges[0].second;

    int startEdge = 0;
    int currentEdge = 0;
    int i = 0;

    do{
        for(int j = 0; j<n ; j++){
            if(currentEdge != j){
                if(genericPoint::coincident(*current, edges[j].first)){
                    current = &edges[j].second;
                    currentEdge = j;

                    break;
                }else if(genericPoint::coincident(*current, edges[j].second)){

                    current = &edges[j].first;
                    currentEdge = j;

                    break;
                }
            }
        }
        i++;
    }while(currentEdge != startEdge && i < n);

    return currentEdge == startEdge;

}

int main() {

    GLcanvas gui(1920,1080);

    //load the mesh
    DrawableTrimesh<> s = DrawableTrimesh("../../../models/two_spheres.stl");
    SurfaceMeshControls<DrawableTrimesh<>> menu(&s, &gui);
    gui.push(&s);
    gui.push(&menu);

    //object to draw the intersections
    DrawableSegmentSoup segments;
    std::vector<Marker> markers;

    std::vector<meshIntersection> global_intersections;
    vector<pair<genericPoint &, genericPoint&>> intersectionEdges;

    std::vector<std::vector<uint>> &tris = s.vector_polys();
    for(int i = 0; i < s.num_polys(); i++){

        std::vector<uint> adj_verts = s.adj_p2v(i);

        for (int j = i+1 ; j < s.num_polys(); j++ ) {

            //check if the two polys are adjacent
            if (s.polys_are_adjacent(i, j)) {
                continue;
            }
            bool poly_adj_vert = false;
            for (auto v: adj_verts) {
                if (s.poly_contains_vert(j, v)) {
                    poly_adj_vert = true;
                    break;
                }
            }
            if (poly_adj_vert) {
                continue;
            }

            //check intersections
            std::vector<uint> t0 = tris[i];
            std::vector<uint> t1 = tris[j];

            intersectionResult result = tri_intersect_tri(s.vert(t0[0])._vec, s.vert(t0[1])._vec, s.vert(t0[2])._vec,
                                                          s.vert(t1[0])._vec, s.vert(t1[1])._vec, s.vert(t1[2])._vec);
            if (result.intersect) {
                meshIntersection intersection;
                intersection.intersction = result;
                intersection.tri0 = i;
                intersection.tri1 = j;

                std::vector<genericPoint *> convertedPoints = to_genericPoint(s.vert(t0[0])._vec, s.vert(t0[1])._vec, s.vert(t0[2])._vec,
                                    s.vert(t1[0])._vec,s.vert(t1[1])._vec, s.vert(t1[2])._vec,
                                    result);

                intersection.points = convertedPoints;
                global_intersections.push_back(intersection);

            }
        }

    }


    //populate the intersection edges
    for(auto g: global_intersections){
        if(g.intersction.getIntersectionEdgesSize() != 0){
            for (int e = 0; e < g.intersction.getIntersectionEdgesSize(); e ++) {
                intersectionEdges.push_back({*g.points[g.intersction.intersectionsEdges[e * 2]], *g.points[g.intersction.intersectionsEdges[e * 2 + 1]]});
            }
        }
    }

    cout << "Number of intersections: " << global_intersections.size() << endl;
    cout << "Number of intersection edges: " << intersectionEdges.size() << endl;

    bool hasChain = hasEdgeChain(intersectionEdges);
    cout << "Has edge chain: " << hasChain << endl;

    for(auto intersection : global_intersections){
       drawIntersections(s, segments, gui, intersection.tri0, intersection.tri1, intersection.intersction);
    }

    gui.push(&segments);
    gui.launch();

}
