#include <iostream>
#undef Split

#include "tri_tri_intersect.h"
#include "printIntersections.h"
#include <implicit_point.h>
#include <cinolib/predicates.h>
#include <cinolib/drawable_triangle_soup.h>


using namespace cinolib;


int main(int argc, char **argv) {


    initFPU();
    typedef enum{
        DOUBLE = 0,
        GMPQ = 1,
        IMPLICIT = 2
    }pointType;

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    //point definition
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    pointType pType = IMPLICIT;


    //DOUBLE
    double t1_0[3] = {6927969884.831744 ,5707912752.136192, 2819502416.855040};
    double t1_1[3] = {9218755653.533695, 5707912752.136192 ,1927106460.123136};
    double t1_2[3] = {6927969884.831744, 5707912752.136192 ,1921459114.999808};

    double t2_0[3] = {6927969884.831744 ,5707912752.136192, 1921459114.999808};
    double t2_1[3] = {9218755653.533695 ,5707912752.136192 ,2819502416.855040};
    double t2_2[3] = {6927969347.960832, 5707912752.136192 ,2520154739.048448};


    //DEBUG
    array<double*,3> dt0 = {t1_0, t1_1, t1_2};
    array<double*,3> dt1 = {t2_0, t2_1, t2_2};



    //GMPQ
    CGAL::Gmpq gt1_0[3] = {CGAL::Gmpq(CGAL::Gmpz(-8), CGAL::Gmpz(1)),
                      CGAL::Gmpq(CGAL::Gmpz(-37), CGAL::Gmpz(1)),
                      CGAL::Gmpq(CGAL::Gmpz(-83886085938), CGAL::Gmpz(1000000))};

    CGAL::Gmpq gt1_1[3] = {CGAL::Gmpq(CGAL::Gmpz(-8015147520), CGAL::Gmpz(1)),
                          CGAL::Gmpq(CGAL::Gmpz(-377789376), CGAL::Gmpz(1)),
                          CGAL::Gmpq(CGAL::Gmpz(1610612750000), CGAL::Gmpz(1000000))};

    CGAL::Gmpq gt1_2[3] = {CGAL::Gmpq(CGAL::Gmpz(-8015164416), CGAL::Gmpz(1)),
                          CGAL::Gmpq(CGAL::Gmpz(-377789376), CGAL::Gmpz(1)),
                          CGAL::Gmpq(CGAL::Gmpz(1644167250000), CGAL::Gmpz(1000000))};

    CGAL::Gmpq gt2_0[3] = {CGAL::Gmpq(CGAL::Gmpz(-8018184192), CGAL::Gmpz(1)),
                          CGAL::Gmpq(CGAL::Gmpz(-377755808), CGAL::Gmpz(1)),
                          CGAL::Gmpq(CGAL::Gmpz(-83886085938), CGAL::Gmpz(1000000))};

    CGAL::Gmpq gt2_1[3] = {CGAL::Gmpq(CGAL::Gmpz(-8021422080), CGAL::Gmpz(1)),
                          CGAL::Gmpq(CGAL::Gmpz(-377789376), CGAL::Gmpz(1)),
                          CGAL::Gmpq(CGAL::Gmpz(1644167250000), CGAL::Gmpz(1000000))};

    CGAL::Gmpq gt2_2[3] = {CGAL::Gmpq(CGAL::Gmpz(-8021438976), CGAL::Gmpz(1)),
                          CGAL::Gmpq(CGAL::Gmpz(-377789376), CGAL::Gmpz(1)),
                          CGAL::Gmpq(CGAL::Gmpz(1610612750000), CGAL::Gmpz(1000000))};




    //IMPLICIT

    const genericPoint * it1_0 = new explicitPoint3D(6927969884.831744, 5707912752.136192, 2819502416.855040);
    const genericPoint * it1_1 = new explicitPoint3D(9218755653.533695, 5707912752.136192, 1927106460.123136);
    const genericPoint * it1_2 = new explicitPoint3D(6927969884.831744, 5707912752.136192, 1921459114.999808 );

    const genericPoint * it2_0 = new explicitPoint3D(8646059077.140480, 5707912752.136192, 1925694758.060032);
    const genericPoint * it2_1 = new explicitPoint3D(8646059077.140480, 5707912752.136192, 601858062.155776);
    const genericPoint * it2_2 = new explicitPoint3D(9218755653.533695, 5707912752.136192, 1927106460.123136);

    //check user input
    if (argc > 18) {
        cout << "caricamento parametri da riga di comando" << endl;
        t1_0[0] = std::stod(argv[1]);
        t1_0[1] = std::stod(argv[2]);
        t1_0[2] = std::stod(argv[3]);

        t1_1[0] = std::stod(argv[4]);
        t1_1[1] = std::stod(argv[5]);
        t1_1[2] = std::stod(argv[6]);

        t1_2[0] = std::stod(argv[7]);
        t1_2[1] = std::stod(argv[8]);
        t1_2[2] = std::stod(argv[9]);

        t2_0[0] = std::stod(argv[10]);
        t2_0[1] = std::stod(argv[11]);
        t2_0[2] = std::stod(argv[12]);

        t2_1[0] = std::stod(argv[13]);
        t2_1[1] = std::stod(argv[14]);
        t2_1[2] = std::stod(argv[15]);

        t2_2[0] = std::stod(argv[16]);
        t2_2[1] = std::stod(argv[17]);
        t2_2[2] = std::stod(argv[18]);

        pType = DOUBLE;


    }

    bigrational b00[3] = {bigrational(1,2,1),
                            bigrational(),
                            bigrational()};
    bigrational b01[3] = {bigrational(1,1,1),
                            bigrational(),
                            bigrational()};
    bigrational b02[3] = {bigrational(1,2,1),
                            bigrational(1,1,1),
                            bigrational()};

    bigrational p[3] = {bigrational(1,2,1),
                        bigrational(1,2,1),
                        bigrational(1,2,-1)};

    int tmp = orient3dT(b00, b01, b02, p);

    cout << "Orient: " << tmp << endl;



//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//triangles definition
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    std::vector<double> t0;
    std::vector<double> t1;

    switch (pType) {
        case DOUBLE:
            t0 = {t1_0[0], t1_0[1], t1_0[2], t1_1[0], t1_1[1], t1_1[2], t1_2[0], t1_2[1], t1_2[2]};
            t1 = {t2_0[0], t2_0[1], t2_0[2], t2_1[0], t2_1[1], t2_1[2], t2_2[0], t2_2[1], t2_2[2]};
            break;
        case GMPQ:
            t0 = tri_to_double(gt1_0, gt1_1, gt1_2);
            t1 = tri_to_double(gt2_0, gt2_1, gt2_2);
            break;
        case IMPLICIT:
            t0 = tri_to_double(*it1_0, *it1_1, *it1_2);
            t1 = tri_to_double(*it2_0, *it2_1, *it2_2);
            break;
    }








    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    //tri_tri_intersect
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    auto start = std::chrono::high_resolution_clock::now();
    intersectionResult result;
    switch(pType){
        case DOUBLE:
            result = tri_intersect_tri(t1_0, t1_1, t1_2, t2_0, t2_1, t2_2);
            break;
        case GMPQ:
            result = tri_intersect_tri(gt1_0, gt1_1, gt1_2, gt2_0, gt2_1, gt2_2);
            break;
        case IMPLICIT:
            result = tri_intersect_tri(it1_0, it1_1, it1_2, it2_0, it2_1, it2_2);
            break;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Time: " << duration.count() << " microseconds" << std::endl;

    result.printIntersectionPoints();
    result.printIntersectionEdges();

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    //printIntersections
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    //gui
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    std::vector<double> coords = {t0[0], t0[1], t0[2],
                                  t0[3], t0[4], t0[5],
                                  t0[6], t0[7], t0[8],
                                  t1[0], t1[1], t1[2],
                                  t1[3], t1[4], t1[5],
                                  t1[6], t1[7], t1[8]};

    std::vector<uint> tris = {0, 1, 2, 3, 4, 5};

    std::vector<cinolib::Color> colors = {cinolib::Color::PASTEL_CYAN(), cinolib::Color::PASTEL_ORANGE()};
    GLcanvas gui(1080,1080);

    DrawableSegmentSoup s;
    s.thickness = 4;
    DrawableTriangleSoup t(coords,tris,colors);
    std::vector<Marker> markers;

    write_OBJ("output.obj",coords,{{0,1,2},{3,4,5}});

    printableIntersections tri_gui = {&markers, &gui, &t, &s};

    float size = 0.1f;

    gui.callback_app_controls = [&]()
    {
        if(ImGui::Checkbox("View triangles", &tri_gui.viewTriangles)){
            if(tri_gui.viewTriangles) gui.push(tri_gui.triangles);
            else gui.pop(tri_gui.triangles);
        }

        if(ImGui::Checkbox("View segments", &tri_gui.viewSegments)){
            if(tri_gui.viewSegments) gui.push(tri_gui.segments);
            else gui.pop(tri_gui.segments);
        }

        if(ImGui::Checkbox("View markers", &tri_gui.viewMarkers)){
            if(tri_gui.viewMarkers) gui.markers = *tri_gui.markers;
            else gui.markers.clear();
        }

        if (ImGui::InputFloat("Marker size", &size)) {}

        if (ImGui::Button("Aggiorna dimensione Marker")) {
            for (auto &m : tri_gui.gui->markers) {
                m.disk_radius = size;
                m.font_size = size;
            }
        }

        ImGui::Text("Intersect: %d", result.intersect);

        if (result.coplanar) ImGui::Text("Coplanar: true"); else ImGui::Text("Coplanar: false");

        ImGui::Text("Intersection Points: %d", result.getIntersectionPointsSize());
        ImGui::Text("Intersection Edges: %d", result.getIntersectionEdgesSize());
        ImGui::Text("Computation Time: %d microseconds", duration.count());





    };

    printTrianglesIntersections(t0,t1,result,tri_gui);

    return gui.launch();
}

