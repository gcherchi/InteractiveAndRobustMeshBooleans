//
// Created by Michele on 06/05/25.
//
#include <iostream>

#include <cinolib/meshes/meshes.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <booleans.h>

using namespace cinolib;
using namespace std;

bool input_mesh_flag = false;
bool is_pushed = false;
int main(int argc, char **argv){

    GLcanvas gui;
    //DrawableTrimesh<> bool_mesh ("output_union.obj");

    std::cout << argv[1] <<std::endl;
    DrawableTrimesh<> bool_mesh (argv[1]);
    SurfaceMeshControls<DrawableTrimesh<>> contr (&bool_mesh, &gui, "bool_mesh");
    DrawableSegmentSoup ray;

    DrawableTrimesh<> input_mesh ("mesh_input.obj");
    SurfaceMeshControls<DrawableTrimesh<>> contr_2 (&input_mesh, &gui, "input_mesh");

    Data data;
    int t_deb;
    int t_deb_inp;
    loadTriangleIDsFromFile(data);

    if(!input_mesh_flag) {
        gui.push(&bool_mesh);
        gui.push(&contr);
        bool_mesh.updateGL();
    }else{
        gui.push(&input_mesh);
        gui.push(&contr_2);
        input_mesh.updateGL();
    }

    std::vector<std::tuple<uint, int, int>> tri_labels = readBelongingFromFile("labels.txt");
    bool trisA = false, trisB = false, trisAB =false;
    bool insideA = false, insideB = false;
    gui.callback_app_controls = [&](){

        if(ImGui::Button("Invisible")){
            for(uint i = 0; i < bool_mesh.num_polys(); ++i){
                bool_mesh.poly_data(i).color =  cinolib::Color(255/255.f,255/255.f,255/255.f,0.0f);
            }
            bool_mesh.updateGL();
        }
    if(ImGui::Button("Intersection parts")){
        for(uint i = 0; i < data.t_ids_intersection.size(); ++i){
                bool_mesh.poly_data(data.t_ids_intersection.at(i)).color = Color::PASTEL_YELLOW();
            }
            bool_mesh.updateGL();
        }
        if(ImGui::Button("Union parts")){
             for(uint i = 0; i < data.t_ids_union.size(); ++i){
                bool_mesh.poly_data(data.t_ids_union.at(i)).color = Color::PASTEL_GREEN();
            }
            bool_mesh.updateGL();
        }
        if(ImGui::Button("Subtraction parts")){
            for(uint i = 0; i < data.t_ids_subtraction.size(); ++i){
                bool_mesh.poly_data(data.t_ids_subtraction.at(i)).color = Color::PASTEL_RED();
            }
        }

        if(ImGui::Button("All parts")){
            for(uint i = 0; i < data.t_ids_intersection.size(); ++i){
                bool_mesh.poly_data(data.t_ids_intersection.at(i)).color = Color::PASTEL_YELLOW();
            }
           for(uint i = 0; i < data.t_ids_union.size(); ++i){
                bool_mesh.poly_data(data.t_ids_union.at(i)).color = Color::PASTEL_GREEN();
            }
            for(uint i = 0; i < data.t_ids_subtraction.size(); ++i){
                bool_mesh.poly_data(data.t_ids_subtraction.at(i)).color = Color::PASTEL_RED();
            }

            bool_mesh.updateGL();
        }

        if(ImGui::Button("Reset")){
            for(uint t_id = 0; t_id < bool_mesh.num_polys(); ++t_id){
                bool_mesh.poly_data(t_id).color = Color::WHITE();
                }
                bool_mesh.updateGL();
        }

        if(ImGui::InputInt("Triangle id", &t_deb, 0, input_mesh.num_verts())){
        }
        if (ImGui::Button("Show Triangle")) {


            if(!input_mesh_flag){
                bool_mesh.poly_data(t_deb).color = cinolib::Color(cinolib::Color::RED());
                cinolib::Marker marker;
                marker.color = cinolib::Color(cinolib::Color::RED());
                marker.text = std::to_string(t_deb);
                marker.disk_radius = 0.3f;
                marker.pos_3d = bool_mesh.poly_centroid(t_deb);
                gui.push(marker);
                bool_mesh.updateGL();
            }else{
                input_mesh.poly_data(t_deb).color = cinolib::Color(cinolib::Color::RED());
                cinolib::Marker marker;
                marker.color = cinolib::Color(cinolib::Color::RED());
                marker.text = std::to_string(t_deb);
                marker.disk_radius = 0.3f;
                marker.pos_3d = input_mesh.poly_centroid(t_deb);
                gui.push(marker);
                input_mesh.updateGL();
                std::cout << "triangle t_id: " << t_deb <<  " v0 " << input_mesh.vert(bool_mesh.poly_vert_id(t_deb, 0)).x() << " " << input_mesh.vert(bool_mesh.poly_vert_id(t_deb, 0)).y() << " " << input_mesh.vert(bool_mesh.poly_vert_id(t_deb, 0)).z() << std::endl;
                std::cout << "triangle t_id: " << t_deb <<  " v1 " << input_mesh.vert(bool_mesh.poly_vert_id(t_deb, 1)).x() << " " << input_mesh.vert(bool_mesh.poly_vert_id(t_deb, 1)).y() << " " << input_mesh.vert(bool_mesh.poly_vert_id(t_deb, 1)).z() << std::endl;
                std::cout << "triangle t_id: " << t_deb <<  " v2 " << input_mesh.vert(bool_mesh.poly_vert_id(t_deb, 2)).x() << " " << input_mesh.vert(bool_mesh.poly_vert_id(t_deb, 2)).y() << " " << input_mesh.vert(bool_mesh.poly_vert_id(t_deb, 2)).z() << std::endl;


            }
        }
        if(ImGui::InputInt("Triangle id input mesh", &t_deb_inp, 0, input_mesh.num_verts())){
        }
        if (ImGui::Button("Show Triangle Input")) {
                if(!is_pushed){
                    gui.push(&input_mesh);
                    gui.push(&contr_2);
                    for(int i = 0; i < input_mesh.num_polys(); ++i){
                        input_mesh.poly_data(i).color =  cinolib::Color(255/255.f,255/255.f,255/255.f,0.0f);
                        input_mesh.show_wireframe(false);
                    }
                    is_pushed =true;
                }


                input_mesh.poly_data(t_deb_inp).color = cinolib::Color(cinolib::Color::PASTEL_CYAN());
                cinolib::Marker marker;
                marker.color = cinolib::Color(cinolib::Color::BLACK());
                marker.text = std::to_string(t_deb_inp);
                marker.disk_radius = 0.3f;
                marker.pos_3d = input_mesh.poly_centroid(t_deb_inp);
                gui.push(marker);
                input_mesh.updateGL();
            }
        if (ImGui::Button("Show Ray")) {
            vec3d ray_start(-7.49611e+09,
                     -3.15328e+08,
                     1.24075e+09
            );
            /*vec3d ray_end(1.11727e+10,
                     -3.15328e+08,
                     1.24075e+09);
            */
            vec3d ray_end(1.61118e+10,
                     bool_mesh.poly_centroid(3794).y(),
                     bool_mesh.poly_centroid(3794).z());

            ray.push_seg(
                    bool_mesh.poly_centroid(3794),
                    //ray_start,
                    ray_end
                    );
            ray.thickness = 5.0f;
            gui.push(&ray);

            std::cout << bool_mesh.poly_area(3794)<< std::endl;
        }

        if (ImGui::Button("Show Multiple Triangles")) {
            std::vector<uint> tris_to_show = { 3794, 3792, 4007, 4015, 4003, 3793, 4005, 3798, 3942, 3797, 3802, 3998, 4006, 4009, 4004, 3796, 4281, 3941, 3945, 3795, 3789, 4013  };

            for(int i = 0 ; i < tris_to_show.size(); ++i){

                if(input_mesh_flag){
                    input_mesh.poly_data(tris_to_show.at(i)).color = cinolib::Color(cinolib::Color::YELLOW());
                    cinolib::Marker marker;
                    marker.color = cinolib::Color(cinolib::Color::RED());
                    marker.text = std::to_string(tris_to_show.at(i));
                    marker.disk_radius = 0.3f;
                    marker.pos_3d = input_mesh.poly_centroid(tris_to_show.at(i));
                    gui.push(marker);
                    input_mesh.updateGL();
                }else{
                    bool_mesh.poly_data(tris_to_show.at(i)).color = cinolib::Color(cinolib::Color::YELLOW());
                    cinolib::Marker marker;
                    marker.color = cinolib::Color(cinolib::Color::RED());
                    marker.text = std::to_string(tris_to_show.at(i));
                    marker.disk_radius = 0.3f;
                    marker.pos_3d = bool_mesh.poly_centroid(tris_to_show.at(i));
                    //gui.push(marker);
                    bool_mesh.updateGL();
                }
            }
        }
        if (ImGui::Button("Show holes")) {

            if(!input_mesh_flag) {
                for(uint eid=0; eid<bool_mesh.num_edges(); ++eid)
                {
                    if(bool_mesh.edge_is_boundary(eid)) {
                        vector<uint> polys_id = bool_mesh.adj_e2p(eid);
                        for (uint k = 0 ; k < polys_id.size(); ++k) {
                            std::cout << polys_id.at(k) << " " << std::endl;
                            Marker marker_app;
                            marker_app.text = std::to_string(polys_id.at(k));
                            marker_app.disk_radius = 1.0f;
                            marker_app.pos_3d = bool_mesh.poly_centroid(polys_id.at(k));
                            gui.push(marker_app);
                        }

                    }
                }
            }
        }

        if (ImGui::Button("Show not manifold")) {

            if(!input_mesh_flag) {
                for(uint e_id=0; e_id<bool_mesh.num_edges(); ++e_id)
                {
                    if( !bool_mesh.edge_is_manifold(e_id))  {
                        vector<uint> polys_id = bool_mesh.adj_e2p(e_id);
                        for (uint k = 0 ; k < polys_id.size(); ++k) {
                            std::cout << polys_id.at(k) << " " << std::endl;
                            Marker marker_app;
                            marker_app.text = std::to_string(polys_id.at(k));
                            marker_app.disk_radius = 1.0f;
                            marker_app.pos_3d = bool_mesh.poly_centroid(polys_id.at(k));
                            gui.push(marker_app);
                        }

                    }
                }
            }
        }
        if (ImGui::Button("Hide all ")) {
            for(uint p = 0 ; p < bool_mesh.num_polys() ; ++p){
                bool_mesh.poly_data(p).flags[HIDDEN] = true;
            }
            bool_mesh.updateGL();
        }
        if (ImGui::Button("Show not manifold triangles")) {

            if(!input_mesh_flag) {

                for(uint e_id=0; e_id<bool_mesh.num_edges(); ++e_id)
                {

                    if( !bool_mesh.edge_is_manifold(e_id)) {
                        vector<uint> polys_id = bool_mesh.adj_e2p(e_id);
                        for (uint k = 0 ; k < polys_id.size(); ++k) {
                            std::cout << "Poly_id : " << polys_id.at(k) << std::endl;
                            uint p_id = polys_id.at(k);
                            auto it = std::find_if(tri_labels.begin(), tri_labels.end(),
                                                   [p_id](const std::tuple<uint, int, int>& t) {
                                                       return std::get<0>(t) == p_id;
                                                   });

                            if (it != tri_labels.end()) {
                                size_t index = std::distance(tri_labels.begin(), it);
                                int mesh_bel = std::get<1>(tri_labels.at(index));
                                int inside_mesh = std::get<2>(tri_labels.at(index));
                                if(mesh_bel == 0 &&  inside_mesh == 0){
                                    bool_mesh.poly_data(p_id).color = cinolib::Color::PASTEL_RED();
                                    bool_mesh.poly_data(p_id).flags[HIDDEN] = false;
                                }else if(mesh_bel == 0 && inside_mesh == 1){
                                    bool_mesh.poly_data(p_id).color = cinolib::Color::RED();
                                    bool_mesh.poly_data(p_id).flags[HIDDEN] = false;
                                }else if(mesh_bel == 1 && inside_mesh == 0){
                                    bool_mesh.poly_data(p_id).color = cinolib::Color::PASTEL_CYAN();
                                    bool_mesh.poly_data(k).flags[HIDDEN] = false;
                                }else if(mesh_bel == 1 && inside_mesh == 1){
                                    bool_mesh.poly_data(p_id).color = cinolib::Color::BLUE();
                                    bool_mesh.poly_data(p_id).flags[HIDDEN] = false;
                                }else if(mesh_bel == 2){
                                    bool_mesh.poly_data(p_id).color = cinolib::Color::PASTEL_GRAY();
                                    bool_mesh.poly_data(p_id).flags[HIDDEN] = false;
                                }else{
                                    std::cout << "Errore" << std::endl;
                                }

                            } else {
                                std::cout << "Valore non trovato." << std::endl;
                            }


                            /*Marker marker_app;
                            marker_app.text = std::to_string(polys_id.at(k));
                            marker_app.disk_radius = 1.0f;
                            marker_app.pos_3d = bool_mesh.poly_centroid(polys_id.at(k));
                            gui.push(marker_app);*/


                        }
                        bool_mesh.updateGL();

                    }
                }

            }
        }

        if(ImGui::Checkbox("A -> RED", &trisA)){
            for(std::tuple<uint, int, int> tri_label : tri_labels){
                if(trisA) {
                    if (std::get<1>(tri_label) == 0){
                        bool_mesh.poly_data(std::get<0>(tri_label)).color = cinolib::Color::PASTEL_RED();
                        //bool_mesh.poly_data(std::get<0>(tri_label)).flags[HIDDEN] = 0;
                    }

                    /*} else {
                        if (!trisB && !trisAB) bool_mesh.poly_data(tri_label.first).color = cinolib::Color(1, 1, 1, 0);
                    }*/
                }else{
                    //bool_mesh.poly_data(tri_label.first).color = cinolib::Color(1, 1, 1, 0);
                    if(trisB && std::get<1>(tri_label) == 1) continue;
                    if(trisAB && std::get<1>(tri_label) == 2) continue;
                    //bool_mesh.poly_data(std::get<0>(tri_label)).flags[HIDDEN] = 1;
                }
            }

            bool_mesh.updateGL();
        }
        if(trisA){
            if(ImGui::Checkbox("A Inside B", &insideB)){
                for(std::tuple<uint, int, int> tri_label : tri_labels){
                    if(std::get<1>(tri_label) == 0 && std::get<2>(tri_label) == 1){
                        if(!insideB){
                            bool_mesh.poly_data(std::get<0>(tri_label)).color = cinolib::Color::PASTEL_RED();
                        }else{
                        bool_mesh.poly_data(std::get<0>(tri_label)).color = cinolib::Color::RED();
                        }
                    }
                }

            }
            bool_mesh.updateGL();
        }

        if(ImGui::Checkbox("B -> CYAN", &trisB)){
            for(std::tuple<uint, int, int> tri_label : tri_labels){
                if(trisB) {
                    if (std::get<1>(tri_label) == 1){
                        bool_mesh.poly_data(std::get<0>(tri_label)).color = cinolib::Color::PASTEL_CYAN();
                        //bool_mesh.poly_data(std::get<0>(tri_label)).flags[HIDDEN] = 0;

                    }
                   /* }else{
                            if(!trisA && !trisAB) bool_mesh.poly_data(tri_label.first).color = cinolib::Color(1,1,1,0);
                    }*/
                }else{
                    //bool_mesh.poly_data(tri_label.first).color = cinolib::Color(1,1,1,0);

                    if(trisA && std::get<1>(tri_label) == 0) continue;
                    if(trisAB && std::get<1>(tri_label) == 2) continue;
                    //bool_mesh.poly_data(std::get<0>(tri_label)).flags[HIDDEN] = 1;

                }
        }

            bool_mesh.updateGL();
        }

        if(trisB){
            if(ImGui::Checkbox("B Inside A", &insideA)){
                for(std::tuple<uint, int, int> tri_label : tri_labels){
                    if(std::get<1>(tri_label) == 1 && std::get<2>(tri_label) == 1){
                        if(!insideA){
                            bool_mesh.poly_data(std::get<0>(tri_label)).color = cinolib::Color::PASTEL_CYAN();
                        }else{
                            bool_mesh.poly_data(std::get<0>(tri_label)).color = cinolib::Color::BLUE();
                           // bool_mesh.poly_data(std::get<0>(tri_label)).flags[HIDDEN] = 0;

                        }
                    }
                }

            }
            bool_mesh.updateGL();
        }

        if(ImGui::Checkbox("AB -> GRAY ", &trisAB)){
            for(std::tuple<uint, int, int> tri_label : tri_labels){
                if(trisAB) {
                    if (std::get<1>(tri_label) == 2){
                        bool_mesh.poly_data(std::get<0>(tri_label)).color = cinolib::Color::PASTEL_GRAY();
                       // bool_mesh.poly_data(std::get<0>(tri_label)).flags[HIDDEN] = 0;
                       }
                   /* }else{
                        if(!trisA && !trisB) bool_mesh.poly_data(tri_label.first).color = cinolib::Color(1,1,1,0);
                        }*/
                    }else{
                        //bool_mesh.poly_data(tri_label.first).color = cinolib::Color(1,1,1,0);
                        if(trisA && std::get<1>(tri_label) == 0) continue;
                        if(trisB && std::get<1>(tri_label) == 1) continue;
                        //bool_mesh.poly_data(std::get<0>(tri_label)).flags[HIDDEN] = 1;
                }
                }
            bool_mesh.updateGL();
        }

    };


    return gui.launch();
}
