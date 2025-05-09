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

bool input_mesh_flag = true;

int main(int argc, char **argv){

    GLcanvas gui;
    DrawableTrimesh<> bool_mesh ("output.obj");
    SurfaceMeshControls<DrawableTrimesh<>> contr (&bool_mesh, &gui, "bool_mesh");
    DrawableSegmentSoup ray;

    DrawableTrimesh<> input_mesh ("mesh_input.obj");
    SurfaceMeshControls<DrawableTrimesh<>> contr_2 (&input_mesh, &gui, "input_mesh");

    Data data;
    int t_deb;
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

    gui.callback_app_controls = [&](){
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
                bool_mesh.poly_data(t_deb).color = cinolib::Color(cinolib::Color::YELLOW());
                cinolib::Marker marker;
                marker.color = cinolib::Color(cinolib::Color::RED());
                marker.text = std::to_string(t_deb);
                marker.disk_radius = 0.3f;
                marker.pos_3d = bool_mesh.poly_centroid(t_deb);
                gui.push(marker);
                bool_mesh.updateGL();
            }else{
                input_mesh.poly_data(t_deb).color = cinolib::Color(cinolib::Color::YELLOW());
                cinolib::Marker marker;
                marker.color = cinolib::Color(cinolib::Color::RED());
                marker.text = std::to_string(t_deb);
                marker.disk_radius = 0.3f;
                marker.pos_3d = input_mesh.poly_centroid(t_deb);
                gui.push(marker);
                input_mesh.updateGL();
            }
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
            std::vector<uint> tris_to_show = {0, 293, 1333, 1560, 1709, 1973, 3042, 3391};
            for(int i = 0 ; i < tris_to_show.size(); ++i){
                input_mesh.poly_data(tris_to_show.at(i)).color = cinolib::Color(cinolib::Color::YELLOW());
                cinolib::Marker marker;
                marker.color = cinolib::Color(cinolib::Color::RED());
                marker.text = std::to_string(tris_to_show.at(i));
                marker.disk_radius = 0.3f;
                marker.pos_3d = input_mesh.poly_centroid(tris_to_show.at(i));
                gui.push(marker);
                input_mesh.updateGL();
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


    };



    return gui.launch();
}
