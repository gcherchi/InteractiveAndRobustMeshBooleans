/*****************************************************************************************
 *              MIT License                                                              *
 *                                                                                       *
 * Copyright (c) 2022 G. Cherchi, F. Pellacini, M. Attene and M. Livesu                  *
 *                                                                                       *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this  *
 * software and associated documentation files (the "Software"), to deal in the Software *
 * without restriction, including without limitation the rights to use, copy, modify,    *
 * merge, publish, distribute, sublicense, and/or sell copies of the Software, and to    *
 * permit persons to whom the Software is furnished to do so, subject to the following   *
 * conditions:                                                                           *
 *                                                                                       *
 * The above copyright notice and this permission notice shall be included in all copies *
 * or substantial portions of the Software.                                              *
 *                                                                                       *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,   *
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A         *
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT    *
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION     *
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE        *
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                *
 *                                                                                       *
 * Authors:                                                                              *
 *      Gianmarco Cherchi (g.cherchi@unica.it)                                           *
 *      https://www.gianmarcocherchi.com                                                 *
 *                                                                                       *
 *      Fabio Pellacini (fabio.pellacini@uniroma1.it)                                    *
 *      https://pellacini.di.uniroma1.it                                                 *
 *                                                                                       *
 *      Marco Attene (marco.attene@ge.imati.cnr.it)                                      *
 *      https://www.cnr.it/en/people/marco.attene/                                       *
 *                                                                                       *
 *      Marco Livesu (marco.livesu@ge.imati.cnr.it)                                      *
 *      http://pers.ge.imati.cnr.it/livesu/                                              *
 *                                                                                       *
 * ***************************************************************************************/

#include <cinolib/meshes/drawable_trimesh.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/drawable_segment_soup.h>
#include <cinolib/gl/surface_mesh_controls.h>

#include "booleans.h"

std::vector<std::string> files;

std::vector<Ray> global_rays;
std::set<uint> global_inters_tris;

bool draw_patches = false;
bool draw_rays = true;
bool draw_inters_tris = true;
bool draw_inside_out = true;

uint curr_label = 0;

int main(int argc, char **argv)
{
    std::vector<double> in_coords, bool_coords;
    std::vector<uint> in_tris, bool_tris;
    std::vector<uint> in_labels;
    std::vector<std::bitset<NBIT>> bool_labels;

    BoolOp op = INTERSECTION;

    files.emplace_back("/Users/gianmarco/Code/InteractiveAndRobustMeshBooleans/data/sphere1.obj");
    files.emplace_back("/Users/gianmarco/Code/InteractiveAndRobustMeshBooleans/data/sphere2.obj");
    std::string file_out = "res.obj";

    loadMultipleFiles(files, in_coords, in_tris, in_labels);

    initFPU();

    point_arena arena;
    std::vector<genericPoint*> arr_verts; // <- it contains the original expl verts + the new_impl verts
    std::vector<uint> arr_in_tris, arr_out_tris;
    std::vector<std::bitset<NBIT>> arr_in_labels;
    std::vector<DuplTriInfo> dupl_triangles;
    Labels labels;
    std::vector<phmap::flat_hash_set<uint>> patches;
    cinolib::Octree octree; // built with arr_in_tris and arr_in_labels

    customArrangementPipeline(in_coords, in_tris, in_labels, arr_in_tris, arr_in_labels, arena, arr_verts,
                              arr_out_tris, labels, octree, dupl_triangles);


    //customBooleanPipeline(arr_verts, arr_in_tris, arr_out_tris, arr_in_labels, dupl_triangles, labels, patches, octree, op, bool_coords, bool_tris, bool_labels);
    //cinolib::write_OBJ(file_out.c_str(), bool_coords, bool_tris, {});
    //return 1;


    std::vector<cinolib::vec3d> drawable_verts(arr_verts.size());

    double x, y, z;
    for(uint vid = 0; vid < arr_verts.size()-4; vid++)
    {
        arr_verts[vid]->getApproxXYZCoordinates(x, y, z);
        drawable_verts[vid] = cinolib::vec3d(x, y, z) ;
    }

    cinolib::GLcanvas gui;
    cinolib::DrawableTrimesh<> arr_mesh(drawable_verts, arr_out_tris);
    cinolib::SurfaceMeshControls<cinolib::DrawableTrimesh<>> menu(&arr_mesh, &gui);
    gui.push(&arr_mesh);
    gui.push(&menu);
    gui.show_side_bar = true;
    gui.side_bar_width = 0.23f;
    gui.side_bar_alpha = 0.6f;

    FastTrimesh tm(arr_verts, arr_out_tris);

    double multiplier = tm.vert(tm.numVerts() - 1)->toExplicit3D().X();
    std::cout << "multiplier: " << multiplier << std::endl;

    computeAllPatches(tm, labels, patches);

    for(uint pid = 0; pid < patches.size(); pid++)
    {
        for(auto ti : patches[pid])
            arr_mesh.poly_data(ti).label = (int)pid;
    }

    if(draw_patches)
    {
        arr_mesh.poly_color_wrt_label();
        arr_mesh.updateGL();
    }


    // the informations about duplicated triangles (removed in arrangements) are restored in the original structures
    addDuplicateTrisInfoInStructures(dupl_triangles, arr_in_tris, arr_in_labels, octree);

    // parse patches with octree and rays
    cinolib::vec3d max_coords(octree.root->bbox.max.x() +0.5, octree.root->bbox.max.y() +0.5, octree.root->bbox.max.z() +0.5);

    //computeInsideOut(tm, patches, octree, arr_verts, arr_in_tris, arr_in_labels, max_coords, labels);


    for(uint pid = 0; pid < patches.size(); pid++)
    {
        auto &patch = patches[pid];
        const std::bitset<NBIT> &patch_surface_label = labels.surface[*patch.begin()]; // label of the first triangle of the patch

        Ray ray;
        findRayEndpoints(tm, patch, max_coords, ray);
        global_rays.push_back(ray);

        phmap::flat_hash_set<uint> tmp_inters;
        cinolib::AABB rayAABB(cinolib::vec3d(ray.v0.X(), ray.v0.Y(), ray.v0.Z()),
                              cinolib::vec3d(ray.v1.X(), ray.v1.Y(), ray.v1.Z()));

        intersects_box(octree, rayAABB, tmp_inters);

        std::vector<uint> sorted_inters;
        pruneIntersectionsAndSortAlongRay(ray, arr_verts, arr_in_tris, arr_in_labels, tmp_inters, patch_surface_label, sorted_inters);

        for(auto ti : sorted_inters) global_inters_tris.insert(ti); // FOR DEBUG

        std::bitset<NBIT> patch_inner_label;
        analyzeSortedIntersections(ray, arr_verts, arr_in_tris, arr_in_labels, sorted_inters, patch_inner_label);

        propagateInnerLabelsOnPatch(patch, patch_inner_label, labels);

    }

    std::vector<cinolib::DrawableSegmentSoup> drawable_rays(global_rays.size());
    for(uint rid = 0; rid < global_rays.size(); rid++)
    {
        double x0, y0, z0, x1, y1, z1;
        global_rays[rid].v0.getApproxXYZCoordinates(x0, y0, z0);
        global_rays[rid].v1.getApproxXYZCoordinates(x1, y1, z1);
        drawable_rays[rid].push_seg(cinolib::vec3d(x0,y0,z0), cinolib::vec3d(x1,y1,z1));
        drawable_rays[rid].set_color(cinolib::Color::BLUE());
        drawable_rays[rid].set_cheap_rendering(true);
        if(draw_rays) gui.push(&drawable_rays[rid], false);
    }


    // DRAW INTERSECTED TRIS
    std::vector<uint> empty;
    cinolib::DrawableTrimesh<> it(in_coords, empty);
    for(uint v = 0; v < it.num_verts(); v++) it.vert(v) = multiplier * it.vert(v);
    for(auto tid : global_inters_tris) it.poly_add(in_tris[3*tid], in_tris[3*tid +1], in_tris[3*tid +2]);
    it.poly_set_color(cinolib::Color::BLUE()); it.updateGL();
    if(draw_inters_tris) gui.push(&it);


    // BOOLEAN OPERATION
    uint num_tris_in_final_solution;
    if(op == INTERSECTION) num_tris_in_final_solution = boolIntersection(tm, labels);
    else if(op == UNION) num_tris_in_final_solution = boolUnion(tm, labels);
    else if(op == SUBTRACTION) num_tris_in_final_solution = boolSubtraction(tm, labels);
    else if(op == XOR) num_tris_in_final_solution = boolXOR(tm, labels);
    else
    {
        std::cerr << "boolean operation not implemented yet" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    computeFinalExplicitResult(tm, labels, num_tris_in_final_solution, bool_coords, bool_tris, bool_labels, true);

    cinolib::GLcanvas gui_out;
    cinolib::DrawableTrimesh<> bool_mesh(bool_coords, bool_tris);
    cinolib::SurfaceMeshControls<cinolib::DrawableTrimesh<>> menu_out(&bool_mesh, &gui_out);
    gui.push(&menu_out);
    gui_out.push(&bool_mesh);

    gui.callback_app_controls = [&]()
    {
        if(ImGui::Checkbox("patches", &draw_patches))
        {
            draw_inside_out = !draw_patches;
            if(draw_patches) arr_mesh.poly_color_wrt_label();
            else arr_mesh.poly_set_color(cinolib::Color::WHITE());
            arr_mesh.updateGL();
        }

        if(ImGui::Checkbox("rays", &draw_rays))
        {
            if(draw_rays) for(auto &r : drawable_rays) gui.push(&r);
            else for(auto &r : drawable_rays) gui.pop(&r);
        }

        if(ImGui::Checkbox("int. tris", &draw_inters_tris))
        {
            if(draw_inters_tris) gui.push(&it);
            else gui.pop(&it);
        }

        if(ImGui::Checkbox("inside/out", &draw_inside_out))
        {
            draw_patches = !draw_patches;
            arr_mesh.poly_set_color(cinolib::Color::WHITE());
            if(draw_inside_out)
            {
                std::cout << "mesh " << curr_label << std::endl;
                for(uint tid = 0; tid < arr_mesh.num_polys(); tid++)
                {
                    if(labels.surface[tid][curr_label] == 1) arr_mesh.poly_data(tid).color = cinolib::Color::PASTEL_YELLOW();
                    if(labels.inside[tid][curr_label] == 1) arr_mesh.poly_data(tid).color = cinolib::Color::PASTEL_ORANGE();
                }
            }

            arr_mesh.updateGL();
        }
    };

    gui.callback_key_pressed = [&](int key, int modifiers) -> bool
    {
        if(key == GLFW_KEY_N) //press "P" to show/hide points
        {
            curr_label = (curr_label + 1) % files.size();
            arr_mesh.poly_set_color(cinolib::Color::WHITE());
            if(draw_inside_out)
            {
                std::cout << "mesh " << curr_label << std::endl;
                for(uint tid = 0; tid < arr_mesh.num_polys(); tid++)
                {
                    if(labels.surface[tid][curr_label] == 1) arr_mesh.poly_data(tid).color = cinolib::Color::PASTEL_YELLOW();
                    if(labels.inside[tid][curr_label] == 1) arr_mesh.poly_data(tid).color = cinolib::Color::PASTEL_ORANGE();
                }
                arr_mesh.updateGL();
            }
        }

        if(key == GLFW_KEY_S)
        {
            cinolib::DrawableTrimesh<> tmp = arr_mesh;
            for(uint vid = 0; vid < arr_mesh.num_verts(); vid++)
                tmp.vert(vid) /= 1;multiplier;

            tmp.poly_set_color(cinolib::Color::WHITE());
            tmp.save("/Users/gianmarco/Desktop/arr_mesh.obj");
        }

        return false;
    };


    return gui.launch({&gui_out});
}