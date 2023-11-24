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

int main(int argc, char **argv)
{
    std::vector<double> in_coords, bool_coords;
    std::vector<uint> in_tris, bool_tris;
    std::vector<uint> in_labels;
    std::vector<std::bitset<NBIT>> bool_labels;

    BoolOp op = INTERSECTION;

    files.emplace_back("../data/sphere1.obj");
    files.emplace_back("../data/sphere2.obj");
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

    computeAllPatches(tm, labels, patches);

    for(uint pid = 0; pid < patches.size(); pid++)
    {
        for(auto ti : patches[pid])
            arr_mesh.poly_data(ti).label = (int)pid;
    }

    arr_mesh.poly_color_wrt_label();
    arr_mesh.updateGL();


    // the informations about duplicated triangles (removed in arrangements) are restored in the original structures
    addDuplicateTrisInfoInStructures(dupl_triangles, arr_in_tris, arr_in_labels, octree);

    // parse patches with octree and rays
    cinolib::vec3d max_coords(octree.root->bbox.max.x() +0.5, octree.root->bbox.max.y() +0.5, octree.root->bbox.max.z() +0.5);

    //computeInsideOut(tm, patches, octree, arr_verts, arr_in_tris, arr_in_labels, max_coords, labels);


    for(auto &patch : patches)
    {
        const std::bitset<NBIT> &patch_surface_label = labels.surface[*patch.begin()]; // label of the first triangle of the patch

        Ray ray;
        findRayEndpoints(tm, patch, max_coords, ray);
        global_rays.push_back(ray);

        phmap::flat_hash_set<uint> tmp_inters;
        cinolib::AABB rayAABB(cinolib::vec3d(ray.v0.X(), ray.v0.Y(), ray.v0.Z()),
                              cinolib::vec3d(ray.v1.X(), ray.v1.Y(), ray.v1.Z()));

        intersects_box(octree, rayAABB, tmp_inters);

        std::vector<uint> sorted_inters;
        //pruneIntersectionsAndSortAlongRay(ray, arr_verts, in_tris, arr_in_labels, tmp_inters, patch_surface_label, sorted_inters);



        for(auto ti : tmp_inters)
            global_inters_tris.insert(ti);

    }

    std::cout << "Num patches: " << patches.size() << " - Num rays: " << global_rays.size() << std::endl;

    std::vector<cinolib::DrawableSegmentSoup> drawable_rays(global_rays.size());
    for(uint rid = 0; rid < global_rays.size(); rid++)
    {
        double x0, y0, z0, x1, y1, z1;
        global_rays[rid].v0.getApproxXYZCoordinates(x0, y0, z0);
        global_rays[rid].v1.getApproxXYZCoordinates(x1, y1, z1);
        drawable_rays[rid].push_seg(cinolib::vec3d(x0,y0,z0), cinolib::vec3d(x1,y1,z1));
        drawable_rays[rid].set_color(cinolib::Color::BLUE());
        drawable_rays[rid].set_cheap_rendering(true);
        gui.push(&drawable_rays[rid], false);
    }


    // DRAW INTERSECTED TRIS
    std::vector<uint> empty;
    cinolib::DrawableTrimesh<> it(in_coords, empty);
    double multiplier = computeMultiplier(in_coords);
    for(uint v = 0; v < it.num_verts(); v++) it.vert(v) = multiplier * it.vert(v);
    for(auto tid : global_inters_tris) it.poly_add(in_tris[3*tid], in_tris[3*tid +1], in_tris[3*tid +2]);
    it.poly_set_color(cinolib::Color::BLUE()); it.updateGL();
    gui.push(&it);


    /*
    // booleans operations
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

     */


    return gui.launch();
}