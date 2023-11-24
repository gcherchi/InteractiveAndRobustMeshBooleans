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
#include <cinolib/gl/surface_mesh_controls.h>

#include "booleans.h"

std::vector<std::string> files;

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

    std::vector<cinolib::vec3d> drawable_verts(arr_verts.size());

    double x, y, z;
    for(uint vid = 0; vid < arr_verts.size()-4; vid++)
    {
        arr_verts[vid]->getApproxXYZCoordinates(x, y, z);
        drawable_verts[vid] = cinolib::vec3d(x, y, z);
    }

    cinolib::GLcanvas gui;
    cinolib::DrawableTrimesh<> arr_mesh(drawable_verts, arr_out_tris);
    cinolib::SurfaceMeshControls<cinolib::DrawableTrimesh<>> menu(&arr_mesh, &gui);
    gui.push(&arr_mesh);
    gui.push(&menu);

    FastTrimesh tm(arr_verts, arr_out_tris);

    computeAllPatches(tm, labels, patches);

    /*

    // the informations about duplicated triangles (removed in arrangements) are restored in the original structures
    addDuplicateTrisInfoInStructures(dupl_triangles, arr_in_tris, arr_in_labels, octree);

    // parse patches with octree and rays
    cinolib::vec3d max_coords(octree.root->bbox.max.x() +0.5, octree.root->bbox.max.y() +0.5, octree.root->bbox.max.z() +0.5);
    computeInsideOut(tm, patches, octree, arr_verts, arr_in_tris, arr_in_labels, max_coords, labels);

    // booleand operations
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

    cinolib::write_OBJ(file_out.c_str(), bool_coords, bool_tris, {});

     */

    return gui.launch();
}