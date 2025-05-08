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
#include "booleans.h"
#include "debug.h"
#include "io_functions.h"
#include <tbb/tbb.h>


#include <cinolib/profiler.h>
bool print_debug = false;
bool profiling = false;


inline void customBooleanPipeline(std::vector<genericPoint*>& arr_verts, std::vector<uint>& arr_in_tris,
                                  std::vector<uint>& arr_out_tris, std::vector<std::bitset<NBIT>>& arr_in_labels,
                                  std::vector<DuplTriInfo>& dupl_triangles, Labels& labels,
                                  std::vector<phmap::flat_hash_set<uint>>& patches, cinolib::Octree& octree,
                                  const BoolOp &op, std::vector<double> &bool_coords, std::vector<uint> &bool_tris,
                                  std::vector< std::bitset<NBIT>> &bool_labels, Data &data)
{
    FastTrimesh tm(arr_verts, arr_out_tris, false);

    computeAllPatches(tm, labels, patches, false);

    // the informations about duplicated triangles (removed in arrangements) are restored in the original structures
    addDuplicateTrisInfoInStructures(dupl_triangles, arr_in_tris, arr_in_labels, octree);

    // parse patches with octree and rays
    cinolib::vec3d max_coords(octree.root->bbox.max.x() +0.5, octree.root->bbox.max.y() +0.5, octree.root->bbox.max.z() +0.5);

    computeInsideOutCustom(tm, patches, octree, arr_verts, arr_in_tris, arr_in_labels, max_coords, labels, data);

    //computeInsideOut(tm, patches, octree, arr_verts, arr_in_tris, arr_in_labels, max_coords, labels);

    // booleand operations
    uint num_tris_in_final_solution;
    if(op == INTERSECTION)
        num_tris_in_final_solution = boolIntersection(tm, labels, data);
    else if(op == UNION)
        num_tris_in_final_solution = boolUnion(tm, labels, data);
    else if(op == SUBTRACTION)
        num_tris_in_final_solution = boolSubtraction(tm, labels, data);
    else if(op == XOR)
        num_tris_in_final_solution = boolXOR(tm, labels);

        else if(op == DEBUG)
        num_tris_in_final_solution = classifyTriangles(tm, labels, data);

    else
    {
        std::cerr << "boolean operation not implemented yet" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    computeFinalExplicitResult(tm, labels, num_tris_in_final_solution, bool_coords, bool_tris, bool_labels, true, data);
}

extern int arr_time;
extern int bool_time;
extern std::vector<std::string> files;

inline void booleanPipeline(const std::vector<double> &in_coords, const std::vector<uint> &in_tris,
                            const std::vector<uint> &in_labels, const BoolOp &op, std::vector<double> &bool_coords,
                            std::vector<uint> &bool_tris, std::vector< std::bitset<NBIT> > &bool_labels, Data &data)
{
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
                              arr_out_tris, labels, octree, dupl_triangles, false);

    customBooleanPipeline(arr_verts, arr_in_tris, arr_out_tris, arr_in_labels, dupl_triangles, labels,
                          patches, octree, op, bool_coords, bool_tris, bool_labels, data);
}



//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

/* a custom arrangement pipeline in witch we can expose the octree used to find the starting intersection list */
inline void customArrangementPipeline(const std::vector<double> &in_coords, const std::vector<uint> &in_tris, const std::vector<uint> &in_labels,
                                      std::vector<uint> &arr_in_tris, std::vector< std::bitset<NBIT>> &arr_in_labels,
                                      point_arena& arena, std::vector<genericPoint *> &vertices, std::vector<uint> &arr_out_tris, Labels &labels,
                                      cinolib::Octree &octree, std::vector<DuplTriInfo> &dupl_triangles, bool parallel)
{
    arr_in_labels.resize(in_labels.size());
    std::bitset<NBIT> mask;

    for(uint i = 0; i < in_labels.size(); i++)
    {
        arr_in_labels[i][in_labels[i]] = true;
        mask[in_labels[i]] = true;
    }

    labels.num = mask.count();

    initFPU();
    double multiplier = computeMultiplier(in_coords);

    mergeDuplicatedVertices(in_coords, in_tris, arena, vertices, arr_in_tris, parallel);

    customRemoveDegenerateAndDuplicatedTriangles(vertices, arr_in_tris, arr_in_labels, dupl_triangles, parallel);

    TriangleSoup ts(arena, vertices, arr_in_tris, arr_in_labels, multiplier, parallel);

    AuxiliaryStructure g;
    customDetectIntersections(ts, g.intersectionList(), octree);

    g.initFromTriangleSoup(ts);

    classifyIntersections(ts, arena, g);

    triangulation(ts, arena, g, arr_out_tris, labels.surface, parallel);
    ts.appendJollyPoints();

    labels.inside.resize(arr_out_tris.size() / 3);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void customRemoveDegenerateAndDuplicatedTriangles(const std::vector<genericPoint *> &verts, std::vector<uint> &tris,
                                                  std::vector<std::bitset<NBIT>> &labels, std::vector<DuplTriInfo> &dupl_triangles,
                                                  bool parallel)
{
    if(parallel)
    {
        using vec3i = std::array<uint, 3>;
        uint num_orig_tris = static_cast<uint>(tris.size() / 3);
        vec3i* data_orig_tris = (vec3i*)tris.data();

        // compute colinear
        auto colinear = vector<bool>(num_orig_tris, false);
        tbb::parallel_for((uint)0, num_orig_tris, [data_orig_tris, &colinear, &verts](uint t_id) {
            auto& t = data_orig_tris[t_id];
            colinear[t_id] = cinolib::points_are_colinear_3d(
                verts[t[0]]->toExplicit3D().ptr(),
                verts[t[1]]->toExplicit3D().ptr(),
                verts[t[2]]->toExplicit3D().ptr());
        });

        // loop as before by use simpler way
        uint t_off = 0, l_off = 0;

        phmap::flat_hash_map < std::array<uint, 3>, std::pair<uint, uint> > tris_map; // tri_vertices -> <l_off, t_off>
        tris_map.reserve(num_orig_tris);

        for(uint t_id = 0; t_id < num_orig_tris; t_id++)
        {
            uint v0_id = tris[(3 * t_id)];
            uint v1_id = tris[(3 * t_id) +1];
            uint v2_id = tris[(3 * t_id) +2];
            std::bitset<NBIT> l = labels[t_id];

            if(!colinear[t_id]) // good triangle
            {
                std::array<uint, 3> tri = {v0_id, v1_id, v2_id};
                std::sort(tri.begin(), tri.end());

                auto ins= tris_map.insert({tri, std::make_pair(l_off, t_off)});

                if(ins.second) // first time for tri v0, v1, v2
                {
                    labels[l_off] = l;
                    l_off++;

                    tris[t_off] = v0_id, tris[t_off +1] = v1_id, tris[t_off +2] = v2_id;
                    t_off += 3;
                }
                else // triangle already present
                {
                    uint pos = ins.first->second.first;
                    labels[pos] |= l; // label for duplicates
                }

                if(!ins.second) // triangle already present -> save info about duplicates
                {
                    uint orig_tri_off = ins.first->second.second;

                    uint mesh_l = bitsetToUint(l);
                    assert(mesh_l >= 0);

                    uint curr_tri_verts[] = {v0_id, v1_id, v2_id};
                    uint orig_tri_verts[] ={tris[orig_tri_off], tris[orig_tri_off +1], tris[orig_tri_off +2]};

                    bool w = consistentWinding(curr_tri_verts, orig_tri_verts);

                    dupl_triangles.push_back({orig_tri_off / 3, // original triangle id
                                            static_cast<uint>(mesh_l), // label of the actual triangle
                                            w}); // winding with respect to the triangle stored in mesh (true -> same, false -> opposite)
                }
            }
        }

        tris.resize(t_off);
        labels.resize(l_off);
    } else {
        uint num_orig_tris = static_cast<uint>(tris.size() / 3);
        uint t_off = 0, l_off = 0;

        phmap::flat_hash_map < std::array<uint, 3>, std::pair<uint, uint> > tris_map; // tri_vertices -> <l_off, t_off>
        tris_map.reserve(num_orig_tris);

        for(uint t_id = 0; t_id < num_orig_tris; t_id++)
        {
            uint v0_id = tris[(3 * t_id)];
            uint v1_id = tris[(3 * t_id) +1];
            uint v2_id = tris[(3 * t_id) +2];
            std::bitset<NBIT> l = labels[t_id];

            if(!cinolib::points_are_colinear_3d(verts[v0_id]->toExplicit3D().ptr(),
                                                verts[v1_id]->toExplicit3D().ptr(),
                                                verts[v2_id]->toExplicit3D().ptr())) // good triangle
            {
                std::array<uint, 3> tri = {v0_id, v1_id, v2_id};
                std::sort(tri.begin(), tri.end());

                auto ins= tris_map.insert({tri, std::make_pair(l_off, t_off)});

                if(ins.second) // first time for tri v0, v1, v2
                {
                    labels[l_off] = l;
                    l_off++;

                    tris[t_off] = v0_id, tris[t_off +1] = v1_id, tris[t_off +2] = v2_id;
                    t_off += 3;
                }
                else // triangle already present
                {
                    uint pos = ins.first->second.first;
                    labels[pos] |= l; // label for duplicates
                }

                if(!ins.second) // triangle already present -> save info about duplicates
                {
                    uint orig_tri_off = ins.first->second.second;

                    uint mesh_l = bitsetToUint(l);
                    assert(mesh_l >= 0);

                    uint curr_tri_verts[] = {v0_id, v1_id, v2_id};
                    uint orig_tri_verts[] ={tris[orig_tri_off], tris[orig_tri_off +1], tris[orig_tri_off +2]};

                    bool w = consistentWinding(curr_tri_verts, orig_tri_verts);

                    dupl_triangles.push_back({orig_tri_off / 3, // original triangle id
                                            static_cast<uint>(mesh_l), // label of the actual triangle
                                            w}); // winding with respect to the triangle stored in mesh (true -> same, false -> opposite)
                }
            }
        }

        tris.resize(t_off);
        labels.resize(l_off);
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void customDetectIntersections(const TriangleSoup &ts, std::vector<std::pair<uint, uint> > &intersection_list, cinolib::Octree &o)
{
    std::vector<cinolib::vec3d> verts(ts.numVerts());

    for(uint v_id = 0; v_id < ts.numVerts(); v_id++)
        verts[v_id] = cinolib::vec3d(ts.vertX(v_id), ts.vertY(v_id), ts.vertZ(v_id));

    o.build_from_vectors(verts, ts.trisVector());

    intersection_list.reserve(ts.numTris());

    tbb::spin_mutex mutex;
    tbb::parallel_for((uint)0, (uint)o.leaves.size(), [&](uint i)
    {
        auto & leaf = o.leaves[i];
        if(leaf->item_indices.empty()) return;
        for(uint j=0;   j<leaf->item_indices.size()-1; ++j)
            for(uint k=j+1; k<leaf->item_indices.size();   ++k)
            {
                uint tid0 = leaf->item_indices[j];
                uint tid1 = leaf->item_indices[k];
                auto T0 = o.items[tid0];
                auto T1 = o.items[tid1];
                if(T0->aabb.intersects_box(T1->aabb)) // early reject based on AABB intersection
                {
                    const cinolib::Triangle *t0 = reinterpret_cast<cinolib::Triangle*>(T0);
                    const cinolib::Triangle *t1 = reinterpret_cast<cinolib::Triangle*>(T1);
                    if(t0->intersects_triangle(t1->v,true)) // precise check (exact if CINOLIB_USES_EXACT_PREDICATES is defined)
                    {
                        std::lock_guard<tbb::spin_mutex> guard(mutex);
                        intersection_list.push_back(cinolib::unique_pair(tid0,tid1));
                    }
                }
            }
    });
    remove_duplicates(intersection_list);
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void addDuplicateTrisInfoInStructures(const std::vector<DuplTriInfo> &dupl_tris, std::vector<uint> &in_tris,
                                             std::vector<std::bitset<NBIT>> &in_labels, cinolib::Octree &octree)
{
    for(auto &item : dupl_tris)
    {
        uint v0_id = in_tris[3 * item.t_id];
        uint v1_id = in_tris[3 * item.t_id + 1];
        uint v2_id = in_tris[3 * item.t_id + 2];

        std::bitset<NBIT> new_label;
        new_label[item.l_id] = true;

        const cinolib::Triangle *orig_tri = dynamic_cast<const cinolib::Triangle*>(octree.items.at(item.t_id));

        uint new_t_id = in_tris.size() / 3;

        if(item.w)
        {
            in_tris.push_back(v0_id);
            in_tris.push_back(v1_id);
            in_tris.push_back(v2_id);
            octree.items.push_back(new cinolib::Triangle(new_t_id, orig_tri->v[0], orig_tri->v[1], orig_tri->v[2]));
        }
        else
        {
            in_tris.push_back(v0_id);
            in_tris.push_back(v2_id);
            in_tris.push_back(v1_id);
            octree.items.push_back(new cinolib::Triangle(new_t_id, orig_tri->v[0], orig_tri->v[1], orig_tri->v[2]));
        }

        in_labels.push_back(new_label); // we add the new_label to the new_triangle
        in_labels[item.t_id][item.l_id] = false; // we remove the dupl label from the orig triangle
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void computeAllPatches(FastTrimesh &tm, const Labels &labels, std::vector<phmap::flat_hash_set<uint>> &patches, bool parallel)
{
    if(parallel) {
        tm.resetVerticesInfo();
        auto adjT2E = tm.adjT2EAll(parallel);

        for(uint t_id = 0; t_id < tm.numTris(); t_id++)
        {
            if(tm.triInfo(t_id) != 1)
            {
                patches.emplace_back();
                computeSinglePatch(tm, t_id, labels, patches.back(), adjT2E);
            }
        }
    } else {
        tm.resetVerticesInfo();

        for(uint t_id = 0; t_id < tm.numTris(); t_id++)
        {
            if(tm.triInfo(t_id) != 1)
            {
                patches.emplace_back();
                computeSinglePatch(tm, t_id, labels, patches.back());
            }
        }
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void computeSinglePatch(FastTrimesh &tm, uint seed_t, const Labels &labels, phmap::flat_hash_set<uint> &patch)
{
    std::bitset<NBIT> ref_l = labels.surface[seed_t];

    std::stack<uint> tris_stack;
    tris_stack.push(seed_t);

    while(!tris_stack.empty())
    {
        uint curr_t = tris_stack.top();
        tris_stack.pop();

        tm.setTriInfo(curr_t, 1); // set triangle as visited
        patch.insert(curr_t);

        for(uint e_id : tm.adjT2E(curr_t))
        {
            if(tm.edgeIsManifold(e_id))
            {
                for(uint t_id : tm.adjE2T(e_id))
                {
                    if(t_id != curr_t && tm.triInfo(t_id) != 1)
                    {
                        assert(labels.surface[t_id] == ref_l);
                        tris_stack.push(t_id);
                    }
                }
            }
            else // e_id is not manifold -> stop flooding
            {
                // we set the vertices in the patch border with 1 (useful for ray computation funcion)
                tm.setVertInfo(tm.edgeVertID(e_id, 0), 1);
                tm.setVertInfo(tm.edgeVertID(e_id, 1), 1);
            }
        }
    }
}

inline void computeSinglePatch(FastTrimesh &tm, uint seed_t, const Labels &labels, phmap::flat_hash_set<uint> &patch, const std::vector<std::array<uint, 3>>& adjT2E)
{
    std::bitset<NBIT> ref_l = labels.surface[seed_t];

    std::stack<uint> tris_stack;
    tris_stack.push(seed_t);

    while(!tris_stack.empty())
    {
        uint curr_t = tris_stack.top();
        tris_stack.pop();

        tm.setTriInfo(curr_t, 1); // set triangle as visited
        patch.insert(curr_t);

        for(uint e_id : adjT2E[curr_t])
        {
            if(tm.edgeIsManifold(e_id))
            {
                for(uint t_id : tm.adjE2T(e_id))
                {
                    if(t_id != curr_t && tm.triInfo(t_id) != 1)
                    {
                        assert(labels.surface[t_id] == ref_l);
                        tris_stack.push(t_id);
                    }
                }
            }
            else // e_id is not manifold -> stop flooding
            {
                // we set the vertices in the patch border with 1 (useful for ray computation funcion)
                tm.setVertInfo(tm.edgeVertID(e_id, 0), 1);
                tm.setVertInfo(tm.edgeVertID(e_id, 1), 1);
            }
        }
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void findRayEndpoints(const FastTrimesh &tm, const phmap::flat_hash_set<uint> &patch, const cinolib::vec3d &max_coords, Ray &ray)
{
    // check for an explicit point (all operations with explicits are faster)
    int v_id = -1;
    for(uint t_id : patch)
    {
        const uint tv[3] = {tm.triVertID(t_id, 0), tm.triVertID(t_id, 1), tm.triVertID(t_id, 2)};

        if (tm.vert(tv[0])->isExplicit3D() && tm.vertInfo(tv[0]) == 0)      v_id = static_cast<int>(tv[0]);
        else if (tm.vert(tv[1])->isExplicit3D() && tm.vertInfo(tv[1]) == 0) v_id = static_cast<int>(tv[1]);
        else if (tm.vert(tv[2])->isExplicit3D() && tm.vertInfo(tv[2]) == 0) v_id = static_cast<int>(tv[2]);

        if (v_id != -1)
        {
            const explicitPoint3D &v = tm.vert(v_id)->toExplicit3D();
            ray.v0  = explicitPoint3D(v.X(), v.Y(), v.Z());
            ray.v1 = explicitPoint3D(max_coords.x(), v.Y(), v.Z());
            return;
        }
    }

    int tri_counter = 0;
    // parse triangles with all implicit points
    for(uint t_id : patch)
    {
        double x0, x1, x2, y0, y1, y2, z0, z1, z2;
        tm.triVert(t_id, 0)->getApproxXYZCoordinates(x0, y0, z0);
        tm.triVert(t_id, 1)->getApproxXYZCoordinates(x1, y1, z1);
        tm.triVert(t_id, 2)->getApproxXYZCoordinates(x2, y2, z2);

        explicitPoint3D tv0(x0, y0, z0), tv1(x1, y1, z1), tv2(x2, y2, z2);
        if(!genericPoint::misaligned(tv0, tv1, tv2)) continue;

        int dir = genericPoint::maxComponentInTriangleNormal(x0, y0, z0, x1, y1, z1, x2, y2, z2);
        if(dir == 0) // dir = X
        {
            ray.v0 = explicitPoint3D(((x0 + x1 + x2) / 3.0) - 0.1, (y0 + y1 + y2) / 3.0, (z0 + z1 + z2) / 3.0);
            ray.v1 = explicitPoint3D(max_coords.x(), ray.v0.Y(), ray.v0.Z());
            ray.dir = 'X';
        }
        else if(dir == 1) // dir = Y
        {
            ray.v0 = explicitPoint3D((x0 + x1 + x2) / 3.0, ((y0 + y1 + y2) / 3.0) -0.1, (z0 + z1 + z2) / 3.0);
            ray.v1 = explicitPoint3D(ray.v0.X(), max_coords.y(), ray.v0.Z());
            ray.dir = 'Y';
        }
        else // dir = Z
        {
            ray.v0 = explicitPoint3D((x0 + x1 + x2) / 3.0, (y0 + y1 + y2) / 3.0, ((z0 + z1 + z2) / 3.0)  -0.1);
            ray.v1 = explicitPoint3D(ray.v0.X(), ray.v0.Y(), max_coords.z());
            ray.dir = 'Z';
        }

        int orf = genericPoint::orient3D(*tm.triVert(t_id, 0), *tm.triVert(t_id, 1), *tm.triVert(t_id, 2), ray.v0);
        int ors = genericPoint::orient3D(*tm.triVert(t_id, 0), *tm.triVert(t_id, 1), *tm.triVert(t_id, 2), ray.v1);

        if((orf < 0 && ors > 0) || (orf > 0 && ors < 0)) // the ray passes through the triangle
        {
            if(checkIntersectionInsideTriangle3DImplPoints(ray, tm.triVert(t_id, 0), tm.triVert(t_id, 1), tm.triVert(t_id, 2))) // the ray passes inside the triangle
            {
                ray.tv[0] = static_cast<int>(tm.triVertID(t_id, 0));
                ray.tv[1] = static_cast<int>(tm.triVertID(t_id, 1));
                ray.tv[2] = static_cast<int>(tm.triVertID(t_id, 2));
                return;
            }
        }
    }

    std::cout << "WARNING: the arrangement contains a fully implicit patch that requires exact rationals for evaluation. This version of the code does not support rationals, therefore the output result may contain open boundaries." << std::endl;

    std::exit(EXIT_FAILURE);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline bool intersects_box(const cinolib::Octree& tree, const cinolib::AABB & b, phmap::flat_hash_set<uint> & ids)
{
    auto root = tree.root;
    auto& items = tree.items;

    std::stack<cinolib::OctreeNode*> lifo;
    if(root && root->bbox.intersects_box(b))
    {
        lifo.push(root);
    }

    while(!lifo.empty())
    {
        cinolib::OctreeNode *node = lifo.top();
        lifo.pop();
        assert(node->bbox.intersects_box(b));

        if(node->is_inner())
        {            
            for(int i = 0; i < 8; ++i)
            {
                if(node->children[i]->bbox.intersects_box(b))
                {
                    lifo.push(node->children[i]);
                }
            }
        }
        else
        {
            for(uint i : node->item_indices)
            {
                if(items.at(i)->aabb.intersects_box(b))
                {
                    ids.insert(items.at(i)->id);
                }
            }
        }
    }

    return !ids.empty();
}

inline void computeInsideOut(const FastTrimesh &tm, const std::vector<phmap::flat_hash_set<uint>> &patches, const cinolib::Octree &octree,
                             const std::vector<genericPoint *> &in_verts, const std::vector<uint> &in_tris,
                             const std::vector<std::bitset<NBIT>> &in_labels, const cinolib::vec3d &max_coords, Labels &labels)
{
    tbb::spin_mutex mutex;
    tbb::parallel_for((uint)0, (uint)patches.size(), [&](uint p_id)
   {
        const phmap::flat_hash_set<uint> &patch_tris = patches[p_id];
        const std::bitset<NBIT> &patch_surface_label = labels.surface[*patch_tris.begin()]; // label of the first triangle of the patch

        Ray ray;
        findRayEndpoints(tm, patch_tris, max_coords, ray);

        // find all the triangles having a bbox intersected by the ray
        phmap::flat_hash_set<uint> tmp_inters;
        cinolib::AABB rayAABB(cinolib::vec3d(ray.v0.X(), ray.v0.Y(), ray.v0.Z()),
                              cinolib::vec3d(ray.v1.X(), ray.v1.Y(), ray.v1.Z()));

        intersects_box(octree, rayAABB, tmp_inters);

        std::vector<uint> sorted_inters;
        pruneIntersectionsAndSortAlongRay(ray, in_verts, in_tris, in_labels, tmp_inters, patch_surface_label,
                                          sorted_inters);

        std::bitset<NBIT> patch_inner_label;
        analyzeSortedIntersections(ray, in_verts, in_tris, in_labels, sorted_inters, patch_inner_label);

        propagateInnerLabelsOnPatch(patch_tris, patch_inner_label, labels);
    });
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void pruneIntersectionsAndSortAlongRay(const Ray &ray, const std::vector<genericPoint*> &in_verts,
                                              const std::vector<uint> &in_tris, const std::vector<std::bitset<NBIT>> &in_labels,
                                              const phmap::flat_hash_set<uint> &tmp_inters, const std::bitset<NBIT> &patch_surface_label,
                                              std::vector<uint> &inters_tris)
{
    phmap::flat_hash_set<uint> visited_tri;
    visited_tri.reserve(tmp_inters.size()/6);
    std::pair<phmap::flat_hash_set<uint>::iterator, bool> ins;

    for(uint t_id : tmp_inters)
    {
        ins = visited_tri.insert(t_id);
        if(!ins.second) continue; // triangle already analyzed or in the one ring of a vert or in the adj of an edge

        const std::bitset<NBIT> tested_tri_label = in_labels[t_id];
        uint uint_tri_label = bitsetToUint(tested_tri_label);
        if(patch_surface_label[uint_tri_label]) continue; // <-- triangle of the same label of the tested patch

        const explicitPoint3D &tv0 = in_verts[in_tris[3 * t_id]]->toExplicit3D();
        const explicitPoint3D &tv1 = in_verts[in_tris[3 * t_id +1]]->toExplicit3D();
        const explicitPoint3D &tv2 = in_verts[in_tris[3 * t_id +2]]->toExplicit3D();

        IntersInfo ii = fast2DCheckIntersectionOnRay(ray, tv0, tv1, tv2);

        if(ii == DISCARD || ii == NO_INT) continue;

        if(ii == INT_IN_TRI)
        {
            inters_tris.push_back(t_id);
        }
        else if(ii == INT_IN_V0 || ii == INT_IN_V1 || ii == INT_IN_V2)
        {
            uint v_id;
            if(ii == INT_IN_V0) v_id = in_tris[3 * t_id];
            else if(ii == INT_IN_V1) v_id = in_tris[3 * t_id +1];
            else v_id = in_tris[3 * t_id +2];

            std::vector<uint> vert_one_ring;
            findVertRingTris(v_id, tested_tri_label, tmp_inters, in_tris, in_labels, vert_one_ring);

            for(uint t : vert_one_ring)
                visited_tri.insert(t); // mark all the one ring as visited

            int winner_tri = -1;
            winner_tri = perturbRayAndFindIntersTri(ray, in_verts, in_tris, vert_one_ring); // the first inters triangle after ray perturbation

            if(winner_tri != -1)
                inters_tris.push_back(winner_tri);
        }
        else if(ii == INT_IN_EDGE01 || ii == INT_IN_EDGE12 || ii == INT_IN_EDGE20)
        {
            uint ev0_id, ev1_id;
            if(ii == INT_IN_EDGE01)
            {
                ev0_id = in_tris[3 * t_id];
                ev1_id = in_tris[3 * t_id +1];
            }
            else if(ii == INT_IN_EDGE12)
            {
                ev0_id = in_tris[3 * t_id +1];
                ev1_id = in_tris[3 * t_id +2];
            }
            else
            {
                ev0_id = in_tris[3 * t_id +2];
                ev1_id = in_tris[3 * t_id];
            }

            std::vector<uint> edge_tris;
            findEdgeTris(ev0_id, ev1_id, tested_tri_label, tmp_inters, in_tris, in_labels, edge_tris);

            for(uint t : edge_tris)
                visited_tri.insert(t); // mark all the one ring as visited

            int winner_tri = -1;
            winner_tri = perturbRayAndFindIntersTri(ray, in_verts, in_tris, edge_tris);

            if(winner_tri != -1)
                inters_tris.push_back(winner_tri);
        }
    }

    if(ray.dir == 'X')
        sortIntersectedTrisAlongX(ray, in_verts, in_tris, inters_tris);
    else if(ray.dir == 'Y')
        sortIntersectedTrisAlongY(ray, in_verts, in_tris, inters_tris);
    else
        sortIntersectedTrisAlongZ(ray, in_verts, in_tris, inters_tris);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void analyzeSortedIntersections(const Ray &ray, const std::vector<genericPoint*> &in_verts, const std::vector<uint> &in_tris,
                                       const std::vector<std::bitset<NBIT>> &in_labels, const std::vector<uint> &sorted_inters,
                                       std::bitset<NBIT> &patch_inner_label)
{
    std::bitset<NBIT> visited_labels;

    for(uint t_id : sorted_inters)
    {
        uint t_label = bitsetToUint(in_labels[t_id]);
        if(visited_labels[t_label]) continue; // already visited patch

        const explicitPoint3D &tv0 = in_verts[in_tris[3 * t_id]]->toExplicit3D();
        const explicitPoint3D &tv1 = in_verts[in_tris[3 * t_id +1]]->toExplicit3D();
        const explicitPoint3D &tv2 = in_verts[in_tris[3 * t_id +2]]->toExplicit3D();

        if(checkTriangleOrientation(ray, tv0, tv1, tv2) == 1) // checkOrientation -> 1 if inside, 0 if outside
            patch_inner_label[t_label] = true;

        visited_labels[t_label] = true;
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline bool triContainsVert(uint t_id, uint v_id, const std::vector<uint> &in_tris)
{
    if(in_tris[3 * t_id]     == v_id) return true;
    if(in_tris[3 * t_id + 1] == v_id) return true;
    if(in_tris[3 * t_id + 2] == v_id) return true;
    return false;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void findVertRingTris(uint v_id, const std::bitset<NBIT> &ref_label, const phmap::flat_hash_set<uint> &inters_tris,
                             const std::vector<uint> &in_tris, const std::vector<std::bitset<NBIT>> &in_labels,
                             std::vector<uint> &one_ring)
{
    for(uint t_id : inters_tris)
    {
        if(in_labels[t_id] == ref_label && triContainsVert(t_id, v_id, in_tris))
            one_ring.push_back(t_id);
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void findEdgeTris(uint ev0_id, uint ev1_id, const std::bitset<NBIT> &ref_label, const phmap::flat_hash_set<uint> &inters_tris,
                         const std::vector<uint> &in_tris, const std::vector<std::bitset<NBIT>> &in_labels,
                         std::vector<uint> &edge_tris)
{
    for(auto t_id : inters_tris)
    {
        if(in_labels[t_id] == ref_label && triContainsVert(t_id, ev0_id, in_tris) && triContainsVert(t_id, ev1_id, in_tris))
            edge_tris.push_back(t_id);
    }

    assert(edge_tris.size() == 2 && "problem in finding edge triangles"); // always true in closed and manifold meshes
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline Ray perturbXRay(const Ray &ray, uint offset) // offset is used to perturb the ray in all the possible directions
{
    Ray new_ray = ray;

    switch (offset)
    {
        case 0: // -> +y
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, ray.v1.Z());
        } break;

        case 1: // -> +y+z
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() + 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, new_z);
        } break;

        case 2: // -> +z
        {
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), ray.v1.Y(), new_z);
        } break;

        case 3: //-> -y+z
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() - 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, new_z);
        } break;

        case 4: // -> -y
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, ray.v1.Z());
        } break;

        case 5: // -> -y-z
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() - 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, new_z);
        } break;

        case 6: // -> -z
        {
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), ray.v1.Y(), new_z);
        } break;

        case 7: // -> +y-z
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() + 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, new_z);
        } break;

        default:
        {
            assert(false && "non-valid offset value");
        } break;
    }

    return new_ray;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline Ray perturbYRay(const Ray &ray, uint offset)
{
    Ray new_ray = ray;

    switch (offset)
    {
        case 0: // -> +x
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() + 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), ray.v1.Z());
        } break;

        case 1: // -> +x+z
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() + 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() + 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), new_z);
        } break;

        case 2: // -> +z
        {
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), ray.v1.Y(), new_z);
        } break;

        case 3: //-> -x+z
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() - 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() + 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), new_z);
        } break;

        case 4: // -> -x
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() - 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), ray.v1.Z());
        } break;

        case 5: // -> -x-z
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() - 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() - 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), new_z);
        } break;

        case 6: // -> -z
        {
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), ray.v1.Y(), new_z);
        } break;

        case 7: // -> +x-z
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() + 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() - 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), new_z);
        } break;

        default:
        {
            assert(false && "non-valid offset value");
        } break;
    }

    return new_ray;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline Ray perturbZRay(const Ray &ray, uint offset)
{
    Ray new_ray = ray;

    switch (offset)
    {
        case 0: // -> +x
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() + 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), ray.v1.Z());
        } break;

        case 1: // -> +x+y
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() + 1.0));
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() + 1.0));
            new_ray.v1 = explicitPoint3D(new_x, new_y, ray.v1.Z());
        } break;

        case 2: // -> +y
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, ray.v1.Z());
        } break;

        case 3: //-> -x+y
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() - 1.0));
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() + 1.0));
            new_ray.v1 = explicitPoint3D(new_x, new_y, ray.v1.Z());
        } break;

        case 4: // -> -x
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() - 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), ray.v1.Z());
        } break;

        case 5: // -> -x-y
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() - 1.0));
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() - 1.0));
            new_ray.v1 = explicitPoint3D(new_x, new_y, ray.v1.Z());
        } break;

        case 6: // -> -y
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, ray.v1.Z());
        } break;

        case 7: // -> +x-y
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() + 1.0));
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() - 1.0));
            new_ray.v1 = explicitPoint3D(new_x, new_y, ray.v1.Z());
        } break;

        default:
        {
            assert(false && "non-valid offset value");
        } break;
    }

    return new_ray;

}

inline int perturbRayAndFindIntersTri(const Ray &ray, const std::vector<genericPoint*> &in_verts, const std::vector<uint> &in_tris,
                                       const std::vector<uint> &tris_to_test)
{
    std::vector<uint> inters_tris;
    Ray p_ray;

    for(uint i = 0; i <= 7; i++)
    {
        if(ray.dir == 'X')      p_ray = perturbXRay(ray, i);
        else if(ray.dir == 'Y') p_ray = perturbYRay(ray, i);
        else if(ray.dir == 'Z') p_ray = perturbZRay(ray, i);

        for(uint t_id : tris_to_test)
        {
            const explicitPoint3D &tv0 = in_verts[in_tris[3 * t_id]]->toExplicit3D();
            const explicitPoint3D &tv1 = in_verts[in_tris[3 * t_id +1]]->toExplicit3D();
            const explicitPoint3D &tv2 = in_verts[in_tris[3 * t_id +2]]->toExplicit3D();

            if(checkIntersectionInsideTriangle3D(p_ray, tv0, tv1, tv2))
                inters_tris.push_back(t_id);

            if(!inters_tris.empty()) break;
        }
    }

    if(inters_tris.empty())
        return -1;

    if(ray.dir == 'X')      sortIntersectedTrisAlongX(p_ray, in_verts, in_tris, inters_tris);
    else if(ray.dir == 'Y') sortIntersectedTrisAlongY(p_ray, in_verts, in_tris, inters_tris);
    else                    sortIntersectedTrisAlongZ(p_ray, in_verts, in_tris, inters_tris);

    return static_cast<int>(inters_tris[0]); // return the first triangle intersected
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline IntersInfo fast2DCheckIntersectionOnRay(const Ray &ray, const explicitPoint3D &tv0, const explicitPoint3D &tv1, const explicitPoint3D &tv2)
{
    double v0[2], v1[2], v2[2], vq[2];

    switch (ray.dir)
    {
        case 'X': // only YZ coordinates
        {
            v0[0] = tv0.Y(); v0[1] = tv0.Z();
            v1[0] = tv1.Y(); v1[1] = tv1.Z();
            v2[0] = tv2.Y(); v2[1] = tv2.Z();
            vq[0] = ray.v1.Y(); vq[1] = ray.v1.Z();
        } break;

        case 'Y': //only XZ coordinates
        {
            v0[0] = tv0.X(); v0[1] = tv0.Z();
            v1[0] = tv1.X(); v1[1] = tv1.Z();
            v2[0] = tv2.X(); v2[1] = tv2.Z();
            vq[0] = ray.v1.X(); vq[1] = ray.v1.Z();
        } break;

        case 'Z': //only XY coordinates
        {
            v0[0] = tv0.X(); v0[1] = tv0.Y();
            v1[0] = tv1.X(); v1[1] = tv1.Y();
            v2[0] = tv2.X(); v2[1] = tv2.Y();
            vq[0] = ray.v1.X(); vq[1] = ray.v1.Y();
        } break;

    }

    double or01 = cinolib::orient2d(v0, v1, vq);
    double or12 = cinolib::orient2d(v1, v2, vq);
    double or20 = cinolib::orient2d(v2, v0, vq);

    if((or01 >= 0 && or12 >= 0 && or20 >= 0) || (or01 <= 0 && or12 <= 0 && or20 <= 0))
    {
        // check if the the ray passes through a vert
        if(v0[0] == vq[0] && v0[1] == vq[1]) return INT_IN_V0;
        if(v1[0] == vq[0] && v1[1] == vq[1]) return INT_IN_V1;
        if(v2[0] == vq[0] && v2[1] == vq[1]) return INT_IN_V2;

        // check if the triangle is coplanar with the ray
        if(or01 == 0 && or12 == 0) return DISCARD;
        if(or12 == 0 && or20 == 0) return DISCARD;
        if(or20 == 0 && or01 == 0) return DISCARD;

        // check if the ray passes through an edge
        if(or01 == 0) return INT_IN_EDGE01;
        if(or12 == 0) return INT_IN_EDGE12;
        if(or20 == 0) return INT_IN_EDGE20;

        return INT_IN_TRI; // so the triangle intersect insede the triangle area
    }

    return NO_INT;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline bool checkIntersectionInsideTriangle3D(const Ray &ray, const explicitPoint3D &tv0, const explicitPoint3D &tv1, const explicitPoint3D &tv2)
{
    // we check the orientation of ray.v1 with respct to the planes v0-v1-ray.v0, v1-v2-ray.v0, v2-v0-ray.v0
    double or01f = cinolib::orient3d(tv0.ptr(), tv1.ptr(), ray.v0.ptr(), ray.v1.ptr());
    double or12f = cinolib::orient3d(tv1.ptr(), tv2.ptr(), ray.v0.ptr(), ray.v1.ptr());
    double or20f = cinolib::orient3d(tv2.ptr(), tv0.ptr(), ray.v0.ptr(), ray.v1.ptr());

    if(or01f > 0 && or12f > 0 && or20f > 0) return true;
    if(or01f < 0 && or12f < 0 && or20f < 0) return true;

    return false;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline bool checkIntersectionInsideTriangle3DImplPoints(const Ray &ray, const genericPoint *tv0, const genericPoint *tv1, const genericPoint *tv2)
{
    // we check the orientation of ray.v1 with respct to the planes v0-v1-ray.v0, v1-v2-ray.v0, v2-v0-ray.v0
    double or01f = genericPoint::orient3D(*tv0, *tv1, ray.v0, ray.v1);
    double or12f = genericPoint::orient3D(*tv1, *tv2, ray.v0, ray.v1);
    double or20f = genericPoint::orient3D(*tv2, *tv0, ray.v0, ray.v1);

    if(or01f > 0 && or12f > 0 && or20f > 0) return true;
    if(or01f < 0 && or12f < 0 && or20f < 0) return true;

    return false;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// sort all intersected triangles from ray.v0 to ray.v1 (intersections before ray.v0 are discarded)
inline void sortIntersectedTrisAlongX(const Ray &ray, const std::vector<genericPoint*> &in_verts,
                                      const std::vector<uint> &in_tris, std::vector<uint> &inters_tris)
{
    phmap::btree_set< std::pair<genericPoint*, uint>, less_than_GP_on_X > inters_set; // <- <t_id, impl_point>
    std::vector<implicitPoint3D_LPI> arena;
    arena.reserve(inters_tris.size());

    for(uint t_id : inters_tris)
    {
        uint v0_id = in_tris[3 * t_id];
        uint v1_id = in_tris[3 * t_id +1];
        uint v2_id = in_tris[3 * t_id +2];

        std::pair<genericPoint*, uint> pair;
        pair.first = &arena.emplace_back(ray.v0, ray.v1,in_verts[v0_id]->toExplicit3D(),
                                         in_verts[v1_id]->toExplicit3D(), in_verts[v2_id]->toExplicit3D());
        pair.second = t_id;

        inters_set.insert(pair);
    }

    inters_tris.clear();
    auto curr_int = inters_set.begin();

    // we discard the intersection before ray.first along X
    if(ray.tv[0] != -1) // the ray is generated
    {
        const genericPoint *tv0 = in_verts[ray.tv[0]];
        const genericPoint *tv1 = in_verts[ray.tv[1]];
        const genericPoint *tv2 = in_verts[ray.tv[2]];

        if(genericPoint::orient3D(*tv0, *tv1, *tv2, ray.v1) > 0)
        {
            while(curr_int != inters_set.end() && genericPoint::orient3D(*tv0, *tv1, *tv2, *curr_int->first) < 0)
                curr_int++;
        }
        else
        {
            while(curr_int != inters_set.end() && genericPoint::orient3D(*tv0, *tv1, *tv2, *curr_int->first) > 0)
                curr_int++;
        }
    }
    else // the ray is composed of 2 real explicit points
    {
        while(curr_int != inters_set.end() && genericPoint::lessThanOnX(*curr_int->first, ray.v0) < 0)
            curr_int++;
    }

    // we save all the intersecting triangles from ray.first to ray.second
    while(curr_int != inters_set.end())
    {
        inters_tris.push_back(curr_int->second);
        curr_int++;
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void sortIntersectedTrisAlongY(const Ray &ray, const std::vector<genericPoint*> &in_verts,
                                      const std::vector<uint> &in_tris, std::vector<uint> &inters_tris)
{
    phmap::btree_set< std::pair<genericPoint*, uint>, less_than_GP_on_Y > inters_set; // <- <t_id, impl_point>
    std::vector<implicitPoint3D_LPI> arena;
    arena.reserve(inters_tris.size());

    for(uint t_id : inters_tris)
    {
        uint v0_id = in_tris[3 * t_id];
        uint v1_id = in_tris[3 * t_id +1];
        uint v2_id = in_tris[3 * t_id +2];

        std::pair<genericPoint*, uint> pair;
        pair.first = &arena.emplace_back(ray.v0, ray.v1,in_verts[v0_id]->toExplicit3D(), in_verts[v1_id]->toExplicit3D(), in_verts[v2_id]->toExplicit3D());
        pair.second = t_id;

        inters_set.insert(pair);
    }

    inters_tris.clear();
    auto curr_int = inters_set.begin();

    // we discard the intersection before ray.first along Y
    if(ray.tv[0] != -1) // the ray is generated
    {
        const genericPoint *tv0 = in_verts[ray.tv[0]];
        const genericPoint *tv1 = in_verts[ray.tv[1]];
        const genericPoint *tv2 = in_verts[ray.tv[2]];

        if(genericPoint::orient3D(*tv0, *tv1, *tv2, ray.v1) > 0)
        {
            while(curr_int != inters_set.end() && genericPoint::orient3D(*tv0, *tv1, *tv2, *curr_int->first) < 0)
                curr_int++;
        }
        else
        {
            while(curr_int != inters_set.end() && genericPoint::orient3D(*tv0, *tv1, *tv2, *curr_int->first) > 0)
                curr_int++;
        }
    }
    else // the ray is composed of 2 real explicit points
    {
        while(curr_int != inters_set.end() && genericPoint::lessThanOnY(*curr_int->first, ray.v0) < 0)
            curr_int++;
    }

    // we save all the intersecting triangles from ray.first to ray.second
    while(curr_int != inters_set.end())
    {
        inters_tris.push_back(curr_int->second);
        curr_int++;
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void sortIntersectedTrisAlongZ(const Ray &ray, const std::vector<genericPoint*> &in_verts,
                                      const std::vector<uint> &in_tris, std::vector<uint> &inters_tris)
{
    phmap::btree_set< std::pair<genericPoint*, uint>, less_than_GP_on_Z> inters_set; // <- <t_id, impl_point>
    std::vector<implicitPoint3D_LPI> arena;
    arena.reserve(inters_tris.size());

    for(uint t_id : inters_tris)
    {
        uint v0_id = in_tris[3 * t_id];
        uint v1_id = in_tris[3 * t_id +1];
        uint v2_id = in_tris[3 * t_id +2];

        std::pair<genericPoint*, uint> pair;
        pair.first = &arena.emplace_back(ray.v0, ray.v1,in_verts[v0_id]->toExplicit3D(), in_verts[v1_id]->toExplicit3D(), in_verts[v2_id]->toExplicit3D());
        pair.second = t_id;

        inters_set.insert(pair);
    }

    inters_tris.clear();
    auto curr_int = inters_set.begin();

    // we discard the intersection before ray.first along Z
    if(ray.tv[0] != -1) // the ray is generated
    {
        const genericPoint *tv0 = in_verts[ray.tv[0]];
        const genericPoint *tv1 = in_verts[ray.tv[1]];
        const genericPoint *tv2 = in_verts[ray.tv[2]];

        if(genericPoint::orient3D(*tv0, *tv1, *tv2, ray.v1) > 0)
        {
            while(curr_int != inters_set.end() && genericPoint::orient3D(*tv0, *tv1, *tv2, *curr_int->first) < 0)
                curr_int++;
        }
        else
        {
            while(curr_int != inters_set.end() && genericPoint::orient3D(*tv0, *tv1, *tv2, *curr_int->first) > 0)
                curr_int++;
        }
    }
    else // the ray is composed of 2 real explicit points
    {
        while(curr_int != inters_set.end() && genericPoint::lessThanOnZ(*curr_int->first, ray.v0) < 0)
            curr_int++;
    }

    // we save all the intersecting triangles from ray.first to ray.second
    while(curr_int != inters_set.end())
    {
        inters_tris.push_back(curr_int->second);
        curr_int++;
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// return 1 if inside, 0 if outside
inline uint checkTriangleOrientation(const Ray &ray, const explicitPoint3D &tv0, const explicitPoint3D &tv1, const explicitPoint3D &tv2)
{
    double res = cinolib::orient3d(tv0.ptr(), tv1.ptr(), tv2.ptr(), ray.v1.ptr());

    assert(res != 0 && "Problem in PointOrientation(...)");

    /* in res we have sign(area(v0, v1, v2, ray.second))
     * if the area is >0 the ray is doing INSIDE -> OUTSIDE, so the patch is INSIDE
     * else the ray is doing OUTSIDE -> INSIDE so the patch is OUTSIDE */
    return (res < 0) ? 1 : 0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void propagateInnerLabelsOnPatch(const phmap::flat_hash_set<uint> &patch_tris, const std::bitset<NBIT> &patch_inner_label, Labels &labels)
{
    for(uint t_id : patch_tris)
        labels.inside[t_id] = patch_inner_label;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void computeFinalExplicitResult(const FastTrimesh &tm, const Labels &labels, uint num_tris_in_final_res,
                                       std::vector<double> &out_coords, std::vector<uint> &out_tris, 
                                       std::vector<std::bitset<NBIT>> &out_label, bool flat_array, Data &data)
{
    if(flat_array)
    {
        // loop over triangles and fix vertex indices
        uint num_vertices = 0;
        std::vector<int>  vertex_index(tm.numVerts(), -1);
        out_tris.resize(3 * num_tris_in_final_res);
        out_label.resize(num_tris_in_final_res);
        uint tri_offset = 0;
        for(uint t_id = 0; t_id < tm.numTris(); t_id++)
        {
            if(tm.triInfo(t_id) == 0){
                continue;
            } // triangle not included in final version
            const uint *triangle = tm.tri(t_id);
            for(uint i = 0; i < 3; i++)
            {
                uint old_vertex = triangle[i];
                if (vertex_index[old_vertex] == -1) {
                    vertex_index[old_vertex] = num_vertices ++;
                }
                out_tris[3 * tri_offset + i] = vertex_index[old_vertex];
            }
            out_label[tri_offset] = labels.surface[t_id];

            if (auto it = std::find(data.t_ids_debug.begin(), data.t_ids_debug.end(), t_id); it != data.t_ids_debug.end()) {
                std::cout <<"Original t_id: " << *it << " new: " << tri_offset << std::endl;}

            tri_offset++;
        }

        // loop over vertices
        out_coords.resize(3 * num_vertices);
        for(uint v_id = 0; v_id < (uint)tm.numVerts(); v_id++) {
            if (vertex_index[v_id] == -1) continue;
            double* v = out_coords.data() + (3 * vertex_index[v_id]);
            tm.vert(v_id)->getApproxXYZCoordinates(v[0], v[1], v[2]);
        }

        // rescale output
        double multiplier = tm.vert(tm.numVerts() - 1)->toExplicit3D().X();
        for(double &c : out_coords) c /= multiplier;
    } else
    {
        out_coords.reserve(3 * 3 * num_tris_in_final_res);
        out_tris.resize(3 * num_tris_in_final_res);
        out_label.resize(num_tris_in_final_res);
        phmap::flat_hash_map<uint, uint> v_map;
        uint tri_offset = 0;

        double multiplier = tm.vert(tm.numVerts() - 1)->toExplicit3D().X();

        for(uint t_id = 0; t_id < tm.numTris(); t_id++)
        {
            if(tm.triInfo(t_id) == 0) continue; // triangle not included in final version

            const uint *v_id = tm.tri(t_id);

            for(uint i = 0; i < 3; i++)
            {
                uint fresh_v_id = (uint)(out_coords.size() / 3);
                auto ins = v_map.insert({v_id[i], fresh_v_id});

                if(ins.second) // vert added
                {
                    double x, y, z;
                    tm.vert(v_id[i])->getApproxXYZCoordinates(x, y, z);
                    out_coords.push_back(x);
                    out_coords.push_back(y);
                    out_coords.push_back(z);
                }

                out_tris[3 * tri_offset + i] = ins.first->second;
            }
            out_label[tri_offset] = labels.surface[t_id];
            tri_offset++;
        }

        for(double &c : out_coords)
            c /= multiplier;
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline uint boolIntersection(FastTrimesh &tm, const Labels &labels, Data &data)
{
    uint num_tris_in_final_solution = 0;
    tm.resetTrianglesInfo();

    for(uint t_id = 0; t_id < tm.numTris(); t_id++)
    {
        if((labels.surface[t_id] ^ labels.inside[t_id]).count() == labels.num) // triangle to keep
        {
            data.t_ids_intersection.push_back(t_id);
            tm.setTriInfo(t_id, 1);
            num_tris_in_final_solution++;
        }
    }

    return num_tris_in_final_solution;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline uint boolUnion(FastTrimesh &tm, const Labels &labels, Data &data)
{
    uint num_tris_in_final_solution = 0;
    tm.resetTrianglesInfo();

    for(uint t_id = 0; t_id < tm.numTris(); t_id++)
    {
        if(labels.inside[t_id].count() == 0) // triangle to keep
        {
            data.t_ids_union.push_back(t_id);
            tm.setTriInfo(t_id, 1);
            num_tris_in_final_solution++;
        }
    }

    return num_tris_in_final_solution;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// if more than 2 models -> model 0 - all the others
inline uint boolSubtraction(FastTrimesh &tm, const Labels &labels, Data &data)
{
    uint num_tris_in_final_solution = 0;
    tm.resetTrianglesInfo();

    for(uint t_id = 0; t_id < tm.numTris(); t_id++)
    {
        if(labels.surface[t_id][0] && labels.inside[t_id].count() == 0) // triangle to keep
        {
            data.t_ids_subtraction.push_back(t_id);
            tm.setTriInfo(t_id, 1);
            num_tris_in_final_solution++;
        }
        else if(!labels.surface[t_id][0] && labels.inside[t_id][0] && labels.inside[t_id].count() == 1)
        {
            data.t_ids_subtraction.push_back(t_id);
            tm.setTriInfo(t_id, 1);
            num_tris_in_final_solution++;
        }
    }

    //fix triangles orientation
    for(uint t_id = 0; t_id < tm.numTris(); t_id++)
    {
        if(tm.triInfo(t_id) == 1 && labels.surface[t_id][0] != 1)
            tm.flipTri(t_id);
    }

    return num_tris_in_final_solution;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline uint boolXOR(FastTrimesh &tm, const Labels &labels)
{
    uint num_tris_in_final_solution = 0;
    tm.resetTrianglesInfo();

    for(uint t_id = 0; t_id < tm.numTris(); t_id++)
    {
        if((labels.inside[t_id].count() == 0) || ((labels.surface[t_id] ^ labels.inside[t_id]).count() == labels.num)) // triangle to keep
        {
            tm.setTriInfo(t_id, 1);
            num_tris_in_final_solution++;
        }
    }

    // fix triangles orientation
    for(uint t_id = 0; t_id < tm.numTris(); t_id++)
    {
        if(tm.triInfo(t_id) == 1 && labels.inside[t_id].count() > 0)
            tm.flipTri(t_id);
    }

    return num_tris_in_final_solution;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline uint bitsetToUint(const bitset<NBIT> &b)
{
    assert(b.count() == 1 && "more than 1 bit set to 1");

    for(uint i = 0; i < NBIT; i++)
        if (b[i]) return i;

    return 0; //warning killer
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool consistentWinding(const uint *t0, const uint *t1) // t0 -> vertex ids of triangle t0, t1 -> vertex ids of triangle t1
{
    int j = 0;

    while(t0[0] != t1[j] && j < 3) j++;
    assert(j < 3 && "not same triangle");

    if (t0[1] == t1[(j+1)%3] && t0[2] == t1[(j+2)%3]) return true;
    else return false;
}

/// :::::::::::::: DEBUG :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline std::string printBitset(const std::bitset<NBIT> &b, uint num_label)
{
    std::string s = b.to_string();
    s = s.substr(b.size()-num_label, b.size());
    std::cerr << s << std::endl;

    return s;
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void saveOutputWithLabels(const std::string &filename, cinolib::Trimesh<> &m, const std::vector<std::bitset<NBIT> > &labels)
{
    std::vector<double> coords(3 * m.num_verts());
    for(uint v_id = 0; v_id < m.num_verts(); v_id++)
    {
        coords[3 * v_id] = m.vert(v_id).x();
        coords[3 * v_id +1] = m.vert(v_id).y();
        coords[3 * v_id +2] = m.vert(v_id).z();
    }

    std::vector<int> int_labels(m.num_polys());
    for(uint t_id = 0; t_id < m.num_polys(); t_id++)
    {
        int l = static_cast<int>(labels[t_id].to_ulong());
        int_labels[t_id] = l;
    }

    cinolib::write_OBJ(filename.c_str(), coords, m.vector_polys(), int_labels);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void loadInputWithLabels(const string &filename, std::vector<double> &coords, std::vector<uint> &tris, std::vector< std::bitset<NBIT> > &labels)
{
    cinolib::Trimesh<> m(filename.c_str());
    m.poly_label_wrt_color();

    coords.reserve(3 * m.num_verts());
    tris.reserve(3 * m.num_polys());
    labels.reserve(m.num_polys());

    for(uint v_id = 0; v_id < m.num_polys(); v_id++)
    {
        coords.push_back(m.vert(v_id).x());
        coords.push_back(m.vert(v_id).y());
        coords.push_back(m.vert(v_id).z());
    }

    for(uint t_id = 0; t_id < m.num_polys(); t_id++)
    {
        tris.push_back(m.poly_vert_id(t_id, 0));
        tris.push_back(m.poly_vert_id(t_id, 1));
        tris.push_back(m.poly_vert_id(t_id, 2));

        std::bitset<NBIT> l(static_cast<unsigned long>(m.poly_data(t_id).label));
        labels.push_back(l);
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void loadInputWithLabels(const string &filename, std::vector<double> &coords, std::vector<uint> &tris, std::vector<uint> &labels)
{
    cinolib::Trimesh<> m(filename.c_str());
    m.poly_label_wrt_color();

    coords.reserve(3 * m.num_verts());
    tris.reserve(3 * m.num_polys());
    labels.reserve(m.num_polys());

    for(uint v_id = 0; v_id < m.num_polys(); v_id++)
    {
        coords.push_back(m.vert(v_id).x());
        coords.push_back(m.vert(v_id).y());
        coords.push_back(m.vert(v_id).z());
    }
    for(uint t_id = 0; t_id < m.num_polys(); t_id++)
    {
        tris.push_back(m.poly_vert_id(t_id, 0));
        tris.push_back(m.poly_vert_id(t_id, 1));
        tris.push_back(m.poly_vert_id(t_id, 2));

        labels.push_back(static_cast<uint>(m.poly_data(t_id).label));
    }
}

///:::::::::::::: RATIONALS FUNCTIONS ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void computeInsideOutCustom(const FastTrimesh &tm, const std::vector<phmap::flat_hash_set<uint>> &patches, const cinolib::Octree &octree,
                                   const std::vector<genericPoint *> &in_verts, const std::vector<uint> &in_tris,
                                   const std::vector<std::bitset<NBIT>> &in_labels, const cinolib::vec3d &max_coords, Labels &labels, Data &data)
{
    tbb::spin_mutex mutex;
    if(print_debug){
        std::cout << "::::INFO PATCHES::::::" << std::endl;
        std::cout << "N° of patches: " << patches.size() << std::endl;
        std::cout << "::::::::::::::::::::::" << std::endl;
    }
    std::bitset <NBIT> labels_patch [patches.size()];


    std::vector<bigrational> in_verts_rational;
    in_verts_rational.resize(in_verts.size()*3);
    bool is_rational = false;

    data.t_ids_debug = {4144,4218,4221,4203,4200,4201,4209,4210,4131};

    //tbb::parallel_for((uint)0, (uint)patches.size(), [&](uint p_id)
    for(uint p_id = 0; p_id < patches.size(); ++p_id) //For each patch
    {
        const phmap::flat_hash_set<uint> &patch_tris = patches[p_id];
        const std::bitset<NBIT> &patch_surface_label = labels.surface[*patch_tris.begin()]; // label of the first triangle of the patch
        std::bitset<NBIT> patch_surface_label_tmp = patch_surface_label;

        Ray ray;
        RationalRay rational_ray;

        //findRayEndpoints(tm, patch_tris, max_coords, ray);
        findRayEndpointsCustom(tm, patch_tris, max_coords, ray, rational_ray, in_verts, in_verts_rational, is_rational,
                               true, data);
        phmap::flat_hash_set<uint> tmp_inters;

        //if rational ray is equal to -1, it means that the ray is defined by floating points
        if(print_debug){
            string full_implicit = is_rational  ? "TRUE" : "FALSE";
            std::cout << "The Patch is full implicit? -> " << full_implicit << std::endl;
        }

        if(rational_ray.tv[0] != -1) {//is defined
           if(print_debug) std::cout << "PROCESSING TRIANGLE IN PATCH N° : " << p_id << std::endl;
            std::vector<IntersectionPointRationals> inter_rat;
            //inter_rat.resize(300);

            findIntersectionsAlongRayRationals(tm, patches, octree, in_verts, in_labels, labels, rational_ray, p_id,
                                               tmp_inters, inter_rat, in_verts_rational, in_tris, data);


            ///profiling and prune intersections and sort along the ray
            std::vector<uint> inters_tris_rat;

            pruneIntersectionsAndSortAlongRayRationals(rational_ray, tm, in_verts, in_tris, in_labels, tmp_inters,
                                                       patch_surface_label, inter_rat, inters_tris_rat,labels,patch_surface_label_tmp,
                                                       patches, in_verts_rational, data);


            ///profiling and analyze sorted intersections
            std::bitset<NBIT> patch_inner_label;

            analyzeSortedIntersectionsRationals(rational_ray, tm, in_verts, inter_rat, patch_inner_label, labels, in_labels, in_verts_rational, in_tris, data);


            ///propagate inner labels on patch
            propagateInnerLabelsOnPatch(patch_tris, patch_inner_label, labels);

        }else {

            cinolib::AABB rayAABB(cinolib::vec3d(ray.v0.X(), ray.v0.Y(), ray.v0.Z()),
                                  cinolib::vec3d(ray.v1.X(), ray.v1.Y(), ray.v1.Z()));

            intersects_box(octree, rayAABB, tmp_inters);
            std::vector<uint> sorted_inters;
            pruneIntersectionsAndSortAlongRay(ray, in_verts, in_tris, in_labels, tmp_inters, patch_surface_label,
                                              sorted_inters);

            std::bitset<NBIT> patch_inner_label;
            analyzeSortedIntersections(ray, in_verts, in_tris, in_labels, sorted_inters, patch_inner_label);
            propagateInnerLabelsOnPatch(patch_tris, patch_inner_label, labels);
        }
    }
    //);
}


inline void setExplicitVertex(const FastTrimesh &tm, std::vector<bigrational> &in_verts_rational, uint vertex_id, bigrational &x, bigrational &y, bigrational &z) {
    if (tm.vert(vertex_id)->isExplicit3D()) {
        x = in_verts_rational[vertex_id * 3];
        y = in_verts_rational[vertex_id * 3 + 1];
        z = in_verts_rational[vertex_id * 3 + 2];
    } else {
        tm.vert(vertex_id)->getExactXYZCoordinates(x, y, z);
   }
}


inline void findRayEndpointsCustom(const FastTrimesh &tm, const phmap::flat_hash_set<uint> &patch, const cinolib::vec3d &max_coords,
                                   Ray &ray, RationalRay &rational_ray, const std::vector<genericPoint *> &in_verts,
                                   std::vector<bigrational> &in_verts_rational, bool &is_rational, bool debug, Data &data)
{
    // check for an explicit point (all operations with explicits are faster)
    int v_id = -1;
    if(!debug){
        for(uint t_id : patch){
            const uint tv[3] = {tm.triVertID(t_id, 0), tm.triVertID(t_id, 1), tm.triVertID(t_id, 2)};

            if (tm.vert(tv[0])->isExplicit3D() && tm.vertInfo(tv[0]) == 0)      v_id = static_cast<int>(tv[0]);
                else if (tm.vert(tv[1])->isExplicit3D() && tm.vertInfo(tv[1]) == 0) v_id = static_cast<int>(tv[1]);
                else if (tm.vert(tv[2])->isExplicit3D() && tm.vertInfo(tv[2]) == 0) v_id = static_cast<int>(tv[2]);

                if (v_id != -1)
                {
                    is_rational = false;
                    const explicitPoint3D &v = tm.vert(v_id)->toExplicit3D();
                    ray.v0  = explicitPoint3D(v.X(), v.Y(), v.Z());
                    ray.v1 = explicitPoint3D(max_coords.x(), v.Y(), v.Z());
                    return;
                }
            }
    }

    //take the patch with triangle with t_id_aux and iterate through the triangles
    in_verts_rational.resize(in_verts.size()*3);

    if(!is_rational) {

        cinolib::Profiler profiler;
        if(profiling)
            profiler.push("Rational conversion");

        for (uint i = 0; i < in_verts.size(); i++) {
            bigrational x, y, z;
            in_verts[i]->getExactXYZCoordinates(x, y, z);
            in_verts_rational[i * 3] = x;
            in_verts_rational[i * 3 + 1] = y;
            in_verts_rational[i * 3 + 2] = z;
        }

        if(profiling)
            profiler.pop(true, "Num Verts: " + std::to_string(in_verts.size()));
    }

    //if the patch contains the triangle with t_id
    for(uint t_id : patch){

        const uint tv[3] = {tm.triVertID(t_id, 0), tm.triVertID(t_id, 1), tm.triVertID(t_id, 2)};

        std::vector<bigrational> x_rat(3), y_rat(3), z_rat(3);

        for (int i = 0; i < 3; ++i) {
            setExplicitVertex(tm, in_verts_rational, tv[i], x_rat[i], y_rat[i], z_rat[i]);
        }

        std::array<std::array<bigrational, 3>, 3> tv_rat;

        tv_rat[0] = {x_rat[0], y_rat[0], z_rat[0]};  // Prima riga
        tv_rat[1] = {x_rat[1], y_rat[1], z_rat[1]};  // Seconda riga
        tv_rat[2] = {x_rat[2], y_rat[2], z_rat[2]};


        int dir = maxComponentInTriangleNormalRationals(x_rat[0], y_rat[0], z_rat[0], x_rat[1], y_rat[1], z_rat[1], x_rat[2], y_rat[2], z_rat[2]);

        bigrational centroid_x = (x_rat[0] + x_rat[1] + x_rat[2]) / bigrational(3);
        bigrational centroid_y = (y_rat[0] + y_rat[1] + y_rat[2]) / bigrational(3);
        bigrational centroid_z = (z_rat[0] + z_rat[1] + z_rat[2]) / bigrational(3);

        rational_ray.v0 = {centroid_x, centroid_y, centroid_z};

        if (dir == 0) { //direction x
            rational_ray.v1 = {bigrational(max_coords.x()), centroid_y, centroid_z};
            rational_ray.dir = 'X';
        } else if (dir == 1) { //direction y
            rational_ray.v1 = {centroid_x, bigrational(max_coords.y()), centroid_z};
            rational_ray.dir = 'Y';
        } else if (dir == 2) { //direction z
            rational_ray.v1 = {centroid_x, centroid_y, bigrational(max_coords.z())};
            rational_ray.dir = 'Z';
        }else{
            std::cout<<"Error in direction"<<std::endl;
            std::exit(EXIT_FAILURE);
        }
        //check if the centroid is inside the triangle
        //if (cinolib::orient3d(&tv_rat[0][0], &tv_rat[1][0], &tv_rat[2][0], &rational_ray.v0[0]).sgn() != 0) {
        if (cinolib::orient3d(&tv_rat[0][0], &tv_rat[1][0], &tv_rat[2][0], &rational_ray.v0[0]) != bigrational()) {
            std::cout << "The centroid is not on the triangle or the orient3d doesn't work properly" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        bigrational e0_rat = cinolib::orient3d(&tv_rat[0][0], &tv_rat[1][0], &rational_ray.v1[0], &rational_ray.v0[0]);
        bigrational e1_rat = cinolib::orient3d(&tv_rat[1][0], &tv_rat[2][0], &rational_ray.v1[0], &rational_ray.v0[0]);
        bigrational e2_rat = cinolib::orient3d(&tv_rat[2][0], &tv_rat[0][0], &rational_ray.v1[0], &rational_ray.v0[0]);

          if((e0_rat > bigrational() && e1_rat > bigrational() && e2_rat > bigrational()) ||
             (e0_rat < bigrational() && e1_rat < bigrational() && e2_rat < bigrational())){

              if(print_debug){
                  std::cout << "Triangle that create the ray: " << t_id <<  " Direction: " << rational_ray.dir << std::endl;
                  //print the coords of the ray
                  std::cout << "Ray v0: " << rational_ray.v0[0] << " " << rational_ray.v0[1] << " " << rational_ray.v0[2] << std::endl;
                  std::cout << "Ray v1: " << rational_ray.v1[0] << " " << rational_ray.v1[1] << " " << rational_ray.v1[2] << std::endl;
                }

            rational_ray.tv[0] = static_cast<int>(tm.triVertID(t_id, 0));
            rational_ray.tv[1] = static_cast<int>(tm.triVertID(t_id, 1));
            rational_ray.tv[2] = static_cast<int>(tm.triVertID(t_id, 2));
            is_rational = true;
            return;
        } else{
            std::exit(EXIT_FAILURE);
        }
    }
    std::cout << "Something went wrong during the creation of the ray with a patch fully composed by implicit points" << std::endl;
    std::exit(EXIT_FAILURE);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline BoundingBox calculateBoundingBox(const std::array<bigrational, 3>& tv0,
                                 const std::array<bigrational, 3>& tv1,
                                 const std::array<bigrational, 3>& tv2)
{
    BoundingBox box;
    box.xmin = std::min({tv0[0], tv1[0], tv2[0]});
    box.xmax = std::max({tv0[0], tv1[0], tv2[0]});
    box.ymin = std::min({tv0[1], tv1[1], tv2[1]});
    box.ymax = std::max({tv0[1], tv1[1], tv2[1]});
    box.zmin = std::min({tv0[2], tv1[2], tv2[2]});
    box.zmax = std::max({tv0[2], tv1[2], tv2[2]});
    return box;
}

inline bool isNormalCorrect(
         bigrational& ov1x,  bigrational& ov1y,  bigrational& ov1z,
         bigrational& ov2x,  bigrational& ov2y,  bigrational& ov2z,
         bigrational& ov3x,  bigrational& ov3y,  bigrational& ov3z,
         bigrational& px,  bigrational& py,  bigrational& pz)
{
    // Calcola la normale al triangolo
    bigrational v3x = ov3x - ov2x;
    bigrational v3y = ov3y - ov2y;
    bigrational v3z = ov3z - ov2z;
    bigrational v2x = ov2x - ov1x;
    bigrational v2y = ov2y - ov1y;
    bigrational v2z = ov2z - ov1z;

    bigrational nvx = v2y * v3z - v2z * v3y;
    bigrational nvy = v3x * v2z - v3z * v2x;
    bigrational nvz = v2x * v3y - v2y * v3x;

    // Punto di riferimento P (può essere il baricentro del triangolo o un punto noto)
    bigrational cx = (ov1x + ov2x + ov3x) / bigrational(3);
    bigrational cy = (ov1y + ov2y + ov3y) / bigrational(3);
    bigrational cz = (ov1z + ov2z + ov3z) / bigrational(3);

    // Vettore dal centroide al punto P
    bigrational px_c = px - cx;
    bigrational py_c = py - cy;
    bigrational pz_c = pz - cz;

    // Prodotto scalare tra la normale e il vettore (centroide -> P)
    bigrational dot_product = nvx * px_c + nvy * py_c + nvz * pz_c;

    // Se il prodotto è positivo, la normale punta verso il punto P (corretta)
    // Se è negativo, la normale è invertita
    return dot_product > bigrational();
}

inline bool rayIntersectAABB(const RationalRay &ray, const BoundingBox &aabb) {

    //Dir X
    if (ray.dir == 'X') {
        //Check if the coordinates Y and Z of the ray are inside the bounding box
        bool yInside = (ray.v0[1] >= aabb.ymin && ray.v0[1] <= aabb.ymax) ||
                       (ray.v1[1] >= aabb.ymin && ray.v1[1] <= aabb.ymax);

        bool zInside = (ray.v0[2] >= aabb.zmin && ray.v0[2] <= aabb.zmax) ||
                       (ray.v1[2] >= aabb.zmin && ray.v1[2] <= aabb.zmax);

        return yInside && zInside;
    }

    //Dir Y
    if (ray.dir == 'Y') {
        //Check if the coordinates X and Z of the ray are inside the bounding box
        bool xInside = (ray.v0[0] >= aabb.xmin && ray.v0[0] <= aabb.xmax) ||
                       (ray.v1[0] >= aabb.xmin && ray.v1[0] <= aabb.xmax);
        bool zInside = (ray.v0[2] >= aabb.zmin && ray.v0[2] <= aabb.zmax) ||
                       (ray.v1[2] >= aabb.zmin && ray.v1[2] <= aabb.zmax);

        // Se entrambe le coordinate X e Z sono dentro, c'è intersezione
        return xInside && zInside;
    }

    //Dir Z
    if (ray.dir == 'Z') {
        //Check if the coordinates X and Y of the ray are inside the bounding box
        bool xInside = (ray.v0[0] >= aabb.xmin && ray.v0[0] <= aabb.xmax) ||
                       (ray.v1[0] >= aabb.xmin && ray.v1[0] <= aabb.xmax);
        bool yInside = (ray.v0[1] >= aabb.ymin && ray.v0[1] <= aabb.ymax) ||
                       (ray.v1[1] >= aabb.ymin && ray.v1[1] <= aabb.ymax);

        return xInside && yInside;
    }

    //The direction is not valid
    return false;
}
inline void findIntersectionsAlongRayRationals(const FastTrimesh &tm,
                                               const std::vector<phmap::flat_hash_set<uint>> &patches,
                                               const cinolib::Octree& tree,
                                               const std::vector<genericPoint *> &in_verts,
                                               const std::vector<std::bitset<NBIT>> &in_labels,
                                               Labels &labels,
                                               const RationalRay &rational_ray,
                                               uint curr_p_id,
                                               phmap::flat_hash_set<uint> &tmp_inters,
                                               std::vector<IntersectionPointRationals> &inter_rat,
                                               std::vector<bigrational> &in_verts_rational,
                                               const std::vector<uint> &in_tris,
                                               Data &data)
{

    tmp_inters.clear();

    for (uint t_id = 0; t_id < in_tris.size() / 3; ++t_id) {

        const uint& id_v0 = in_tris[3 * t_id];
        const uint& id_v1 = in_tris[3 * t_id + 1];
        const uint& id_v2 = in_tris[3 * t_id + 2];

        const bigrational& x0 = in_verts_rational[3 * id_v0];
        const bigrational& y0 = in_verts_rational[3 * id_v0 + 1];
        const bigrational& z0 = in_verts_rational[3 * id_v0 + 2];

        const bigrational& x1 = in_verts_rational[3 * id_v1];
        const bigrational& y1 = in_verts_rational[3 * id_v1 + 1];
        const bigrational& z1 = in_verts_rational[3 * id_v1 + 2];

        const bigrational& x2 = in_verts_rational[3 * id_v2];
        const bigrational& y2 = in_verts_rational[3 * id_v2 + 1];
        const bigrational& z2 = in_verts_rational[3 * id_v2 + 2];

        std::array<bigrational, 3> tv0 = {x0, y0, z0};
        std::array<bigrational, 3> tv1 = {x1, y1, z1};
        std::array<bigrational, 3> tv2 = {x2, y2, z2};

        // Calcola la bounding box del triangolo
        BoundingBox box = calculateBoundingBox(tv0, tv1, tv2);

        cinolib::Profiler p;
        //p.push("::: Time of one test ray triangle intersection --> ");

         if (rayIntersectAABB(rational_ray, box)) {
             int intersection = segment_triangle_intersect_3d(&rational_ray.v0[0], &rational_ray.v1[0], &tv0[0], &tv1[0], &tv2[0]);
             if (intersection) {

                if(print_debug){
                    string type = intersection == 1 ? "Simplicial Complex" : intersection == 2 ? "Intersect" : "Overlap";
                    std::cout << "t_id of the triangle that is intersected by the ray: " << t_id <<  " type: " << type << std::endl;
                    data.t_ids_inters_ray.push_back(t_id);
                }

                IntersectionPointRationals p_int;
                p_int.setTriId(t_id);

                plane_line_intersection(&tv0[0] ,&tv1[0], &tv2[0], &rational_ray.v0[0], &rational_ray.v1[0], &p_int[0]);

                tmp_inters.insert(t_id);

                inter_rat.emplace_back(p_int[0], p_int[1], p_int[2], t_id);
                inter_rat.emplace_back(p_int);

                 if(print_debug){
                    p_int.printIntersectionPoint();
                 }

            }
        }
       // p.pop();
    }

    if(print_debug){
        std::cout << "\n:::: Triangles that are intersected by the ray: \n";
        for (const uint t_id : tmp_inters) {
            std::cout << t_id << std::endl;
        }
        std::cout << ":::::::::::::::::::::::::::::::::::::::::::::::::\n";
    }
}

bool isIntersectionValid(const std::vector<bigrational>& inter, const RationalRay& rational_ray) {

    int coord_index = (rational_ray.dir == 'X') ? 0 : (rational_ray.dir == 'Y') ? 1 : 2;

    const bigrational& ray_origin_coord = rational_ray.v0[coord_index];
    const bigrational& ray_end_coord = rational_ray.v1[coord_index];

    bool positive_direction = ray_end_coord > ray_origin_coord;

    const bigrational& inter_coord = (coord_index == 0) ? inter.at(0) : (coord_index == 1) ? inter.at(1) : inter.at(2);
    return positive_direction ? (inter_coord >= ray_origin_coord) : (inter_coord <= ray_origin_coord);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline IntersInfo fast2DCheckIntersectionOnRayRationals(const RationalRay &ray, const std::vector<bigrational> &tv0, const std::vector<bigrational> &tv1, const std::vector<bigrational> &tv2)
{
    bigrational v0_rat[2], v1_rat[2], v2_rat[2], vq_rat[2];

    bigrational ray_v1_x, ray_v1_y, ray_v1_z;

    ray_v1_x = ray.v1[0];
    ray_v1_y = ray.v1[1];
    ray_v1_z = ray.v1[2];

    switch (ray.dir)
    {
        case 'X': // only YZ coordinates
            v0_rat[0] = tv0[1]; v0_rat[1] = tv0[2];
            v1_rat[0] = tv1[1]; v1_rat[1] = tv1[2];
            v2_rat[0] = tv2[1]; v2_rat[1] = tv2[2];
            vq_rat[0] = ray_v1_y; vq_rat[1] = ray_v1_z;
            break;

        case 'Y': // only XZ coordinates
            v0_rat[0] = tv0[0]; v0_rat[1] = tv0[2];
            v1_rat[0] = tv1[0]; v1_rat[1] = tv1[2];
            v2_rat[0] = tv2[0]; v2_rat[1] = tv2[2];
            vq_rat[0] = ray_v1_x; vq_rat[1] = ray_v1_z;
            break;

        case 'Z': // only XY coordinates
            v0_rat[0] = tv0[0]; v0_rat[1] = tv0[1];
            v1_rat[0] = tv1[0]; v1_rat[1] = tv1[1];
            v2_rat[0] = tv2[0]; v2_rat[1] = tv2[1];
            vq_rat[0] = ray_v1_x; vq_rat[1] = ray_v1_y;
            break;
    }

    bigrational or01_rat = cinolib::orient2d(&v0_rat[0], &v1_rat[0], &vq_rat[0]);
    bigrational or12_rat = cinolib::orient2d(&v1_rat[0], &v2_rat[0], &vq_rat[0]);
    bigrational or20_rat = cinolib::orient2d(&v2_rat[0], &v0_rat[0], &vq_rat[0]);


    //int or01_rat = orient2dT(&v0_rat[0], &v1_rat[0], &vq_rat[0]);
    //int or12_rat = orient2dT(&v1_rat[0], &v2_rat[0], &vq_rat[0]);
    //int or20_rat = orient2dT(&v2_rat[0], &v0_rat[0], &vq_rat[0]);

    //if (or01_rat.sgn() == 0 && or12_rat.sgn() == 0 && or20_rat.sgn() == 0)
    if (or01_rat == bigrational() && or12_rat == bigrational() && or20_rat == bigrational())
        //if (or01_rat == 0 && or12_rat == 0 && or20_rat == 0)
    {
        switch (ray.dir)
        {
            case 'X':
                if (tv0[0] == ray_v1_x && tv1[0] == ray_v1_x && tv2[0] == ray_v1_x)
                    return DISCARD;
                break;
            case 'Y':
                if (tv0[1] == ray_v1_y && tv1[1] == ray_v1_y && tv2[1] == ray_v1_y)
                    return DISCARD;
                break;
            case 'Z':
                if (tv0[2] == ray_v1_z && tv1[2] == ray_v1_z && tv2[2] == ray_v1_z)
                    return DISCARD;
                break;
        }
    }

    if (((or01_rat < bigrational() || or01_rat.sgn() == 0) && (or12_rat < bigrational() || or12_rat.sgn() == 0) && (or20_rat < bigrational() || or20_rat.sgn() == 0)) ||
       ((or01_rat > bigrational() || or01_rat.sgn() == 0) && (or12_rat > bigrational() || or12_rat.sgn() == 0) && (or20_rat > bigrational() || or20_rat.sgn() == 0)))
      //if (((or01_rat <= 0) && (or12_rat <= 0) && (or20_rat <= 0)) ||
        //((or01_rat >= 0) && (or12_rat >= 0) && (or20_rat >= 0)))
    {
        if (v0_rat[0] == vq_rat[0] && v0_rat[1] == vq_rat[1]) return INT_IN_V0;
        if (v1_rat[0] == vq_rat[0] && v1_rat[1] == vq_rat[1]) return INT_IN_V1;
        if (v2_rat[0] == vq_rat[0] && v2_rat[1] == vq_rat[1]) return INT_IN_V2;

        if (or01_rat == bigrational()) return INT_IN_EDGE01;
        if (or12_rat == bigrational()) return INT_IN_EDGE12;
        if (or20_rat == bigrational()) return INT_IN_EDGE20;

        //if (or01_rat == 0) return INT_IN_EDGE01;
        //if (or12_rat == 0) return INT_IN_EDGE12;
        //if (or20_rat == 0) return INT_IN_EDGE20;


        return INT_IN_TRI;
    }

    return NO_INT;
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline bool checkIntersectionInsideTriangle3DRationals(const RationalRay &ray, const std::array<bigrational,3> &tv0, const std::array<bigrational,3> &tv1, const std::array<bigrational,3> &tv2)
{
    bigrational or01f_rat = cinolib::orient3d(&tv0[0], &tv1[0], &ray.v0[0], &ray.v1[0]);
    bigrational or12f_rat = cinolib::orient3d(&tv1[0], &tv2[0], &ray.v0[0], &ray.v1[0]);
    bigrational or20f_rat = cinolib::orient3d(&tv2[0], &tv0[0], &ray.v0[0], &ray.v1[0]);
    //int or01f_rat = orient3dT(&tv0[0], &tv1[0], &ray.v0[0], &ray.v1[0]);
    //int or12f_rat = orient3dT(&tv1[0], &tv2[0], &ray.v0[0], &ray.v1[0]);
    //int or20f_rat = orient3dT(&tv2[0], &tv0[0], &ray.v0[0], &ray.v1[0]);

    if(or01f_rat > bigrational() && or12f_rat > bigrational() && or20f_rat > bigrational()) return true;
    if(or01f_rat < bigrational() && or12f_rat < bigrational() && or20f_rat < bigrational()) return true;

    //if(or01f_rat > 0 && or12f_rat > 0 && or20f_rat > 0) return true;
    //if(or01f_rat < 0 && or12f_rat < 0 && or20f_rat < 0) return true;

    return false;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// return 1 if inside, 0 if outside
inline uint checkTriangleOrientationRationals(const RationalRay &ray, const std::vector<bigrational> &tv0, const std::vector<bigrational> &tv1, const std::vector<bigrational> &tv2)
{

    bigrational res = cinolib::orient3d(&tv0[0],&tv1[0],&tv2[0], &ray.v1[0]);
    //int res = orient3dT(&tv0[0],&tv1[0],&tv2[0], &ray.v1[0]);

    assert(res != bigrational() && "Problem in PointOrientation(...)");

    /* in res we have sign(area(v0, v1, v2, ray.second))
     * if the area is >0 the ray is doing INSIDE -> OUTSIDE, so the patch is INSIDE
     * else the ray is doing OUTSIDE -> INSIDE so the patch is OUTSIDE */

    // Determina se il raggio è crescente o decrescente
    bool increasingOrder = (ray.dir == 'X') ? (ray.v0[0] < ray.v1[0]) :
                           (ray.dir == 'Y') ? (ray.v0[1] < ray.v1[1]) :
                           (ray.v0[2] < ray.v1[2]);

    return (increasingOrder ? (res < bigrational()) : (res > bigrational())) ? 1 : 0;
    //return (increasingOrder ? (res < 0) : (res > 0)) ? 1 : 0;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void pruneIntersectionsAndSortAlongRayRationals(const RationalRay &ray, const FastTrimesh &tm, const std::vector<genericPoint*> &in_verts,
                                                       const std::vector<uint> &in_tris, const std::vector<std::bitset<NBIT>> &in_labels,
                                                       const phmap::flat_hash_set<uint> &tmp_inters, const std::bitset<NBIT> &patch_surface_label,
                                                       std::vector<IntersectionPointRationals> &inter_rat, std::vector<uint> &inters_tris_rat, Labels &labels, std::bitset<NBIT> &patch_surface_label_tmp,
                                                       const std::vector<phmap::flat_hash_set<uint>> &patches, std::vector<bigrational> &in_verts_rational, Data &data){

    phmap::flat_hash_set<uint> visited_tri;
    visited_tri.reserve(inter_rat.size() / 6);
    std::pair<phmap::flat_hash_set<uint>::iterator, bool> ins;

    //print tmp_inters
    if(print_debug && !tmp_inters.empty()) {

        std::cout << ":::: Triangles that are intersecated by the ray: " << std::endl;
        for (uint t_id: tmp_inters) {
            std::cout << t_id << std::endl;
        }
        std::cout << ":::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
        std::cout << "inter rat size : " << inter_rat.size() << std::endl;
    }
    //print inter_rat

    if(print_debug && !inter_rat.empty()) {
        std::cout << ":::: Intersections points: " << std::endl;
        for (IntersectionPointRationals ipr: inter_rat) {
            std::cout << "Intersection point of t_id: " << ipr.getTriId() << " with coords: " << ipr.getX().get_d() << " "
                      << ipr.getY().get_d() << " " << ipr.getZ().get_d() << std::endl;

            std::cout << ":::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
        }
    }
    std::vector<IntersectionPointRationals> inter_rat_tmp;
    inter_rat_tmp.reserve(inter_rat.size());

    for (uint t_id_int: tmp_inters)
    {
        ins = visited_tri.insert(t_id_int);
        if (!ins.second){
            continue; }// triangle already analyzed or in the one ring of a vert or in the adj of an edge

        const std::bitset<NBIT> tested_tri_label = in_labels.at(t_id_int);
        uint uint_tri_label = bitsetToUint(tested_tri_label);

        if (patch_surface_label_tmp[uint_tri_label]) continue; // <-- triangle of the same label of the tested patch

        bigrational tv0_x, tv0_y, tv0_z,
                    tv1_x, tv1_y, tv1_z,
                    tv2_x, tv2_y, tv2_z;

        const int t_base = 3 * t_id_int;

        const int i0 = 3 * in_tris[t_base];
        const int i1 = 3 * in_tris[t_base + 1];
        const int i2 = 3 * in_tris[t_base + 2];

        tv0_x = in_verts_rational[i0];
        tv0_y = in_verts_rational[i0 + 1];
        tv0_z = in_verts_rational[i0 + 2];

        tv1_x = in_verts_rational[i1];
        tv1_y = in_verts_rational[i1 + 1];
        tv1_z = in_verts_rational[i1 + 2];

        tv2_x = in_verts_rational[i2];
        tv2_y = in_verts_rational[i2 + 1];
        tv2_z = in_verts_rational[i2 + 2];

        const std::vector<bigrational> tv0_exact = {tv0_x, tv0_y, tv0_z};
        const std::vector<bigrational> tv1_exact = {tv1_x, tv1_y, tv1_z};
        const std::vector<bigrational> tv2_exact = {tv2_x, tv2_y, tv2_z};

        IntersInfo ii = fast2DCheckIntersectionOnRayRationals(ray, tv0_exact, tv1_exact, tv2_exact);

        //std::cout << std::endl;
        if (ii == DISCARD){
            if(print_debug) std::cout<<"DISCARD the intersection detected on tri: "<<t_id_int << std::endl;
            continue;
        }
        if (ii == NO_INT){
            if(print_debug) std::cout<<"NO Intersection detected on tri: " << t_id_int <<std::endl;
            continue;
        }

        if (ii == INT_IN_TRI) {
            if(print_debug || data.t_id_debug == 4218) std::cout << "Intersection detected IN TRIANGLE with t_id: " << t_id_int << std::endl;
            inters_tris_rat.push_back(t_id_int);

            bool copy_found = copyIntersectionPoint(inter_rat, inter_rat_tmp, t_id_int);

            if(!copy_found){
             IntersectionPointRationals p_int;
             plane_line_intersection(&tv0_exact[0], &tv1_exact[0], &tv2_exact[0], &ray.v0[0], &ray.v1[0], &p_int[0]);
             inter_rat_tmp.emplace_back(p_int);
            }

        } else if (ii == INT_IN_V0 || ii == INT_IN_V1 || ii == INT_IN_V2) {

            uint v_id;
            if (ii == INT_IN_V0){
                if(print_debug || data.t_id_debug == 4218) std::cout << "Intersection detected IN VERTEX v0 on triangle with t_id: " << t_id_int  << " on vertex with id: " <<  tm.triVertID(t_id_int,0 ) << std::endl;
                v_id = in_tris[3 * t_id_int];}
            else if (ii == INT_IN_V1){
                if(print_debug || data.t_id_debug == 4218) std::cout << "Intersection detected IN VERTEX v1 on triangle with t_id: "<< t_id_int << " on vertex with id: "<< tm.triVertID(t_id_int,1)<< std::endl;
                v_id = in_tris[3 * t_id_int + 1];}
            else{
                if(print_debug || data.t_id_debug == 4218) std::cout << "Intersection detected IN VERTEX v2 on triangle with t_id: " << t_id_int << " on vertex with id: " << tm.triVertID(t_id_int,2) << std::endl;
                v_id = in_tris[3 * t_id_int + 2];}

            std::vector<uint> vert_one_ring;
            findVertRingTris(v_id, tested_tri_label, tmp_inters, in_tris, in_labels, vert_one_ring);

            for (uint t: vert_one_ring)
                visited_tri.insert(t); // mark all the one ring as visited

            int winner_tri = -1;// the first inters triangle after ray perturbation
            winner_tri = perturbRayAndFindIntersTriRationals(ray, in_verts, in_tris,
                                                    vert_one_ring, in_verts_rational); // the first inters triangle after ray perturbation

            if (winner_tri != -1){
                if(print_debug) std::cout << "Winner triangle: " << winner_tri << std::endl;

                inters_tris_rat.push_back(winner_tri);

                bool copy_found = copyIntersectionPoint(inter_rat, inter_rat_tmp, winner_tri);
                if(!copy_found){
                    IntersectionPointRationals p_int;
                    plane_line_intersection(&tv0_exact[0], &tv1_exact[0], &tv2_exact[0], &ray.v0[0], &ray.v1[0], &p_int[0]);
                    inter_rat_tmp.emplace_back(p_int);
                }
            }

        } else if (ii == INT_IN_EDGE01 || ii == INT_IN_EDGE12 || ii == INT_IN_EDGE20) {

            int ev0_id, ev1_id;
            if (ii == INT_IN_EDGE01) {
                if(print_debug) std::cout << "Intersection detected IN EDGE e01 on triangle with t_id: " << t_id_int << std::endl;
                ev0_id = in_tris[3 * t_id_int];
                ev1_id = in_tris[3 * t_id_int + 1];
            } else if (ii == INT_IN_EDGE12) {
                if(print_debug) std::cout << "Intersection detected IN EDGE e12 on triangle with t_id: " << t_id_int << std::endl;
                ev0_id = in_tris[3 * t_id_int + 1];
                ev1_id = in_tris[3 * t_id_int + 2];
            } else {
                if(print_debug) std::cout << "Intersection detected IN EDGE e20 on triangle with t_id: " << t_id_int << std::endl;
                ev0_id = in_tris[3 * t_id_int + 2];
                ev1_id = in_tris[3 * t_id_int];
            }

            std::vector<uint> edge_tris;
            findEdgeTris(ev0_id, ev1_id, tested_tri_label, tmp_inters, in_tris, in_labels, edge_tris);

            for (uint t: edge_tris)
                visited_tri.insert(t); // mark all the one ring as visited

            int winner_tri = -1;
            winner_tri = perturbRayAndFindIntersTriRationals(ray, in_verts, in_tris, edge_tris, in_verts_rational);

            if(print_debug) std::cout << "Winner triangle: " << winner_tri << std::endl;

            if (winner_tri != -1){
                inters_tris_rat.push_back(winner_tri);

                bool copy_found = copyIntersectionPoint(inter_rat, inter_rat_tmp, winner_tri);
                if(!copy_found){
                    IntersectionPointRationals p_int;
                    plane_line_intersection(&tv0_exact[0], &tv1_exact[0], &tv2_exact[0], &ray.v0[0], &ray.v1[0], &p_int[0]);
                    inter_rat_tmp.emplace_back(p_int);
                }
            }
        }
        //std::cout << std::endl;
    }
    if((print_debug) && !inter_rat.empty()){
        std::cout << ":::: Triangles that are intersecated AFTER PRUNING by the ray: " << std::endl;
        for (uint t_id: inters_tris_rat){
            std::cout << t_id << std::endl;
        }
        std::cout << ":::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
    }

    if((print_debug) && !inter_rat.empty()){
        std::cout << ":::: Intersections points BEFORE COPY AFTER PRUNING: " << std::endl;
        for (IntersectionPointRationals ipr: inter_rat){
            std::cout << "Intersection point of t_id: " <<ipr.getTriId() << " with coords: " << ipr.getX().get_d() << " " << ipr.getY().get_d() << " " << ipr.getZ().get_d() << std::endl;
        }
        std::cout << ":::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
    }
    inter_rat = inter_rat_tmp;

    //print inter_rat
    if((print_debug) && !inter_rat.empty()){
        std::cout << ":::: Intersections points AFTER COPY AFTER PRUNING: " << std::endl;
        for (IntersectionPointRationals ipr: inter_rat){
            std::cout << "Intersection point of t_id: " <<ipr.getTriId() << " with coords: " << ipr.getX().get_d() << " " << ipr.getY().get_d() << " " << ipr.getZ().get_d() << std::endl;
        }
        std::cout << ":::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
    }
    bool increasingOrder = (ray.dir == 'X') ? (ray.v0[0] < ray.v1[0]) :
                           (ray.dir == 'Y') ? (ray.v0[1] < ray.v1[1]) :
                           (ray.v0[2] < ray.v1[2]);

    auto sortComparator = [&ray, increasingOrder](const IntersectionPointRationals &a, const IntersectionPointRationals &b) {
        if (ray.dir == 'X') return increasingOrder ? a.lessThanX(b) : b.lessThanX(a);
        if (ray.dir == 'Y') return increasingOrder ? a.lessThanY(b) : b.lessThanY(a);
        return increasingOrder ? a.lessThanZ(b) : b.lessThanZ(a);
    };

    std::sort(inter_rat.begin(), inter_rat.end(), sortComparator);

    if((print_debug) && !inter_rat.empty()){
        std::cout << ":::: Intersections points AFTER PRUNING AND SORTING: " << std::endl;
        for (IntersectionPointRationals ipr: inter_rat){
            std::cout << "Intersection point of t_id: " <<ipr.getTriId() << " with coords: " << ipr.getX() << " " << ipr.getY() << " " << ipr.getZ() << std::endl;
        }
        std::cout << ":::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void analyzeSortedIntersectionsRationals(const RationalRay &rational_ray, const FastTrimesh &tm, const std::vector<genericPoint*> &in_verts,
                                                std::vector<IntersectionPointRationals> &inter_rat, std::bitset<NBIT> &patch_inner_label, Labels &labels,
                                                const std::vector<std::bitset<NBIT>> &in_labels, std::vector<bigrational> &in_verts_rational,
                                                const std::vector<uint> &in_tris, Data &data){


    //print intersection point ids
    if(print_debug) {
        std::cout << ":::: Intersections points INTO ANALYZING: " << std::endl;
        for (IntersectionPointRationals ipr: inter_rat){
            std::cout << "Intersection point of t_id: " <<ipr.getTriId() << " with coords: " << ipr.getX().get_d() << " " << ipr.getY().get_d() << " " << ipr.getZ().get_d() << std::endl;
        }
        std::cout << ":::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
    }
    std::bitset<NBIT> visited_labels;
    for (const IntersectionPointRationals &intersection : inter_rat) {
        uint t_id = intersection.getTriId();
        if(print_debug){
            std::cout << t_id << " analysizing..." << std::endl;
        }
        uint t_label = bitsetToUint(in_labels[t_id]);
        if (visited_labels[t_label]){
            if(print_debug){
                std::cout << "t_label: " << t_label << " already analyzed" << std::endl;
                std::cout << "visited_labels: " << visited_labels[t_label] << std::endl;
                std::cout << "Triangle with t_id: " << t_id << " already analyzed" << std::endl;}
            continue;}

        bigrational x0, x1, x2, y0, y1, y2, z0, z1, z2;

        const uint i0 = in_tris[3 * t_id];
        const uint i1 = in_tris[3 * t_id + 1];
        const uint i2 = in_tris[3 * t_id + 2];

        x0 = in_verts_rational[3 * i0];
        y0 = in_verts_rational[3 * i0 + 1];
        z0 = in_verts_rational[3 * i0 + 2];

        x1 = in_verts_rational[3 * i1];
        y1 = in_verts_rational[3 * i1 + 1];
        z1 = in_verts_rational[3 * i1 + 2];

        x2 = in_verts_rational[3 * i2];
        y2 = in_verts_rational[3 * i2 + 1];
        z2 = in_verts_rational[3 * i2 + 2];

        std::vector<bigrational> tv0 = {x0, y0, z0};
        std::vector<bigrational> tv1 = {x1, y1, z1};
        std::vector<bigrational> tv2 = {x2, y2, z2};

       bool check = checkTriangleOrientationRationals(rational_ray, tv0, tv1, tv2);

       if(print_debug) {
           string check_str = (check == 1) ? "INSIDE" : "OUTSIDE";
           std::cout << "Analysis of intersection on t_id: " << t_id << " check: " << check_str << std::endl;
       }

       if (check == 1) {// checkOrientation -> 1 if inside, 0 if outside
            patch_inner_label[t_label] = true;
        }
        visited_labels[t_label] = true;
    }
}
/*
inline bigrational next_after(const bigrational& x, const bigrational& target) {
    // If x equals target, return x
    if (x.get_num() == target.get_num() && x.get_den() == target.get_den() && x.sgn() == target.sgn()) {
        return x;
    }

    bigrational result = x; // Start with x

    // Compare x and target
    bool moving_up = false;
    if (x.sgn() < target.sgn()) {  // If signs are different
        moving_up = (target.sgn() > 0);
    } else if (x.sgn() > target.sgn()) {
        moving_up = false;
    } else {  // If same signs, compare magnitudes
        bignatural left = x.get_num() * target.get_den();
        bignatural right = target.get_num() * x.get_den();

        moving_up = (x.sgn() > 0) ? (left < right) : (left > right);
    }

    // Compute the next fraction
    if (moving_up) {
        bigrational tmp(result.get_num() * 2 + result.get_den(), result.get_den() * 2, result.sgn());
        // Increase fraction slightly
        result = tmp;
    } else {
        // Decrease fraction slightly
        if (result.get_num() > result.get_den()) {
            bigrational tmp(result.get_num() * 2 - result.get_den(), result.get_den(), result.sgn());
            result = tmp;
        } else {
            bigrational tmp(result.get_den() - (result.get_den() * 2 - result.get_num()), result.get_den(), result.sgn());
            result = tmp;
        }
        bigrational tmp_den(result.get_num(), result.get_den() * 2, result.sgn());
        result = tmp_den;
    }

    return result;
}*/

// Function to calculate a small perturbation (epsilon) based on a given rational value
inline bigrational getEpsilon(const bigrational& value) {
    // Return epsilon as a small rational number based on the denominator of the input value
    // Epsilon is calculated as 1 / (denominator of the value + 1), with the same sign as the original value
    /*return bigrational(static_cast<uint32_t>(1),
                        value.get_den() + static_cast<uint32_t>(1),
                        value.sgn()); */
    return bigrational(static_cast<uint64_t>(1), std::numeric_limits<uint64_t>::max(), value.sgn());
}

// Function to perturb a ray along the X direction with small displacements
inline RationalRay perturbXRayRationals(const RationalRay &ray, uint offset) {
    RationalRay new_ray = ray;
    bigrational epsilon = getEpsilon(ray.v1[1]);  // Epsilon based on the y value

    switch (offset) {
        case 0:  // Slight perturbation in the Y direction (positive epsilon)
            new_ray.v1[1] = ray.v1[1] + epsilon;
            break;
        case 1:  // Perturbation in both Y and Z directions (positive epsilon)
            new_ray.v1[1] = ray.v1[1] + epsilon;
            new_ray.v1[2] = ray.v1[2] + epsilon;
            break;
        case 2:  // Perturbation in the Z direction (positive epsilon)
            new_ray.v1[2] = ray.v1[2] + epsilon;
            break;
        case 3:  // Perturbation in the Y direction (negative epsilon) and Z direction (positive epsilon)
            new_ray.v1[1] = ray.v1[1] - epsilon;
            new_ray.v1[2] = ray.v1[2] + epsilon;
            break;
        case 4:  // Perturbation in the Y direction (negative epsilon)
            new_ray.v1[1] = ray.v1[1] - epsilon;
            break;
        case 5:  // Perturbation in both Y and Z directions (negative epsilon)
            new_ray.v1[1] = ray.v1[1] - epsilon;
            new_ray.v1[2] = ray.v1[2] - epsilon;
            break;
        case 6:  // Perturbation in the Z direction (negative epsilon)
            new_ray.v1[2] = ray.v1[2] - epsilon;
            break;
        case 7:  // Perturbation in the Y direction (positive epsilon) and Z direction (negative epsilon)
            new_ray.v1[1] = ray.v1[1] + epsilon;
            new_ray.v1[2] = ray.v1[2] - epsilon;
            break;
        default:  // Invalid offset
            assert(false && "Invalid offset value");
            break;
    }
    return new_ray;
}

// Function to perturb a ray along the Y direction with small displacements
inline RationalRay perturbYRayRationals(const RationalRay &ray, uint offset) {
    RationalRay new_ray = ray;
    bigrational epsilon = getEpsilon(ray.v1[0]);  // Epsilon based on the x value

    switch (offset) {
        case 0:  // Slight perturbation in the X direction (positive epsilon)
            new_ray.v1[0] = ray.v1[0] + epsilon;
            break;
        case 1:  // Perturbation in both X and Z directions (positive epsilon)
            new_ray.v1[0] = ray.v1[0] + epsilon;
            new_ray.v1[2] = ray.v1[2] + epsilon;
            break;
        case 2:  // Perturbation in the Z direction (positive epsilon)
            new_ray.v1[2] = ray.v1[2] + epsilon;
            break;
        case 3:  // Perturbation in the X direction (negative epsilon) and Z direction (positive epsilon)
            new_ray.v1[0] = ray.v1[0] - epsilon;
            new_ray.v1[2] = ray.v1[2] + epsilon;
            break;
        case 4:  // Perturbation in the X direction (negative epsilon)
            new_ray.v1[0] = ray.v1[0] - epsilon;
            break;
        case 5:  // Perturbation in both X and Z directions (negative epsilon)
            new_ray.v1[0] = ray.v1[0] - epsilon;
            new_ray.v1[2] = ray.v1[2] - epsilon;
            break;
        case 6:  // Perturbation in the Z direction (negative epsilon)
            new_ray.v1[2] = ray.v1[2] - epsilon;
            break;
        case 7:  // Perturbation in the X direction (positive epsilon) and Z direction (negative epsilon)
            new_ray.v1[0] = ray.v1[0] + epsilon;
            new_ray.v1[2] = ray.v1[2] - epsilon;
            break;
        default:  // Invalid offset
            assert(false && "Invalid offset value");
            break;
    }
    return new_ray;
}

// Function to perturb a ray along the Z direction with small displacements
inline RationalRay perturbZRayRationals(const RationalRay &ray, uint offset) {
    RationalRay new_ray = ray;
    bigrational epsilon = getEpsilon(ray.v1[0]);  // Epsilon based on the x value

    switch (offset) {
        case 0:  // Slight perturbation in the X direction (positive epsilon)
            new_ray.v1[0] = ray.v1[0] + epsilon;
            break;
        case 1:  // Perturbation in both X and Y directions (positive epsilon)
            new_ray.v1[0] = ray.v1[0] + epsilon;
            new_ray.v1[1] = ray.v1[1] + epsilon;
            break;
        case 2:  // Perturbation in the Y direction (positive epsilon)
            new_ray.v1[1] = ray.v1[1] + epsilon;
            break;
        case 3:  // Perturbation in the X direction (negative epsilon) and Y direction (positive epsilon)
            new_ray.v1[0] = ray.v1[0] - epsilon;
            new_ray.v1[1] = ray.v1[1] + epsilon;
            break;
        case 4:  // Perturbation in the X direction (negative epsilon)
            new_ray.v1[0] = ray.v1[0] - epsilon;
            break;
        case 5:  // Perturbation in both X and Y directions (negative epsilon)
            new_ray.v1[0] = ray.v1[0] - epsilon;
            new_ray.v1[1] = ray.v1[1] - epsilon;
            break;
        case 6:  // Perturbation in the Y direction (negative epsilon)
            new_ray.v1[1] = ray.v1[1] - epsilon;
            break;
        case 7:  // Perturbation in the X direction (positive epsilon) and Y direction (negative epsilon)
            new_ray.v1[0] = ray.v1[0] + epsilon;
            new_ray.v1[1] = ray.v1[1] - epsilon;
            break;
        default:  // Invalid offset
            assert(false && "Invalid offset value");
            break;
    }
    return new_ray;
}



inline int perturbRayAndFindIntersTriRationals(const RationalRay &ray, const std::vector<genericPoint*> &in_verts, const std::vector<uint> &in_tris,
                                               const std::vector<uint> &tris_to_test,  std::vector<bigrational> &in_verts_rational)
{
    std::vector<uint> inters_tris;
    std::vector<IntersectionPointRationals> inter_rat;
    std::unordered_set<uint> unique_tris;

    for(uint i = 0; i <= 7; i++)
    {
        RationalRay p_ray;
        if(ray.dir == 'X')      p_ray = perturbXRayRationals(ray, i);
        else if(ray.dir == 'Y') p_ray = perturbYRayRationals(ray, i);
        else if(ray.dir == 'Z') p_ray = perturbZRayRationals(ray, i);

        for(uint t_id : tris_to_test)
        {
            if (unique_tris.find(t_id) != unique_tris.end()) continue;

            const uint idx0 = 3 * in_tris[3 * t_id];
            const uint idx1 = 3 * in_tris[3 * t_id + 1];
            const uint idx2 = 3 * in_tris[3 * t_id + 2];

            std::array<bigrational, 3> tv0 = {in_verts_rational[idx0], in_verts_rational[idx0 + 1], in_verts_rational[idx0 + 2]};
            std::array<bigrational, 3> tv1 = {in_verts_rational[idx1], in_verts_rational[idx1 + 1], in_verts_rational[idx1 + 2]};
            std::array<bigrational, 3> tv2 = {in_verts_rational[idx2], in_verts_rational[idx2 + 1], in_verts_rational[idx2 + 2]};

            if (checkIntersectionInsideTriangle3DRationals(p_ray, tv0, tv1, tv2))
            {
                inters_tris.push_back(t_id);
                unique_tris.insert(t_id);
                int intersection = segment_triangle_intersect_3d(&ray.v0[0], &ray.v1[0], &tv0[0], &tv1[0], &tv2[0]);
                if (intersection)
                {
                    std::array<bigrational, 3> p_int;
                    plane_line_intersection(&tv0[0], &tv1[0], &tv2[0], &ray.v0[0], &ray.v1[0], &p_int[0]);
                    inter_rat.emplace_back(p_int[0], p_int[1], p_int[2], t_id);
                }
            }
        }
    }

    if (inter_rat.empty()) return -1;

    auto compare = [&ray](const IntersectionPointRationals &a, const IntersectionPointRationals &b) {
        if (ray.dir == 'X') return a.lessThanX(b);
        else if (ray.dir == 'Y') return a.lessThanY(b);
        else return a.lessThanZ(b);
    };

    std::sort(inter_rat.begin(), inter_rat.end(), compare);

    return static_cast<int>(inter_rat.front().getTriId());
}
inline uint findPatchIdByTriId(const std::vector<phmap::flat_hash_set<uint>> &patches, uint t_id){

    for(const auto& patch : patches){
        for(uint t_id_patch : patch){
            if (t_id == t_id_patch) return t_id_patch;
        }
    }
    assert(false && "Triangle not found in any patch");
}

inline void eraseIntersectionPoints(std::vector<IntersectionPointRationals>& inter_rat, uint t_id_int) {
    inter_rat.erase(
            std::remove_if(inter_rat.begin(), inter_rat.end(),
                           [t_id_int](const IntersectionPointRationals &ipr) {
                               return ipr.getTriId() == t_id_int;
                           }),
            inter_rat.end());
}

inline bool copyIntersectionPoint(const std::vector<IntersectionPointRationals>& inter_rat, std::vector<IntersectionPointRationals>& inter_rat_tmp, uint t_id_int) {
    for (const IntersectionPointRationals &ipr : inter_rat) {
        if (ipr.getTriId() == t_id_int) {

            //if it'a already in the vector, don't add it
            bool found = false;
            for (const IntersectionPointRationals &ipr_tmp : inter_rat_tmp) {
                if (ipr_tmp.getTriId() == ipr.getTriId()) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                inter_rat_tmp.push_back(ipr);
                return true;
            }
        }
    }
    return false;
}

///:::::::::::::: DEBUG RATIONALS FUNCTIONS ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


inline void printInfoTriangleRationals(RationalRay &rational_ray, std::vector <bigrational> &tv0_aux, std::vector <bigrational> &tv1_aux, std::vector <bigrational> &tv2_aux,
                                       uint *tv_aux, uint &t_id, bool &print_ray){

    if(print_ray){
        std::cout << "Info about ray:" << std::endl;
        std::cout << "Start point: " << "v0: " << rational_ray.v0[0].get_d() << " v1: " << rational_ray.v0[1].get_d() << " v2: " << rational_ray.v0[2].get_d() << std::endl;
        std::cout << "End point: "   << "v0: " << rational_ray.v1[0].get_d() << " v1: " << rational_ray.v1[1].get_d() << " v2: " << rational_ray.v1[2].get_d() << std::endl;
        print_ray = false;
    }

    std::cout << "Info about triangle whith t_id: " << t_id << std::endl;
    std::cout << "v_id: " <<tv_aux[0]<< " v0: " << tv0_aux[0].get_d() << " v1: " << tv0_aux[1].get_d() << " v2: " << tv0_aux[2].get_d() << std::endl;
    std::cout << "v_id: " <<tv_aux[1]<< " v0: " << tv1_aux[0].get_d() << " v1: " << tv1_aux[1].get_d() << " v2: " << tv1_aux[2].get_d() << std::endl;
    std::cout << "v_id: " <<tv_aux[2]<< " v0: " << tv2_aux[0].get_d() << " v1: " << tv2_aux[1].get_d() << " v2: " << tv2_aux[2].get_d() << std::endl;


}

inline int maxComponentInTriangleNormalRationals(bigrational &ov1x, bigrational &ov1y, bigrational &ov1z, bigrational &ov2x, bigrational &ov2y, bigrational &ov2z, bigrational &ov3x, bigrational &ov3y, bigrational &ov3z)
{

    bigrational v3x = ov3x - ov2x;
    bigrational v3y = ov3y - ov2y;
    bigrational v3z = ov3z - ov2z;
    bigrational v2x = ov2x - ov1x;
    bigrational v2y = ov2y - ov1y;
    bigrational v2z = ov2z - ov1z;

    bigrational nvx = v2y * v3z - v2z * v3y;
    bigrational nvy = v3x * v2z - v3z * v2x;
    bigrational nvz = v2x * v3y - v2y * v3x;

    bigrational nvxc = fabs(nvx);
    bigrational nvyc = fabs(nvy);
    bigrational nvzc = fabs(nvz);

    if (nvxc >= nvyc && nvxc >= nvzc) return 0;
    if (nvyc >= nvxc && nvyc >= nvzc) return 1;
    return 2;
}


inline bigrational fabs(bigrational x)
{
   return (x < bigrational()) ? x.negation() : x;
}



///::::::::::::::DEBUG PARSER DIFF :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::///

inline void savePartsToFile(const std::vector<std::vector<std::vector<unsigned int>>>& parts_to_color,
                     const std::string& filename,
                     bool debug_impl, bool patch_view) {
    // Open a file for writing
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Define the color labels
    std::vector<std::string> labels = {"GRAY Parts", "RED Parts", "GREEN Parts", "WHITE Parts", "PATCH Parts"};

    // Adjust the number of parts to save
    size_t num_parts = debug_impl ? 4 : 3;
    num_parts = patch_view ? (num_parts + 1) : num_parts;

    // Iterate through the parts
    for (size_t i = 0; i < num_parts; ++i) {
        // Write the label
        if (i < labels.size()) {
            outFile << labels[i] << "\n";
        } else {
            outFile << "UNKNOWN Parts\n";
        }

        // Write the indices of vertices
        for (const auto& vertex : parts_to_color[i]) {
            for (unsigned int index : vertex) {
                outFile << index << " ";
            }
            outFile << "\n";
        }
    }

    // Close the file
    outFile.close();
}

inline bool parseFileToParts(const std::string& filename,
                      std::vector<std::vector<std::vector<unsigned int>>>& parts_to_color,
                      bool debug_impl, bool patch_view) {
    // Open the file for reading
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }

    // Define labels and corresponding indices
    std::unordered_map<std::string, size_t> labelToIndex = {
            {"GRAY Parts", 0},
            {"RED Parts", 1},
            {"GREEN Parts", 2},
            {"WHITE Parts", 3},
            {"PATCH Parts", 4}// New label for debug_impl
    };

    // Adjust the parts_to_color size
    size_t num_parts = debug_impl ? 4 : 3;
    num_parts = patch_view ? (num_parts + 1) : num_parts;
    parts_to_color.resize(num_parts);

    // Parsing state
    std::string line;
    size_t currentPart = 0;

    // Parse the file line by line
    while (std::getline(inFile, line)) {
        // Check if the line is a label
        if (labelToIndex.find(line) != labelToIndex.end()) {
            currentPart = labelToIndex[line];
        } else {
            // Parse indices in the current line
            std::istringstream iss(line);
            std::vector<unsigned int> vertices;
            unsigned int vertex;
            while (iss >> vertex) {
                vertices.push_back(vertex);
            }
            // Add parsed vertices to the correct part
            if (currentPart < parts_to_color.size()) {
                parts_to_color[currentPart].push_back(vertices);
            }
        }
    }

    inFile.close();
    return true;
}


inline void savePatchesTriangles(const std::string& filename, const std::vector<int>& p_ids, const std::vector<phmap::flat_hash_set<uint>>& patches) {
    // Open the file to write
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    // Create a hash set for quick lookup of p_ids
    std::unordered_set<int> p_id_set(p_ids.begin(), p_ids.end());

    // Iterate through patches
    for (size_t i = 0; i < patches.size(); ++i) {
         // Check if the current patch ID is in p_ids
            // Write patch header
            outfile << "Patch " << i << "\n";

            // Write all triangle IDs in the current patch
            for (const auto& triangle_id : patches[i]) {
                outfile << triangle_id << '\n';
            }

    }

    outfile.close();
    std::cout << "Triangles saved to " << filename << " successfully.\n";
}

inline bool parsePatches(const std::string& filename, std::vector<phmap::flat_hash_set<uint>>& patches) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for reading.\n";
        return false;
    }

    patches.clear(); // Ensure the patches vector is empty before populating
    std::string line;
    phmap::flat_hash_set<uint> currentPatch;
    size_t currentPatchIndex = 0;

    // Regular expression to identify a "Patch X" line
    std::regex patchHeaderRegex(R"(Patch\s+(\d+))");

    while (std::getline(infile, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        if (line.empty()) continue;

        std::smatch match;
        if (std::regex_match(line, match, patchHeaderRegex)) {
            // If we encounter a new patch header, store the previous patch
            if (!currentPatch.empty()) {
                patches.push_back(std::move(currentPatch));
                currentPatch.clear();
            }
            // Get the patch index (though we use the order of appearance, not the index explicitly)
            currentPatchIndex = std::stoul(match[1].str());
        } else {
            // Parse triangle ID
            try {
                uint triangleID = static_cast<uint>(std::stoul(line));
                currentPatch.insert(triangleID);
            } catch (const std::exception& e) {
                std::cerr << "Error parsing triangle ID: " << line << "\n";
                return false;
            }
        }
    }

    // Push the last patch if it's not empty
    if (!currentPatch.empty()) {
        patches.push_back(std::move(currentPatch));
    }

    infile.close();
    return true;
}




inline uint classifyTriangles(FastTrimesh &tm, const Labels &labels, Data &data)
{
    tm.resetTrianglesInfo();
    data.t_ids_intersection.clear();
    data.t_ids_union.clear();
    data.t_ids_subtraction.clear();

    uint num_tris_in_final_solution = 0;

    for(uint t_id = 0; t_id < tm.numTris(); ++t_id)
    {
        const auto &surf = labels.surface[t_id];
        const auto &in = labels.inside[t_id];

        bool is_intersection = (surf ^ in).count() == labels.num;
        bool is_union        = in.count() == 0;
        bool is_subtraction  = (surf[0] && in.count() == 0) ||
                               (!surf[0] && in[0] && in.count() == 1);

        bool counted = false;

        if(is_intersection)
        {
            data.t_ids_intersection.push_back(t_id);
            counted = true;
        }

        if(is_union)
        {
            data.t_ids_union.push_back(t_id);
            counted = true;
        }

        if(is_subtraction)
        {
            data.t_ids_subtraction.push_back(t_id);
            counted = true;
        }

        if(counted)
        {
            if(tm.triInfo(t_id) == 0) // only count if not already included
                num_tris_in_final_solution++;

            tm.setTriInfo(t_id, 1);
        }
    }

    // Fix orientation for subtraction triangles not from model 0
    for(uint t_id : data.t_ids_subtraction)
    {
        if(labels.surface[t_id][0] != 1)
            tm.flipTri(t_id);
    }

    return num_tris_in_final_solution;
}

inline void saveTriangleIDsToFile(const Data &data)
{
    std::ofstream f_intersection("intersection_ids.txt");
    for(auto id : data.t_ids_intersection)
        f_intersection << id << "\n";

    std::ofstream f_union("union_ids.txt");
    for(auto id : data.t_ids_union)
        f_union << id << "\n";

    std::ofstream f_subtraction("subtraction_ids.txt");
    for(auto id : data.t_ids_subtraction)
        f_subtraction << id << "\n";
}

inline void loadTriangleIDsFromFile(Data &data)
{
    data.t_ids_intersection.clear();
    data.t_ids_union.clear();
    data.t_ids_subtraction.clear();

    std::ifstream f_intersection("intersection_ids.txt");
    std::ifstream f_union("union_ids.txt");
    std::ifstream f_subtraction("subtraction_ids.txt");

    uint id;
    while(f_intersection >> id)
        data.t_ids_intersection.push_back(id);

    while(f_union >> id)
        data.t_ids_union.push_back(id);

    while(f_subtraction >> id)
        data.t_ids_subtraction.push_back(id);
}

