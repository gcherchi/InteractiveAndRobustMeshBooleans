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

#ifndef EXACT_BOOLEANS_BOOLEANS_H
#define EXACT_BOOLEANS_BOOLEANS_H

#include "processing.h"
#include "aux_structure.h"
#include "triangle_soup.h"
#include "intersection_classification.h"
#include "triangulation.h"
#include <cinolib/octree.h>

#include "intersect_custom.h"
#include "rationals_code/intersect_point_rationals.h"
#include <regex>

#include <bitset>

struct Data {
    std::vector<uint> t_ids_intersection ;
    std::vector<uint> t_ids_union ;
    std::vector<uint> t_ids_subtraction ;
    std::vector<uint> t_ids_debug;
    std::vector <uint> t_ids_inters_ray;
    uint t_id_debug = -1;
};

struct Labels
{
    std::vector< std::bitset<NBIT> > surface;
    std::vector< std::bitset<NBIT> > inside;
    uint num;
};

struct Ray
{
    explicitPoint3D v0;
    explicitPoint3D v1;
    char dir = 'X';
    int tv[3] = {-1, -1, -1};
};

struct DuplTriInfo
{
    uint t_id;
    uint l_id;
    bool w;
};

enum BoolOp {UNION, INTERSECTION, SUBTRACTION, XOR, DEBUG, NONE};

enum IntersInfo {DISCARD, NO_INT, INT_IN_V0, INT_IN_V1, INT_IN_V2, INT_IN_EDGE01, INT_IN_EDGE12, INT_IN_EDGE20, INT_IN_TRI};

struct less_than_GP_on_X // lessThan GenericPoint along X
{
    bool operator() (const std::pair<genericPoint*, uint> &p0, const std::pair<genericPoint*, uint> &p1) const
    {
        return (genericPoint::lessThanOnX(*p0.first, *p1.first) <= 0);
    }
};

struct less_than_GP_on_Y // lessThan GenericPoint along Y
{
    bool operator() (const std::pair<genericPoint*, uint> &p0, const std::pair<genericPoint*, uint> &p1) const
    {
        return (genericPoint::lessThanOnY(*p0.first, *p1.first) <= 0);
    }
};

struct less_than_GP_on_Z // lessThan GenericPoint along Z
{
    bool operator() (const std::pair<genericPoint*, uint> &p0, const std::pair<genericPoint*, uint> &p1) const
    {
        return (genericPoint::lessThanOnZ(*p0.first, *p1.first) <= 0);
    }
};

inline void customBooleanPipeline(std::vector<genericPoint*>& arr_verts, std::vector<uint>& arr_in_tris,
                                  std::vector<uint>& arr_out_tris, std::vector<std::bitset<NBIT>>& arr_in_labels,
                                  std::vector<DuplTriInfo>& dupl_triangles, Labels& labels,
                                  std::vector<phmap::flat_hash_set<uint>>& patches, cinolib::Octree& octree,
                                  const BoolOp &op, std::vector<double> &bool_coords, std::vector<uint> &bool_tris,
                                  std::vector< std::bitset<NBIT>> &bool_labels, Data &data);

inline void booleanPipeline(const std::vector<double> &in_coords, const std::vector<uint> &in_tris,
                            const std::vector<uint> &in_labels, const BoolOp &op, std::vector<double> &bool_coords,
                            std::vector<uint> &bool_tris, std::vector< std::bitset<NBIT> > &bool_labels, Data &data);


inline void customArrangementPipeline(const std::vector<double> &in_coords, const std::vector<uint> &in_tris, const std::vector<uint> &in_labels,
                                      std::vector<uint> &arr_in_tris, std::vector< std::bitset<NBIT>> &arr_in_labels,
                                      point_arena& arena, std::vector<genericPoint *> &vertices, std::vector<uint> &arr_out_tris, Labels &labels,
                                      cinolib::Octree &octree, std::vector<DuplTriInfo> &dupl_triangles, bool parallel);

inline void customRemoveDegenerateAndDuplicatedTriangles(const std::vector<genericPoint*> &verts, std::vector<uint> &tris,
                                                         std::vector< std::bitset<NBIT> > &labels, std::vector<DuplTriInfo> &dupl_triangles,
                                                         bool parallel);

inline void customDetectIntersections(const TriangleSoup &ts, std::vector<std::pair<uint, uint> > &intersection_list, cinolib::Octree &o);

inline void addDuplicateTrisInfoInStructures(const std::vector<DuplTriInfo> &dupl_tris, std::vector<uint> &in_tris,
                                             std::vector<std::bitset<NBIT>> &in_labels, cinolib::Octree &octree);

inline void computeAllPatches(FastTrimesh &tm, const Labels &labels, std::vector<phmap::flat_hash_set<uint>> &patches, bool parallel);

inline void computeSinglePatch(FastTrimesh &tm, uint seed_t, const Labels &labels, phmap::flat_hash_set<uint> &patch);
inline void computeSinglePatch(FastTrimesh &tm, uint seed_t, const Labels &labels, phmap::flat_hash_set<uint> &patch, const std::vector<std::array<uint, 3>>& adjT2E);

inline void findRayEndpoints(const FastTrimesh &tm, const phmap::flat_hash_set<uint> &patch, const cinolib::vec3d &max_coords, Ray &ray);

inline bool intersects_box(const cinolib::Octree& tree, const cinolib::AABB & b, phmap::flat_hash_set<uint> & ids);

inline void computeInsideOut(const FastTrimesh &tm, const std::vector<phmap::flat_hash_set<uint>> &patches, const cinolib::Octree &octree,
                             const std::vector<genericPoint *> &in_verts, const std::vector<uint> &in_tris,
                             const std::vector<std::bitset<NBIT>> &in_labels, const cinolib::vec3d &max_coords, Labels &labels);

inline void pruneIntersectionsAndSortAlongRay(const Ray &ray, const std::vector<genericPoint*> &in_verts,
                                              const std::vector<uint> &in_tris, const std::vector<std::bitset<NBIT>> &in_labels,
                                              const phmap::flat_hash_set<uint> &tmp_inters, const std::bitset<NBIT> &patch_surface_label,
                                              std::vector<uint> &inters_tris);

inline void analyzeSortedIntersections(const Ray &ray, const std::vector<genericPoint*> &in_verts, const std::vector<uint> &in_tris,
                                       const std::vector<std::bitset<NBIT>> &in_labels, const std::vector<uint> &sorted_inters,
                                       std::bitset<NBIT> &patch_inner_label);

inline bool triContainsVert(uint t_id, uint v_id, const std::vector<uint> &in_tris);

inline void findVertRingTris(uint v_id, const std::bitset<NBIT> &ref_label, const phmap::flat_hash_set<uint> &inters_tris,
                             const std::vector<uint> &in_tris, const std::vector<std::bitset<NBIT>> &in_labels,
                             std::vector<uint> &one_ring);

inline void findEdgeTris(uint ev0_id, uint ev1_id, const std::bitset<NBIT> &ref_label, const phmap::flat_hash_set<uint> &inters_tris,
                             const std::vector<uint> &in_tris, const std::vector<std::bitset<NBIT>> &in_labels,
                             std::vector<uint> &edge_tris);

inline Ray perturbXRay(const Ray &ray, uint offset);

inline Ray perturbYRay(const Ray &ray, uint offset);

inline Ray perturbZRay(const Ray &ray, uint offset);

inline int perturbRayAndFindIntersTri(const Ray &ray, const std::vector<genericPoint*> &in_verts, const std::vector<uint> &in_tris,
                                       const std::vector<uint> &tris_to_test);

inline IntersInfo fast2DCheckIntersectionOnRay(const Ray &ray, const explicitPoint3D &tv0, const explicitPoint3D &tv1, const explicitPoint3D &tv2);

inline bool checkIntersectionInsideTriangle3D(const Ray &ray, const explicitPoint3D &tv0, const explicitPoint3D &tv1, const explicitPoint3D &tv2);

inline bool checkIntersectionInsideTriangle3DImplPoints(const Ray &ray, const genericPoint *tv0, const genericPoint *tv1, const genericPoint *tv2);

inline void sortIntersectedTrisAlongX(const Ray &ray, const std::vector<genericPoint*> &in_verts,
                                      const std::vector<uint> &in_tris, std::vector<uint> &inters_tris);

inline void sortIntersectedTrisAlongY(const Ray &ray, const std::vector<genericPoint*> &in_verts,
                                      const std::vector<uint> &in_tris, std::vector<uint> &inters_tris);

inline void sortIntersectedTrisAlongZ(const Ray &ray, const std::vector<genericPoint*> &in_verts,
                                      const std::vector<uint> &in_tris, std::vector<uint> &inters_tris);

inline uint checkTriangleOrientation(const Ray &ray, const explicitPoint3D &tv0, const explicitPoint3D &tv1, const explicitPoint3D &tv2);

inline void propagateInnerLabelsOnPatch(const phmap::flat_hash_set<uint> &patch_tris, const std::bitset<NBIT> &patch_inner_label, Labels &labels);

inline void computeFinalExplicitResult(const FastTrimesh &tm, const Labels &labels, uint num_tris_in_final_res,
                                       std::vector<double> &out_coords, std::vector<uint> &out_tris, std::vector<std::bitset<NBIT>> &out_label, bool flat_array, Data &data);

inline uint boolIntersection(FastTrimesh &tm, const Labels &labels, Data &data);

inline uint boolUnion(FastTrimesh &tm, const Labels &labels, Data &data);

inline uint boolSubtraction(FastTrimesh &tm, const Labels &labels, Data &data);

inline uint boolXOR(FastTrimesh &tm, const Labels &labels);

inline uint bitsetToUint(const std::bitset<NBIT> &b);

inline bool consistentWinding(const uint *t0, const uint *t1);

/// :::::::::::::: DEBUG :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline std::string printBitset(const std::bitset<NBIT> &b, uint num_label); // just for debug

inline void saveOutputWithLabels(const std::string &filename, cinolib::Trimesh<> &m, const std::vector<std::bitset<NBIT>> &labels);

inline void loadInputWithLabels(const std::string &filename, std::vector<double> &coords, std::vector<uint> &tris, std::vector<std::bitset<NBIT> > &labels);

inline void loadInputWithLabels(const std::string &filename, std::vector<double> &coords, std::vector<uint> &tris, std::vector<uint> &labels);

///_:::::::::::::::::: RATIONALS STRUCTS ::::::::::::::::::::::::::::::::::::::::::::
struct RationalRay{
    std::array<bigrational,3> v0;
    std::array<bigrational,3> v1;
    char dir = 'X';
    int tv[3] = {-1, -1, -1};

};

struct BoundingBox {
    bigrational xmin, xmax, ymin, ymax, zmin, zmax;
};




///::::::::::::::::::: RATIONALS FUNCTIONS ::::::::::::::::::::::::::::::::::::::::::

///::::::::::::::::::: RATIONALS FUNCTIONS ::::::::::::::::::::::::::::::::::::::::::

inline void setExplicitVertex(const FastTrimesh &tm, std::vector<bigrational> &in_verts_rational, uint vertex_id, bigrational &x, bigrational &y, bigrational &z);
inline void computeInsideOutCustom(const FastTrimesh &tm, const std::vector<phmap::flat_hash_set<uint>> &patches, const cinolib::Octree &octree,
                                   const std::vector<genericPoint *> &in_verts, const std::vector<uint> &in_tris,
                                   const std::vector<std::bitset<NBIT>> &in_labels, const cinolib::vec3d &max_coords, Labels &labels, Data &data);
inline void findRayEndpointsCustom(const FastTrimesh &tm, const phmap::flat_hash_set<uint> &patch, const cinolib::vec3d &max_coords, Ray &ray, RationalRay &rational_ray, const std::vector<genericPoint *> &in_verts, std::vector<bigrational> &in_verts_rational, bool &is_rational, bool debug, Data &data);

inline void findIntersectionsAlongRayRationals(const FastTrimesh &tm, const std::vector<phmap::flat_hash_set<uint>> &patches, const cinolib::Octree& tree, const std::vector<genericPoint *> &in_verts,
                                               const std::vector<std::bitset<NBIT>> &in_labels, Labels &labels, const RationalRay &rational_ray,
                                               uint curr_p_id, phmap::flat_hash_set<uint> &tmp_inters, std::vector<IntersectionPointRationals> &inter_rat, std::vector<bigrational> &in_verts_rational, const std::vector<uint> &in_tris, Data &data);

inline IntersInfo fast2DCheckIntersectionOnRayRationals(const RationalRay &ray, const std::vector<bigrational> &tv0, const std::vector<bigrational> &tv1, const std::vector<bigrational> &tv2);

inline bool checkIntersectionInsideTriangle3DRationals(const RationalRay &ray, const std::array<bigrational,3> &tv0, const std::array<bigrational,3> &tv1, const std::array<bigrational,3> &tv2);
inline uint checkTriangleOrientationRationals(const RationalRay &ray, const std::vector<bigrational> &tv0, const std::vector<bigrational> &tv1, const std::vector<bigrational> &tv2);

inline void pruneIntersectionsAndSortAlongRayRationals(const RationalRay &ray, const FastTrimesh &tm, const std::vector<genericPoint*> &in_verts,
                                                       const std::vector<uint> &in_tris, const std::vector<std::bitset<NBIT>> &in_labels,
                                                       const phmap::flat_hash_set<uint> &tmp_inters, const std::bitset<NBIT> &patch_surface_label,
                                                       std::vector<IntersectionPointRationals> &inter_rat, std::vector<uint> &inters_tris_rat, Labels &labels, std::bitset<NBIT> &patch_surface_label_tmp,
                                                       const std::vector<phmap::flat_hash_set<uint>> &patches, std::vector<bigrational> &in_verts_rational, Data &data);
inline void analyzeSortedIntersectionsRationals(const RationalRay &rational_ray, const FastTrimesh &tm, const std::vector<genericPoint*> &in_verts,
                                                std::vector<IntersectionPointRationals> &inter_rat, std::bitset<NBIT> &patch_inner_label, Labels &labels, const std::vector<std::bitset<NBIT>> &in_labels,
                                                std::vector<bigrational> &in_verts_rational, const std::vector<uint> &in_tris, Data &data);

inline int perturbRayAndFindIntersTriRationals(const RationalRay &ray, const std::vector<genericPoint*> &in_verts, const std::vector<uint> &in_tris,
                                                const std::vector<uint> &tris_to_test, std::vector<bigrational> &in_verts_rational);

inline bigrational next_after(const bigrational& x, const bigrational& target);

inline uint findPatchIdByTriId(const std::vector<phmap::flat_hash_set<uint>> &patches, uint t_id);

inline RationalRay perturbXRayRationals(const RationalRay &ray, uint offset);
inline RationalRay perturbYRayRationals(const RationalRay &ray, uint offset);
inline RationalRay perturbZRayRationals(const RationalRay &ray, uint offset);

inline void eraseIntersectionPoints(std::vector<IntersectionPointRationals>& inter_rat, uint t_id_int);

inline BoundingBox calculateBoundingBox(const std::array<bigrational, 3>& tv0,
                                 const std::array<bigrational, 3>& tv1,
                                 const std::array<bigrational, 3>& tv2);
inline bool rayIntersectAABB(const RationalRay &ray, const BoundingBox & aabb);

inline bool copyIntersectionPoint(const std::vector<IntersectionPointRationals>& inter_rat, std::vector<IntersectionPointRationals>& inter_rat_tmp, uint t_id_int);

inline bool isNormalCorrect(
         bigrational& ov1x,  bigrational& ov1y,  bigrational& ov1z,
         bigrational& ov2x,  bigrational& ov2y,  bigrational& ov2z,
         bigrational& ov3x,  bigrational& ov3y,  bigrational& ov3z,
         bigrational& px,  bigrational& py,  bigrational& pz);

bool isIntersectionValid(const std::vector<bigrational>& inter, const RationalRay& rational_ray);

inline int maxComponentInTriangleNormalRationals(bigrational &ov1x, bigrational &ov1y, bigrational &ov1z, bigrational &ov2x, bigrational &ov2y, bigrational &ov2z, bigrational &ov3x, bigrational &ov3y, bigrational &ov3z);
inline bigrational fabs(bigrational x);

inline bigrational getEpsilon(const bigrational& value);

////::::::::::: DEBUG CUSTOM ::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void printInfoTriangleRationals(RationalRay &rational_ray, std::vector <bigrational> &tv0_aux, std::vector <bigrational> &tv1_aux, std::vector <bigrational> &tv2_aux,
                                       uint *tv_aux, uint &t_id, bool &print_ray);


///::::::::: DEUBUG PARSER DIFF :::::::::::::::::::::::::::::::::::::::::::::::::::
inline bool parseFileToParts(const std::string& filename,
                      std::vector<std::vector<std::vector<unsigned int>>>& parts_to_color,
                      bool debug_impl);

inline void savePartsToFile(const std::vector<std::vector<std::vector<unsigned int>>>& parts_to_color,
                     const std::string& filename,
                     bool debug_impl, bool patch_view);

inline void savePatchesTriangles(const std::string& filename, const std::vector<int>& p_ids, const std::vector<phmap::flat_hash_set<uint>>& patches);
inline bool parsePatches(const std::string& filename, std::vector<phmap::flat_hash_set<uint>>& patches);
inline uint classifyTriangles(FastTrimesh &tm, const Labels &labels, Data &data);


#include <fstream>
#include <sstream>
inline void saveTriangleIDsToFile(const Data &data);
inline void loadTriangleIDsFromFile(Data &data);





#include "booleans.cpp"

#endif //EXACT_BOOLEANS_BOOLEANS_H
