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

#ifdef _MSC_VER // Workaround for known bugs and issues on MSVC
#define _HAS_STD_BYTE 0  // https://developercommunity.visualstudio.com/t/error-c2872-byte-ambiguous-symbol/93889
#define NOMINMAX // https://stackoverflow.com/questions/1825904/error-c2589-on-stdnumeric-limitsdoublemin
#endif

#include "intersect_custom.h"
#include "numerics.h"
#include <vector>
#include <random>
#include <cinolib/meshes/meshes.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <booleans.h>
#include <tbb/parallel_for.h>
#include <iostream>

#include <tbb/tbb.h>
#include <mutex>

void debugIntersectionTestParallel(const std::vector<bigrational> &verts,
                                   const std::vector<uint> &tris,
                                   const RationalRay &ray)
{
    std::mutex print_mutex;  // per sincronizzare cout

    tbb::parallel_for(tbb::blocked_range<uint>(0, tris.size() / 3),
                      [&](const tbb::blocked_range<uint>& r) {
                          for (uint tid = r.begin(); tid != r.end(); ++tid) {
                              const uint id0 = tris[3 * tid];
                              const uint id1 = tris[3 * tid + 1];
                              const uint id2 = tris[3 * tid + 2];

                              std::array<bigrational, 3> v0, v1, v2;
                              for (int i = 0; i < 3; ++i) {
                                  v0[i] = verts[3 * id0 + i];
                                  v1[i] = verts[3 * id1 + i];
                                  v2[i] = verts[3 * id2 + i];
                              }

                              bigrational n[3], l[3];
                              triangle_normal(&v0[0], &v1[0], &v2[0], &n[0]);

                              for (int i = 0; i < 3; ++i)
                                  l[i] = ray.v1[i] - ray.v0[i];

                              if (dot(&l[0], &n[0]) == bigrational()) {
                                  std::lock_guard<std::mutex> lock(print_mutex);
                                  std::cout << "Ray and triangle " << tid << " are coplanar" << std::endl;
                                  continue;
                              }

                              bigrational p[3];
                              plane_line_intersection(&v0[0], &v1[0], &v2[0], &ray.v0[0], &ray.v1[0], &p[0]);

                              {
                                  std::lock_guard<std::mutex> lock(print_mutex);
                                  std::cout << "Intersection at triangle " << tid << ": ("
                                            << p[0] << ", " << p[1] << ", " << p[2] << ")" << std::endl;
                              }
                          }
                      });
}


void generateGridMeshInput(std::vector<bigrational> &verts_rational,
                           std::vector<uint> &tris,
                           RationalRay &ray,
                           int grid_size = 20) {
    verts_rational.clear();
    tris.clear();

    // Create vertices in a 3D grid
    for (int z = 0; z <= grid_size; ++z) {
        for (int y = 0; y <= grid_size; ++y) {
            for (int x = 0; x <= grid_size; ++x) {
                verts_rational.push_back(bigrational(x) + bigrational(1/(x + 1)));
                verts_rational.push_back(bigrational(y) + bigrational(1/(x + 1)));
                verts_rational.push_back(bigrational(z) + bigrational(1/(x + 1)));
            }
        }
    }

    auto index = [grid_size](int x, int y, int z) {
        return x + (grid_size + 1) * (y + (grid_size + 1) * z);
    };

    // Build triangles (two per cube face)
    for (int z = 0; z < grid_size; ++z) {
        for (int y = 0; y < grid_size; ++y) {
            for (int x = 0; x < grid_size; ++x) {
                uint v0 = index(x, y, z);
                uint v1 = index(x + 1, y, z);
                uint v2 = index(x, y + 1, z);
                uint v3 = index(x + 1, y + 1, z);
                uint v4 = index(x, y, z + 1);
                uint v5 = index(x + 1, y, z + 1);
                uint v6 = index(x, y + 1, z + 1);
                uint v7 = index(x + 1, y + 1, z + 1);

                // Bottom face
                tris.push_back(v0); tris.push_back(v1); tris.push_back(v2);
                tris.push_back(v1); tris.push_back(v3); tris.push_back(v2);

                // Top face
                tris.push_back(v4); tris.push_back(v5); tris.push_back(v6);
                tris.push_back(v5); tris.push_back(v7); tris.push_back(v6);

                // Front face
                tris.push_back(v0); tris.push_back(v1); tris.push_back(v4);
                tris.push_back(v1); tris.push_back(v5); tris.push_back(v4);
            }
        }
    }

    // Define a ray that cuts through the grid diagonally
    ray.v0 = { bigrational(-10), bigrational(-10), bigrational(-10) };
    ray.v1 = { bigrational(grid_size * 2), bigrational(grid_size * 2), bigrational(grid_size * 2) };
}


int main(int argc, char *argv[]) {

    std::vector<bigrational> verts_rational;
    std::vector<uint> tris;
    RationalRay ray;
    generateGridMeshInput(verts_rational, tris, ray, 20);
    std::cout << "Generated grid mesh with " << verts_rational.size() / 3 << " vertices and "
              << tris.size() / 3 << " triangles." << std::endl;
    // Debug intersection test
    debugIntersectionTestParallel(verts_rational, tris, ray);


    return 0;

}
