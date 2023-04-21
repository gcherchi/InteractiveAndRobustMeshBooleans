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

#include <cinolib/meshes/trimesh.h>
#include <cinolib/find_intersections.h>

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool manifold(const Trimesh<> & m)
{
    for(uint vid=0; vid<m.num_verts(); ++vid)
    {
        if(!m.vert_is_manifold(vid)) return false;
    }
    return true;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool watertight(const Trimesh<> & m)
{
    for(uint eid=0; eid<m.num_edges(); ++eid)
    {
        if(m.edge_is_boundary(eid)) return false;
    }
    return true;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool local_orientation(const Trimesh<> & m)
{
    std::vector<uint> visited(m.num_polys(),false);
    for(uint pid=0; pid<m.num_polys(); ++pid)
    {
        if(visited[pid]) continue;
        std::queue<uint> q;
        q.push(pid);
        visited[pid] = true;
        while(!q.empty())
        {
            uint pid = q.front();
            q.pop();
            for(uint nbr : m.adj_p2p(pid))
            {
                uint eid = m.edge_shared(pid, nbr);
                if(m.edge_is_CCW(eid,pid)==m.edge_is_CCW(eid,nbr)) return false;
                if(!visited[nbr])
                {
                    visited[nbr] = true;
                    q.push(nbr);
                }
            }
        }
    }
    return true;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// WARNING: this is correct only if the mesh has a single connected component
bool global_orientation(const Trimesh<> & m)
{
    vec3d  O   = m.bbox().min - vec3d(1,0,0); // surely external
    double vol = 0;
    for(uint pid=0; pid<m.num_polys(); ++pid)
    {
        // sum signed volumes
        vol += tet_volume(m.poly_vert(pid,0),
                          m.poly_vert(pid,1),
                          m.poly_vert(pid,2), O);
    }
    // if faces are CCW volume is negative
    return (vol<0);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool intersection_free(const Trimesh<> & m)
{
    std::set<ipair> intersections;
    find_intersections(m, intersections);
    return (intersections.empty());
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main(int argc, char **argv)
{
    Trimesh<> m(argv[1]);
    std::cout << "Manifold check: "           << (manifold(m)          ?"passed":"failed") << std::endl;
    std::cout << "Watertight check: "         << (watertight(m)        ?"passed":"failed") << std::endl;
    std::cout << "Local  Orientation check: " << (local_orientation(m) ?"passed":"failed") << std::endl;
    std::cout << "Global Orientation check: " << (global_orientation(m)?"passed":"failed") << std::endl;
    std::cout << "Intersection check: "       << (intersection_free(m) ?"passed":"failed") << std::endl;
    return EXIT_SUCCESS;
}
