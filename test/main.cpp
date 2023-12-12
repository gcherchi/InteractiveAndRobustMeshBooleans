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

#include <cinolib/meshes/trimesh.h>
#include <cinolib/find_intersections.h>
#include <cinolib/connected_components.h>

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
    std::string file = argv[1];
    Trimesh<> m(file.c_str());

    bool is_manifold = manifold(m);
    bool is_watertight = watertight(m);
    bool is_local_well_oriented = local_orientation(m);
    bool is_global_well_oriented = global_orientation(m);
    bool is_intersection_free = intersection_free(m);

    if(is_manifold && is_watertight && is_local_well_oriented && is_global_well_oriented && is_intersection_free)
    {
        std::cout << "file " << file << " OK" << std::endl;

        //std::string base_name = file.substr(0, file.size()-6);
        //m.rotate(vec3d(0.0, 0.0, 1.0), 3.14/3.0);
        //m.save((base_name + "_c.obj").c_str());
        return EXIT_SUCCESS;
    }

    if(!is_manifold) std::cout << file << " is not manifold" << std::endl;
    if(!is_watertight) std::cout << file << " is not watertight" << std::endl;
    if(!is_local_well_oriented) std::cout << file << " is not local well oriented" << std::endl;
    if(!is_global_well_oriented) std::cout << file << " is not global well oriented" << std::endl;
    if(!is_intersection_free) std::cout << file << " is not intersection free" << std::endl;

    return EXIT_SUCCESS;
}


/*
int main(int argc, char **argv)
{
    std::string f1 = argv[1];
    std::string f2 = argv[2];

    Trimesh<> m1(f1.c_str());
    Trimesh<> m2(f2.c_str());

    bool ok = true;

    if(manifold(m1) != manifold(m2))                            ok = false;
    if(watertight(m1) != watertight(m2))                        ok = false;
    if(local_orientation(m1) != local_orientation(m2))          ok = false;
    if(global_orientation(m1) != local_orientation(m2))         ok = false;
    if(connected_components(m1) != connected_components(m2))    ok = false;
    if(m1.Euler_characteristic() != m2.Euler_characteristic())  ok = false;

    std::string base_name = f1.substr(0, f2.size()-6);

    if(!ok)
    {
        std::cout << base_name << " PROBLEM!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << base_name << " OK" << std::endl;
    return EXIT_SUCCESS;
}
*/