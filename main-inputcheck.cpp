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
#include "booleans.h"

int main(int argc, char **argv)
{
    if(argc != 2)
    {
        std::cout << "syntax error!" << std::endl;
        std::cout << "./exact_boolean_inputcheck <file-to-check>" << std::endl;
        return -1;
    }

    cinolib::Trimesh<> m(argv[1]);

    // CHECK if mesh is MANIFOLD and WATERTIGHT
    for(uint v_id = 0; v_id < m.num_verts(); v_id++)
    {
        if(!m.vert_is_manifold(v_id))
        {
            std::cout << "- your input is NOT manifold!" << std::endl;
            std::cout << "- CHECK FAILED!" << std::endl;
            return EXIT_FAILURE;
        }

        if(m.vert_is_boundary(v_id))
        {
            std::cout << "- your input is NOT watertight!" << std::endl;
            std::cout << "- CHECK FAILED!" << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::cout << "- your input is manifold" << std::endl;
    std::cout << "- your input is watertight" << std::endl;

    // CHECK triangle ORIENTATION
    uint min_v = 0;
    double min_v_x = m.vert(0).x();

    for(uint v_id = 0; v_id < m.num_verts(); v_id++)
        if(m.vert(v_id).x() < min_v_x)
        {
            min_v_x = m.vert(v_id).x();
            min_v = v_id;
        }

    int highest_edge = -1;
    double max_height = 0.0;

    for(uint e_id : m.adj_v2e(min_v))
    {
        double edge_height = std::abs(m.edge_vert(e_id, 0).y() - m.edge_vert(e_id, 1).y());

        if(edge_height > max_height)
        {
            max_height = edge_height;
            highest_edge = (int)e_id;
        }
    }

    if(highest_edge == -1)
    {
        std::cout << "this should not be happening" << std::endl;
        return EXIT_FAILURE;
    }

    uint other_endpoint = m.vert_opposite_to(highest_edge, min_v);

    uint t_id0 = m.adj_e2p(highest_edge)[0];
    uint opp_vert0 = m.vert_opposite_to(t_id0, min_v, other_endpoint);
    double v0_z = m.vert(opp_vert0).z();

    uint t_id1 = m.adj_e2p(highest_edge)[1];
    uint opp_vert1 = m.vert_opposite_to(t_id1, min_v, other_endpoint);
    double v1_z = m.vert(opp_vert1).z();

    double min_v_z = m.vert(min_v).z();
    int seed_t = -1;

    if((v0_z < min_v_z && v1_z > min_v_z) || (v0_z > min_v_z && v1_z < min_v_z))
        seed_t = (int)t_id0;
    else if((v0_z < min_v_z && v1_z == min_v_z) || (v0_z > min_v_z && v1_z == min_v_z))
        seed_t = (int)t_id0;
    else if((v0_z == min_v_z && v1_z < min_v_z) || (v0_z == min_v_z && v1_z > min_v_z))
        seed_t = (int)t_id1;
    else if((v0_z > min_v_z && v1_z > min_v_z) || (v0_z < min_v_z && v1_z < min_v_z))
    {
        if(m.vert(opp_vert0).x() < m.vert(opp_vert1).x())
            seed_t = (int)t_id0;
        else
            seed_t = (int)t_id1;
    }

    if(seed_t == -1)
    {
        std::cout << "this should not be happening" << std::endl;
        return EXIT_FAILURE;
    }

    // check orientation of the first triangle
    cinolib::vec3d external_v(m.vert(min_v).x() -0.5, m.vert(min_v).y(), m.vert(min_v).z());
    double orient = cinolib::orient3d(m.poly_vert(seed_t, 0), m.poly_vert(seed_t, 1), m.poly_vert(seed_t, 2), external_v);

    if(orient == 0)
    {
        std::cout << "this should not be happening" << std::endl;
        return EXIT_FAILURE;
    }

    if(orient > 0)
    {
        std::cout << "- your input is NOT well oriented!" << std::endl;
        std::cout << "- CHECK FAILED!" << std::endl;
        return EXIT_FAILURE;
    }

    // flooding and triangle winding check
    std::vector<bool> visited_t(m.num_polys(), false);
    std::stack<uint> tris_stack;
    tris_stack.push(seed_t);

    while(!tris_stack.empty())
    {
        uint curr_t = tris_stack.top();
        tris_stack.pop();
        visited_t[curr_t] = true;

        for(uint adj_t : m.adj_p2p(curr_t))
        {
            if(!visited_t[adj_t])
            {
                uint shared_edge = m.edge_shared(curr_t, adj_t);
                if(m.edge_is_CCW(shared_edge, curr_t) != m.edge_is_CCW(shared_edge, adj_t))
                {
                    tris_stack.push(adj_t);
                }
                else // stop flooding
                {
                    std::cout << "- your input is NOT well oriented!" << std::endl;
                    std::cout << "- CHECK FAILED!" << std::endl;
                    return EXIT_FAILURE;
                }
            }
        }
    }

    std::cout << "- your input is well oriented" << std::endl;

    // CHECK self INTERSECTIONS
    std::set<cinolib::ipair> intersections;
    cinolib::find_intersections(m, intersections);

    if(!intersections.empty())
    {
        std::cout << "- your input is NOT self-intersections free!" << std::endl;
        std::cout << "- CHECK FAILED!" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "- your input is self-intersections free" << std::endl;

    std::cout << "- CHECK PASSED!" << std::endl;
    return EXIT_SUCCESS;
}
