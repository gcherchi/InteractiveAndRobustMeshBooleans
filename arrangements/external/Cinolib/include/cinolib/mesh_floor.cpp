/********************************************************************************
*  This file is part of CinoLib                                                 *
*  Copyright(C) 2022: Marco Livesu                                              *
*                                                                               *
*  The MIT License                                                              *
*                                                                               *
*  Permission is hereby granted, free of charge, to any person obtaining a      *
*  copy of this software and associated documentation files (the "Software"),   *
*  to deal in the Software without restriction, including without limitation    *
*  the rights to use, copy, modify, merge, publish, distribute, sublicense,     *
*  and/or sell copies of the Software, and to permit persons to whom the        *
*  Software is furnished to do so, subject to the following conditions:         *
*                                                                               *
*  The above copyright notice and this permission notice shall be included in   *
*  all copies or substantial portions of the Software.                          *
*                                                                               *
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR   *
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
*  FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE *
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      *
*  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS *
*  IN THE SOFTWARE.                                                             *
*                                                                               *
*  Author(s):                                                                   *
*                                                                               *
*     Marco Livesu (marco.livesu@gmail.com)                                     *
*     http://pers.ge.imati.cnr.it/livesu/                                       *
*                                                                               *
*     Italian National Research Council (CNR)                                   *
*     Institute for Applied Mathematics and Information Technologies (IMATI)    *
*     Via de Marini, 6                                                          *
*     16149 Genoa,                                                              *
*     Italy                                                                     *
*********************************************************************************/
#include <cinolib/mesh_floor.h>

namespace cinolib
{

template<class M, class V, class E, class P>
CINO_INLINE
void mesh_floor(const AbstractMesh<M,V,E,P> & m,
                      std::vector<vec3d>    & verts,
                      std::vector<uint>     & tris,
                const int                     plane)
{
    AABB bb  = m.bbox();
    auto tmp = bb.corners();

    switch(plane)
    {
        case XY_MIN: verts = std::vector<vec3d>({ tmp[0], tmp[1], tmp[2], tmp[3] }); break;
        case XZ_MIN: verts = std::vector<vec3d>({ tmp[0], tmp[1], tmp[5], tmp[4] }); break;
        case YZ_MIN: verts = std::vector<vec3d>({ tmp[0], tmp[3], tmp[7], tmp[4] }); break;
        case XY_MAX: verts = std::vector<vec3d>({ tmp[6], tmp[7], tmp[4], tmp[5] }); break;
        case XZ_MAX: verts = std::vector<vec3d>({ tmp[6], tmp[7], tmp[3], tmp[2] }); break;
        case YZ_MAX: verts = std::vector<vec3d>({ tmp[6], tmp[5], tmp[1], tmp[2] }); break;
        default : assert(false && "Unknown floor plane!");
    }

    tris  = std::vector<uint>({ 0, 1, 2, 0, 2, 3 });

    // re-center and scale floor to make it big enough
    vec3d c = (verts[0] + verts[1] + verts[2] + verts[3])*0.25;
    for(auto & v : verts) v -= c;
    for(auto & v : verts) v *= 6;
    for(auto & v : verts) v += c;
}

}
