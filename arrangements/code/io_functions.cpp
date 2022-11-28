/*****************************************************************************************
 *              MIT License                                                              *
 *                                                                                       *
 * Copyright (c) 2020 Gianmarco Cherchi, Marco Livesu, Riccardo Scateni e Marco Attene   *
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
 *      Marco Livesu (marco.livesu@ge.imati.cnr.it)                                      *
 *      http://pers.ge.imati.cnr.it/livesu/                                              *
 *                                                                                       *
 *      Riccardo Scateni (riccardo@unica.it)                                             *
 *      https://people.unica.it/riccardoscateni/                                         *
 *                                                                                       *
 *      Marco Attene (marco.attene@ge.imati.cnr.it)                                      *
 *      https://www.cnr.it/en/people/marco.attene/                                       *
 *                                                                                       *
 * ***************************************************************************************/

#include "io_functions.h"

#include <cinolib/octree.h>

inline void load(const std::string &filename, std::vector<double> &coords, std::vector<uint> &tris)
{
    std::vector<cinolib::vec3d> tmp_verts;

    std::string filetype = filename.substr(filename.size() - 4, 4);

    if (filetype.compare(".off") == 0 || filetype.compare(".OFF") == 0)
    {
        std::vector< std::vector<uint> > tmp_tris;
        cinolib::read_OFF(filename.c_str(), tmp_verts, tmp_tris);
        tris = cinolib::serialized_vids_from_polys(tmp_tris);
    }
    else if (filetype.compare(".obj") == 0 || filetype.compare(".OBJ") == 0)
    {
        std::vector< std::vector<uint> > tmp_tris;
        cinolib::read_OBJ(filename.c_str(), tmp_verts, tmp_tris);
        tris = cinolib::serialized_vids_from_polys(tmp_tris);
    }
    else if (filetype.compare(".stl") == 0 || filetype.compare(".STL") == 0)
    {
        cinolib::read_STL(filename.c_str(), tmp_verts, tris, false);
    }
    else
    {
        std::cerr << "ERROR: file format not supported yet " << std::endl;
    }

    coords = cinolib::serialized_xyz_from_vec3d(tmp_verts);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void loadMultipleFiles(const std::vector<std::string> &files, std::vector<double> &coords, std::vector<uint> &tris, std::vector<uint> &labels)
{
    for(uint f_id = 0; f_id < files.size(); f_id++)
    {
        std::vector<double> tmp_coords;
        std::vector<uint> tmp_tris;

        load(files[f_id], tmp_coords, tmp_tris);

        uint off = static_cast<uint>(coords.size() / 3); // prev num verts

        coords.insert(coords.end(), tmp_coords.begin(), tmp_coords.end());

        for(auto &i : tmp_tris) tris.push_back(i + off);

        for(uint i = 0; i < tmp_tris.size() / 3; i++)
            labels.push_back(f_id);
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void loadMultipleFilesWithVertFix(const std::vector<std::string> &files, std::vector<double> &coords, std::vector<uint> &tris, std::vector<uint> &labels)
{
    for(uint f_id = 0; f_id < files.size(); f_id++)
    {
        cinolib::Trimesh<> tm(files[f_id].c_str());

        bool fix = fixCoincidentVertices(tm);

        auto end_f = chrono::steady_clock::now();

        uint off = static_cast<uint>(coords.size() / 3); // prev num verts

        for(cinolib::vec3d &v : tm.vector_verts())
        {
            coords.push_back(v.x());
            coords.push_back(v.y());
            coords.push_back(v.z());
        }

        for(auto &t : tm.vector_polys())
        {
            tris.push_back(t[0] + off);
            tris.push_back(t[1] + off);
            tris.push_back(t[2] + off);
        }

        for(uint t_id = 0; t_id < tm.num_polys(); t_id++)
            labels.push_back(f_id);
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline bool fixCoincidentVertices(cinolib::Trimesh<> &m) // return true if the fix is applied
{
    bool res = false;
    m.vert_set_flag(cinolib::MARKED, false);

    cinolib::Octree o;
    double r = m.edge_avg_length() * 0.01;

    for(uint v_id = 0; v_id < m.num_verts(); ++v_id)
        o.push_sphere(v_id, m.vert(v_id), r);

    o.build();

    for(uint v_id = 0; v_id < m.num_verts(); ++v_id)
    {
        cinolib::vec3d p = m.vert(v_id);
        std::unordered_set<uint> ids;
        o.contains(p, false, ids);
        if(ids.size() > 1)
        {
            for(uint id : ids)
                m.vert_data(id).flags[cinolib::MARKED] = true;
        }
    }

    for(uint v_id = 0; v_id < m.num_verts(); ++v_id)
    {
        if(m.vert_data(v_id).flags[cinolib::MARKED])
        {
            m.vert(v_id) += m.vert_data(v_id).normal * r * 2.0;
            res = true;
        }
    }

    return res;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void save(const std::string &filename, std::vector<double> &coords, std::vector<uint> &tris)
{
    std::string filetype = filename.substr(filename.size() - 4, 4);

    std::vector< std::vector<uint> > tmp_tris = cinolib::polys_from_serialized_vids(tris, 3);

    if (filetype.compare(".off") == 0 || filetype.compare(".OFF") == 0)
    {
        cinolib::write_OFF(filename.c_str(), coords, tmp_tris);
    }
    else if (filetype.compare(".obj") == 0 || filetype.compare(".OBJ") == 0)
    {
        cinolib::write_OBJ(filename.c_str(), coords, tmp_tris);
    }
    else
    {
        std::cerr << "ERROR: file format not supported yet " << std::endl;
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void writeIMPL(const string &filename, const std::vector<genericPoint *> &verts, const std::vector<uint> &tris, const std::vector<std::bitset<NBIT> > &labels)
{
    setlocale(LC_NUMERIC, "en_US.UTF-8"); // makes sure "." is the decimal separator

    FILE *fp = fopen(filename.c_str(), "w");

    if(!fp)
    {
        std::cerr << "ERROR : " << __FILE__ << ", line " << __LINE__ << " : writeIMPL() : couldn't open input file " << filename << std::endl;
        exit(-1);
    }

    fprintf(fp, "nv %d\n", static_cast<uint>(verts.size()));    // num verts
    fprintf(fp, "nt %d\n", static_cast<uint>(tris.size() /3)); // num triangles

    std::map<std::tuple<double, double, double>, uint> v_map;

    // explicit points map creation
    for(uint v_id = 0; v_id < verts.size(); v_id++)
    {
        if(verts[v_id]->isExplicit3D())
        {
            const explicitPoint3D &ev = verts[v_id]->toExplicit3D();
            std::tuple<double,double,double> v(ev.X(), ev.Y(), ev.Z());
            v_map[v] = v_id;
        }
    }

    for(uint v_id = 0; v_id < verts.size(); v_id++)
    {
        if(verts[v_id]->isExplicit3D())
        {
            const explicitPoint3D &ev = verts[v_id]->toExplicit3D();
            fprintf(fp, "e %.17g %.17g %.17g\n", ev.X(), ev.Y(), ev.Z());
        }
        else if(verts[v_id]->isLPI())
        {
            const explicitPoint3D &ep = verts[v_id]->toLPI().P();   std::tuple<double,double,double> tp(ep.X(), ep.Y(), ep.Z());
            const explicitPoint3D &eq = verts[v_id]->toLPI().Q();   std::tuple<double,double,double> tq(eq.X(), eq.Y(), eq.Z());
            const explicitPoint3D &er = verts[v_id]->toLPI().R();   std::tuple<double,double,double> tr(er.X(), er.Y(), er.Z());
            const explicitPoint3D &es = verts[v_id]->toLPI().S();   std::tuple<double,double,double> ts(es.X(), es.Y(), es.Z());
            const explicitPoint3D &et = verts[v_id]->toLPI().T();   std::tuple<double,double,double> tt(et.X(), et.Y(), et.Z());

            uint id_p = v_map.find(tp)->second;     uint id_q = v_map.find(tq)->second;
            uint id_r = v_map.find(tr)->second;     uint id_s = v_map.find(ts)->second;     uint id_t = v_map.find(tt)->second;

            fprintf(fp, "l %d %d %d %d %d\n", id_p, id_q, id_r, id_s, id_t);
        }
        else if(verts[v_id]->isTPI())
        {
            const explicitPoint3D &eu1 = verts[v_id]->toTPI().U1();   std::tuple<double,double,double> tu1(eu1.X(), eu1.Y(), eu1.Z());
            const explicitPoint3D &eu2 = verts[v_id]->toTPI().U2();   std::tuple<double,double,double> tu2(eu2.X(), eu2.Y(), eu2.Z());
            const explicitPoint3D &eu3 = verts[v_id]->toTPI().U3();   std::tuple<double,double,double> tu3(eu3.X(), eu3.Y(), eu3.Z());
            const explicitPoint3D &ev1 = verts[v_id]->toTPI().V1();   std::tuple<double,double,double> tv1(ev1.X(), ev1.Y(), ev1.Z());
            const explicitPoint3D &ev2 = verts[v_id]->toTPI().V2();   std::tuple<double,double,double> tv2(ev2.X(), ev2.Y(), ev2.Z());
            const explicitPoint3D &ev3 = verts[v_id]->toTPI().V3();   std::tuple<double,double,double> tv3(ev3.X(), ev3.Y(), ev3.Z());
            const explicitPoint3D &ew1 = verts[v_id]->toTPI().W1();   std::tuple<double,double,double> tw1(ew1.X(), ew1.Y(), ew1.Z());
            const explicitPoint3D &ew2 = verts[v_id]->toTPI().W2();   std::tuple<double,double,double> tw2(ew2.X(), ew2.Y(), ew2.Z());
            const explicitPoint3D &ew3 = verts[v_id]->toTPI().W3();   std::tuple<double,double,double> tw3(ew3.X(), ew3.Y(), ew3.Z());

            uint id_u1 = v_map.find(tu1)->second;   uint id_u2 = v_map.find(tu2)->second;   uint id_u3 = v_map.find(tu3)->second;
            uint id_v1 = v_map.find(tv1)->second;   uint id_v2 = v_map.find(tv2)->second;   uint id_v3 = v_map.find(tv3)->second;
            uint id_w1 = v_map.find(tw1)->second;   uint id_w2 = v_map.find(tw2)->second;   uint id_w3 = v_map.find(tw3)->second;

            fprintf(fp, "t %d %d %d %d %d %d %d %d %d\n", id_u1, id_u2, id_u3, id_v1, id_v2, id_v3, id_w1, id_w2, id_w3);
        }
    }

    for(uint t_id = 0; t_id < tris.size()/3; t_id++)
    {
        fprintf(fp, "f %d %d %d %lu\n", tris[3 * t_id], tris[3 * t_id +1], tris[3 * t_id +2], labels[t_id].to_ulong());
    }

    fclose(fp);
}

inline void readIMPL(const std::string &filename, std::vector<genericPoint*> &verts, std::vector<uint> &tris, std::vector<std::bitset<NBIT>> &labels)
{
    std::ifstream fp(filename);
    if(!fp.is_open())
    {
        std::cerr << "ERROR : " << __FILE__ << ", line " << __LINE__ << " : readIMPL() : couldn't open input file " << filename << std::endl;
        exit(-1);
    }

    uint curr_v_id = 0;
    uint curr_t_id = 0;

    std::vector<std::vector<uint>> impl_points; // v[0] -> pos vert, v[1]..v[n] verts id

    std::string line;
    while(std::getline(fp, line))
    {
        switch(line[0])
        {
            case 'e': // Explicit point
            {
                double x, y, z;
                sscanf(line.data(), "e %lf %lf %lf", &x, &y, &z);
                verts[curr_v_id++] = new explicitPoint3D(x, y, z);
            } break;

            case 'l': // LPI point
            {
                uint p, q, r, s, t;
                sscanf(line.data(), "l %d %d %d %d %d", &p, &q, &r, &s, &t);
                impl_points.push_back({curr_v_id++, p, q, r, s, t});
            } break;

            case 't': // TPI point
            {
                uint u1, u2, u3, v1, v2, v3, w1, w2, w3;
                sscanf(line.data(), "t %d %d %d %d %d %d %d %d %d", &u1, &u2, &u3, &v1, &v2, &v3, &w1, &w2, &w3);
                impl_points.push_back({curr_v_id++, u1, u2, u3, v1, v2, v3, w1, w2, w3});
            } break;

            case 'f': // Triangle
            {
                uint v0, v1, v2;
                unsigned long int l;
                sscanf(line.data(), "f %d %d %d %lu", &v0, &v1, &v2, &l);
                tris[3 * curr_t_id] = v0;
                tris[3 * curr_t_id +1] = v1;
                tris[3 * curr_t_id +2] = v2;
                labels[curr_t_id] = static_cast<unsigned long>(l);
                curr_t_id++;
            } break;

            case 'n': // Preliminary info
            {
                uint size;
                if(line[1] == 'v')
                {
                    sscanf(line.data(), "nv %d", &size);
                    verts.clear();
                    verts.resize(size);
                }
                else if(line[1] == 't')
                {
                    sscanf(line.data(), "nt %d", &size);
                    tris.clear();
                    labels.clear();
                    tris.resize(3 * size);
                    labels.resize(3 * size);
                }
            } break;
        }
    }

    fp.close();

    // parse implicit_points
    for(auto &point : impl_points)
    {
        if(point.size() == 6) // LPI
        {
            verts[point[0]] = new implicitPoint3D_LPI(verts[point[1]]->toExplicit3D(), verts[point[2]]->toExplicit3D(),
                                                      verts[point[3]]->toExplicit3D(), verts[point[4]]->toExplicit3D(), verts[point[5]]->toExplicit3D());
        }
        else // TPI
        {
            verts[point[0]] = new implicitPoint3D_TPI(verts[point[1]]->toExplicit3D(), verts[point[2]]->toExplicit3D(), verts[point[3]]->toExplicit3D(),
                                                      verts[point[4]]->toExplicit3D(), verts[point[5]]->toExplicit3D(), verts[point[6]]->toExplicit3D(),
                                                      verts[point[7]]->toExplicit3D(), verts[point[8]]->toExplicit3D(), verts[point[9]]->toExplicit3D());
        }
    }
}
