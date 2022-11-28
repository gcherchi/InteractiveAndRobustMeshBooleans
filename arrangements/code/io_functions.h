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

#ifndef IO_FUNCTIONS_H
#define IO_FUNCTIONS_H

#include<implicit_point.h>

#include <cinolib/vector_serialization.h>

#include <cinolib/meshes/trimesh.h>

#include <cinolib/io/read_OFF.h>
#include <cinolib/io/read_OBJ.h>
#include <cinolib/io/read_STL.h>

#include <cinolib/io/write_OFF.h>
#include <cinolib/io/write_OBJ.h>
#include <cinolib/io/write_STL.h>

#include <common.h>


inline void load(const std::string &filename, std::vector<double> &coords, std::vector<uint> &tris);

inline void loadMultipleFiles(const std::vector<std::string> &files, std::vector<double> &coords, std::vector<uint> &tris, std::vector<uint> &labels);

inline void loadMultipleFiles(const std::vector<std::string> &files, std::vector<double> &coords, std::vector<uint> &tris, std::vector<uint> &labels, int &vert_offset);

inline void loadMultipleFilesWithVertFix(const std::vector<std::string> &files, std::vector<double> &coords, std::vector<uint> &tris, std::vector<uint> &labels);

inline bool fixCoincidentVertices(cinolib::Trimesh<> &m);

inline void save(const std::string &filename, std::vector<double> &coords, std::vector<uint> &tris);

inline void writeIMPL(const std::string &filename, const std::vector<genericPoint*> &verts, const std::vector<uint> &tris, const std::vector<std::bitset<NBIT>> &labels);

inline void readIMPL(const std::string &filename, std::vector<genericPoint*> &verts, std::vector<uint> &tris, std::vector<std::bitset<NBIT>> &labels);

#include "io_functions.cpp"

#endif // IO_FUNCTIONS_H
