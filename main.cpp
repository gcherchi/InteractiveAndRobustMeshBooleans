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

#include "booleans.h"
#include "filesystem"
#include "cinolib/meshes/trimesh.h"

std::vector<std::string> files;
namespace fs = std::filesystem;
bool test = true;
bool debug = true;
bool test_multiple = true;

std::string replaceKeyword(const std::string& input, const std::string& from, const std::string& to) {
    std::string result = input;
    size_t pos = result.find(from);
    if (pos != std::string::npos) {
        result.replace(pos, from.length(), to);
    }
    return result;
}

int main(int argc, char **argv)
{
    BoolOp op;
    std::string file_out;

    if(debug) {
        std::cout << "Debug mode enabled" << std::endl;
        op = UNION;
        //mesh non manifold individuato ma non risolto
        //files.emplace_back("../modelli_filtrati/Tinghi10K/39929_sf_a.obj");
        //files.emplace_back("../modelli_filtrati/mesh_rotated/39929_sf_a.obj");
        files.emplace_back("../modelli_filtrati/Tinghi10K/47568_sf_a.obj");
        files.emplace_back("../modelli_filtrati/mesh_rotated/47568_sf_a.obj");
        //files.emplace_back("../mostro5.obj");
        //files.emplace_back("../mostro4.obj");

        file_out = "output_union.obj";
    }
    if(!debug){
        if(argc < 5)
        {
            std::cout << "syntax error!" << std::endl;
            std::cout << "./exact_boolean BOOL_OPERATION (intersection OR union OR subtraction) input1.obj input2.obj output.obj" << std::endl;
            return -1;
        }
        else
        {
            if (strcmp(argv[1], "intersection") == 0)       op = INTERSECTION;
            else if (strcmp(argv[1], "union") == 0)         op = UNION;
            else if (strcmp(argv[1], "subtraction") == 0)   op = SUBTRACTION;
            else if (strcmp(argv[1], "xor") == 0)           op = XOR;
        }
        for(int i = 2; i < (argc -1); i++)
            files.emplace_back(argv[i]);

        file_out = argv[argc-1];
    }


    std::vector<double> in_coords, bool_coords;
    std::vector<uint> in_tris, bool_tris;
    std::vector<uint> in_labels;
    std::vector<std::bitset<NBIT>> bool_labels;

    loadMultipleFiles(files, in_coords, in_tris, in_labels);
    Data data;

    data.test_multiple = test_multiple;
    data.num_poly_input = in_tris.size() / 3;
    data.num_vert_input = in_coords.size() / 3;
    if (debug) {
        cinolib::write_OBJ("mesh_input.obj", in_coords, in_tris, {});

        //data.t_id_debug = 3794;
    }

    cinolib::Profiler p;
    p.push("Boolean Pipeline time: ");
    booleanPipeline(in_coords, in_tris, in_labels, op, bool_coords, bool_tris, bool_labels, data);
    p.pop();
    std::cout << "Dimension input: n poly - " << data.num_poly_input << std::endl;
    std::cout << "Dimension arrangement: n poly - " << data.num_poly_arrang << std::endl;


    if(debug) {
        //create a parser of the previous file
        saveTriangleIDsToFile(data);
        Data data_tmp;
        loadTriangleIDsFromFile(data_tmp);
    }

    if(test_multiple){
        fs::path script_dir = fs::absolute(fs::path(argv[0])).parent_path();
        fs::path exePath = fs::absolute(script_dir / "mesh_booleans_inputcheck");
        const std::string exe = exePath.string();

        //union
        cinolib::write_OBJ(file_out.c_str(), data.bool_coords_union, data.bool_tris_union, {});
        std::string command_union = std::string(exe) + " " + file_out.c_str() + " union";


        std::string intersection_output = replaceKeyword(file_out, "union", "intersection");
        cinolib::write_OBJ(intersection_output.c_str(), data.bool_coords_intersection, data.bool_tris_intersection, {});
        std::string command_intersection = std::string(exe) + " " + intersection_output.c_str() + " intersection";

        std::string subtraction_output = replaceKeyword(file_out, "union", "subtraction");
        cinolib::write_OBJ(subtraction_output.c_str(), data.bool_coords_subtraction, data.bool_tris_subtraction, {});
        std::string command_subtraction = std::string(exe) + " " + subtraction_output.c_str() + " subtraction";

        int result = system(command_union.c_str());
        if (result != 0) {
            std::cerr << "Error in the execution of the command - union" << std::endl;
        }

        result = system(command_intersection.c_str());
        if (result != 0) {
            std::cerr << "Error in the execution of the command - intersection" << std::endl;
        }

        result = system(command_subtraction.c_str());
        if (result != 0) {
            std::cerr << "Error in the execution of the command - subtraction" << std::endl;
        }

    }
    cinolib::write_OBJ(file_out.c_str(), bool_coords, bool_tris, {});


    if(!test_multiple && test) {
        fs::path script_dir = fs::absolute(fs::path(argv[0])).parent_path();
        fs::path exePath = fs::absolute(script_dir / "mesh_booleans_inputcheck");
        const std::string exe = exePath.string();

        // Costruzione del comando da eseguire
        std::string command = std::string(exe) + " " + file_out.c_str();
        // Esecuzione del comando
        std::cout<<"mesh_name: " << file_out.c_str() << std::endl;
        int result = system(command.c_str());
        if (result != 0) {
            std::cerr << "Error in the execution of the command" << std::endl;
        }
    }
    return 0;
}