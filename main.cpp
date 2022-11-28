#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include "booleans.h"

std::vector<std::string> files;

int main(int argc, char **argv)
{
    BoolOp op;
    std::string file_out;

    if(argc < 5)
    {
        std::cout << "syntax error!" << std::endl;
        std::cout << "./exact_boolean BOOL_OPERATION (intersection OR union OR subtraction) input1.obj input2.obj output.obj" << std::endl;
        return -1;
    }
    else
    {
        if (strcmp(argv[1], "intersection") == 0) op = INTERSECTION;
        else if (strcmp(argv[1], "union") == 0) op = UNION;
        else if (strcmp(argv[1], "subtraction") == 0) op = SUBTRACTION;
    }

    for(int i = 2; i < (argc -1); i++)
        files.emplace_back(argv[i]);

    file_out = argv[argc-1];

    std::vector<double> in_coords, bool_coords;
    std::vector<uint> in_tris, bool_tris;
    std::vector<uint> in_labels;
    std::vector<std::bitset<NBIT>> bool_labels;

    loadMultipleFiles(files, in_coords, in_tris, in_labels);

    booleanPipeline(in_coords, in_tris, in_labels, op, bool_coords, bool_tris, bool_labels);

    cinolib::write_OBJ(file_out.c_str(), bool_coords, bool_tris, {});

    return 0;
}