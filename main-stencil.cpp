#include <cinolib/gl/glcanvas.h>
#include <cinolib/drawable_triangle_soup.h>
#include <cinolib/meshes/drawable_trimesh.h>
#include <thread>
#include "booleans.h"

int main(int argc, char **argv)
{
    std::vector<std::string> in_files;
    in_files.emplace_back("../data/fertility.obj"); // base model
    std::string stencil_path = "/Users/gianmarco/Code/InteractiveAndRobustMeshBooleans/data/spheres/";

    std::vector<double> in_coords, bool_coords;
    std::vector<uint> in_tris;
    std::vector<uint> in_labels;
    BoolOp op = SUBTRACTION;

    int num_sub = 30;

    if(NBIT < num_sub)
    {
        std::cerr << "Set NBIT to a higher number!\n";
        return 0;
    }

    for(int i = 0; i < num_sub; i++)
        in_files.push_back(stencil_path + std::to_string(i) + ".obj");

    std::vector<uint> bool_tris;
    std::vector<std::bitset<NBIT>> bool_labels;

    loadMultipleFiles(in_files, in_coords, in_tris, in_labels);

    booleanPipeline(in_coords, in_tris, in_labels, op, bool_coords, bool_tris, bool_labels);

    cinolib::write_OBJ("output.obj", bool_coords, bool_tris, {});

    return 0;
}