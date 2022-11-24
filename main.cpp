#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include "booleans.h"
//#include "interactive_rot.h"
//#include "interactive_arap.h"
//#include "stencil_demo.h"

int  vert_offset = 0; // at what index do the coordinates of the second mesh start?

std::vector<std::string> files;
uint int_info;


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

/*

void print_stats(bool verbose) {
    std::cout << "boolean operations complete" << std::endl;
    std::cout << "arr time: " << arr_time << " ms"<< std::endl;
    std::cout << "bool time: " << bool_time << " ms"<< std::endl;
    std::cout << "TOT time: " << arr_time + bool_time << " ms"<< std::endl;

    if(!verbose) return;

    std::cout << std::endl << "arr times ----------" << std::endl;
    for(auto time : arr_times) {
        std::cout << "arr times: " << time.elapsed << " ms " << time.percent << "% " << time.title << std::endl;
    }
    std::cout << std::endl << "bool times ----------" << std::endl;
    for(auto time : bool_times) {
        std::cout << "bool times: " << time.elapsed << " ms " << time.percent << "% " << time.title << std::endl;
    }
}

void translateMesh(std::vector<double>& in_coords, const std::array<double, 3>& v = {0.1,0,0}) {
    for(auto idx = 0; idx < in_coords.size(); idx += 3) {
        in_coords[idx + 0] += v[0];
        in_coords[idx + 1] += v[1];
        in_coords[idx + 2] += v[2];
    }
}

int normalize_meshes(int argc, char **argv) {
    for(int idx = 2; idx < argc; idx++) {
        string filename(argv[idx]);
        std::cout << "from " << filename << "\n";
        cinolib::Trimesh mesh(filename.c_str());
        mesh.center_bbox();
        mesh.normalize_bbox();
        string outname = filename.substr(0, filename.size()-4) + "-normalized.obj";
        std::cout << "to " << outname << "\n";
        mesh.save(outname.c_str());
    }
    return 0;
}

int align_meshes(int argc, char **argv) {
    for(int idx = 2; idx < argc; idx++) {
        string filename(argv[idx]);
        std::cout << "from " << filename << "\n";
        cinolib::Trimesh mesh(filename.c_str());
        double min = mesh.bbox().min[1];
        for(auto idx = 0; idx < mesh.num_verts(); idx++) mesh.vert(idx)[1] -= min;
        string outname = filename.substr(0, filename.size()-4) + "-aligned.obj";
        std::cout << "to " << outname << "\n";
        mesh.save(outname.c_str());
    }
    return 0;
}

int poisson_sample(int argc, char **argv)
{
    cinolib::DrawableTrimesh<> m("./data/pumpkin.obj");
    cinolib::GLcanvas gui;
    gui.push(&m);
    gui.push(new cinolib::SurfaceMeshControls<cinolib::DrawableTrimesh<>>(&m,&gui));

    cinolib::Octree o;
    o.build_from_mesh_polys(m);

    std::vector<cinolib::vec3d> cloud;
    std::vector<std::vector<uint>> dummy;
    cinolib::read_OBJ("./data/pumpkin-samples.obj", cloud, dummy);
    std::cout << cloud.size() << std::endl;

    cinolib::DrawableTrimesh<> all_spheres;
    cinolib::DrawableTrimesh<> ref("./data/MAKIES_Jack_O_Lantern_Lid.STL");
    ref.center_bbox();
    // ref.scale(0.25); // THIS IS A BUG
    for(auto idx = 0; idx < ref.num_verts(); idx++) ref.vert(idx) *= 0.5;
    ref.update_bbox();
    double offset = - ref.bbox().min[2];
    std::vector<cinolib::DrawableTrimesh<>> spheres(cloud.size(),ref);
    for(uint i=0; i<cloud.size(); ++i)
    {
        uint pid;
        cinolib::vec3d pos;
        double dist;
        o.closest_point(cloud[i], pid, pos, dist);

        cinolib::vec3d n = m.poly_data(pid).normal;
        cinolib::vec3d Z(0,0,1);
        cinolib::vec3d axis = Z.cross(n);
        axis.normalize();
        spheres[i].rotate(axis,Z.angle_rad(n) + M_PI);
        spheres[i].translate(cloud[i] - n * offset * 0.9);
        spheres[i].updateGL();
        spheres[i].save(("./data/pumpkin-deco" + std::to_string(i) + ".obj").c_str());
        all_spheres += spheres[i];
        gui.push(&spheres[i],false);
    }
    all_spheres.save("./data/pumpkin-alldeco.obj");

    return gui.launch();
}

int main(int argc, char **argv)
{
    if(strcmp(argv[1], "normalize") == 0)
        return normalize_meshes(argc, argv);
    if(strcmp(argv[1], "sample") == 0)
        return poisson_sample(argc, argv);
    if(strcmp(argv[1], "align") == 0)
        return align_meshes(argc, argv);

    BoolOp op;
    std::string file_out;
    bool interactive = false;
    int runs = 1;
    bool verbose = false;
    bool inter_rot_demo = false;
    bool inter_ARAP_demo = false;
    bool normalize = false;
    bool stencil_demo = false;

    if(argc < 5)
    {
        //std::cout << "syntax error!" << std::endl;
        //std::cout << "./exact_boolean BOOL_OPERATION (intersection OR union OR subtraction) input1.obj input2.obj output.obj" << std::endl;
        //return -1;
    }
    else
    {
        if(strcmp(argv[1], "intersection") == 0) op = INTERSECTION;
        else if(strcmp(argv[1], "union") == 0)   op = UNION;
        else if(strcmp(argv[1], "subtraction") == 0)   op = SUBTRACTION;
        else if(strcmp(argv[1], "xor") == 0)   op = XOR;
        else if(strcmp(argv[1], "interactive-intersection") == 0) { op = INTERSECTION; interactive = true; }
        else if(strcmp(argv[1], "interactive-union") == 0) { op = UNION; interactive = true; }
        else if(strcmp(argv[1], "interactive-subtraction") == 0) { op = SUBTRACTION; interactive = true; }
        else if(strcmp(argv[1], "interactive-xor") == 0) { op = XOR; interactive = true; }
        else if(strcmp(argv[1], "verbose-intersection") == 0) { op = INTERSECTION; verbose = true; }
        else if(strcmp(argv[1], "verbose-union") == 0) { op = UNION; verbose = true; }
        else if(strcmp(argv[1], "verbose-subtraction") == 0) { op = SUBTRACTION; verbose = true; }
        else if(strcmp(argv[1], "verbose-xor") == 0) { op = XOR; verbose = true; }
        else if(strcmp(argv[1], "fast-intersection") == 0) { op = INTERSECTION; verbose = true; runs = 1; }
        else if(strcmp(argv[1], "fast-union") == 0) { op = UNION; verbose = true; runs = 1; }
        else if(strcmp(argv[1], "fast-subtraction") == 0) { op = SUBTRACTION; verbose = true; runs = 1; }
        else if(strcmp(argv[1], "fast-xor") == 0) { op = XOR; verbose = true; runs = 1; }
        else if(strcmp(argv[1], "inter-rot-demo") == 0) { inter_rot_demo = true; }
        else if(strcmp(argv[1], "inter-arap-demo") == 0) { inter_ARAP_demo = true; }
        else if(strcmp(argv[1], "stencil") == 0) {stencil_demo = true; }
        else
        {
            std::cout << "syntax error!" << std::endl;
            std::cout << "./exact_boolean BOOL_OPERATION (intersection OR union OR subtraction) input1.obj input2.obj output.obj" << std::endl;
            return -1;
        }

        for(int i = 2; i < (argc -1); i++)
            files.emplace_back(argv[i]);

        file_out = argv[argc-1];
    }

    if(inter_rot_demo)
    {
        std::vector<double> in_coords;
        std::vector<uint> in_tris;
        std::vector<uint> in_labels;
        loadMultipleFiles(files, in_coords, in_tris, in_labels);
        return interactive_rot(in_coords, in_tris, in_labels);
    }

    if(inter_ARAP_demo)
    {
        std::vector<double> in_coords;
        std::vector<uint> in_tris;
        std::vector<uint> in_labels;
        files = std::vector<std::string>({"../data/armadillo.obj","../data/armadillo.obj"});
        loadMultipleFiles(files, in_coords, in_tris, in_labels);
        cinolib::DrawableTrimesh<> m1(files.begin()->c_str());
        cinolib::DrawableTrimesh<> m2 = m1;
        m2.translate(cinolib::vec3d(m1.bbox().delta_y()*0.8,0,0));
        m2.update_bbox();
        m2.updateGL();
        return interactive_arap(in_coords, in_tris, in_labels, m1, m2);
    }

    if(stencil_demo)
    {
        std::vector<std::string> in_files;
        std::string out_file;

        std::string stencil_path = argv[2];
        in_files.emplace_back(stencil_path + "/base.obj");
        uint num_sub = stoi(argv[3]);
        out_file = argv[4];

        if(NBIT < num_sub)
        {
            std::cerr << "Set NBIT to a higher number!\n";
            return 0;
        }

        for(int i = 0; i < num_sub; i++)
            in_files.push_back(stencil_path + "/rock" + std::to_string(i) + ".obj");

        return stencil(in_files, out_file);
    }

    std::vector<int> more_arr_time;
    std::vector<int> more_bool_time;

    for(uint i = 0; i < runs; i++)
    {
        arr_time = 0;
        bool_time = 0;

        std::vector<double> in_coords, bool_coords;
        std::vector<uint> in_tris, bool_tris;
        std::vector<uint> in_labels;
        std::vector<std::bitset<NBIT>> bool_labels;

        loadMultipleFiles(files, in_coords, in_tris, in_labels);

        //////////// BOOLEANS PIPELINE START //////////////////////////////////////////////////////////////////////////////

        booleanPipeline(in_coords, in_tris, in_labels, op, bool_coords, bool_tris, bool_labels);

        //saveStatisticsOnFile(files[0], files[1], arr_time, bool_time, int_info);
        //saveStatisticsOnFile(files[0], files[1], arr_time, bool_time, int_info, arr_times, bool_times);

        /////////// BOOLEANS PIPELINE END //////////////////////////////////////////////////////////////////////////////

        if(i==runs-1)
        {
            std::vector<cinolib::Color> colors(bool_tris.size()/3);
            for(uint i=0; i<colors.size(); ++i)
            {
                colors[i] = (bool_labels[i][1]) ? cinolib::Color::PASTEL_ORANGE() : cinolib::Color::PASTEL_CYAN();
            }
            dump_colored_mesh(bool_coords, bool_tris,colors);
            cinolib::write_OBJ(file_out.c_str(), bool_coords, bool_tris, {});
        }

        more_arr_time.push_back(arr_time);
        more_bool_time.push_back(bool_time);

        if(verbose) {

            std::cout << "arr  time: " << std::setw(4) << (int)std::round(arr_time / 1000.0) << " ms " << std::setw(8) << arr_time << " us " << std::endl;
            std::cout << "bool time: " << std::setw(4) << (int)std::round(bool_time / 1000.0) << " ms " << std::setw(8) << bool_time << " us " << std::endl;
            std::cout << "\n";
            for(auto& time : arr_times) {
                std::cout << "arr  time: " << std::setw(4) << (int)std::round(time.elapsed / 1000.0) << " ms " << std::setw(8) << time.elapsed << " us " << std::setw(2) << time.percent << "% " << time.title << std::endl;
            }
            for(auto& time : bool_times) {
                std::cout << "bool time: " << std::setw(4) << (int)std::round(time.elapsed / 1000.0) << " ms " << std::setw(8) << time.elapsed << " us " << std::setw(2) << time.percent << "% " << time.title << std::endl;
            }
            std::cout << "\n";
        }
    }

    int best_arr_time = *std::min_element(more_arr_time.begin(), more_arr_time.end());
    int best_bool_time = *std::min_element(more_bool_time.begin(), more_bool_time.end());

    saveStatisticsOnFile(files[0], files[1], best_arr_time, best_bool_time, 0);

    std::cout << "best arr  time: " << std::setw(4) << (int)std::round(best_arr_time / 1000.0) << " ms " << std::setw(8) << arr_time << " us " << std::endl;
    std::cout << "best bool time: " << std::setw(4) << (int)std::round(best_bool_time / 1000.0) << " ms " << std::setw(8) << bool_time << " us " << std::endl;
    std::cout << "best TOT  time: " << std::setw(4) << (int)std::round((best_bool_time+best_arr_time)/1000.0) << " ms " << std::endl;



    return 0;
}

 */
