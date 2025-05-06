//
// Created by Michele on 06/05/25.
//
#include <iostream>

#include <cinolib/meshes/meshes.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>

using namespace cinolib;
using namespace std;

int main(int argc, char **argv){

    GLcanvas gui;
    DrawableTrimesh<> bool_mesh ;
    SurfaceMeshControls<DrawableTrimesh<>> contr (&bool_mesh, &gui, "bool_mesh");

    gui.push(&bool_mesh);
    gui.push(&contr);
    bool_mesh.updateGL();

    return gui.launch();
}
