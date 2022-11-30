#include <cinolib/gl/glcanvas.h>
#include <cinolib/drawable_triangle_soup.h>
#include <cinolib/meshes/drawable_trimesh.h>
#include <cinolib/ARAP.h>
#include <thread>
#include "booleans.h"

std::vector<uint> handles;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

uint define_handle(const cinolib::Trimesh<> & m1, const cinolib::Trimesh<> & m2, const cinolib::vec3d & click_3d)
{
    uint v0 = m1.pick_vert(click_3d);
    uint v1 = m2.pick_vert(click_3d);
    if(click_3d.dist_sqrd(m1.vert(v0)) < click_3d.dist_sqrd(m2.vert(v1))) return v0;
    return m1.num_verts()+v1;
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


uint pick_handle(const std::vector<cinolib::vec3d> & m1, const std::vector<cinolib::vec3d> & m2, const cinolib::vec3d & click_3d)
{
    return *std::min_element(handles.begin(),
                             handles.end(),
                             [&](const uint a, const uint b)
                             {
                                 cinolib::vec3d pa = (a<m1.size()) ? m1[a] : m2[a-m1.size()];
                                 cinolib::vec3d pb = (b<m1.size()) ? m1[b] : m2[b-m1.size()];
                                 return click_3d.dist(pa) < click_3d.dist(pb);
                             });
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void update_input_coords(std::vector<double> & in_coords, std::vector<cinolib::vec3d> & m, const int offset)
{
    int count = offset;
    for(uint i = 0; i < m.size(); ++i)
    {
        in_coords[count++] = m[i][0];
        in_coords[count++] = m[i][1];
        in_coords[count++] = m[i][2];
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void update_handle_helpers(cinolib::GLcanvas & gui,
                           cinolib::DrawableTrimesh<> & m1,
                           cinolib::DrawableTrimesh<> & m2)
{
    gui.pop_all_markers();
    for(uint vid : handles)
    {
        if(vid<m1.num_verts())
            gui.push_marker(m1.vert(vid), "", cinolib::Color::PASTEL_ORANGE(), 5);
        else
            gui.push_marker(m2.vert(vid-m1.num_verts()), "", cinolib::Color::PASTEL_CYAN(), 5);
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main(int argc, char **argv)
{
    BoolOp op = UNION;
    std::vector<std::string> files = {"../data/armadillo.obj", "../data/armadillo.obj"};

    std::vector<double> in_coords;
    std::vector<uint> in_tris;
    std::vector<uint> in_labels;

    loadMultipleFiles(files, in_coords, in_tris, in_labels);

    cinolib::DrawableTrimesh<> m1(files.begin()->c_str());
    cinolib::DrawableTrimesh<> m2 = m1;

    m2.translate(cinolib::vec3d(m1.bbox().delta_y()*0.8,0,0));
    m2.update_bbox();
    m2.updateGL();

    std::cout << "Commands:" << std::endl;
    std::cout << "- press shift + left click to add handles" << std::endl;
    std::cout << "- press I for Intersection" << std::endl;
    std::cout << "- press U for Union" << std::endl;
    std::cout << "- press S for Subtraction" << std::endl;
    std::cout << "- press Space for the ARAP factorization" << std::endl;
    std::cout << "- press shift + left click to drag handles" << std::endl;

    // ARAP pipeline
    bool mode_init = true;

    cinolib::ARAP_data arap_m1;
    cinolib::ARAP_data arap_m2;
    update_input_coords(in_coords,m1.vector_verts(),0);
    update_input_coords(in_coords,m2.vector_verts(), m1.num_verts() * 3);

    std::vector<double>            bool_coords;
    std::vector<uint>              bool_tris;
    std::vector<std::bitset<NBIT>> bool_labels;
    booleanPipeline(in_coords, in_tris, in_labels, op, bool_coords, bool_tris, bool_labels);
    const cinolib::Color & c0 = cinolib::Color::PASTEL_ORANGE();
    const cinolib::Color & c1 = cinolib::Color::PASTEL_CYAN();
    uint n_tri = bool_tris.size()/3;
    std::vector<cinolib::Color> tri_colors(n_tri, c0);

    for(uint id=0; id<n_tri; ++id)
        if(bool_labels[id][1])
            tri_colors[id] = c1;

    cinolib::DrawableTriangleSoup soup(bool_coords, bool_tris, tri_colors);

    cinolib::GLcanvas gui;
    gui.push(&soup);
    gui.depth_cull_markers = false;
    update_handle_helpers(gui,m1,m2);

    uint curr_handle = 0;//*handles.begin();
    GLdouble zbuf = 0;
    gui.callback_mouse_left_click = [&](int mod) -> bool
    {
        cinolib::vec3d click_3d;
        cinolib::vec2d click_2d = gui.cursor_pos();
        if(mod & GLFW_MOD_SHIFT)
        {
            if(gui.unproject(click_2d, click_3d))
            {
                if(mode_init)
                {
                    handles.push_back(define_handle(m1,m2, click_3d));
                    update_handle_helpers(gui,m1,m2);
                    if(handles.back()<m1.num_verts()) arap_m1.bcs[handles.back()] = m1.vert(handles.back());
                    else                              arap_m2.bcs[handles.back()-m1.num_verts()] = m2.vert(handles.back()-m1.num_verts());
                }
                else
                {
                    curr_handle = pick_handle(arap_m1.xyz_out, arap_m2.xyz_out, click_3d);
                    zbuf = gui.query_Z_buffer(click_2d);
                }
            }
            return true;
        }
        return false;
    };

    /* FPS should count draw calls, not booleans completed!
    */

    std::atomic<int> count = 0;

    gui.callback_mouse_moved = [&](double x_pos, double y_pos) -> bool
    {
        if(mode_init) return false;

        bool left  = glfwGetMouseButton(gui.window,GLFW_MOUSE_BUTTON_LEFT)==GLFW_PRESS;
        bool shift = glfwGetKey(gui.window,GLFW_KEY_LEFT_SHIFT )==GLFW_PRESS ||
                     glfwGetKey(gui.window,GLFW_KEY_RIGHT_SHIFT)==GLFW_PRESS;
        if(left && shift)
        {
            // link the proper mesh depending on the current handle
            bool is_m1 = curr_handle < m1.num_verts();
            auto *m_ptr    = (is_m1) ? &m1 : & m2;
            auto *data_ptr = (is_m1) ? &arap_m1 : &arap_m2;
            int   handle   = (is_m1) ? curr_handle : curr_handle-m1.num_verts();

            cinolib::vec3d p;
            gui.unproject(cinolib::vec2d(x_pos,y_pos),zbuf, p);
            cinolib::vec3d delta = p - m_ptr->vert(handle);
            data_ptr->bcs[handle] += delta;
            ARAP(*m_ptr,*data_ptr);
            //
            bool_coords.clear();
            bool_tris.clear();
            if(is_m1) update_input_coords(in_coords,arap_m1.xyz_out,0);
            else      update_input_coords(in_coords,arap_m2.xyz_out,m1.num_verts()*3);
            booleanPipeline(in_coords, in_tris, in_labels, op, bool_coords, bool_tris, bool_labels);
            count++;
            n_tri = bool_tris.size()/3;
            tri_colors.resize(n_tri, c0);
            for(uint id=0; id<n_tri; ++id)
                tri_colors[id] = (bool_labels[id][0]) ? c0 : c1;
            soup = cinolib::DrawableTriangleSoup(bool_coords, bool_tris, tri_colors);

            update_handle_helpers(gui,m1,m2);
            gui.draw();

            return true;
        }
        return false;
    };

    gui.callback_key_pressed = [&](int key, int mod) -> bool
    {
        if(key==GLFW_KEY_I) op = INTERSECTION; else
        if(key==GLFW_KEY_U) op = UNION;        else
        if(key==GLFW_KEY_S) op = SUBTRACTION;  else
        if(key==GLFW_KEY_SPACE)
        {
            mode_init = false;
            ARAP(m1,arap_m1);
            ARAP(m2,arap_m2);
            std::cout << "ARAP has been initialized. Interactive Session Starts" << std::endl;
        }
        else if(key==GLFW_KEY_W)
        {
            std::vector<cinolib::Color> colors(bool_tris.size()/3);
            for(uint i=0; i<colors.size(); ++i)
            {
                colors[i] = (bool_labels[i][0]) ? cinolib::Color::PASTEL_YELLOW() : cinolib::Color::PASTEL_ORANGE();
            }
        }
        else return false;
        bool_coords.clear();
        bool_tris.clear();
        booleanPipeline(in_coords, in_tris, in_labels, op, bool_coords, bool_tris, bool_labels);
        n_tri = bool_tris.size()/3;
        tri_colors.resize(n_tri, c0);
        for(uint id=0; id<n_tri; ++id)
            tri_colors[id] = (bool_labels[id][0]) ? c0 : c1;
        soup = cinolib::DrawableTriangleSoup(bool_coords, bool_tris, tri_colors);
        return true;
    };

    std::promise<void> fps_exit_signal;
    std::thread fps_thread([&]()
    {
       std::future<void> exit_condition = fps_exit_signal.get_future();
       while(exit_condition.wait_for(std::chrono::milliseconds(1000)) == std::future_status::timeout)
       {
           count = 0;
       }
    });

    return gui.launch();
}

