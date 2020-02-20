#include "driver_state.h"
#include <cstring>
#include <iostream>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color = new pixel[width * height];
    state.image_depth=0;
    for(int i = 0; i < width*height; i++) {
        state.image_color[i] = make_pixel(0,0,0);
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    switch (type) {
        case render_type::triangle:
            for(int i = 0; i < state.num_vertices; i += 3) {
                data_geometry* geoData[3];
                geoData[0] = new data_geometry;
                geoData[1] = new data_geometry;
                geoData[2] = new data_geometry;
                for(int j = 0; j < 3; j++) {//set geo data for single triangle here
                    data_vertex currVertex;
                    currVertex.data = state.vertex_data + (j+i) * state.floats_per_vertex;
                    state.vertex_shader(currVertex, *geoData[j], state.uniform_data);
                }
                rasterize_triangle(state, (const data_geometry **)geoData);
                
                delete geoData[0];
                delete geoData[1];
                delete geoData[2];
            }
            
            break;
        case render_type::indexed:
            break;
        case render_type::fan:
            break;
        case render_type::strip:
            break;
        case render_type::invalid:
            break;
    }
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    float ax = in[0]->gl_Position[0]/in[0]->gl_Position[3];
    float ay = in[0]->gl_Position[1]/in[0]->gl_Position[3];
    
    float bx = in[1]->gl_Position[0]/in[1]->gl_Position[3];
    float by = in[1]->gl_Position[1]/in[1]->gl_Position[3];
    
    float cx = in[2]->gl_Position[0]/in[2]->gl_Position[3];
    float cy = in[2]->gl_Position[1]/in[2]->gl_Position[3];
    
    vec2 a((ax + 1)/2 * state.image_width - 0.5f, (ay + 1)/2 * state.image_height -0.5f);
    vec2 b((bx + 1)/2 * state.image_width - 0.5f, (by + 1)/2 * state.image_height - 0.5f);
    vec2 c((cx + 1)/2 * state.image_width - 0.5f, (cy + 1)/2 * state.image_height -0.5f);
    
    float tri = (b[0]*c[1] - c[0]*b[1]) - (a[0]*c[1] - c[0]*a[1]) + (a[0]*b[1] - b[0]*a[1]);
    for(unsigned i = 0; i < state.image_width; i++) {
        for(unsigned j = 0; j < state.image_height; j++) {
            float pbc = (b[0]*c[1] - c[0]*b[1]) - (i*c[1] - c[0]*j) + (i*b[1] - b[0]*j);
            float apc = (i*c[1] - c[0]*j) - (a[0]*c[1] - c[0]*a[1]) + (a[0]*j - i*a[1]);
            float abp = (b[0]*j - i*b[1]) - (a[0]*j - i*a[1]) + (a[0]*b[1] - b[0]*a[1]);
            
            float alpha = pbc / tri;
            float beta = apc / tri;
            float phi = abp / tri;
           
            if(alpha >= 0 && beta >= 0 && phi >= 0) {
                state.image_color[i+j*state.image_width] = make_pixel(255,255,255);
            }
        }
    }
}

