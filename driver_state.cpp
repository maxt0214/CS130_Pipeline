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
    state.image_depth = new float[width * height];
    for(int i = 0; i < width*height; i++) {
        state.image_color[i] = make_pixel(0,0,0);
        state.image_depth[i] = 99999.0f;
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
                    geoData[j]->data = currVertex.data;
                    state.vertex_shader(currVertex, *geoData[j], state.uniform_data);
                }
                clip_triangle(state, (const data_geometry **)geoData,0);
                
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
    int axis, bound = 2 * (face % 2) - 1;// -1 for -w, 1 for w
    if(face < 2) {
        axis = 0; //x-axis
    } else if(face < 4) {
        axis = 1; //y-axis
    } else {
        axis = 2; //z-axis
    }
    
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    
    bool isInside[3];
    for(int i = 0; i < 3; i++) {
        if(bound < 0) {
            isInside[i] = in[i]->gl_Position[axis] >= bound * in[i]->gl_Position[3];
        } else {
            isInside[i] = in[i]->gl_Position[axis] <= bound * in[i]->gl_Position[3];
        }
    }
    
    if(!isInside[0] && !isInside[1] && !isInside[2]) {//none inside
        return;
    }
    
    if(isInside[0] && isInside[1] && isInside[2]) {//all inside
        clip_triangle(state,in,face+1);
        return;
    }
    
    float alpha, beta;//bary to get new vertex
    data_geometry* tri[3];
    
    tri[0] = new data_geometry;
    tri[1] = new data_geometry;
    tri[2] = new data_geometry;
    
    tri[0]->data = new float[state.floats_per_vertex];
    tri[1]->data = new float[state.floats_per_vertex];
    tri[2]->data = new float[state.floats_per_vertex];
    //only one vertex inside
    if(isInside[0] && !isInside[1] && !isInside[2]) {
        alpha = BaryOfTwoPoints(in[0]->gl_Position, in[1]->gl_Position);//A B
        beta = BaryOfTwoPoints(in[0]->gl_Position, in[2]->gl_Position);//A C
        
        SetNewTri(state, tri[0], in[0], nullptr, 1); //A
        SetNewTri(state, tri[1], in[0], in[1], alpha);//new B
        SetNewTri(state, tri[2], in[0], in[2], beta);//new C
        
        clip_triangle(state,(const data_geometry **)tri,face+1);
    }
    else if(!isInside[0] && isInside[1] && !isInside[2]) {
        alpha = BaryOfTwoPoints(in[1]->gl_Position, in[0]->gl_Position);//B A
        beta = BaryOfTwoPoints(in[1]->gl_Position, in[2]->gl_Position);//B C
        
        SetNewTri(state, tri[0], in[1], nullptr, 1);
        SetNewTri(state, tri[1], in[1], in[0], alpha);
        SetNewTri(state, tri[2], in[1], in[2], beta);
        
        clip_triangle(state,(const data_geometry **)tri,face+1);
    }
    else if(!isInside[0] && !isInside[1] && isInside[2]) {
        alpha = BaryOfTwoPoints(in[2]->gl_Position, in[0]->gl_Position);//C A
        beta = BaryOfTwoPoints(in[2]->gl_Position, in[1]->gl_Position);//C B
        
        SetNewTri(state, tri[0], in[2], nullptr, 1);
        SetNewTri(state, tri[1], in[2], in[0], alpha);
        SetNewTri(state, tri[2], in[2], in[1], beta);
        
        clip_triangle(state,(const data_geometry **)tri,face+1);
    }
    
    //two vertices inside
    else if(isInside[0] && isInside[1] && !isInside[2]) {
        alpha = BaryOfTwoPoints(in[0]->gl_Position, in[2]->gl_Position);//A C
        beta = BaryOfTwoPoints(in[1]->gl_Position, in[2]->gl_Position);// B C
        
        SetNewTri(state, tri[0], in[0], nullptr, 1);
        SetNewTri(state, tri[1], in[1], nullptr, 1);
        SetNewTri(state, tri[2], in[0], in[2], alpha);
        clip_triangle(state,(const data_geometry **)tri,face+1);
        
        SetNewTri(state, tri[0], in[1], nullptr, 1);
        SetNewTri(state, tri[1], in[0], in[2], alpha);
        SetNewTri(state, tri[2], in[1], in[2], beta);
        clip_triangle(state,(const data_geometry **)tri,face+1);
    }
    else if(isInside[0] && !isInside[1] && isInside[2]) {
        alpha = BaryOfTwoPoints(in[0]->gl_Position, in[1]->gl_Position);//A B
        beta = BaryOfTwoPoints(in[2]->gl_Position, in[1]->gl_Position);// C B
        
        SetNewTri(state, tri[0], in[0], nullptr, 1);
        SetNewTri(state, tri[1], in[2], nullptr, 1);
        SetNewTri(state, tri[2], in[0], in[1], alpha);
        clip_triangle(state,(const data_geometry **)tri,face+1);
        
        SetNewTri(state, tri[0], in[2], nullptr, 1);
        SetNewTri(state, tri[1], in[0], in[1], alpha);
        SetNewTri(state, tri[2], in[2], in[1], beta);
        clip_triangle(state,(const data_geometry **)tri,face+1);
    }
    else {
        alpha = BaryOfTwoPoints(in[1]->gl_Position, in[0]->gl_Position);//B A
        beta = BaryOfTwoPoints(in[2]->gl_Position, in[0]->gl_Position);// C A
        
        SetNewTri(state, tri[0], in[1], nullptr, 1);
        SetNewTri(state, tri[1], in[2], nullptr, 1);
        SetNewTri(state, tri[2], in[1], in[0], alpha);
        clip_triangle(state,(const data_geometry **)tri,face+1);
        
        SetNewTri(state, tri[0], in[2], nullptr, 1);
        SetNewTri(state, tri[1], in[1], in[0], alpha);
        SetNewTri(state, tri[2], in[2], in[0], beta);
        clip_triangle(state,(const data_geometry **)tri,face+1);
    }
    
    delete[] tri[0]->data;
    delete[] tri[1]->data;
    delete[] tri[2]->data;
    delete tri[0];
    delete tri[1];
    delete tri[2];
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    //declare frag data for later use
    data_fragment frag;
    frag.data = new float[state.floats_per_vertex];
    
    //calculate each vertex of the tri on the screen
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
    //an array of vertex bary for each pixel
    float* triBary = new float[3];
    
    //calculate the z value of each vertex
    float z[3];
    float depth;
    z[0] = in[0]->gl_Position[2]/in[0]->gl_Position[3]; //z value for A
    z[1] = in[1]->gl_Position[2]/in[1]->gl_Position[3]; //z value for A
    z[2] = in[2]->gl_Position[2]/in[2]->gl_Position[3]; //z value for A
    
    for(unsigned i = 0; i < state.image_width; i++) {
        for(unsigned j = 0; j < state.image_height; j++) {
            float pbc = (b[0]*c[1] - c[0]*b[1]) - (i*c[1] - c[0]*j) + (i*b[1] - b[0]*j);
            float apc = (i*c[1] - c[0]*j) - (a[0]*c[1] - c[0]*a[1]) + (a[0]*j - i*a[1]);
            float abp = (b[0]*j - i*b[1]) - (a[0]*j - i*a[1]) + (a[0]*b[1] - b[0]*a[1]);
            
            triBary[0] = pbc / tri;
            triBary[1] = apc / tri;
            triBary[2] = abp / tri;
            //find the depth of current pixel for z-buffering
            depth = triBary[0] * z[0] + triBary[1] * z[1] + triBary[2] * z[2];
            
            if(triBary[0] >= 0 && triBary[1] >= 0 && triBary[2] >= 0 && depth < state.image_depth[i+j*state.image_width]) {
                state.image_color[i+j*state.image_width] = PixelColor(state,in,frag,triBary);
                state.image_depth[i+j*state.image_width] = depth;
            }
        }
    }
    delete[] triBary;
    delete[] frag.data;
}

pixel PixelColor(driver_state& state, const data_geometry* in[3], data_fragment& frag, const float* triBary) {
    data_output out;
    
    for(int i = 0; i < state.floats_per_vertex; i++) {
        switch (state.interp_rules[i]) {
            case interp_type::flat:
                frag.data[i] = in[0]->data[i];
                break;
            case interp_type::smooth:
                break;
            case interp_type::noperspective:
                frag.data[i] = in[0]->data[i] * triBary[0] + in[1]->data[i] * triBary[1] + in[2]->data[i] * triBary[2];
                break;
            default:
                break;
        }
    }
    state.fragment_shader(frag, out, state.uniform_data);
    
    return make_pixel(out.output_color[0] * 255, out.output_color[1] * 255, out.output_color[2] * 255);
}

void SetNewTri(driver_state& state, data_geometry* newVer,const data_geometry* in1,const data_geometry* in2,float bary) {
    if(in2 == nullptr) {
        newVer->gl_Position = in1->gl_Position;
        for(int i = 0; i < state.floats_per_vertex; i++) {
            newVer->data[i] = in1->data[i];
        }
        return;
    }
    
    newVer->gl_Position = bary * in1->gl_Position + (1-bary) * in2->gl_Position;
    
    for(int i = 0; i < state.floats_per_vertex; i++) {
        newVer->data[i] = bary * in1->data[i] + (1-bary) * in2->data[i];
    }
}

float BaryOfTwoPoints(vec4 A, vec4 B) {
    return (B[3] - B[0]) / (A[0] - A[3] + B[3] - B[0]);
}
