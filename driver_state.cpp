#include "driver_state.h"
#include <cstring>
#include <limits>

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
    state.image_color= new pixel[height * width];
    pixel temp_pix = make_pixel(0,0,0);
    for (int i = 0; i < (width*height);i++){
        state.image_color[i] = temp_pix;
    }
    state.image_depth= new float[height * width];
    for (int i = 0; i < (height * width); i++){
        //state.image_depth[i] = std::numeric_limits<float>::max();
        state.image_depth[i] = 1;
    }
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
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
    //std::cout<<"TODO: implement rendering."<<std::endl;
    switch(type){
        case render_type::triangle:
        {
            const data_geometry * out[3];
            data_geometry d_geo[3];
            data_vertex v[3];

            int numTriangles = state.num_vertices  / 3.0;
            int p = 0;
           
            for(int c = 0; c < numTriangles; c++){
                for(int i = 0; i < 3; i++){
                    v[i].data = &state.vertex_data[p];
    
                    d_geo[i].data = v[i].data;

                    state.vertex_shader(v[i], d_geo[i], state.uniform_data);

                    out[i] = &d_geo[i];

                    p += state.floats_per_vertex; 
  
                        
                }
                //std::cout << c << std::endl;
                //rasterize_triangle(state,out);
                clip_triangle(state,out,0);
            }
            break;
        }

        case render_type::indexed:
        {
            // const data_geometry * out[3];
            // data_geometry d_geo[3];
            // data_vertex v[3];        

            // for(int i = 0,j = 0; i < 3*state.num_triangles ; i++, j++){
            //     v[i].data = (float*)state.index_data[i];
            //     std::cout << v[i].data << std::endl;
            // }
            // break;
        }
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
    const data_geometry * new_data[3];
    data_geometry temp_dg[3];

    new_data[0] = in[0];
    new_data[1] = in[1];
    new_data[2] = in[2];

    temp_dg[0] = *in[0];
    temp_dg[1] = *in[1];
    temp_dg[2] = *in[2];

    vec4 a = in[0]->gl_Position;
    vec4 b = in[1]->gl_Position;
    vec4 c = in[2]->gl_Position;


    switch(face){
        // case 0: //right plane
        // {
        //     if(a[0] <= a[3] && b[0] <= b[3] && c[0] <= c[3]){
        //         new_data[0] = in[0];
        //         new_data[1] = in[1];
        //         new_data[2] = in[2];
        //         //std::cout << "penis1\n";
        //     }
        //     if(c[0] <= c[3] && b[0] > b[3] && a[0] > a[3]){
        //         //std::cout << "penis2\n";
        //         new_data[0] = in[2]; // new position will be the c position
        //         alpha = (c[3] - c[0]) / (a[0] - a[3] + c[3] - c[0]); // a -> c
        //         p = alpha * a + (1-alpha) * c;
        //         temp_dg[1].gl_Position = p;
        //         new_data[1] = &temp_dg[1];
        //         alpha = (b[3] - b[0]) / (c[0] - c[3] + b[3] - b[0]); // c -> b
        //         p = alpha * c + (1-alpha) * b;
        //         temp_dg[2].gl_Position = p;
        //         new_data[2] = &temp_dg[2]; 
        //     }

        //     if(b[0] <= b[3] && c[0] > c[3] && a[0] > a[3]){
        //         //std::cout << "penis2\n";
        //         new_data[0] = in[1]; // new position will be the c position
        //         alpha = (b[3] - b[0]) / (a[0] - a[3] + b[3] - b[0]); // a -> c
        //         p = alpha * a + (1-alpha) * b;
        //         temp_dg[1].gl_Position = p;
        //         new_data[1] = &temp_dg[1];
        //         alpha = (c[3] - c[0]) / (b[0] - b[3] + c[3] - c[0]); // c -> b
        //         p = alpha * b + (1-alpha) * c;
        //         temp_dg[2].gl_Position = p;
        //         new_data[2] = &temp_dg[2]; 
        //     }
        //     if(a[0] <= a[3] && b[0] > b[3] && c[0] > c[3]){
        //         // std::cout << "penis2\n";
        //         new_data[0] = in[0]; // new position will be the c position
        //         alpha = (a[3] - a[0]) / (b[0] - b[3] + a[3] - a[0]); // a -> c
        //         p = alpha * b + (1-alpha) * a;
        //         temp_dg[1].gl_Position = p;
        //         new_data[1] = &temp_dg[1];
        //         alpha = (c[3] - c[0]) / (a[0] - a[3] + c[3] - c[0]); // c -> b
        //         p = alpha * a + (1-alpha) * c;
        //         temp_dg[2].gl_Position = p;
        //         new_data[2] = &temp_dg[2]; 
        //     }

        //     if(b[0] <= b[3] && c[0] <= c[3] && a[0] > a[3]){
        //         vec4 p1;
        //         vec4 p2;

        //         alpha = (a[3] - a[0]) / (c[0] - c[3] + a[3] - a[0]);
        //         p1 = alpha * c + (1-alpha) * a;
        //         alpha = (b[3] - b[0]) / (a[0] - a[3] + b[3] - b[0]);
        //         p2 = alpha * a + (1-alpha) * b;


        //     }

        //     break;
        // }
        case 1:
        {
            if(c[0] >= -c[3] && a[0] < -a[3] && b[0] < -b[3]){
                //std::cout << "c in\n";
            }

            if(b[0] >= -b[3] && a[0] < -a[3] && c[0] < -c[3]){
                //std::cout << "b in\n";
            }

            if(a[0] >= -c[3] && c[0] < -c[3] && b[0] < -b[3]){
                //std::cout << "a in\n";
            }

            if(b[0] >= -b[3] && c[0] >= -c[3] && a[0] < -a[3]){
                //std::cout << "b and c in\n";
            }

            if(a[0] >= -a[3] && c[0] >= -c[3] && b[0] < -b[3]){
                //std::cout << "a and c in\n";
            }

            if(a[0] >= -a[3] && b[0] >= -b[3] && c[0] < -c[3]){
                //std::cout << "a and b in\n";
            }
            if(a[0] >= -a[3] && c[0] >= -c[3] && b[0] >= -b[3]){
                //std::cout << "a,b,c in\n";
            }
            break;
        }

        case 3: //bottom case
        {
            if(a[1] >= -a[3] && b[1] < -b[3] && c[1] < -c[3]){
               // std::cout << "a in near \n";
                float a1;
                float a2;
                a1 = (-1*a[3] - a[1]) / (c[1] + c[3] - a[3] - a[1]);
                vec4 c_p = a1 * c + (1-a1) * a; // c to a
                a2 = (-1*b[3] - b[1]) / (a[1] + a[3] - b[3] - b[1]);
                vec4 a_p = a2 * a + (1-a2) * b; // a to b

                temp_dg[1].gl_Position = c_p;
                temp_dg[2].gl_Position = a_p;
                new_data[0] = in[0];
                new_data[1] = &temp_dg[1];
                new_data[2] = &temp_dg[2];
            }

            if(b[1] >= -b[3] && a[1] < -a[3] && c[1] < -c[3]){
                //std::cout << "b in near \n";
                float a1;
                float a2;
                a1 = (-1*b[3] - b[1]) / (a[1] + a[3] - b[3] - b[1]);
                vec4 a_p = a1 * a + (1-a1) * b; // a to b
                a2 = (-1*c[3] - c[1]) / (b[1] + b[3] - c[3] - c[1]);
                vec4 b_p = a2 * b + (1-a2) * c; // b to c

                temp_dg[2].gl_Position = b_p;
                temp_dg[0].gl_Position = a_p;
                new_data[0] = &temp_dg[0];
                new_data[1] = in[1];
                new_data[2] = &temp_dg[2];
                
            }

            if(c[1] >= -c[3] && b[1] < -b[3] && a[1] < -a[3]){
                //std::cout << "c in near \n";
                float a1;
                float a2;
                a1 = (-1*b[3] - b[1]) / (c[1] + c[3] - b[3] - b[1]);
                vec4 b_p = a1 * b + (1-a1) * c; // b -> c
                a2 = (-1*a[3] - a[1]) / (c[1] + c[3] - a[3] - a[1]);
                vec4 c_p = a2 * c + (1-a2) * a; // c->a

                temp_dg[0].gl_Position = c_p;
                temp_dg[1].gl_Position = b_p;
                new_data[0] = &temp_dg[0];
                new_data[1] = &temp_dg[1];
                new_data[2] = in[2];
                
            }

            if(a[1] >= -a[3] && b[1] >= -b[3] && c[1] < -c[3]){ // a and b in
                //std::cout << "a, b in near \n";
                float a1;
                float a2;
                a1 = (-1*a[3] - a[1]) / (c[1] + c[3] - a[3] - a[1]);
                vec4 c_p = a1 * c + (1-a1) * a; // c -> a
                a2 = (-1*c[3] - c[1]) / (b[1] + b[3] - c[3] - c[1]);
                vec4 b_p = a2 * b + (1-a2) * c; // b to c

                temp_dg[2].gl_Position = c_p;
                new_data[0] = in[0];
                new_data[1] = in[1];
                new_data[2] =  &temp_dg[2];
               
                clip_triangle(state, new_data, face+1);
                temp_dg[2].gl_Position = b_p;
                temp_dg[0].gl_Position = c_p;
                new_data[1] = in[1];
                new_data[2] = &temp_dg[2];
                new_data[0] = &temp_dg[0];
            } 

            if(a[1] >= -a[3] && c[1] >= -c[3] && b[1] < -b[3]){
                //std::cout << "a, c in near \n";
                float a1;
                float a2;
                a1 = (-1*c[3] - c[1]) / (b[1] + b[3] - c[3] - c[1]);
                vec4 b_p = a1 * b + (1-a1) * c; // b -> c
                a2 = (-1*b[3] - b[1]) / (a[1] + a[3] - b[3] - b[1]);
                vec4 a_p = a2 * a + (1-a2) * b; // a->b

                temp_dg[1].gl_Position = b_p;
                new_data[0] = in[0];
                new_data[1] =  &temp_dg[1];
                new_data[2] = in[2];
                clip_triangle(state, new_data, face+1);
                temp_dg[1].gl_Position = a_p;
                temp_dg[2].gl_Position = b_p;
                new_data[0] = in[0];
                new_data[1] = &temp_dg[1];
                new_data[2] = &temp_dg[2];

            }

            // if(b[1] >= -b[3] && c[1] >= -c[3] && a[1] < -a[3]){
            //     //std::cout << "b, c  in  \n";
            //     float a1;
            //     float a2;
            //     a1 = (-1*b[3] - b[1]) / (a[1] + a[3] - b[3] - b[1]);
            //     vec4 a_p = a1 * a + (1-a1) * b; // a to b
            //     a2 = (-1*c[3] - c[2]) / (a[2] + a[3] - c[3] - c[2]);
            //     vec4 p2 = a2 * a + (1-a2) * c;

            //     new_data[1] = in[1];
            //     new_data[2] = in[2];
            //     temp_dg[0].gl_Position = p1;
            //     new_data[0] = &temp_dg[0];
            //     clip_triangle(state, new_data, face+1);
            //     temp_dg[0].gl_Position = p1;
            //     temp_dg[1].gl_Position = p2;
            //     new_data[2] = in[2];
            //     new_data[1] = &temp_dg[1];
            //     new_data[0] = &temp_dg[0];

            // }
            if(b[1] >= -b[3] && c[1] >= -c[3] && a[1] >= -a[3]){
                new_data[0] = in[0];
                new_data[1] = in[1];
                new_data[2] = in[2];
            }
            break;
        }

        case 4: //far 
        {
            if(a[2] <= a[3] && c[2] > c[3] && b[2] > b[3]){
                //std::cout << "a in\n";
                float a1;
                float a2;
                a1 = (b[3] - b[2]) / (a[2] - a[3] + b[3] - b[2]);
                vec4 a_p = a1 * a + (1-a1) * b; //a to b
                a2 = (a[3] - a[2]) / (c[2] - c[3] + a[3] - a[2]);
                vec4 b_p = a2 * c + (1-a2) * a; // c to a

                temp_dg[0].gl_Position = a_p;
                temp_dg[1].gl_Position = b_p;
                new_data[0] = in[0];
                new_data[1] = &temp_dg[0];
                new_data[2] = &temp_dg[1];

            }

            if(c[2] <= b[3] && b[2] > b[3] && a[2] > a[3]){
                //std::cout << "b in\n";
                float a1;
                float a2;
                a1 = (a[3] - a[2]) / (c[2] - c[3] + a[3] - a[2]);
                vec4 c_p = a1 * c + (1-a1) * a; //c to a
                a2 = (c[3] - c[2]) / (b[2] - b[3] + c[3] - c[2]);
                vec4 b_p = a2 * b + (1-a2) * c; // c to a

                temp_dg[0].gl_Position = c_p;
                temp_dg[1].gl_Position = b_p;
                new_data[2] = in[2];
                new_data[0] = &temp_dg[0];
                new_data[1] = &temp_dg[1];
            }

            if(b[2] <= c[3] && a[2] > a[3] && c[2] > c[3]){
                //std::cout << "c in\n";
                float a1;
                float a2;
                a1 = (c[3] - c[2]) / (b[2] - b[3] + c[3] - c[2]);
                vec4 b_p = a1 * b + (1-a1) * c; //c to a
                a2 = (b[3] - b[2]) / (a[2] - a[3] + b[3] - b[2]);
                vec4 a_p = a2 * a + (1-a2) * b; // c to a

                temp_dg[0].gl_Position = b_p;
                temp_dg[1].gl_Position = a_p;
                new_data[1] = in[1];
                new_data[0] = &temp_dg[0];
                new_data[2] = &temp_dg[1];

            }

            if(a[2] <= a[3] && b[2] <= c[3] && c[2] > c[3]){
                //std::cout << "b and c in\n";
                float a1;
                float a2;
                a1 = (c[3] - c[2]) / (b[2] - b[3] + c[3] - c[2]);
                vec4 b_p = a1 * b + (1-a1) * c; //c to a
                a2 = (a[3] - a[2]) / (c[2] - c[3] + a[3] - a[2]);
                vec4 c_p = a2 * c + (1-a2) * a; // c to a

                new_data[0] = in[0];
                new_data[1] = in[1];
                temp_dg[2].gl_Position =  c_p;
                new_data[2] = &temp_dg[2];
                clip_triangle(state, new_data, face+1);
                new_data[1] = in[1];
                temp_dg[2].gl_Position = b_p;
                temp_dg[0].gl_Position = c_p;
                new_data[2] = &temp_dg[2];
                new_data[0] = &temp_dg[0];

            }

            if(a[2] <= a[3] && c[2] <= c[3] && b[2] > b[3]){
                //std::cout << "a and c in\n";
                float a1;
                float a2;
                a1 = (b[3] - b[2]) / (a[2] - a[3] + b[3] - b[2]);
                vec4 a_p = a1 * a + (1-a1) * b; //a to b
                a2 = (c[3] - c[2]) / (b[2] - b[3] + c[3] - c[2]);
                vec4 b_p = a2 * b + (1-a2) * c; // b to c

                new_data[0] = in[0];
                new_data[2] = in[2];
                temp_dg[1].gl_Position =  b_p;
                new_data[1] = &temp_dg[1];
                clip_triangle(state, new_data, face+1);
                new_data[0] = in[0];
                temp_dg[1].gl_Position = a_p;
                temp_dg[2].gl_Position = b_p;
                new_data[1] = &temp_dg[1];
                new_data[2] = &temp_dg[2];

            }

            if(b[2] <= b[3] && c[2] <= c[3] && a[2] > a[3]){
                //std::cout << "a and b in\n";
                float a1;
                float a2;
                a1 = (a[3] - a[2]) / (c[2] - c[3] + a[3] - a[2]);
                vec4 c_p = a1 * c + (1-a1) * a; //c to a
                a2 = (b[3] - b[2]) / (a[2] - a[3] + b[3] - b[2]);
                vec4 a_p = a2 * a + (1-a2) * b; // b to a

                new_data[1] = in[1];
                new_data[2] = in[2];
                temp_dg[0].gl_Position =  a_p;
                new_data[0] = &temp_dg[0];
                clip_triangle(state, new_data, face+1);
                new_data[2] = in[2];
                temp_dg[0].gl_Position = c_p;
                temp_dg[1].gl_Position = a_p;
                new_data[0] = &temp_dg[0];
                new_data[1] = &temp_dg[1];
            }
            if(a[2] <= a[3] && c[2] <= c[3] && b[2] <= b[3]){
                //std::cout << "a,b,c in\n";
                new_data[0] = in[0];
                new_data[1] = in[1];
                new_data[2] = in[2];
            }
            break;
        }

        case 5:
        {
            if(a[2] >= -a[3] && b[2] < -b[3] && c[2] < -c[3]){
               // std::cout << "a in near \n";
                float a1;
                float a2;
                a1 = (-1*b[3] - b[2]) / (a[2] + a[3] - b[3] - b[2]);
                vec4 a_p = a1 * a + (1-a1) * b; // a -> b
                a2 = (-1*a[3] - a[2]) / (c[2] + c[3] - a[3] - a[2]);
                vec4 c_p = a2 * c + (1-a2) * a; // c -> a

                temp_dg[0].gl_Position = a_p;
                temp_dg[1].gl_Position = c_p;
                new_data[0] = in[0];
                new_data[1] = &temp_dg[0];
                new_data[2] = &temp_dg[1];
            }

            if(b[2] >= -b[3] && a[2] < -a[3] && c[2] < -c[3]){
                //std::cout << "b in near \n";
                float a1;
                float a2;
                a1 = (-1*a[3] - a[2]) / (b[2] + b[3] - a[3] - a[2]);
                vec4 a_p = a1 * b + (1-a1) * a; // b -> a
                a2 = (-1*b[3] - b[2]) / (c[2] + c[3] - b[3] - b[2]);
                vec4 b_p = a2 * c + (1-a2) * b; // c -> b

                temp_dg[0].gl_Position = b_p;
                temp_dg[2].gl_Position = a_p;
                new_data[0] = &temp_dg[0];
                new_data[1] = in[1];
                new_data[2] = &temp_dg[2];
                
            }

            if(c[2] >= -c[3] && b[2] < -b[3] && a[2] < -a[3]){
                //std::cout << "c in near \n";
                float a1;
                float a2;
                a1 = (-1*c[3] - c[2]) / (b[2] + b[3] - c[3] - c[2]);
                vec4 b_p = a1 * b + (1-a1) * c; // b -> c
                a2 = (-1*a[3] - a[2]) / (c[2] + c[3] - a[3] - a[2]);
                vec4 a_p = a2 * c + (1-a2) * a; // c->a

                temp_dg[0].gl_Position = a_p;
                temp_dg[1].gl_Position = b_p;
                new_data[0] = &temp_dg[0];
                new_data[1] = &temp_dg[1];
                new_data[2] = in[2];
                
            }

            if(a[2] >= -a[3] && b[2] >= -b[3] && c[2] < -c[3]){ // a and b in
                //std::cout << "a, b in near \n";
                float a1;
                float a2;
                a1 = (-1*c[3] - c[2]) / (a[2] + a[3] - c[3] - c[2]);
                vec4 a_p = a1 * a + (1-a1) * c; // a -> c
                a2 = (-1*b[3] - b[2]) / (c[2] + c[3] - b[3] - b[2]);
                vec4 b_p = a2 * c + (1-a2) * b; // c->b

                temp_dg[2].gl_Position = a_p;
                new_data[0] = in[0];
                new_data[1] = in[1];
                new_data[2] =  &temp_dg[2];
               
                clip_triangle(state, new_data, face+1);
                temp_dg[0].gl_Position = a_p;
                temp_dg[1].gl_Position = b_p;
                new_data[0] = &temp_dg[0];
                new_data[1] = &temp_dg[1];
                new_data[2] = in[1];
            }

            if(a[2] >= -a[3] && c[2] >= -c[3] && b[2] < -b[3]){
                //std::cout << "a, c in near \n";
                float a1;
                float a2;
                a1 = (-1*c[3] - c[2]) / (b[2] + b[3] - c[3] - c[2]);
                vec4 c_p = a1 * b + (1-a1) * c; // b -> c
                a2 = (-1*b[3] - b[2]) / (a[2] + a[3] - b[3] - b[2]);
                vec4 a_p = a2 * a + (1-a2) * b; // a->b

                temp_dg[0].gl_Position = a_p;
                new_data[0] = in[0];
                new_data[1] =  &temp_dg[0];
                new_data[2] = in[2];
                clip_triangle(state, new_data, face+1);
                temp_dg[0].gl_Position = a_p;
                temp_dg[1].gl_Position = c_p;
                new_data[0] = &temp_dg[0];
                new_data[1] = &temp_dg[1];
                new_data[2] = in[2];
            }

            if(b[2] >= -b[3] && c[2] >= -c[3] && a[2] < -a[3]){
                //std::cout << "b, c  in near \n";
                float a1;
                float a2;
                a1 = (-1*b[3] - b[2]) / (a[2] + a[3] - b[3] - b[2]);
                vec4 a_p = a1 * a + (1-a1) * b;
                a2 = (-1*a[3] - a[2]) / (c[2] + c[3] - a[3] - a[2]);
                vec4 c_p = a2 * c + (1-a2) * a;

                temp_dg[0].gl_Position = a_p;
                new_data[0] = &temp_dg[0];
                new_data[1] = in[1];
                new_data[2] = in[2];
                clip_triangle(state, new_data, face+1);
                temp_dg[0].gl_Position = c_p;
                temp_dg[1].gl_Position = a_p;
                new_data[0] = &temp_dg[0];
                new_data[1] = &temp_dg[1];
                new_data[2] = in[2];

            }
            if(b[2] >= -b[3] && c[2] >= -c[3] && a[2] >= -a[3]){
                new_data[0] = in[0];
                new_data[1] = in[1];
                new_data[2] = in[2];
            }
            break;
        }
    } 

    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    //clip_triangle(state,in,face+1);
    clip_triangle(state,new_data,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    //for converting to NDC
    float w0 = in[0]->gl_Position[3];
    float w1 = in[1]->gl_Position[3];
    float w2 = in[2]->gl_Position[3];
  
    float x0 = in[0]->gl_Position[0] / w0;
    float x1 = in[1]->gl_Position[0] / w1;
    float x2 = in[2]->gl_Position[0] / w2;

    float y0 = in[0]->gl_Position[1] / w0;
    float y1 = in[1]->gl_Position[1] / w1;
    float y2 = in[2]->gl_Position[1] / w2;

    float z0 = in[0]->gl_Position[2] / w0;
    float z1 = in[1]->gl_Position[2] / w1;
    float z2 = in[2]->gl_Position[2] / w2;

    // std::cout << "vertex a:" << std::endl;
    // for (int i = 0; i < state.floats_per_vertex; i++){
    //     std::cout << in[0]->gl_Position[i] << " ";
    // }
    // std::cout << std::endl;

    vec3 Z = {z0, z1, z2};

   

    float Ax = x0 * (state.image_width/2) + ((state.image_width/2) - 0.5);
    float Ay = y0 * (state.image_height/2) + ((state.image_height/2) - 0.5);
   

    float Bx = x1 * (state.image_width/2) + ((state.image_width/2) - 0.5);
    float By = y1 * (state.image_height/2) + ((state.image_height/2) - 0.5);
   

    float Cx = x2 * (state.image_width/2) + ((state.image_width/2) - 0.5);
    float Cy = y2 * (state.image_height/2) + ((state.image_height/2) - 0.5);


    vec2 A = {Ax,Ay};
    vec2 B = {Bx,By};
    vec2 C = {Cx,Cy};


    float area_abc = (C[0] - A[0]) * (B[1] - A[1]) - (C[1] - A[1]) * (B[0] - A[0]);


    for(int i = 0; i < state.image_width; i++){
        for (int j = 0; j < state.image_height; j++){

            vec2 P = {(float)i,(float)j};
            int index = i + j * state.image_width;

            float area_pbc = (P[0] - B[0]) * (C[1] - B[1]) - (P[1] - B[1]) * (C[0] - B[0]); // cross((p - b), (c-b))
            //std::cout << "area pbc: " << area_pbc << " ";
            float area_pca = (P[0] - C[0]) * (A[1] - C[1]) - (P[1] - C[1]) * (A[0] - C[0]); // cross ((P-C),(A-C))
            //std::cout << "area pca: " << area_pca << " ";
            float area_pba = (P[0] - A[0]) * (B[1] - A[1]) - (P[1] - A[1]) * (B[0] - A[0]); // cross ((P-A), (B-A))
            //std::cout << "area pba: " << area_pba << " " << std::endl;

            float alpha_p = area_pbc / area_abc;
            float beta_p = area_pca / area_abc;
            float gamma_p = area_pba / area_abc;


            float alpha,beta,gamma;

            if(alpha_p >= 0 && beta_p >= 0 && gamma_p >= 0){
                data_output out;
                data_fragment frag;
                frag.data = new float[MAX_FLOATS_PER_VERTEX];
                float depth = (alpha_p * Z[0]) + (beta_p * Z[1]) + (gamma_p * Z[2]);
                if(state.image_depth[index] > depth){
                    for(int c = 0; c < state.floats_per_vertex; c++){
                        if (state.interp_rules[c] == interp_type::flat){
                            frag.data[c] = in[0]->data[c];
                        }
                        if (state.interp_rules[c] == interp_type::noperspective){
                            frag.data[c] = (alpha_p * in[0]->data[c]) + (beta_p * in[1]->data[c]) + (gamma_p * in[2]->data[c]);
                        }
                        if (state.interp_rules[c] == interp_type::smooth){
                            float k = (alpha_p / w0) + (beta_p / w1) + (gamma_p / w2);
                            alpha = alpha_p / w0 / k; 
                            beta = beta_p / w1 / k; 
                            gamma = gamma_p / w2 / k; 
                            frag.data[c] = (alpha * in[0]->data[c]) + (beta * in[1]->data[c]) + (gamma * in[2]->data[c]);

                        }

                    }
                    state.fragment_shader(frag, out, state.uniform_data);
                    float r = out.output_color[0] * 255;
                    float g = out.output_color[1] * 255;
                    float b = out.output_color[2] * 255;
                    state.image_color[index] = make_pixel(r,g,b);
                    state.image_depth[index] = depth;
                }
            }
        }
    }


    //std::cout<<"TODO: implement rasterization"<<std::endl;
}

