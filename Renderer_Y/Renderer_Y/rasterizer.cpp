// clang-format off
//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include <vector>
#include "rasterizer.hpp"
//#include <opencv2/opencv.hpp>
#include <math.h>
#include <eigen3/Eigen/Geometry>
#include <iostream>

rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f> &positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return {id};
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i> &indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return {id};
}

rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f> &cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return {id};
}

auto to_vec4(const Eigen::Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}


static bool insideTriangle(int x, int y, const Vector3f* _v)
{   
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]
    return false;
}

static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector3f* v)
{
    float c1 = (x*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*y + v[1].x()*v[2].y() - v[2].x()*v[1].y()) / (v[0].x()*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*v[0].y() + v[1].x()*v[2].y() - v[2].x()*v[1].y());
    float c2 = (x*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*y + v[2].x()*v[0].y() - v[0].x()*v[2].y()) / (v[1].x()*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*v[1].y() + v[2].x()*v[0].y() - v[0].x()*v[2].y());
    float c3 = (x*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*y + v[0].x()*v[1].y() - v[1].x()*v[0].y()) / (v[2].x()*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*v[2].y() + v[0].x()*v[1].y() - v[1].x()*v[0].y());
    return {c1,c2,c3};
}

void rst::rasterizer::draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type)
{
    auto& buf = pos_buf[pos_buffer.pos_id];
    auto& ind = ind_buf[ind_buffer.ind_id];
    auto& col = col_buf[col_buffer.col_id];

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model;
    for (auto& i : ind)
    {
        Triangle t;
        Eigen::Vector4f v[] = {
                mvp * to_vec4(buf[i[0]], 1.0f),
                mvp * to_vec4(buf[i[1]], 1.0f),
                mvp * to_vec4(buf[i[2]], 1.0f)
        };
        //Homogeneous division
        for (auto& vec : v) {
            vec /= vec.w();
        }
        //Viewport transformation
        for (auto & vert : v)
        {
            vert.x() = 0.5*width*(vert.x()+1.0);
            vert.y() = 0.5*height*(vert.y()+1.0);
            vert.z() = vert.z() * f1 + f2;
        }

        for (int i = 0; i < 3; ++i)
        {
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
        }

        auto col_x = col[i[0]];
        auto col_y = col[i[1]];
        auto col_z = col[i[2]];
        
        t.setColor(0, col_x[0], col_x[1], col_x[2]);
        t.setColor(1, col_y[0], col_y[1], col_y[2]);
        t.setColor(2, col_z[0], col_z[1], col_z[2]);

        rasterize_triangle(t);
        
        // settings x
//        Triangle test;
//        Eigen::Vector3f ta(300,600,1);
//        Eigen::Vector3f tb(670,30,1);
//        Eigen::Vector3f tc(60,200,1);
//
//        test.setVertex(0, ta);
//        test.setVertex(1, tb);
//        test.setVertex(2, tc);
//
//        test.setColor(0, col_x[0], col_x[1], col_x[2]);
//        test.setColor(1, col_y[0], col_y[1], col_y[2]);
//        test.setColor(2, col_z[0], col_z[1], col_z[2]);
//
//        rasterize_triangle(test);
    }
}

//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle& t) {
//    auto v = t.toVector4();
    //
/*
 {2, 0, -2},
                    {0, 2, -2},
                    {-2, 0, -2},
                    {3.5, -1, -5},
                    {2.5, 1.5, -5},
                    {-1, 0.5, -5}
 */
    auto xmin = std::min(t.v[0].x(), t.v[1].x());
    xmin = std::min(xmin, t.v[2].x());
    auto xmax = std::max(t.v[0].x(), t.v[1].x());
    xmax = std::max(xmax, t.v[2].x());

    auto ymin = std::min(t.v[0].y(), t.v[1].y());
    ymin = std::min(xmin, t.v[2].y());
    auto ymax = std::max(t.v[0].y(), t.v[1].y());
    ymax = std::max(xmax, t.v[2].y());
    
    Eigen::Vector3f a(t.v[0].x(),t.v[0].y(),1);
    Eigen::Vector3f b(t.v[1].x(),t.v[1].y(),1);
    Eigen::Vector3f c(t.v[2].x(),t.v[2].y(),1);
    
    std::cout << "=========== one raster ============\n\n\n\n";
    int tick = 0;
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
//            if (i < xmin || i > xmax || j < ymin || j > ymax) {
//                Eigen::Vector3f pt(i,j,1);
//                Vector3f tempColor(0,0,0);
//                set_pixel(pt, tempColor);
//                continue;
//            }
            
//            Eigen::Vector3f p(i + 0.5f,j + 0.5f,1);
            Eigen::Vector3f p(i,j,1);
            Eigen::Vector3f ap = p - a;
            Eigen::Vector3f bp = p - b;
            Eigen::Vector3f cp = p - c;
            
            float ret1 = ap.cross(b - a).z();
            float ret2 = bp.cross(c - b).z();
            float ret3 = cp.cross(a - c).z();
            
            Eigen::Vector3f pt(i,j,1);
            if ((ret1 > 0 && ret2 > 0 && ret3 > 0) || (ret1 < 0 && ret2 < 0 && ret3 < 0))  {
                auto color = t.getColor();
                set_pixel(pt, color);
                tick++;
                if (tick > 50) {
                    tick = 0;
                    printf("x:%f,y:%f,r:%f,g:%f,b:%f ",pt.x(),pt.y(),color.x(),color.y(),color.z());
                }
                
            }else{
                Vector3f tmpColor(0,0,0);
//                set_pixel(pt, tmpColor);
            }
        }
    }

    
    // TODO : Find out the bou.nding box of current triangle.
    // iterate through the pixel and find if the current pixel is inside the triangle

    // If so, use the following code to get the interpolated z value.
    //auto[alpha, beta, gamma] = computeBarycentric2D(x, y, t.v);
    //float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
    //float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
    //z_interpolated *= w_reciprocal;

    // TODO : set the current pixel (use the set_pixel function) to the color of the triangle (use getColor function) if it should be painted.
}

void rst::rasterizer::set_model(const Eigen::Matrix4f& m)
{
    model = m;
}

void rst::rasterizer::set_view(const Eigen::Matrix4f& v)
{
    view = v;
}

void rst::rasterizer::set_projection(const Eigen::Matrix4f& p)
{
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{0, 0, 0});
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), std::numeric_limits<float>::infinity());
    }
}

rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    frame_buf.resize(w * h);
    depth_buf.resize(w * h);
}

int rst::rasterizer::get_index(int x, int y)
{
    return (height-1-y)*width + x;
}

void rst::rasterizer::set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    //old index: auto ind = point.y() + point.x() * width;
    auto ind = (height-1-point.y())*width + point.x();
    frame_buf[ind] = color;

}

// clang-format on
