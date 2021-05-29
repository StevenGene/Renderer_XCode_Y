//
//  main.cpp
//  Renderer_Y
//
//  Created by St on 2021/5/25.
//  Copyright © 2021 St. All rights reserved.
//

#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include "rasterizer.hpp"
#include "global.hpp"
#include "Triangle.hpp"


#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

#define INT2CHAR_BIT(num, bit) (unsigned char)(((num) >> (bit)) & 0xff)
#define INT2CHAR(num) INT2CHAR_BIT((num),0), INT2CHAR_BIT((num),8), INT2CHAR_BIT((num),16), INT2CHAR_BIT((num),24)

using namespace Eigen;

const int PIC_WIDTH = 700;
const int PIC_HEIGHT = 700;

void SaveBmp (std::vector<Eigen::Vector3f> &frameBuffer, int width, int height, std::string file) {
    unsigned char buf[54] = {
        'B', 'M', INT2CHAR(54 + width*height * 32), INT2CHAR (0), INT2CHAR (54), INT2CHAR (40), INT2CHAR (width), INT2CHAR (height), 1, 0, 32, 0
    };
    
    std::ofstream ofs (file, std::ios_base::out | std::ios_base::binary);
    ofs.write ((char *)buf, sizeof (buf));
    for (auto &color : frameBuffer) {
        buf[0] = (unsigned char)std::min (255, (int)(color.z() * 1));
        buf[1] = (unsigned char)std::min (255, (int)(color.y() * 1));
        buf[2] = (unsigned char)std::min (255, (int)(color.x() * 1));
//        buf[3] = (unsigned char)std::min (255, (int)(color.w() * 255));
        ofs.write ((char *)buf, 4);
    }
    
} // bmp's color is bgra order

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1,0,0,-eye_pos[0],
                 0,1,0,-eye_pos[1],
                 0,0,1,-eye_pos[2],
                 0,0,0,1;

    view = translate*view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();
    return model;
}

Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio, float zNear, float zFar)
{
    float n = zNear;
    float f = zFar;
    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f P2O = Eigen::Matrix4f::Identity();
    P2O<<n, 0, 0, 0,
         0, n, 0, 0,
         0, 0, n+f,(-1)*f*n,
         0, 0, 1, 0;

    float halfRadian = eye_fov / 2.0 / 180.0 * 3.1415;
    float t = tan(halfRadian) * n;
    float b = -t;
    float r = aspect_ratio * t;
    float l = -r;




//     float halfEyeAngelRadian = eye_fov/2.0/180.0*3.1415;
//     float t = n*std::tan(halfEyeAngelRadian);
//     float r=t*aspect_ratio;
//     float l=(-1)*r;
//     float b=(-1)*t;

    Eigen::Matrix4f ortho1=Eigen::Matrix4f::Identity();
    ortho1<<2/(r-l),0,0,0,
        0,2/(t-b),0,0,
        0,0,2/(n-f),0,
        0,0,0,1;

    Eigen::Matrix4f ortho2 = Eigen::Matrix4f::Identity();
    ortho2<<1,0,0,(-1)*(r+l)/2,
        0,1,0,(-1)*(t+b)/2,
        0,0,1,(-1)*(n+f)/2,
        0,0,0,1;
    Eigen::Matrix4f Matrix_ortho = ortho1 * ortho2;
    projection = Matrix_ortho * P2O;
    return projection;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc == 2)
    {
        command_line = true;
        filename = std::string(argv[1]);
    }

    rst::rasterizer r(700, 700);

    Eigen::Vector3f eye_pos = {0,0,5};


    std::vector<Eigen::Vector3f> poses
            {
                    {2, 0, -2},
                    {0, 2, -2},
                    {-2, 0, -2},
                    {3.5, -1, -5},
                    {2.5, 1.5, -5},
                    {-1, 0.5, -5}
            };

    std::vector<Eigen::Vector3i> indexes
            {
                    {0, 1, 2},
                    {3, 4, 5}
            };

    std::vector<Eigen::Vector3f> colors
            {
                    {217.0, 238.0, 185.0},
                    {217.0, 238.0, 185.0},
                    {217.0, 238.0, 185.0},
                    {185.0, 217.0, 238.0},
                    {185.0, 217.0, 238.0},
                    {185.0, 217.0, 238.0}
            };

    auto pos_id = r.load_positions(poses);
    auto ind_id = r.load_indices(indexes);
    auto col_id = r.load_colors(colors);

    int key = 0;
    int frame_count = 0;

    if (command_line)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);
//        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
//        image.convertTo(image, CV_8UC3, 1.0f);
//        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);
//
//        cv::imwrite(filename, image);

        return 0;
    }

    while(key != 27)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);

//        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
//        image.convertTo(image, CV_8UC3, 1.0f);
//        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);
//        cv::imshow("image", image);
//        key = cv::waitKey(10);
        
//        std::vector<Eigen::Vector3f> frameBuffer(PIC_WIDTH * PIC_HEIGHT,Eigen::Vector3f(0.5,0.5,0.5));
//            SaveBmp (frameBuffer , PIC_WIDTH, PIC_HEIGHT, "output1.bmp");

        SaveBmp (r.frame_buffer() , PIC_WIDTH, PIC_HEIGHT, "output2.png");
        
        std::cout << "frame count: " << frame_count++ << '\n';
        break;
    }

    return 0;
}

