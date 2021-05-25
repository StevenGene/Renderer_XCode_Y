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

void SaveBmp (std::vector<Eigen::Vector4f> &frameBuffer, int width, int height, std::string file) {
    unsigned char buf[54] = {
        'B', 'M', INT2CHAR(54 + width*height * 32), INT2CHAR (0), INT2CHAR (54), INT2CHAR (40), INT2CHAR (width), INT2CHAR (height), 1, 0, 32, 0
    };
    
    std::ofstream ofs (file, std::ios_base::out | std::ios_base::binary);
    ofs.write ((char *)buf, sizeof (buf));
    for (auto &color : frameBuffer) {
        buf[0] = (unsigned char)std::min (255, (int)(color.z() * 255));
        buf[1] = (unsigned char)std::min (255, (int)(color.y() * 255));
        buf[2] = (unsigned char)std::min (255, (int)(color.x() * 255));
        buf[3] = (unsigned char)std::min (255, (int)(color.w() * 255));
        ofs.write ((char *)buf, 4);
    }
    
} // bmp's color is bgra order



int main(int argc, const char * argv[]) {
    
    Vector3f a(1123,451,0);
    Vector3f b(31313,4563,0);
    auto c = b - a;
    auto d = a.cross(b);
    
    printf("%f,%f \n",c.x(),c.y());
    printf("%f,%f,%f \n",d.x(),d.y(),d.z());
    
    d = b.cross(a);
    printf("%f,%f,%f \n",d.x(),d.y(),d.z());
    
    printf("\n\n");

    std::vector<Eigen::Vector4f> frameBuffer(PIC_WIDTH * PIC_HEIGHT,Eigen::Vector4f(0.5,0.5,0.5,1));
    SaveBmp (frameBuffer, PIC_WIDTH, PIC_HEIGHT, "output1.bmp");

    return 0;
}

