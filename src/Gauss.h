#pragma once
#include <vector>

namespace Gauss
{
    const std::vector<std::vector<float>>* get_gauss_integration_points(unsigned short ndim, unsigned short ngp);
    const std::vector<float>* get_gauss_integration_weights(unsigned short ndim, unsigned short ngp);
    namespace _3D
    {
        // location and weights for 3D 3-by-3-by-3 points gauss integration, used in C3D20
        extern const std::vector<std::vector<float>> integration_points_3_by_3_by_3;    
        extern const std::vector<float> gauss_weights_3_by_3_by_3;
        // location and weights for 3D 2-by-2-by-8 points gauss integration, used in C3D8
        extern const std::vector<std::vector<float>> integration_points_2_by_2_by_2;    
        extern const std::vector<float> gauss_weights_2_by_2_by_2;
        // location and weights for 3D 4 points gauss integration, used in C3D10
        extern const std::vector<std::vector<float>> integration_points_4;    
        extern const std::vector<float> gauss_weights_4;
    }
    namespace _2D
    {
        // location and weights for 2D 4 gauss points integration, used in CPS4
        extern const std::vector<std::vector<float>> integration_points_2_by_2;    
        extern const std::vector<float> gauss_weights_2_by_2;
    } 
} 
