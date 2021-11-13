#pragma once
#include <array>


namespace Gauss
{
    namespace _3D
    {
        // location and weights for 3-by-3-by-3 gauss integration, used in C3D20
        extern const std::array<std::array<float, 3>,27> integration_points_3_by_3_by_3;    
        extern const std::array<float, 27> gauss_weights_3_by_3_by_3;
        // location and weights for 2-by-2-by-8 gauss integration, used in C3D8
        extern const std::array<std::array<float, 3>,8> integration_points_2_by_2_by_2;    
        extern const std::array<float, 8> gauss_weights_2_by_2_by_2;
        // location and weights for 4 gauss integration, used in C3D10
        extern const std::array<std::array<float, 3>,4> integration_points_4;    
        extern const std::array<float, 4> gauss_weights_4;

    } 
    namespace _2D
    {
        // location and weights for 4 gauss integration, used in CPS4
        extern const std::array<std::array<float, 2>,4> integration_points_2_by_2;    
        extern const std::array<float, 4> gauss_weights_2_by_2;
    } 
    
    

    
    
    

} 
