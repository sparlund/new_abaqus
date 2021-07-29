#include "Gauss.h"


const std::vector<std::vector<float>>* Gauss::get_gauss_integration_points(unsigned short ndim, unsigned short ngp){
    if(ndim == 3){
        if(ngp==27){
            return &_3D::integration_points_3_by_3_by_3;
        }
        if(ngp == 8){
            return &_3D::integration_points_2_by_2_by_2;
        }
        if(ngp == 4){
            return &_3D::integration_points_4;
        }
    }
    else if (ndim == 2)
    {
        if(ngp == 4){
            return &_2D::integration_points_2_by_2;
        }
    }
    else{ 
        // dim == 1, not implemented yet
    }
}

const std::vector<float>* Gauss::get_gauss_integration_weights(unsigned short ndim, unsigned short ngp){
    if(ndim == 3){
        if(ngp==27){
            return &_3D::gauss_weights_3_by_3_by_3;
        }
        if(ngp == 8){
            return &_3D::gauss_weights_2_by_2_by_2;
        }
        if(ngp == 4){
            return &_3D::gauss_weights_4;
        }
    }
    else if (ndim == 2)
    {
        if(ngp == 4){
            return &_2D::gauss_weights_2_by_2;
        }
    }
    else{ 
        // dim == 1, not implemented yet
    }
}



// for i in [-sqrt(3/5), 0, -sqrt(3/5)]:
//  for j in [-sqrt(3/5), 0, -sqrt(3/5)]:
//      for k in [-sqrt(3/5), 0, -sqrt(3/5)]:
//          gp = [i,j,k]
const std::vector<std::vector<float>> Gauss::_3D::integration_points_3_by_3_by_3 = {{
                                            {-0.774596669241483f,  -0.774596669241483f,  -0.774596669241483f},
                                            {0.000000000000000f,  -0.774596669241483f,  -0.774596669241483f},
                                            {0.774596669241483f,  -0.774596669241483f,  -0.774596669241483f},
                                            {-0.774596669241483f,   0.000000000000000f,  -0.774596669241483f},
                                            {0.000000000000000f,   0.000000000000000f,  -0.774596669241483f},
                                            {0.774596669241483f,   0.000000000000000f,  -0.774596669241483f},
                                            {-0.774596669241483f,   0.774596669241483f,  -0.774596669241483f},
                                            {0.000000000000000f,   0.774596669241483f,  -0.774596669241483f},
                                            {0.774596669241483f,   0.774596669241483f,  -0.774596669241483f},
                                            {-0.774596669241483f,  -0.774596669241483f,   0.000000000000000f},
                                            {0.000000000000000f,  -0.774596669241483f,   0.000000000000000f},
                                            {0.774596669241483f,  -0.774596669241483f,   0.000000000000000f},
                                            {-0.774596669241483f,   0.000000000000000f,   0.000000000000000f},
                                            {0.000000000000000f,   0.000000000000000f,   0.000000000000000f},
                                            {0.774596669241483f,   0.000000000000000f,   0.000000000000000f},
                                            {-0.774596669241483f,   0.774596669241483f,   0.000000000000000f},
                                            {0.000000000000000f,   0.774596669241483f,   0.000000000000000f},
                                            {0.774596669241483f,   0.774596669241483f,   0.000000000000000f},
                                            {-0.774596669241483f,  -0.774596669241483f,   0.774596669241483f},
                                            {0.000000000000000f,  -0.774596669241483f,   0.774596669241483f},
                                            {0.774596669241483f,  -0.774596669241483f,   0.774596669241483f},
                                            {-0.774596669241483f,   0.000000000000000f,   0.774596669241483f},
                                            {0.000000000000000f,   0.000000000000000f,   0.774596669241483f},
                                            {0.774596669241483f,   0.000000000000000f,   0.774596669241483f},
                                            {-0.774596669241483f,   0.774596669241483f,   0.774596669241483f},
                                            {0.000000000000000f,   0.774596669241483f,   0.774596669241483f},
                                            {0.774596669241483f,   0.774596669241483f,   0.774596669241483f}   }};   
// for i in [5/9, 8/9, 5/9]:
//  for j in [5/9, 8/9, 5/9]:
//      for k in [5/9, 8/9, 5/9]:
//          w = i*j*k
const std::vector<float> Gauss::_3D::gauss_weights_3_by_3_by_3 = {0.171467764060357f,   0.274348422496571f,   0.171467764060357f,
                                                                     0.274348422496571f,   0.438957475994513f,   0.274348422496571f,
                                                                     0.171467764060357f,   0.274348422496571f,   0.171467764060357f,
                                                                     0.274348422496571f,   0.438957475994513f,   0.274348422496571f,
                                                                     0.438957475994513f,   0.702331961591221f,   0.438957475994513f,
                                                                     0.274348422496571f,   0.438957475994513f,   0.274348422496571f,
                                                                     0.171467764060357f,   0.274348422496571f,   0.171467764060357f,
                                                                     0.274348422496571f,   0.438957475994513f,   0.274348422496571f,
                                                                     0.171467764060357f,   0.274348422496571f,   0.171467764060357f};



// for i in [-np.sqrt(1/3), np.sqrt(1/3)]:
//  for j in [-np.sqrt(1/3), np.sqrt(1/3)]
//      for j in [-np.sqrt(1/3), np.sqrt(1/3)]
//          gp = [i,j,k]
const std::vector<std::vector<float>> Gauss::_3D::integration_points_2_by_2_by_2 = {{
                                                                      {-0.5773502691896257f, -0.5773502691896257f, -0.5773502691896257f},
                                                                      {-0.5773502691896257f, -0.5773502691896257f, 0.5773502691896257f},
                                                                      {-0.5773502691896257f, 0.5773502691896257f, -0.5773502691896257f},
                                                                      {-0.5773502691896257f, 0.5773502691896257f, 0.5773502691896257f},
                                                                      {0.5773502691896257f, -0.5773502691896257f, -0.5773502691896257f},
                                                                      {0.5773502691896257f, -0.5773502691896257f, 0.5773502691896257f},
                                                                      {0.5773502691896257f, 0.5773502691896257f, -0.5773502691896257f},
                                                                      {0.5773502691896257f, 0.5773502691896257f, 0.5773502691896257f}}};    
// w same for all points
const std::vector<float> Gauss::_3D::gauss_weights_2_by_2_by_2 = {1,1,1,1,1,1,1,1};

// a = (5-sqrt(5))/20
// b = (5 + 3*sqrt(5))/20
// gp = [a a a;
//       a a b;
//       a b a;
//       b a a]
// w = [1 1 1 1]/24
const std::vector<std::vector<float>> Gauss::_3D::integration_points_4 = {{
                                                                  {0.138196601125011f,   0.138196601125011f,   0.138196601125011f},
                                                                  {0.138196601125011f,   0.138196601125011f,   0.585410196624969f},
                                                                  {0.138196601125011f,   0.585410196624969f,   0.138196601125011f},
                                                                  {0.585410196624969f,   0.138196601125011f,   0.138196601125011f}}};
const std::vector<float> Gauss::_3D::gauss_weights_4 = {0.0416666666666667f,   0.0416666666666667f,   0.0416666666666667f,   0.0416666666666667f};

// for i in [-sqrt(1/3), sqrt(1/3)]:
//  for j in [-sqrt(1/3), sqrt(1/3)]:
//      gp = [i,j]
// w = [1 1 1 1]
const std::vector<std::vector<float>> Gauss::_2D::integration_points_2_by_2 = {{{-0.5773502691896257f, -0.5773502691896257f},
                                                                                   {-0.5773502691896257f, 0.5773502691896257f},
                                                                                   {0.5773502691896257f, -0.5773502691896257f},
                                                                                   {0.5773502691896257f, 0.5773502691896257f}}};
 
const std::vector<float> Gauss::_2D::gauss_weights_2_by_2 = {1.0f, 1.0f, 1.0f, 1.0f};



