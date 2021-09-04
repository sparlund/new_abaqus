#include "CPS3.h"
#include "../dof.h"
#include <iostream>


void CPS3::calculate_Ke(){
    Eigen::Matrix<float,2,2> J;
    Eigen::Matrix<float,3,6> B;
    // Want to find elements addition to the stiffness- & load matrix
    float x1 = connectivity.at(0)->x;
    float x2 = connectivity.at(1)->x;
    float x3 = connectivity.at(2)->x;
    float y1 = connectivity.at(0)->y;
    float y2 = connectivity.at(1)->y;
    float y3 = connectivity.at(2)->y;
    // create coord matrix needed to find area
    coord << 1, x1, y1,
                1, x2, y2,
                1, x3, y3;
    area = 0.5*coord.determinant();
    J << x1-x3, y1 - y3,
            x2-x3, y2 - y3;
    B << J(1,1), 0, -J(0,1), 0, -J(1,1)+J(0,1), 0,
            0, -J(1,0), 0, J(0,0), 0, J(1,0)-J(0,0),
            -J(1,0), J(1,1), J(0,0), -J(0,1), J(1,0)-J(0,0), -J(1,1)+J(0,1);
    B = 1/J.determinant()*B; 
    Eigen::Matrix<float,2,6> global_N_matrix;
    // global_N_matrix << 1 x1
    // TODO: plane strain
    Mid* mid = pid->get_mid();
    float v = mid->get_v();
    float E = mid->get_E();        
    float t = 1.0;
    Eigen::Matrix<float,3,3> D_temp;
    // plane stress:
    D_temp << 1, v,       0,
                v, 1,       0,
                0, 0, (1-v)/2;
    Eigen::Matrix<float,3,3> D = (E/( 1-(v*v) )) * D_temp;
    // // Finally compute elements contribution to stiffness matrix and load vector:
    Ke = t*0.5*J.determinant()*B.transpose()*D*B;
}
void CPS3::calculate_Me(){
    // for this simple element there exists analytical expression
    Me << 2,0,1,0,1,0,
          0,2,0,1,0,1,
          1,0,2,0,1,0,
          0,1,0,2,0,1,
          1,0,1,0,2,0,
          0,1,0,1,0,2;   
    Me *= (pid->get_mid()->get_density()*area/12.0f);
}

CPS3::CPS3(unsigned int                        id,
           std::vector<Node*>  connectivity,
           Pid*                pid,
           const unsigned short                nnodes,
           const unsigned short                ndofs,
           const unsigned short                vtk_identifier,
           const unsigned short                ngp,
           const unsigned short                dimensions,
           std::string                         element_type):
Element{id,connectivity,pid,nnodes,ndofs,vtk_identifier,ngp,dimensions,element_type}{}