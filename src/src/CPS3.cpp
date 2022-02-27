#include "../include/CPS3.h"
#include "../include/dof.h"
#include <iostream>
#include <iomanip>


void CPS3::calculate_Ke(){
    Eigen::Matrix<double,2,2> J;
    Eigen::Matrix<double,3,6> B;
    auto D = pid->get_mid()->D_2D_linear_continuum_mechanics;
    // Want to find elements addition to the stiffness- & load matrix
    double x1 = connectivity.at(0)->x;
    double x2 = connectivity.at(1)->x;
    double x3 = connectivity.at(2)->x;
    double y1 = connectivity.at(0)->y;
    double y2 = connectivity.at(1)->y;
    double y3 = connectivity.at(2)->y;
    J << x1-x3, y1 - y3,
            x2-x3, y2 - y3;
    B << J(1,1), 0, -J(0,1), 0, -J(1,1)+J(0,1), 0,
            0, -J(1,0), 0, J(0,0), 0, J(1,0)-J(0,0),
            -J(1,0), J(1,1), J(0,0), -J(0,1), J(1,0)-J(0,0), -J(1,1)+J(0,1);
    B = 1/J.determinant()*B; 
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
    Me *= (pid->get_mid()->get_density()*area/12.d);
}

CPS3::CPS3(unsigned int                        id,
           std::vector<Node*>                  connectivity,
           Pid*                                pid):
Element{id,connectivity,pid,ElementType::CPS3,3,3*2,6,1,2}{}