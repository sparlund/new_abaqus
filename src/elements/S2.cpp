#include "S2.h"
#include <iostream>



S2::~S2(){}
S2::S2(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid):
    id(id),connectivity(connectivity),pid(pid){
        // Want to find elements addition to the stiffness- & load matrix
        float x1 = connectivity.at(0)->x;
        float x2 = connectivity.at(1)->x;
        float x3 = connectivity.at(2)->x;
        float y1 = connectivity.at(0)->y;
        float y2 = connectivity.at(1)->y;
        float y3 = connectivity.at(2)->y;

        C <<  1, x1, y1,  0,  0,  0,
              0,  0,  0,  1, x1, y1,
              1, x2, y2,  0,  0,  0,
              0,  0,  0,  1, x2, y2, 
              1, x3, y3,  0,  0,  0,  
              0,  0,  0,  1, x3, y3;
        Eigen::Matrix<float,3,3> area_temp;
        area_temp << 1, x1, y1,
                     1, x2, y2,
                     1, x3, y3;
        A = 0.5*area_temp.determinant();
        // // plane stress
        B << 0, 1, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 1,
             0, 0, 1, 0, 1, 0;
        B = B * C.inverse();
        // // // TODO: plane strain
        float t = 1.0;
        float E = 100;
        float v = 0.33;
        Eigen::Matrix<float,3,3> D_temp;
        // plane strain:
        // D_temp << 1-v,   v,   v,          0,
        //             v, 1-v,   v,          0,
        //             v,   v, 1-v,          0,
        //             0,   0,   0,    (1-2*v)/2;
        // plane stress:
        D_temp << 1, v,       0,
                  v, 1,       0,
                  0, 0, (1-v)/2;
        Eigen::Matrix<float,3,3> D = (E/( 1-(v*v) )) * D_temp;
        std::cout << D << "\n";
        // Finally compute elements contribution to stiffness matrix and load vector:
        Ke = B.transpose()*D*B*A*t;
        fe << 0,0,0,0,0,0;
              
    }
