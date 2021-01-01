#include "CPS3.h"
#include "../dof.h"
#include <iostream>

const std::string CPS3::element_type = "CPS3";
CPS3::~CPS3(){}
CPS3::CPS3(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid):
    id(id),connectivity(connectivity),pid(pid){
        // add dofs to each node. can be done first now because now we know how many dofs each node should have
        for (unsigned int i = 0; i < connectivity.size(); i++)
        {
            // Check if current node already has dofs or if we need to create
            if (connectivity.at(i)->dofs.size() != ndofs/nnodes)
            {
                // create 2 dofs
                Dof x = Dof();
                Dof y = Dof();
                connectivity.at(i)->dofs.push_back(x);
                connectivity.at(i)->dofs.push_back(y);
                dofs_id.push_back(x.id);
                dofs_id.push_back(y.id);
            }
            else
            {
                // find dofs from node and add to dofs_id vector
                dofs_id.push_back(connectivity.at(i)->dofs.at(0).id);
                dofs_id.push_back(connectivity.at(i)->dofs.at(1).id);
            }
            
        }
        
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
        std::shared_ptr<Mid> mid = pid->get_mid();
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
        Eigen::Matrix<float,6,6> Me;
        // for this simple element there exists analytical expression
        Me << 2,0,1,0,1,0,
              0,2,0,1,0,1,
              1,0,2,0,1,0,
              0,1,0,2,0,1,
              1,0,1,0,2,0,
              0,1,0,1,0,2;
        
        Me *= (mid->get_density()*area/12.0);
        Element::print_element_info_to_log();
}