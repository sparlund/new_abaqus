#include "S2.h"
#include "../dof.h"
#include <iostream>

const std::string S2::element_type = "S2";

S2::~S2(){}
S2::S2(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid):
    id(id),connectivity(connectivity),pid(pid){
        // add dofs to each node. can be done first now because now we know how many dofs each node should have
        for (unsigned int i = 0; i < connectivity.size(); i++)
        {
            // Check if current node already has dofs or if we need to create
            if (connectivity.at(i)->dofs.size() != ndofs/nnodes)
            {
                // create 2 dofs
                std::cout << "creating 2 dofs in S2 element..\n";
                std::cout << connectivity.at(i)->dofs.size() << "\n";
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
        J << x1-x3, y1 - y3,
             x2-x3, y2 - y3;
        B << J(1,1), 0, -J(0,1), 0, -J(1,1)+J(0,1), 0,
             0, -J(1,0), 0, J(0,0), 0, J(1,0)-J(0,0),
             -J(1,0), J(1,1), J(0,0), -J(0,1), J(1,0)-J(0,0), -J(1,1)+J(0,1);
        B = 1/J.determinant()*B; 
        // TODO: plane strain
        
        float t = 1.0;
        float E = 210e3;
        float v = 0.33;
        Eigen::Matrix<float,3,3> D_temp;
        // plane stress:
        D_temp << 1, v,       0,
                  v, 1,       0,
                  0, 0, (1-v)/2;
        Eigen::Matrix<float,3,3> D = (E/( 1-(v*v) )) * D_temp;
        std::cout << "element  " << id << " created" << "\n"; 
        // // Finally compute elements contribution to stiffness matrix and load vector:
        Ke = t*0.5*J.determinant()*B.transpose()*D*B;
}
        // std::cout << "Ke:\n";
        // std::cout << Ke << "\n";
// float detJ = J.determinant();
// std::cout << "---\n";
// std::cout << "det(J):\n";
// std::cout << detJ << "\n";
// std::cout << "J:\n";
// std::cout << J << "\n";
// std::cout << "B:\n";
// std::cout << B << "\n";

// std::cout << "---\n";
// plane strain:
// D_temp << 1-v,   v,   v,          0,
//             v, 1-v,   v,          0,
//             v,   v, 1-v,          0,
//             0,   0,   0,    (1-2*v)/2;
