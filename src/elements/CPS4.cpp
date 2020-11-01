#include <iostream>
#include <cmath>
#include "CPS4.h"

const std::string CPS4::element_type = "CPS4";

CPS4::~CPS4(){}

CPS4::CPS4(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid):id(id),connectivity(connectivity),pid(pid){
    // add dofs to each node. can be done first now because now we know how many dofs each node should have    
    for (unsigned int i = 0; i < connectivity.size(); i++)
    {
        // Check if current node already has dofs or if we need to create
        if (connectivity.at(i)->dofs.size() != ndofs/nnodes)
        {
            // create 2 dofs
            std::cout << "creating 2 dofs in CPS4 element with id=" << id << std::endl;
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
    float v = 0.33;
    float E = 210e3;
    float t = 1.0f;
    D << 1, v, 0,
         v, 1, 0,
         0, 0, 0.5*(1-v);
    D *= E/(1-(v*v));
    xhi << -std::sqrt(0.33),  std::sqrt(0.33), -std::sqrt(0.33), std::sqrt(0.33);
    eta << -std::sqrt(0.33), -std::sqrt(0.33),  std::sqrt(0.33), std::sqrt(0.33);
    for (unsigned char i = 0; i < ngp; i++)
    {            
        // shape functions
        N(0,i) = 0.25*(1-xhi(i))*(1-eta(i)),
        N(1,i) = 0.25*(1+xhi(i))*(1-eta(i)),
        N(2,i) = 0.25*(1+xhi(i))*(1+eta(i)),
        N(3,i) = 0.25*(1-xhi(i))*(1+eta(i));
        // shape functions derivatives wrt xhi
        dNdXhi(0,i) = -0.25*(1-eta(i));
        dNdXhi(1,i) =  0.25*(1-eta(i));
        dNdXhi(2,i) =  0.25*(1+eta(i));
        dNdXhi(3,i) = -0.25*(1+eta(i));
        // shape functions derivatives wrt eta
        dNdEta(0,i) = -0.25*(1-xhi(i));
        dNdEta(1,i) = -0.25*(1+xhi(i)); 
        dNdEta(2,i) =  0.25*(1+xhi(i)); 
        dNdEta(3,i) =  0.25*(1-xhi(i)); 
    }

    // create coord matrix needed to find Jacobian
    for (unsigned char i = 0; i < nnodes; i++)
    {
        coord(i,0) = connectivity.at(i)->x;
        coord(i,1) = connectivity.at(i)->y;
    }
    // Find jacobian J and det(J) for each gauss point
    float w = 1;
    for (unsigned int i = 0; i < ngp; i++)
    {
        Eigen::Matrix<float,2,4> dNdXhidEta;
        dNdXhidEta.row(0) = dNdXhi.col(i).transpose();
        dNdXhidEta.row(1) = dNdEta.col(i).transpose();
        J << dNdXhidEta*coord;
        detJ.push_back(J.determinant());
        dNdxdy = J.inverse()*dNdXhidEta;
        B <<    dNdxdy(0,0),            0, dNdxdy(0,1),           0, dNdxdy(0,2),           0, dNdxdy(0,3),           0,
                        0, dNdxdy(1,0),           0, dNdxdy(1,1),           0, dNdxdy(1,2),           0, dNdxdy(1,3),
                dNdxdy(1,0), dNdxdy(0,0),dNdxdy(1,1), dNdxdy(0,1),dNdxdy(1,2), dNdxdy(0,2),dNdxdy(1,3), dNdxdy(0,3); 
        Ke += w*detJ.at(i)*t*B.transpose()*D*B;
    }
}



