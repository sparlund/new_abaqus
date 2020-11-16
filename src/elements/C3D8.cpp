#include <iostream>
#include <cmath>
#include "C3D8.h"

const std::string C3D8::element_type = "C3D8";

C3D8::~C3D8(){}

C3D8::C3D8(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid):id(id),connectivity(connectivity),pid(pid){
    // add dofs to each node. can be done first now because now we know how many dofs each node should have    
    for (unsigned int i = 0; i < connectivity.size(); i++)
    {
        // Check if current node already has dofs or if we need to create
        if (connectivity.at(i)->dofs.size() != ndofs/nnodes)
        {
            // create 3 dofs
            Dof x = Dof();
            Dof y = Dof();
            Dof z = Dof();
            connectivity.at(i)->dofs.push_back(x);
            connectivity.at(i)->dofs.push_back(y);
            connectivity.at(i)->dofs.push_back(z);
            dofs_id.push_back(x.id);
            dofs_id.push_back(y.id);
            dofs_id.push_back(z.id);
        }
        else
        {
            // find dofs from node and add to dofs_id vector
            dofs_id.push_back(connectivity.at(i)->dofs.at(0).id);
            dofs_id.push_back(connectivity.at(i)->dofs.at(1).id);
            dofs_id.push_back(connectivity.at(i)->dofs.at(2).id);
        }
        
    }
    std::shared_ptr<Mid> mid = pid->get_mid();
    float v = mid->get_v();
    float E = mid->get_E();
    Eigen::Matrix<float,6,6> D;
    // Constitutive matrix (linear continuum mechanics)
    D << 1-v,     v,     v,           0,          0,          0,
           v,   1-v,     v,           0,          0,          0,
           v,     v,   1-v,           0,          0,          0,
           0,     0,     0,   (1-2*v)/2,          0,          0,
           0,     0,     0,           0,  (1-2*v)/2,          0,
           0,     0,     0,           0,          0,  (1-2*v)/2;
    D *= E/((1+v)*(1-2*v));
    // shape functions:
    // N1=1/8*(1−xhi)(1−eta)(1−my)
    // N2=1/8*(1+xhi)(1−eta)(1−my)
    // N3=1/8*(1+xhi)(1+eta)(1−my)
    // N4=1/8*(1−xhi)(1+eta)(1−my)
    // N5=1/8*(1−xhi)(1−eta)(1+my)
    // N6=1/8*(1+xhi)(1−eta)(1+my)
    // N7=1/8*(1+xhi)(1+eta)(1+my)
    // N8=1/8*(1−xhi)(1+eta)(1+my)
    
    // coord
    for (unsigned short i = 0; i < nnodes; i++)
    {
        coord(i,0) = connectivity.at(i)->x;
        coord(i,1) = connectivity.at(i)->y;
        coord(i,2) = connectivity.at(i)->z;
    }
    // gauss points same in all directions!
    xhi << -std::sqrt(0.33), std::sqrt(0.33);
    eta << -std::sqrt(0.33), std::sqrt(0.33);
    my  << -std::sqrt(0.33), std::sqrt(0.33);

    // set Ke zero
    Ke.setZero();
    // dim(J) = ndim x ndim
    Eigen::Matrix<float,3,3> J;
    Eigen::Matrix<float,6,24> B;
    
    // all weights are 1, dont need w
    // float w = 1;
    Eigen::Matrix<float,1,8> q_x;
    Eigen::Matrix<float,1,8> q_y;
    Eigen::Matrix<float,1,8> q_z;
    Eigen::Matrix<float,1,8> zero;
    zero.setZero();
    for (unsigned short i = 0; i < ngp_per_dim; i++)
    {
        for (unsigned short j = 0; j < ngp_per_dim; j++)
        {
            for (unsigned short k = 0; k < ngp_per_dim; k++)
            {
                dNdXhidEtadMy << -(1-eta(0,j))*(1-my(0,k)),(1-eta(0,j))*(1-my(0,k)),(1+eta(0,j))*(1-my(0,k)),-(1+eta(0,j))*(1-my(0,k)),-(1-eta(0,j))*(1+my(0,k)),(1-eta(0,j))*(1+my(0,k)),(1+eta(0,j))*(1+my(0,k)),-(1+eta(0,j))*(1+my(0,k)),
                                 -(1-xhi(0,i))*(1-my(0,k)),-(1+xhi(0,i))*(1-my(0,k)),(1+xhi(0,i))*(1-my(0,k)),(1-xhi(0,i))*(1-my(0,k)),-(1-xhi(0,i))*(1+my(0,k)),-(1+xhi(0,i))*(1+my(0,k)),(1+xhi(0,i))*(1+my(0,k)),(1-xhi(0,i))*(1+my(0,k)),
                                 -(1-xhi(0,i))*(1-eta(0,j)),-(1+xhi(0,i))*(1-eta(0,j)),-(1+xhi(0,i))*(1+eta(0,j)),-(1-xhi(0,i))*(1+eta(0,j)),(1-xhi(0,i))*(1-eta(0,j)),(1+xhi(0,i))*(1-eta(0,j)),(1+xhi(0,i))*(1+eta(0,j)),(1-xhi(0,i))*(1+eta(0,j));
                dNdXhidEtadMy *= 1/8;
                J = dNdXhidEtadMy*coord;
                dNdxdydz = J.inverse()*dNdXhidEtadMy;
                // make sure all values are zero
                q_x << dNdxdydz(0,0),dNdxdydz(0,1),dNdxdydz(0,2),dNdxdydz(0,3),dNdxdydz(0,4),dNdxdydz(0,5),dNdxdydz(0,6),dNdxdydz(0,7);
                q_y << dNdxdydz(1,0),dNdxdydz(1,1),dNdxdydz(1,2),dNdxdydz(1,3),dNdxdydz(1,4),dNdxdydz(1,5),dNdxdydz(1,6),dNdxdydz(1,7);
                q_z << dNdxdydz(2,0),dNdxdydz(2,1),dNdxdydz(2,2),dNdxdydz(2,3),dNdxdydz(2,4),dNdxdydz(2,5),dNdxdydz(2,6),dNdxdydz(2,7);
                // Find B
                B <<  q_x, zero, zero,
                     zero,  q_y, zero,
                     zero, zero,  q_z,
                      q_y,  q_x, zero,
                     zero,  q_z,  q_y,
                      q_z, zero, q_x;
                Ke = Ke + B.transpose()*D*B*J.determinant();
            }
        }
    }
    print_element_info_to_log();
}



