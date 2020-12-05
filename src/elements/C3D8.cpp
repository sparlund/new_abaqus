#include <iostream>
#include "C3D8.h"
#include <stdlib.h>

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
    // std::cout << D << std::endl;
    
    // coord system is located in the middle of 
    // the element, as weights are from -1 to 1
    for (unsigned short i = 0; i < nnodes; i++)
    {
        coord(i,0) = connectivity.at(i)->x;
        coord(i,1) = connectivity.at(i)->y;
        coord(i,2) = connectivity.at(i)->z;
    }
    // std::cout << "coord=" << coord << std::endl;

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
    Eigen::Matrix<float,1,8> N;
    // size(global_N_matrix) = element ndof*element ndof
    Eigen::Matrix<float,3,24> global_N_matrix;
    Eigen::Matrix<float,1,8> dNdXhi;
    Eigen::Matrix<float,1,8> dNdEta;
    Eigen::Matrix<float,1,8> dNdMy;
    // dim(dNdxdydz) = dim(dNdXhidEtadMy) = ndim x nnodes
    Eigen::Matrix<float,3,8> dNdxdydz;
    Eigen::Matrix<float,3,8> dNdXhidEtadMy;

    // init Ke & Me to zero
    Ke.setZero();
    Me.setZero();
    // gauss points same in all directions!
    // 1/sqrt(3) = 0.57735026919
    xhi << -0.57735,-0.57735,-0.57735,-0.57735, 0.57735, 0.57735,0.57735,0.57735;
    eta << -0.57735,-0.57735, 0.57735,0.57735,-0.57735,-0.57735,0.57735,0.57735;
    my  << -0.57735, 0.57735,-0.57735,0.57735,-0.57735, 0.57735,-0.57735,0.57735;
    for (unsigned int i = 0; i < ngp; i++)
    {
        // shape functions, 1/8=0.125
        N(0) = (1-xhi(i))*(1-eta(i))*(1-my(i))*0.125;
        N(1) = (1+xhi(i))*(1-eta(i))*(1-my(i))*0.125; 
        N(2) = (1+xhi(i))*(1+eta(i))*(1-my(i))*0.125; 
        N(3) = (1-xhi(i))*(1+eta(i))*(1-my(i))*0.125; 
        N(4) = (1-xhi(i))*(1-eta(i))*(1+my(i))*0.125;
        N(5) = (1+xhi(i))*(1-eta(i))*(1+my(i))*0.125;
        N(6) = (1+xhi(i))*(1+eta(i))*(1+my(i))*0.125;
        N(7) = (1-xhi(i))*(1+eta(i))*(1+my(i))*0.125;
        // derive shape functions wrt xhi
        dNdXhi(0) = -(1-eta(i))*(1-my(i))*0.125;
        dNdXhi(1) =  (1-eta(i))*(1-my(i))*0.125;
        dNdXhi(2) =  (1+eta(i))*(1-my(i))*0.125; 
        dNdXhi(3) = -(1+eta(i))*(1-my(i))*0.125; 
        dNdXhi(4) = -(1-eta(i))*(1+my(i))*0.125; 
        dNdXhi(5) =  (1-eta(i))*(1+my(i))*0.125; 
        dNdXhi(6) =  (1+eta(i))*(1+my(i))*0.125; 
        dNdXhi(7) = -(1+eta(i))*(1+my(i))*0.125; 
        // derive shape functions wrt eta
        dNdEta(0) = -(1-xhi(i))*(1-my(i))*0.125;
        dNdEta(1) = -(1+xhi(i))*(1-my(i))*0.125;
        dNdEta(2) =  (1+xhi(i))*(1-my(i))*0.125;
        dNdEta(3) =  (1-xhi(i))*(1-my(i))*0.125;
        dNdEta(4) = -(1-xhi(i))*(1+my(i))*0.125;
        dNdEta(5) = -(1+xhi(i))*(1+my(i))*0.125;
        dNdEta(6) =  (1+xhi(i))*(1+my(i))*0.125;
        dNdEta(7) =  (1-xhi(i))*(1+my(i))*0.125;
        // derive shape functions wrt my
        dNdMy(0)  = -(1-xhi(i))*(1-eta(i))*0.125;
        dNdMy(1)  = -(1+xhi(i))*(1-eta(i))*0.125;
        dNdMy(2)  = -(1+xhi(i))*(1+eta(i))*0.125;
        dNdMy(3)  = -(1-xhi(i))*(1+eta(i))*0.125;
        dNdMy(4)  =  (1-xhi(i))*(1-eta(i))*0.125;
        dNdMy(5)  =  (1+xhi(i))*(1-eta(i))*0.125;
        dNdMy(6)  =  (1+xhi(i))*(1+eta(i))*0.125;
        dNdMy(7)  =  (1-xhi(i))*(1+eta(i))*0.125;
        // find Jacobian for each Gauss point
        dNdXhidEtadMy.row(0) = dNdXhi;
        dNdXhidEtadMy.row(1) = dNdEta;
        dNdXhidEtadMy.row(2) = dNdMy;
        J = dNdXhidEtadMy*coord;
        // find shape functions derivate matrix wrt x, y & z
        dNdxdydz = J.inverse()*dNdXhidEtadMy;
        // construct B matrix
        B << dNdxdydz(0,0),0,0,dNdxdydz(0,1),0,0,dNdxdydz(0,2),0,0,dNdxdydz(0,3),0,0,dNdxdydz(0,4),0,0,dNdxdydz(0,5),0,0,dNdxdydz(0,6),0,0,dNdxdydz(0,7),0,0,
             0,dNdxdydz(1,0),0,0,dNdxdydz(1,1),0,0,dNdxdydz(1,2),0,0,dNdxdydz(1,3),0,0,dNdxdydz(1,4),0,0,dNdxdydz(1,5),0,0,dNdxdydz(1,6),0,0,dNdxdydz(1,7),0,
             0,0,dNdxdydz(2,0),0,0,dNdxdydz(2,1),0,0,dNdxdydz(2,2),0,0,dNdxdydz(2,3),0,0,dNdxdydz(2,4),0,0,dNdxdydz(2,5),0,0,dNdxdydz(2,6),0,0,dNdxdydz(2,7),
             dNdxdydz(1,0),dNdxdydz(0,0),0,dNdxdydz(1,1),dNdxdydz(0,1),0,dNdxdydz(1,2),dNdxdydz(0,2),0,dNdxdydz(1,3),dNdxdydz(0,3),0,dNdxdydz(1,4),dNdxdydz(0,4),0,dNdxdydz(1,5),dNdxdydz(0,5),0,dNdxdydz(1,6),dNdxdydz(0,6),0,dNdxdydz(1,7),dNdxdydz(0,7),0,
             0,dNdxdydz(2,0),dNdxdydz(0,0),0,dNdxdydz(2,1),dNdxdydz(0,1),0,dNdxdydz(2,2),dNdxdydz(0,2),0,dNdxdydz(2,3),dNdxdydz(0,3),0,dNdxdydz(2,4),dNdxdydz(0,4),0,dNdxdydz(2,5),dNdxdydz(0,5),0,dNdxdydz(2,6),dNdxdydz(0,6),0,dNdxdydz(2,7),dNdxdydz(0,7),
             dNdxdydz(2,0),0,dNdxdydz(0,0),dNdxdydz(2,1),0,dNdxdydz(0,1),dNdxdydz(2,2),0,dNdxdydz(0,2),dNdxdydz(2,3),0,dNdxdydz(0,3),dNdxdydz(2,4),0,dNdxdydz(0,4),dNdxdydz(2,5),0,dNdxdydz(0,5),dNdxdydz(2,6),0,dNdxdydz(0,6),dNdxdydz(2,7),0,dNdxdydz(0,7);
        // compute Gauss point contribution to element Ke matrix
        // § 8.5. The Mass Matrix, chapter 8
        global_N_matrix << N(0),0,0,N(1),0,0,N(2),0,0,N(3),0,0,N(4),0,0,N(5),0,0,N(6),0,0,N(7),0,0,
                           0,N(0),0,0,N(1),0,0,N(2),0,0,N(3),0,0,N(4),0,0,N(5),0,0,N(6),0,0,N(7),0,
                           0,0,N(0),0,0,N(1),0,0,N(2),0,0,N(3),0,0,N(4),0,0,N(5),0,0,N(6),0,0,N(7);
        Ke = Ke + (B.transpose()*D*B*J.determinant());
        Me = Me + mid->get_density()*global_N_matrix.transpose()*global_N_matrix;
    }
    print_element_info_to_log();
}
