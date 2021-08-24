#include <iostream>
#include "C3D8.h"
#include <stdlib.h>

const std::string C3D8::element_type = "C3D8";

void C3D8::calculate_Ke(){
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
    
    Eigen::Matrix<float,1,8> dNdXhi;
    Eigen::Matrix<float,1,8> dNdEta;
    Eigen::Matrix<float,1,8> dNdMy;
    // dim(dNdxdydz) = dim(dNdXhidEtadMy) = ndim x nnodes
    Eigen::Matrix<float,3,8> dNdxdydz;
    Eigen::Matrix<float,3,8> dNdXhidEtadMy;

    // init Ke & Me to zero
    Ke.setZero();
    float dN1dXhi,dN2dXhi,dN3dXhi,dN4dXhi,dN5dXhi,dN6dXhi,dN7dXhi,dN8dXhi;
    float dN1dEta,dN2dEta,dN3dEta,dN4dEta,dN5dEta,dN6dEta,dN7dEta,dN8dEta;
    float dN1dMy,dN2dMy,dN3dMy,dN4dMy,dN5dMy,dN6dMy,dN7dMy,dN8dMy;
    float xhi,eta,my,w;
    for (unsigned int i = 0; i < gauss_points->size(); i++)
    {
        xhi = gauss_points->at(i).at(0);
        eta = gauss_points->at(i).at(1);
        my  = gauss_points->at(i).at(2);
        w   = gauss_weights->at(i);
        // derivates of shape functions wrt xhi
        dN1dXhi = -(1-eta)*(1-my)*0.125;
        dN2dXhi =  (1-eta)*(1-my)*0.125;
        dN3dXhi =  (1+eta)*(1-my)*0.125; 
        dN4dXhi = -(1+eta)*(1-my)*0.125; 
        dN5dXhi = -(1-eta)*(1+my)*0.125; 
        dN6dXhi =  (1-eta)*(1+my)*0.125; 
        dN7dXhi =  (1+eta)*(1+my)*0.125; 
        dN8dXhi = -(1+eta)*(1+my)*0.125; 
        // derivates of shape functions wrt eta
        dN1dEta = -(1-xhi)*(1-my)*0.125;
        dN2dEta = -(1+xhi)*(1-my)*0.125;
        dN3dEta =  (1+xhi)*(1-my)*0.125;
        dN4dEta =  (1-xhi)*(1-my)*0.125;
        dN5dEta = -(1-xhi)*(1+my)*0.125;
        dN6dEta = -(1+xhi)*(1+my)*0.125;
        dN7dEta =  (1+xhi)*(1+my)*0.125;
        dN8dEta =  (1-xhi)*(1+my)*0.125;
        // derivates of shape functions wrt my
        dN1dMy  = -(1-xhi)*(1-eta)*0.125;
        dN2dMy  = -(1+xhi)*(1-eta)*0.125;
        dN3dMy  = -(1+xhi)*(1+eta)*0.125;
        dN4dMy  = -(1-xhi)*(1+eta)*0.125;
        dN5dMy  =  (1-xhi)*(1-eta)*0.125;
        dN6dMy  =  (1+xhi)*(1-eta)*0.125;
        dN7dMy  =  (1+xhi)*(1+eta)*0.125;
        dN8dMy  =  (1-xhi)*(1+eta)*0.125;
        // find Jacobian for each Gauss point
        dNdXhidEtadMy.row(0) << dN1dXhi,dN2dXhi,dN3dXhi,dN4dXhi,dN5dXhi,dN6dXhi,dN7dXhi,dN8dXhi;
        dNdXhidEtadMy.row(1) << dN1dEta,dN2dEta,dN3dEta,dN4dEta,dN5dEta,dN6dEta,dN7dEta,dN8dEta;
        dNdXhidEtadMy.row(2) << dN1dMy, dN2dMy, dN3dMy, dN4dMy, dN5dMy, dN6dMy, dN7dMy, dN8dMy;
        J = dNdXhidEtadMy*coord;
        detJ.push_back(J.determinant());
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
        Ke = Ke + w*(B.transpose()*D*B*detJ.back());
    }
}
void C3D8::calculate_Me(){
    // size(global_N_matrix) = element ndof*element ndof
    Eigen::Matrix<float,3,24> N;
    Me.setZero();
    float N1,N2,N3,N4,N5,N6,N7,N8;
    float xhi,eta,my,w;
    for (unsigned int i = 0; i < gauss_points->size(); i++)
    {
        xhi = gauss_points->at(i).at(0);
        eta = gauss_points->at(i).at(1);
        my  = gauss_points->at(i).at(2);
        w   = gauss_weights->at(i);
        // shape functions, 1/8=0.125
        N1 = (1-xhi)*(1-eta)*(1-my)*0.125;
        N2 = (1+xhi)*(1-eta)*(1-my)*0.125; 
        N3 = (1+xhi)*(1+eta)*(1-my)*0.125; 
        N4 = (1-xhi)*(1+eta)*(1-my)*0.125; 
        N5 = (1-xhi)*(1-eta)*(1+my)*0.125;
        N6 = (1+xhi)*(1-eta)*(1+my)*0.125;
        N7 = (1+xhi)*(1+eta)*(1+my)*0.125;
        N8 = (1-xhi)*(1+eta)*(1+my)*0.125;
        // § 8.5. The Mass Matrix, chapter 8
        N << N1,0.0f,0.0f,N2,0.0f,0.0f,N3,0.0f,0.0f,N4,0.0f,0.0f,N5,0.0f,0.0f,N6,0.0f,0.0f,N7,0.0f,0.0f,N8,0.0f,0.0f,
             0.0f,N1,0.0f,0.0f,N2,0.0f,0.0f,N3,0.0f,0.0f,N4,0.0f,0.0f,N5,0.0f,0.0f,N6,0.0f,0.0f,N7,0.0f,0.0f,N8,0.0f,
             0.0f,0.0f,N1,0.0f,0.0f,N2,0.0f,0.0f,N3,0.0f,0.0f,N4,0.0f,0.0f,N5,0.0f,0.0f,N6,0.0f,0.0f,N7,0.0f,0.0f,N8;
        Me = Me + w*(pid->get_mid()->get_density()*N.transpose()*N*detJ.at(i));
    }

}

C3D8::C3D8(unsigned int                        id,
           std::vector<std::shared_ptr<Node>>  connectivity,
           std::shared_ptr<Pid>                pid,
           const unsigned short                nnodes,
           const unsigned short                ndofs,
           const unsigned short                vtk_identifier,
           const unsigned short                ngp,
           const unsigned short                dimensions):
Element{id,connectivity,pid,nnodes,ndofs,vtk_identifier,ngp,dimensions}{
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
            // put Dof object itself in list of Dofs for element
            connectivity.at(i)->dofs.push_back(x);
            connectivity.at(i)->dofs.push_back(y);
            connectivity.at(i)->dofs.push_back(z);
            // put indiviual dof id's in a list for easy access?
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
    print_element_info_to_log();
}
