#include "../include/C3D8.h"
#include <iostream>
#include <stdlib.h>

void C3D8::calculate_Ke(){    
    setup_coord();
    // dim(J) = ndim x ndim
    Eigen::Matrix<double,3,3> J;
    Eigen::Matrix<double,6,24> B;
    
    Eigen::Matrix<double,1,8> dNdXhi;
    Eigen::Matrix<double,1,8> dNdEta;
    Eigen::Matrix<double,1,8> dNdMy;
    // dim(dNdxdydz) = dim(dNdXhidEtadMy) = ndim x nnodes
    Eigen::Matrix<double,3,8> dNdxdydz;
    Eigen::Matrix<double,3,8> dNdXhidEtadMy;

    double dN1dXhi,dN2dXhi,dN3dXhi,dN4dXhi,dN5dXhi,dN6dXhi,dN7dXhi,dN8dXhi;
    double dN1dEta,dN2dEta,dN3dEta,dN4dEta,dN5dEta,dN6dEta,dN7dEta,dN8dEta;
    double dN1dMy,dN2dMy,dN3dMy,dN4dMy,dN5dMy,dN6dMy,dN7dMy,dN8dMy;
    double xhi,eta,my,w;
    auto D = pid->get_mid()->D_3D_linear_continuum_mechanics;
    for (unsigned int i = 0; i < gauss_points->size(); i++)
    {
        xhi     = gauss_points->at(i).at(0);
        eta     = gauss_points->at(i).at(1);
        my      = gauss_points->at(i).at(2);
        w       = gauss_weights->at(i);
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
        // detJ is not used but can be good for fault tracing later
        detJ.push_back(J.determinant());
        // find shape functions derivate matrix wrt x, y & z
        dNdxdydz = J.inverse()*dNdXhidEtadMy;
        // construct B matrix
        B << dNdxdydz(0,0),            0,0,dNdxdydz(0,1),0,0,dNdxdydz(0,2),0,0,dNdxdydz(0,3),0,0,dNdxdydz(0,4),0,0,dNdxdydz(0,5),0,0,dNdxdydz(0,6),0,0,dNdxdydz(0,7),0,0,
                         0,dNdxdydz(1,0),0,0,dNdxdydz(1,1),0,0,dNdxdydz(1,2),0,0,dNdxdydz(1,3),0,0,dNdxdydz(1,4),0,0,dNdxdydz(1,5),0,0,dNdxdydz(1,6),0,0,dNdxdydz(1,7),0,
                         0,            0,dNdxdydz(2,0),0,0,dNdxdydz(2,1),0,0,dNdxdydz(2,2),0,0,dNdxdydz(2,3),0,0,dNdxdydz(2,4),0,0,dNdxdydz(2,5),0,0,dNdxdydz(2,6),0,0,dNdxdydz(2,7),
             dNdxdydz(1,0),dNdxdydz(0,0),0,dNdxdydz(1,1),dNdxdydz(0,1),0,dNdxdydz(1,2),dNdxdydz(0,2),0,dNdxdydz(1,3),dNdxdydz(0,3),0,dNdxdydz(1,4),dNdxdydz(0,4),0,dNdxdydz(1,5),dNdxdydz(0,5),0,dNdxdydz(1,6),dNdxdydz(0,6),0,dNdxdydz(1,7),dNdxdydz(0,7),0,
                         0,dNdxdydz(2,0),dNdxdydz(0,0),0,dNdxdydz(2,1),dNdxdydz(0,1),0,dNdxdydz(2,2),dNdxdydz(0,2),0,dNdxdydz(2,3),dNdxdydz(0,3),0,dNdxdydz(2,4),dNdxdydz(0,4),0,dNdxdydz(2,5),dNdxdydz(0,5),0,dNdxdydz(2,6),dNdxdydz(0,6),0,dNdxdydz(2,7),dNdxdydz(0,7),
             dNdxdydz(2,0),            0,dNdxdydz(0,0),dNdxdydz(2,1),0,dNdxdydz(0,1),dNdxdydz(2,2),0,dNdxdydz(0,2),dNdxdydz(2,3),0,dNdxdydz(0,3),dNdxdydz(2,4),0,dNdxdydz(0,4),dNdxdydz(2,5),0,dNdxdydz(0,5),dNdxdydz(2,6),0,dNdxdydz(0,6),dNdxdydz(2,7),0,dNdxdydz(0,7);
        // compute Gauss point contribution to element Ke matrix
        Ke = Ke + w*(B.transpose()*D*B*detJ.back());
    }
}
void C3D8::calculate_Me(){
    setup_coord();
    // size(global_N_matrix) = element ndof*element ndof
    Eigen::Matrix<double,3,24> N;
    double N1,N2,N3,N4,N5,N6,N7,N8;
    double xhi,eta,my,w;
    for (unsigned int i = 0; i < gauss_points->size(); i++)
    {
        xhi = gauss_points->at(i).at(0);
        eta = gauss_points->at(i).at(1);
        my  = gauss_points->at(i).at(2);
        w   = gauss_weights->at(i);
        // shape functions, 1/8=0.125
        N1  = (1-xhi)*(1-eta)*(1-my)*0.125;
        N2  = (1+xhi)*(1-eta)*(1-my)*0.125; 
        N3  = (1+xhi)*(1+eta)*(1-my)*0.125; 
        N4  = (1-xhi)*(1+eta)*(1-my)*0.125; 
        N5  = (1-xhi)*(1-eta)*(1+my)*0.125;
        N6  = (1+xhi)*(1-eta)*(1+my)*0.125;
        N7  = (1+xhi)*(1+eta)*(1+my)*0.125;
        N8  = (1-xhi)*(1+eta)*(1+my)*0.125;
        // § 8.5. The Mass Matrix, chapter 8
        N <<   N1,0.d,0.d,N2,0.d,0.d,N3,0.d,0.d,N4,0.d,0.d,N5,0.d,0.d,N6,0.d,0.d,N7,0.d,0.d,N8,0.d,0.d,
             0.d,  N1,0.d,0.d,N2,0.d,0.d,N3,0.d,0.d,N4,0.d,0.d,N5,0.d,0.d,N6,0.d,0.d,N7,0.d,0.d,N8,0.d,
             0.d,0.d,N1,0.d,0.d,N2,0.d,0.d,N3,0.d,0.d,N4,0.d,0.d,N5,0.d,0.d,N6,0.d,0.d,N7,0.d,0.d,N8;
        Me = Me + w*(pid->get_mid()->get_density()*N.transpose()*N*detJ.at(i));
    }

}

C3D8::C3D8(unsigned int                        id,
           std::vector<Node*>                  connectivity,
           Pid*                                pid):
Element{id,connectivity,pid,ElementType::C3D8,8,24,12,8,3}{}
