#include <iostream>
#include <cmath>
#include "CPS4.h"


void CPS4::calculate_Ke(){
    Mid* mid = pid->get_mid();
    float v = mid->get_v();
    float E = mid->get_E();
    Eigen::Matrix<float,3,3> D;
    D << 1, v, 0,
         v, 1, 0,
         0, 0, 0.5*(1-v);
    D *= E/(1-(v*v));
    // add dofs to each node. can be done first now because now we know how many dofs each node should have    
    setup_dofs();
     // create coord matrix needed to find Jacobian
    setup_coord();
    Eigen::Matrix<float,2,4> dNdXhidEta;
    Eigen::Matrix<float,2,4> dNdxdy;
    Eigen::Matrix<float,2,2> J;
    Eigen::Matrix<float,3,8> B;
    float xhi,eta,w;
    float dN1dXhi,dN2dXhi,dN3dXhi,dN4dXhi;
    float dN1dEta,dN2dEta,dN3dEta,dN4dEta;
    for (unsigned char i = 0; i < gauss_points->size(); i++)
    {            
        xhi = gauss_points->at(i).at(0);
        eta = gauss_points->at(i).at(1);
        w = gauss_weights->at(i);
        // shape functions derivatives wrt xhi
        dN1dXhi = -0.25*(1-eta);
        dN2dXhi =  0.25*(1-eta);
        dN3dXhi =  0.25*(1+eta);
        dN4dXhi = -0.25*(1+eta);
        // shape functions derivatives wrt eta
        dN1dEta = -0.25*(1-xhi);
        dN2dEta = -0.25*(1+xhi); 
        dN3dEta =  0.25*(1+xhi); 
        dN4dEta =  0.25*(1-xhi); 
        // 
        dNdXhidEta.row(0) << dN1dXhi,dN2dXhi,dN3dXhi,dN4dXhi;
        dNdXhidEta.row(1) << dN1dEta,dN2dEta,dN3dEta,dN4dEta;
        J = dNdXhidEta*coord;
        detJ.push_back(J.determinant());
        dNdxdy = J.inverse()*dNdXhidEta;
        B <<    dNdxdy(0,0),            0, dNdxdy(0,1),           0, dNdxdy(0,2),           0, dNdxdy(0,3),           0,
                        0, dNdxdy(1,0),           0, dNdxdy(1,1),           0, dNdxdy(1,2),           0, dNdxdy(1,3),
                dNdxdy(1,0), dNdxdy(0,0),dNdxdy(1,1), dNdxdy(0,1),dNdxdy(1,2), dNdxdy(0,2),dNdxdy(1,3), dNdxdy(0,3); 
        
        Ke = Ke + w*(B.transpose()*D*B*detJ.back());
    }
};
void CPS4::calculate_Me(){
    Eigen::Matrix<float,2,8> N;
    float xhi,eta,w;
    float N1,N2,N3,N4;
    for (unsigned char i = 0; i < gauss_points->size(); i++)
    {            
        xhi = gauss_points->at(i).at(0);
        eta = gauss_points->at(i).at(1);
        w = gauss_weights->at(i);
        // shape functions
        N1 = 0.25*(1-xhi)*(1-eta),
        N2 = 0.25*(1+xhi)*(1-eta),
        N3 = 0.25*(1+xhi)*(1+eta),
        N4 = 0.25*(1-xhi)*(1+eta);
        N <<  N1,0.0f,0.0f,N2,0.0f,0.0f,N3,0.0f,0.0f,N4,0.0f,0.0f,
            0.0f,N1,0.0f,0.0f,N2,0.0f,0.0f,N3,0.0f,0.0f,N4,0.0f;
        Me = Me + w*pid->get_mid()->get_density()*N.transpose()*N*detJ.at(i);
    }
};


CPS4::CPS4(unsigned int                        id,
           std::vector<Node*>  connectivity,
           Pid*                pid,
           const unsigned short                nnodes,
           const unsigned short                ndofs,
           const unsigned short                vtk_identifier,
           const unsigned short                ngp,
           const unsigned short                dimensions,
           std::string                         element_type):
Element{id,connectivity,pid,nnodes,ndofs,vtk_identifier,ngp,dimensions,element_type}{}



