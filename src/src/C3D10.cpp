#include "../include/C3D10.h"
#include <iostream>
#include <stdlib.h>

void C3D10::calculate_Ke(){
    setup_coord();
    Eigen::Matrix<double,3,3> J;
    Eigen::Matrix<double,6,30> B;
    
    // size(dNdXhi) = ngp*ngp, 4 vectors of (1x4)
    Eigen::Matrix<double,3,10> dNdXhidEtadMy;
    Eigen::Matrix<double,3,10> dNdxdydz;

    // init Ke zero
    Ke.setZero();
    double dN1dXhi, dN2dXhi, dN3dXhi, dN4dXhi, dN5dXhi, dN6dXhi, dN7dXhi, dN8dXhi, dN9dXhi, dN10dXhi;
    double dN1dEta, dN2dEta, dN3dEta, dN4dEta, dN5dEta, dN6dEta, dN7dEta, dN8dEta, dN9dEta, dN10dEta;
    double dN1dMy, dN2dMy, dN3dMy, dN4dMy, dN5dMy, dN6dMy, dN7dMy, dN8dMy, dN9dMy, dN10dMy;
    double xhi,eta,my,w;
    auto D = pid->get_mid()->D_3D_linear_continuum_mechanics;
    for (unsigned int i = 0; i < gauss_points->size(); i++)
    {
        xhi = gauss_points->at(i).at(0);
        eta = gauss_points->at(i).at(1);
        my  = gauss_points->at(i).at(2);
        w   = gauss_weights->at(i);    
        // derive shape functions wrt xhi
        dN1dXhi = 0;
        dN2dXhi = 0;
        dN3dXhi = 4*eta + 4*my + 4*xhi - 3;
        dN4dXhi = 4*xhi - 1;
        dN5dXhi = 0;
        dN6dXhi = -4*my;
        dN7dXhi = -4*eta;
        dN8dXhi = 4*eta;
        dN9dXhi = 4*my;
        dN10dXhi = -4*eta - 4*my - 8*xhi + 4;
        // derive shape functions wrt eta
        dN1dEta = 4*eta - 1;
        dN2dEta = 0;
        dN3dEta = 4*eta + 4*my + 4*xhi - 3;
        dN4dEta = 0;
        dN5dEta = 4*my;
        dN6dEta = -4*my;
        dN7dEta = -8*eta - 4*my - 4*xhi + 4;
        dN8dEta = 4*xhi;
        dN9dEta = 0;
        dN10dEta = -4*xhi;
        // derive shape functions wrt my
        dN1dMy = 0;
        dN2dMy = 4*my - 1;
        dN3dMy = 4*eta + 4*my + 4*xhi - 3;
        dN4dMy = 0;
        dN5dMy = 4*eta;
        dN6dMy = -4*eta - 8*my - 4*xhi + 4;
        dN7dMy = -4*eta;
        dN8dMy = 0;
        dN9dMy = 4*xhi;
        dN10dMy = -4*xhi;
        // 
        dNdXhidEtadMy.row(0) << dN1dXhi,dN2dXhi,dN3dXhi,dN4dXhi,dN5dXhi,dN6dXhi,dN7dXhi,dN8dXhi,dN9dXhi,dN10dXhi;
        dNdXhidEtadMy.row(1) << dN1dEta,dN2dEta,dN3dEta,dN4dEta,dN5dEta,dN6dEta,dN7dEta,dN8dEta,dN9dEta,dN10dEta;
        dNdXhidEtadMy.row(2) << dN1dMy, dN2dMy, dN3dMy, dN4dMy, dN5dMy, dN6dMy, dN7dMy, dN8dMy, dN9dMy, dN10dMy;
        // find Jacobian for each Gauss point
        J = dNdXhidEtadMy*coord;
        detJ.push_back(J.determinant());
        // find shape functions derivate matrix wrt x, y & z
        dNdxdydz = J.inverse()*dNdXhidEtadMy;
        // construct B matrix
        B << dNdxdydz(0,0),0.d,0.d,dNdxdydz(0,1),0.d,0.d,dNdxdydz(0,2),0.d,0.d,dNdxdydz(0,3),0.d,0.d,dNdxdydz(0,4),0.d,0.d,dNdxdydz(0,5),0.d,0.d,dNdxdydz(0,6),0.d,0.d,dNdxdydz(0,7),0.d,0.d,dNdxdydz(0,8),0.d,0.d,dNdxdydz(0,9),0.d,0.d,
             0.d,dNdxdydz(1,0),0.d,0.d,dNdxdydz(1,1),0.d,0.d,dNdxdydz(1,2),0.d,0.d,dNdxdydz(1,3),0.d,0.d,dNdxdydz(1,4),0.d,0.d,dNdxdydz(1,5),0.d,0.d,dNdxdydz(1,6),0.d,0.d,dNdxdydz(1,7),0.d,0.d,dNdxdydz(1,8),0.d,0.d,dNdxdydz(1,9),0.d,
             0.d,0.d,dNdxdydz(2,0),0.d,0.d,dNdxdydz(2,1),0.d,0.d,dNdxdydz(2,2),0.d,0.d,dNdxdydz(2,3),0.d,0.d,dNdxdydz(2,4),0.d,0.d,dNdxdydz(2,5),0.d,0.d,dNdxdydz(2,6),0.d,0.d,dNdxdydz(2,7),0.d,0.d,dNdxdydz(2,8),0.d,0.d,dNdxdydz(2,9),
             dNdxdydz(1,0),dNdxdydz(0,0),0.d,dNdxdydz(1,1),dNdxdydz(0,1),0.d,dNdxdydz(1,2),dNdxdydz(0,2),0.d,dNdxdydz(1,3),dNdxdydz(0,3),0.d,dNdxdydz(1,4),dNdxdydz(0,4),0.d,dNdxdydz(1,5),dNdxdydz(0,5),0.d,dNdxdydz(1,6),dNdxdydz(0,6),0.d,dNdxdydz(1,7),dNdxdydz(0,7),0.d,dNdxdydz(1,8),dNdxdydz(0,8),0.d,dNdxdydz(1,9),dNdxdydz(0,9),0.d,
             0.d,dNdxdydz(2,0),dNdxdydz(1,0),0.d,dNdxdydz(2,1),dNdxdydz(1,1),0.d,dNdxdydz(2,2),dNdxdydz(1,2),0.d,dNdxdydz(2,3),dNdxdydz(1,3),0.d,dNdxdydz(2,4),dNdxdydz(1,4),0.d,dNdxdydz(2,5),dNdxdydz(1,5),0.d,dNdxdydz(2,6),dNdxdydz(1,6),0.d,dNdxdydz(2,7),dNdxdydz(1,7),0.d,dNdxdydz(2,8),dNdxdydz(1,8),0.d,dNdxdydz(2,9),dNdxdydz(1,9),
             dNdxdydz(2,0),0.d,dNdxdydz(0,0),dNdxdydz(2,1),0.d,dNdxdydz(0,1),dNdxdydz(2,2),0.d,dNdxdydz(0,2),dNdxdydz(2,3),0.d,dNdxdydz(0,3),dNdxdydz(2,4),0.d,dNdxdydz(0,4),dNdxdydz(2,5),0.d,dNdxdydz(0,5),dNdxdydz(2,6),0.d,dNdxdydz(0,6),dNdxdydz(2,7),0.d,dNdxdydz(0,7),dNdxdydz(2,8),0.d,dNdxdydz(0,8),dNdxdydz(2,9),0.d,dNdxdydz(0,9);
        // compute Gauss point contribution to element Ke matrix
        Ke = Ke + w*(B.transpose()*D*B*detJ.back());
    }

}
void C3D10::calculate_Me(){
    setup_coord();
    Eigen::Matrix<double,3,30> N;
    // init Ke zero
    Me.setZero();
    double N1,N2,N3,N4,N5,N6,N7,N8, N9, N10;
    double xhi,eta,my,w;
    for (unsigned int i = 0; i < gauss_points->size(); i++)
    {
        xhi = gauss_points->at(i).at(0);
        eta = gauss_points->at(i).at(1);
        my  = gauss_points->at(i).at(2);
        w   = gauss_weights->at(i);
        // shape functions
        N1 = eta*(2*eta - 1.d);
        N2 = my*(2*my - 1.d);
        N3 = (1.d - xhi - eta - my)*(1 - 2*xhi - 2*eta - 2*my);
        N4 = xhi*(2*xhi - 1.d);
        N5 = 4*eta*my;
        N6 = 4*my*(1.d - xhi - eta - my);
        N7 = 4*eta*(1.d - xhi - eta - my);
        N8 = 4*xhi*eta;
        N9 = 4*xhi*my;
        N10= 4*xhi*(1.d - xhi - eta - my);
        N << N1,0.d,0.d,N2,0.d,0.d,N3,0.d,0.d,N4,0.d,0.d,N5,0.d,0.d,N6,0.d,0.d,N7,0.d,0.d,N8,0.d,0.d,N9,0.d,0.d,N10,0.d,0.d,
             0.d,N1,0.d,0.d,N2,0.d,0.d,N3,0.d,0.d,N4,0.d,0.d,N5,0.d,0.d,N6,0.d,0.d,N7,0.d,0.d,N8,0.d,0.d,N9,0.d,0.d,N10,0.d,
             0.d,0.d,N1,0.d,0.d,N2,0.d,0.d,N3,0.d,0.d,N4,0.d,0.d,N5,0.d,0.d,N6,0.d,0.d,N7,0.d,0.d,N8,0.d,0.d,N9,0.d,0.d,N10;
        Me = Me + w*(pid->get_mid()->get_density()*N.transpose()*N*detJ.at(i));
    }

}
// C3D10 is 10 node tetrahedron element
C3D10::C3D10(unsigned int                        id,
             std::vector<Node*>                  connectivity,
             Pid*                                pid):
Element{id,connectivity,pid,ElementType::C3D10,10,30,24,4,3}{
    print_element_info_to_log();
}


