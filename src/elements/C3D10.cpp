#include <iostream>
#include "C3D10.h"
#include <stdlib.h>
// Eigen::Matrix<float,3,10> C3D10::shape_functions(Eigen::Matrix<float,1,4> evaluation_points){
//     return;
// }

const std::string C3D10::element_type = "C3D10";


// C3D10 is 10 node tetrahedron element
C3D10::C3D10(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid){
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
    std::cout << D << std::endl;
    for (unsigned short i = 0; i < nnodes; i++)
    {
        coord(i,0) = connectivity.at(i)->x;
        coord(i,1) = connectivity.at(i)->y;
        coord(i,2) = connectivity.at(i)->z;
    }
    std::cout << "coord=" << coord << std::endl;



    Eigen::Matrix<float,4,4> J;
    Eigen::Matrix<float,6,30> B;
    
    // size(N) = ngp*nnodes?
    Eigen::Matrix<float,4,10> N;
    // size(dNdXhi) = ngp*ngp, 4 vectors of (1x4)
    Eigen::Matrix<float,1,10> dNdt;
    Eigen::Matrix<float,1,10> dNdXhi;
    Eigen::Matrix<float,1,10> dNdEta;
    Eigen::Matrix<float,1,10> dNdMy;
    Eigen::Matrix<float,4,10> dNdtdXhidEtadMy;
    Eigen::Matrix<float,4,10> dNdxdydz;
    Eigen::Matrix<float,1,4> t;
    Eigen::Matrix<float,1,4> xhi;
    Eigen::Matrix<float,1,4> eta;
    Eigen::Matrix<float,1,4> my;

    // init Ke zero
    Ke.setZero();
    // weights same for all Gauss points
    unsigned short W = 0.041666667; // 0.25/6 = 0.041666667
    // gauss points same in all directions!
    // 1/sqrt(3) = 0.57735026919
    t   <<  0.58541020,0.13819660,0.13819660,0.13819660;
    xhi <<  0.13819660,0.58541020,0.13819660,0.13819660; 
    eta <<  0.13819660,0.13819660,0.58541020,0.13819660; 
    my  <<  0.13819660,0.13819660,0.13819660,0.58541020;
// alfa=(5.0+3.0*sqrt(5.0))/20 = 0.58
// beta =(5.0-sqrt(5.0))/20 = 0.13
// GaussPoints=[alfa beta beta beta;beta alfa beta beta;beta beta alfa  beta;beta beta beta alfa];

    for (unsigned int i = 0; i < ngp; i++)
    {
        // shape functions, 1/8=0.125
        N(0) = t(i)*(2*t(i)-1);
        N(1) = xhi(i)*(2*xhi(i)-1);
        N(2) = eta(i)*(2*eta(i)-1); 
        N(3) = my(i)*(2*my(i)-1); 
        N(4) = 4*t(i)*xhi(i);
        N(5) = 4*xhi(i)*eta(i);
        N(6) = 4*eta(i)*t(i);
        N(7) = 4*t(i)*my(i);
        N(8) = 4*xhi(i)*my(i);
        N(9) = 4*eta(i)*my(i);
        // derive shape functions wrt t
        dNdt(0) = 4*t(i)-1;
        dNdt(1) = 0;
        dNdt(2) = 0;
        dNdt(3) = 0;
        dNdt(4) = 4*xhi(i);
        dNdt(5) = 0;
        dNdt(6) = 4*eta(i);
        dNdt(7) = 4*my(i);
        dNdt(8) = 0;
        dNdt(9) = 0;
        // derive shape functions wrt xhi
        dNdXhi(0) = 0;
        dNdXhi(1) = 4*xhi(i)-1;
        dNdXhi(2) = 0;
        dNdXhi(3) = 0;
        dNdXhi(4) = 4*t(i);
        dNdXhi(5) = 4*eta(i);
        dNdXhi(6) = 0;
        dNdXhi(7) = 0;
        dNdXhi(8) = 4*my(i);
        dNdXhi(9) = 0;
        // derive shape functions wrt eta
        dNdEta(0) = 0; 
        dNdEta(1) = 0;
        dNdEta(2) = 4*eta(i)-1;
        dNdEta(3) = 0;
        dNdEta(4) = 0;
        dNdEta(5) = 4*xhi(i);
        dNdEta(6) = 4*t(i);
        dNdEta(7) = 0;
        dNdEta(8) = 0;
        dNdEta(9) = 4*my(i);
        // derive shape functions wrt my
        dNdMy(0)  = 0;
        dNdMy(1)  = 0;
        dNdMy(2)  = 0;
        dNdMy(3)  = 4*my(i)-1;
        dNdMy(4)  = 0;
        dNdMy(5)  = 0;
        dNdMy(6)  = 0;
        dNdMy(7)  = 4*t(i);
        dNdMy(8)  = 4*xhi(i);
        dNdMy(9)  = 4*eta(i);
        // find Jacobian for each Gauss point
        dNdtdXhidEtadMy.row(0) = dNdt;
        dNdtdXhidEtadMy.row(1) = dNdXhi;
        dNdtdXhidEtadMy.row(2) = dNdEta;
        std::cout << ":)" << std::endl;
        std::cout << "coord(4,0)*t(i)+coord(1,0)*(xhi(i)-1)+coord(5,0)*(eta(i)-1)+coord(8,0)*my(i)=" << coord(4,0)*t(i)       +coord(1,0)*(xhi(i)-1)+ coord(5,0)*(eta(i)-1) + coord(8,0)*my(i) << std::endl;
        std::cout << coord(4,0)*t(i) << std::endl;
        std::cout << coord(1,0)*(xhi(i)-1) << std::endl;
        std::cout << coord(5,0)*(eta(i)-1) << std::endl;
        std::cout << coord(8,0)*my(i)    << std::endl;
        
        J << 0.25, coord(0,0)*(t(i)-0.25)+coord(4,0)*xhi(i)         + coord(6,0)*eta(i)         + coord(7,0)*my(i)       , coord(0,1)*(t(i)-0.25)+ coord(4,1)*xhi(i)        + coord(6,1)*eta(i)         +coord(7,1)*my(i)       , coord(0,2)*(t(i)-0.25) + coord(4,2)*xhi(i)        + coord(6,2)*eta(i)         + coord(7,2)*my(i),
             0.25, coord(4,0)*t(i)       +coord(1,0)*(xhi(i)-0.25)  + coord(5,0)*(eta(i)-0.25)  + coord(8,0)*my(i)       , coord(4,1)*t(i)       + coord(1,1)*(xhi(i)-0.25) + coord(5,1)*(eta(i)-0.25)  +coord(8,1)*my(i)       , coord(4,2)*t(i)        + coord(1,2)*(xhi(i)-0.25) + coord(5,2)*(eta(i)-0.25)  + coord(8,2)*my(i),
             0.25, coord(6,0)*t(i)       +coord(5,0)*xhi(i)         + coord(2,0)*eta(i)         + coord(9,0)*(my(i)-0.25), coord(6,1)*t(i)       + coord(5,1)*xhi(i)        + coord(2,1)*eta(i)         +coord(9,1)*(my(i)-0.25), coord(6,2)*t(i)        + coord(5,2)*xhi(i)        + coord(2,2)*eta(i)         + coord(9,2)*(my(i)-0.25),
             0.25, coord(7,0)*t(i)       +coord(8,0)*xhi(i)         + coord(9,0)*eta(i)         + coord(3,0)*my(i)       , coord(7,1)*t(i)       + coord(8,1)*xhi(i)        + coord(9,1)*eta(i)         +coord(3,1)*my(i)       , coord(7,2)*t(i)        + coord(8,2)*xhi(i)        + coord(9,2)*eta(i)         + coord(3,2)*my(i);  
        J *= 4;
        J.transposeInPlace();
        std::cout << "J=[" << J << std::endl;

        // find shape functions derivate matrix wrt x, y & z
        std::cout << "inv(J)=[" << J.inverse() << std::endl;
        dNdxdydz = J.inverse()*dNdtdXhidEtadMy;
        std::cout << "dNdxdydz=[" << dNdxdydz << std::endl;
        // find strain-displacement matrix B
        // B << dNdxdydz(0,0), 0, 0, dNdxdydz(0,1), 0, 0,dNdxdydz(0,2), 0, 0,dNdxdydz(0,3), 0, 0,dNdxdydz(0,4), 0, 0,dNdxdydz(0,5), 0, 0,dNdxdydz(0,6), 0, 0,dNdxdydz(0,7), 0, 0,dNdxdydz(0,8), 0, 0,dNdxdydz(0,9), 0, 0, 
        //      0,dNdxdydz(1,0),0,0,dNdxdydz(1,1),0,0,dNdxdydz(1,2),0,0,dNdxdydz(1,3),0,0,dNdxdydz(1,4),0,0,dNdxdydz(1,5),0,0,dNdxdydz(1,6),0,0,dNdxdydz(1,7),0,0,dNdxdydz(1,8),0,0,dNdxdydz(1,9),0,
        //      0,0,dNdxdydz(2,0),0,0,dNdxdydz(2,1),0,0,dNdxdydz(2,2),0,0,dNdxdydz(2,3),0,0,dNdxdydz(2,4),0,0,dNdxdydz(2,5),0,0,dNdxdydz(2,6),0,0,dNdxdydz(2,7),0,0,dNdxdydz(2,8),0,0,dNdxdydz(2,9),
        //      dNdxdydz(1,0),dNdxdydz(0,0),0,dNdxdydz(1,1),dNdxdydz(0,1),0,dNdxdydz(1,2),dNdxdydz(0,2),0,dNdxdydz(1,3),dNdxdydz(0,3),0,dNdxdydz(1,4),dNdxdydz(0,4),0,dNdxdydz(1,5),dNdxdydz(0,5),0,dNdxdydz(1,6),dNdxdydz(0,6),0,dNdxdydz(1,7),dNdxdydz(0,7),0,dNdxdydz(1,8),dNdxdydz(0,8),0,dNdxdydz(1,9),dNdxdydz(0,9),0,
        //      0,dNdxdydz(2,0),dNdxdydz(1,0),0,dNdxdydz(2,1),dNdxdydz(1,1),0,dNdxdydz(2,2),dNdxdydz(1,2),0,dNdxdydz(2,3),dNdxdydz(1,3),0,dNdxdydz(2,4),dNdxdydz(1,4),0,dNdxdydz(2,5),dNdxdydz(1,5),0,dNdxdydz(2,6),dNdxdydz(1,6),0,dNdxdydz(2,7),dNdxdydz(1,7),0,dNdxdydz(2,8),dNdxdydz(1,8),0,dNdxdydz(2,9),dNdxdydz(1,9),
        //      dNdxdydz(2,0),0,dNdxdydz(0,0),dNdxdydz(2,1),0,dNdxdydz(0,1),dNdxdydz(2,2),0,dNdxdydz(0,2),dNdxdydz(2,3),0,dNdxdydz(0,3),dNdxdydz(2,4),0,dNdxdydz(0,4),dNdxdydz(2,5),0,dNdxdydz(0,5),dNdxdydz(2,6),0,dNdxdydz(0,6),dNdxdydz(2,7),0,dNdxdydz(0,7),dNdxdydz(2,8),0,dNdxdydz(0,8),dNdxdydz(2,9),0,dNdxdydz(0,9);
        Ke = Ke + W*(B.transpose()*D*B*J.determinant());
        std::cout << "Ke=[" << Ke << std::endl;
        // exit(0);
    }
    

}
C3D10::~C3D10(){

}


