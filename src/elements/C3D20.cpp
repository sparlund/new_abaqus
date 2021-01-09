#include <iostream>
#include "C3D20.h"
#include <stdlib.h>

const std::string C3D20::element_type = "C3D20";


C3D20::~C3D20(){}

C3D20::C3D20(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid):id(id),connectivity(connectivity),pid(pid){
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
    // Weight = volume*density
    // Formula for volume of generic hexahedron found here
    // Don't account for mid nodes...
    volume = ((coord.row(6) - coord.row(0)).dot( (coord.row(1) - coord.row(0)).cross( (coord.row(3) - coord.row(5)) ))  + 
              (coord.row(6) - coord.row(0)).dot( (coord.row(4) - coord.row(0)).cross( (coord.row(5) - coord.row(3)) ))  +
              (coord.row(6) - coord.row(0)).dot( (coord.row(3) - coord.row(0)).cross( (coord.row(7) - coord.row(2)) )))*(0.16667);
    weight = volume*mid->get_density();
    // abaqus node number -> Felipe node numbering
    // 17 -> 13
    // 18 -> 14
    // 19 -> 15
    // 20 -> 16
    // 13 -> 17
    // 14 -> 18
    // 15 -> 19
    // 16 -> 20

    // dim(J) = ndim x ndim
    Eigen::Matrix<float,3,3> J;
    Eigen::Matrix<float,3,3> invJ;
    // dim(B) = 6 x nnodes*3
    Eigen::Matrix<float,6,60> B;
    
    // all weights are 1, dont need w
    // float w = 1;
    Eigen::Matrix<float,1,8> N;
    // size(global_N_matrix) = dof per node x (nnodes*dof per node)
    Eigen::Matrix<float,3,60> global_N_matrix;
    // dim(dNdxdydz) = dim(dNdXhidEtadMy) = ndim x nnodes
    Eigen::Matrix<float,3,20> dNdxdydz;
    Eigen::Matrix<float,3,20> dNdXhidEtadMy;

    // init dense Ke & Me to zero
    Ke.setZero();
    Me.setZero();
    float N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,N14,N15,N16,N17,N18,N19,N20;
    float dN1dXhi,dN2dXhi,dN3dXhi,dN4dXhi,dN5dXhi,dN6dXhi,dN7dXhi,dN8dXhi,dN9dXhi,dN10dXhi,dN11dXhi,dN12dXhi,dN13dXhi,dN14dXhi,dN15dXhi,dN16dXhi,dN17dXhi,dN18dXhi,dN19dXhi,dN20dXhi;
    float dN1dEta,dN2dEta,dN3dEta,dN4dEta,dN5dEta,dN6dEta,dN7dEta,dN8dEta,dN9dEta,dN10dEta,dN11dEta,dN12dEta,dN13dEta,dN14dEta,dN15dEta,dN16dEta,dN17dEta,dN18dEta,dN19dEta,dN20dEta;   
    float dN1dMy,dN2dMy,dN3dMy,dN4dMy,dN5dMy,dN6dMy,dN7dMy,dN8dMy,dN9dMy,dN10dMy,dN11dMy,dN12dMy,dN13dMy,dN14dMy,dN15dMy,dN16dMy,dN17dMy,dN18dMy,dN19dMy,dN20dMy;
    float w,xhi,eta,my;
    std::cout.precision(10);
    for (unsigned int i = 0; i < ngp; i++) {
                xhi = gauss_points.at(i).at(0);
                eta = gauss_points.at(i).at(1);
                my = gauss_points.at(i).at(2);
                w = gauss_weights.at(i);
                // https://www.code-aster.org/V2/doc/v11/en/man_r/r3/r3.01.01.pdf
                // shape functions
                N1      = (0.125 - 0.125*xhi)*(1 - eta)*(1 - my)*(-eta - my - xhi - 2);
                dN1dXhi = -(0.125 - 0.125*xhi)*(1 - eta)*(1 - my) - 0.125*(1 - eta)*(1 - my)*(-eta - my - xhi - 2);
                dN1dEta = -(0.125 - 0.125*xhi)*(1 - eta)*(1 - my) + (0.125 - 0.125*xhi)*(my - 1)*(-eta - my - xhi - 2);
                dN1dMy  = -(0.125 - 0.125*xhi)*(1 - eta)*(1 - my) + (0.125 - 0.125*xhi)*(eta - 1)*(-eta - my - xhi - 2);
                //
                N2      = (1 - eta)*(1 - my)*(0.125*xhi + 0.125)*(-eta - my + xhi - 2);
                dN2dXhi = (1 - eta)*(1 - my)*(0.125*xhi + 0.125) + 0.125*(1 - eta)*(1 - my)*(-eta - my + xhi - 2);
                dN2dEta = -(1 - eta)*(1 - my)*(0.125*xhi + 0.125) + (my - 1)*(0.125*xhi + 0.125)*(-eta - my + xhi - 2);
                dN2dMy  = -(1 - eta)*(1 - my)*(0.125*xhi + 0.125) + (eta - 1)*(0.125*xhi + 0.125)*(-eta - my + xhi - 2);
                //
                N3      = (1 - my)*(eta + 1)*(0.125*xhi + 0.125)*(eta - my + xhi - 2);
                dN3dXhi = (1 - my)*(eta + 1)*(0.125*xhi + 0.125) + 0.125*(1 - my)*(eta + 1)*(eta - my + xhi - 2);
                dN3dEta = (1 - my)*(eta + 1)*(0.125*xhi + 0.125) + (1 - my)*(0.125*xhi + 0.125)*(eta - my + xhi - 2);
                dN3dMy  = -(1 - my)*(eta + 1)*(0.125*xhi + 0.125) + (-eta - 1)*(0.125*xhi + 0.125)*(eta - my + xhi - 2);
                //
                N4      = (0.125 - 0.125*xhi)*(1 - my)*(eta + 1)*(eta - my - xhi - 2);
                dN4dXhi = -(0.125 - 0.125*xhi)*(1 - my)*(eta + 1) - 0.125*(1 - my)*(eta + 1)*(eta - my - xhi - 2);
                dN4dEta = (0.125 - 0.125*xhi)*(1 - my)*(eta + 1) + (0.125 - 0.125*xhi)*(1 - my)*(eta - my - xhi - 2);
                dN4dMy  = -(0.125 - 0.125*xhi)*(1 - my)*(eta + 1) + (0.125 - 0.125*xhi)*(-eta - 1)*(eta - my - xhi - 2);
                //
                N5      = (0.125 - 0.125*xhi)*(1 - eta)*(my + 1)*(-eta + my - xhi - 2);
                dN5dXhi = -(0.125 - 0.125*xhi)*(1 - eta)*(my + 1) - 0.125*(1 - eta)*(my + 1)*(-eta + my - xhi - 2);
                dN5dEta = -(0.125 - 0.125*xhi)*(1 - eta)*(my + 1) + (0.125 - 0.125*xhi)*(-my - 1)*(-eta + my - xhi - 2);
                dN5dMy  = (0.125 - 0.125*xhi)*(1 - eta)*(my + 1) + (0.125 - 0.125*xhi)*(1 - eta)*(-eta + my - xhi - 2);
                //
                N6      = (1 - eta)*(my + 1)*(0.125*xhi + 0.125)*(-eta + my + xhi - 2);
                dN6dXhi = (1 - eta)*(my + 1)*(0.125*xhi + 0.125) + 0.125*(1 - eta)*(my + 1)*(-eta + my + xhi - 2);
                dN6dEta = -(1 - eta)*(my + 1)*(0.125*xhi + 0.125) + (-my - 1)*(0.125*xhi + 0.125)*(-eta + my + xhi - 2);
                dN6dMy  = (1 - eta)*(my + 1)*(0.125*xhi + 0.125) + (1 - eta)*(0.125*xhi + 0.125)*(-eta + my + xhi - 2);
                //
                N7      = (eta + 1)*(my + 1)*(0.125*xhi + 0.125)*(eta + my + xhi - 2);
                dN7dXhi = (eta + 1)*(my + 1)*(0.125*xhi + 0.125) + 0.125*(eta + 1)*(my + 1)*(eta + my + xhi - 2);
                dN7dEta = (eta + 1)*(my + 1)*(0.125*xhi + 0.125) + (my + 1)*(0.125*xhi + 0.125)*(eta + my + xhi - 2);
                dN7dMy  = (eta + 1)*(my + 1)*(0.125*xhi + 0.125) + (eta + 1)*(0.125*xhi + 0.125)*(eta + my + xhi - 2);
                //
                N8      = (0.125 - 0.125*xhi)*(eta + 1)*(my + 1)*(eta + my - xhi - 2);
                dN8dXhi = -(0.125 - 0.125*xhi)*(eta + 1)*(my + 1) - 0.125*(eta + 1)*(my + 1)*(eta + my - xhi - 2);
                dN8dEta = (0.125 - 0.125*xhi)*(eta + 1)*(my + 1) + (0.125 - 0.125*xhi)*(my + 1)*(eta + my - xhi - 2);
                dN8dMy  = (0.125 - 0.125*xhi)*(eta + 1)*(my + 1) + (0.125 - 0.125*xhi)*(eta + 1)*(eta + my - xhi - 2);
                //
                N9      = (0.25 - 0.25*std::pow(xhi,2.0))*(1 - eta)*(1 - my);
                dN9dXhi = -0.5*xhi*(1 - eta)*(1 - my);
                dN9dEta = (0.25 - 0.25*std::pow(xhi,2.0))*(my - 1);
                dN9dMy  = (0.25 - 0.25*std::pow(xhi,2.0))*(eta - 1);
                //
                N10      = (0.25 - 0.25*std::pow(eta,2.0))*(1 - my)*(xhi + 1);
                dN10dXhi = (0.25 - 0.25*std::pow(eta,2.0))*(1 - my);
                dN10dEta = -0.5*eta*(1 - my)*(xhi + 1);
                dN10dMy  = (0.25 - 0.25*std::pow(eta,2.0))*(-xhi - 1);
                //
                N11      = (0.25 - 0.25*std::pow(xhi,2.0))*(1 - my)*(eta + 1);
                dN11dXhi = -0.5*xhi*(1 - my)*(eta + 1);
                dN11dEta = (0.25 - 0.25*std::pow(xhi,2.0))*(1 - my);
                dN11dMy  = (0.25 - 0.25*std::pow(xhi,2.0))*(-eta - 1);
                //
                N12      = (0.25 - 0.25*std::pow(eta,2.0))*(1 - my)*(1 - xhi);
                dN12dXhi = (0.25 - 0.25*std::pow(eta,2.0))*(my - 1);
                dN12dEta = -0.5*eta*(1 - my)*(1 - xhi);
                dN12dMy  = (0.25 - 0.25*std::pow(eta,2.0))*(xhi - 1);
                //
                N17      = (0.25 - 0.25*std::pow(my,2.0))*(1 - eta)*(1 - xhi);
                dN17dXhi = (0.25 - 0.25*std::pow(my,2.0))*(eta - 1);
                dN17dEta = (0.25 - 0.25*std::pow(my,2.0))*(xhi - 1);
                dN17dMy  = -0.5*my*(1 - eta)*(1 - xhi);
                //
                N18      = (0.25 - 0.25*std::pow(my,2.0))*(1 - eta)*(xhi + 1);
                dN18dXhi = (0.25 - 0.25*std::pow(my,2.0))*(1 - eta);
                dN18dEta = (0.25 - 0.25*std::pow(my,2.0))*(-xhi - 1);
                dN18dMy  = -0.5*my*(1 - eta)*(xhi + 1);
                //
                N19      = (0.25 - 0.25*std::pow(my,2.0))*(eta + 1)*(xhi + 1);
                dN19dXhi = (0.25 - 0.25*std::pow(my,2.0))*(eta + 1);
                dN19dEta = (0.25 - 0.25*std::pow(my,2.0))*(xhi + 1);
                dN19dMy  = -0.5*my*(eta + 1)*(xhi + 1);
                //
                N20      = (0.25 - 0.25*std::pow(my,2.0))*(1 - xhi)*(eta + 1);
                dN20dXhi = (0.25 - 0.25*std::pow(my,2.0))*(-eta - 1);
                dN20dEta = (0.25 - 0.25*std::pow(my,2.0))*(1 - xhi);
                dN20dMy  = -0.5*my*(1 - xhi)*(eta + 1);
                //
                N13      = (0.25 - 0.25*std::pow(xhi,2.0))*(1 - eta)*(my + 1);
                dN13dXhi = -0.5*xhi*(1 - eta)*(my + 1);
                dN13dEta = (0.25 - 0.25*std::pow(xhi,2.0))*(-my - 1);
                dN13dMy  = (0.25 - 0.25*std::pow(xhi,2.0))*(1 - eta);
                //
                N14      = (0.25 - 0.25*std::pow(eta,2.0))*(my + 1)*(xhi + 1);
                dN14dXhi = (0.25 - 0.25*std::pow(eta,2.0))*(my + 1);
                dN14dEta = -0.5*eta*(my + 1)*(xhi + 1);
                dN14dMy  = (0.25 - 0.25*std::pow(eta,2.0))*(xhi + 1);
                //
                N15      = (0.25 - 0.25*std::pow(xhi,2.0))*(eta + 1)*(my + 1);
                dN15dXhi = -0.5*xhi*(eta + 1)*(my + 1);
                dN15dEta = (0.25 - 0.25*std::pow(xhi,2.0))*(my + 1);
                dN15dMy  = (0.25 - 0.25*std::pow(xhi,2.0))*(eta + 1);
                //
                N16      = (0.25 - 0.25*std::pow(eta,2.0))*(1 - xhi)*(my + 1);
                dN16dXhi = (0.25 - 0.25*std::pow(eta,2.0))*(-my - 1);
                dN16dEta = -0.5*eta*(1 - xhi)*(my + 1);
                dN16dMy  = (0.25 - 0.25*std::pow(eta,2.0))*(1 - xhi);
                //
                // find Jacobian for each Gauss point
                dNdXhidEtadMy.row(0) << dN1dXhi,dN2dXhi,dN3dXhi,dN4dXhi,dN5dXhi,dN6dXhi,dN7dXhi,dN8dXhi,dN9dXhi,dN10dXhi,dN11dXhi,dN12dXhi,dN13dXhi,dN14dXhi,dN15dXhi,dN16dXhi,dN17dXhi,dN18dXhi,dN19dXhi,dN20dXhi;
                dNdXhidEtadMy.row(1) << dN1dEta,dN2dEta,dN3dEta,dN4dEta,dN5dEta,dN6dEta,dN7dEta,dN8dEta,dN9dEta,dN10dEta,dN11dEta,dN12dEta,dN13dEta,dN14dEta,dN15dEta,dN16dEta,dN17dEta,dN18dEta,dN19dEta,dN20dEta;   
                dNdXhidEtadMy.row(2) << dN1dMy, dN2dMy, dN3dMy, dN4dMy, dN5dMy, dN6dMy, dN7dMy, dN8dMy, dN9dMy, dN10dMy, dN11dMy, dN12dMy, dN13dMy, dN14dMy, dN15dMy, dN16dMy, dN17dMy, dN18dMy, dN19dMy, dN20dMy;
                J = dNdXhidEtadMy*coord;
                // find shape functions derivate matrix wrt x, y & z
                // not always safe to just use member inverse(), instead use singlue value decomposition (SVD)
                if (J.determinant() < 0.1f)
                {
                    std::cout << "WARNING: Jacobian determinant less than 0.1 for element #" << this->get_id() << std::endl;
                }
                // This stupid shit doesnt work lol, misunderstood?
                // Eigen::JacobiSVD<Eigen::Matrix<float,3,3>> svd(J, Eigen::ComputeFullV | Eigen::ComputeFullU);
                // invJ = svd.matrixV()*(svd.singularValues().unaryExpr(&inv_div_by1).asDiagonal().inverse())*svd.matrixU().transpose();
                // std::cout << "J=[" << J << "];" << std::endl;
                // std::cout << "invJ=[" << invJ << "];" << std::endl;
                // exit(0);
                dNdxdydz = J.inverse()*dNdXhidEtadMy;
                // construct B matrix
                B << dNdxdydz(0,0),0.0f,0.0f,dNdxdydz(0,1),0.0f,0.0f,dNdxdydz(0,2),0.0f,0.0f,dNdxdydz(0,3),0.0f,0.0f,dNdxdydz(0,4),0.0f,0.0f,dNdxdydz(0,5),0.0f,0.0f,dNdxdydz(0,6),0.0f,0.0f,dNdxdydz(0,7),0.0f,0.0f,dNdxdydz(0,8),0.0f,0.0f,dNdxdydz(0,9),0.0f,0.0f,dNdxdydz(0,10),0.0f,0.0f,dNdxdydz(0,11),0.0f,0.0f,dNdxdydz(0,12),0.0f,0.0f,dNdxdydz(0,13),0.0f,0.0f,dNdxdydz(0,14),0.0f,0.0f,dNdxdydz(0,15),0.0f,0.0f,dNdxdydz(0,16),0.0f,0.0f,dNdxdydz(0,17),0.0f,0.0f,dNdxdydz(0,18),0.0f,0.0f,dNdxdydz(0,19),0.0f,0.0f,
                     0.0f,dNdxdydz(1,0),0.0f,0.0f,dNdxdydz(1,1),0.0f,0.0f,dNdxdydz(1,2),0.0f,0.0f,dNdxdydz(1,3),0.0f,0.0f,dNdxdydz(1,4),0.0f,0.0f,dNdxdydz(1,5),0.0f,0.0f,dNdxdydz(1,6),0.0f,0.0f,dNdxdydz(1,7),0.0f,0.0f,dNdxdydz(1,8),0.0f,0.0f,dNdxdydz(1,9),0.0f,0.0f,dNdxdydz(1,10),0.0f,0.0f,dNdxdydz(1,11),0.0f,0.0f,dNdxdydz(1,12),0.0f,0.0f,dNdxdydz(1,13),0.0f,0.0f,dNdxdydz(1,14),0.0f,0.0f,dNdxdydz(1,15),0.0f,0.0f,dNdxdydz(1,16),0.0f,0.0f,dNdxdydz(1,17),0.0f,0.0f,dNdxdydz(1,18),0.0f,0.0f,dNdxdydz(1,19),0.0f,
                     0.0f,0.0f,dNdxdydz(2,0),0.0f,0.0f,dNdxdydz(2,1),0.0f,0.0f,dNdxdydz(2,2),0.0f,0.0f,dNdxdydz(2,3),0.0f,0.0f,dNdxdydz(2,4),0.0f,0.0f,dNdxdydz(2,5),0.0f,0.0f,dNdxdydz(2,6),0.0f,0.0f,dNdxdydz(2,7),0.0f,0.0f,dNdxdydz(2,8),0.0f,0.0f,dNdxdydz(2,9),0.0f,0.0f,dNdxdydz(2,10),0.0f,0.0f,dNdxdydz(2,11),0.0f,0.0f,dNdxdydz(2,12),0.0f,0.0f,dNdxdydz(2,13),0.0f,0.0f,dNdxdydz(2,14),0.0f,0.0f,dNdxdydz(2,15),0.0f,0.0f,dNdxdydz(2,16),0.0f,0.0f,dNdxdydz(2,17),0.0f,0.0f,dNdxdydz(2,18),0.0f,0.0f,dNdxdydz(2,19),
                     dNdxdydz(1,0),dNdxdydz(0,0),0.0f,dNdxdydz(1,1),dNdxdydz(0,1),0.0f,dNdxdydz(1,2),dNdxdydz(0,2),0.0f,dNdxdydz(1,3),dNdxdydz(0,3),0.0f,dNdxdydz(1,4),dNdxdydz(0,4),0.0f,dNdxdydz(1,5),dNdxdydz(0,5),0.0f,dNdxdydz(1,6),dNdxdydz(0,6),0.0f,dNdxdydz(1,7),dNdxdydz(0,7),0.0f,dNdxdydz(1,8),dNdxdydz(0,8),0.0f,dNdxdydz(1,9),dNdxdydz(0,9),0.0f,dNdxdydz(1,10),dNdxdydz(0,10),0.0f,dNdxdydz(1,11),dNdxdydz(0,11),0.0f,dNdxdydz(1,12),dNdxdydz(0,12),0.0f,dNdxdydz(1,13),dNdxdydz(0,13),0.0f,dNdxdydz(1,14),dNdxdydz(0,14),0.0f,dNdxdydz(1,15),dNdxdydz(0,15),0.0f,dNdxdydz(1,16),dNdxdydz(0,16),0.0f,dNdxdydz(1,17),dNdxdydz(0,17),0.0f,dNdxdydz(1,18),dNdxdydz(0,18),0.0f,dNdxdydz(1,19),dNdxdydz(0,19),0.0f,
                     0.0f,dNdxdydz(2,0),dNdxdydz(1,0),0.0f,dNdxdydz(2,1),dNdxdydz(1,1),0.0f,dNdxdydz(2,2),dNdxdydz(1,2),0.0f,dNdxdydz(2,3),dNdxdydz(1,3),0.0f,dNdxdydz(2,4),dNdxdydz(1,4),0.0f,dNdxdydz(2,5),dNdxdydz(1,5),0.0f,dNdxdydz(2,6),dNdxdydz(1,6),0.0f,dNdxdydz(2,7),dNdxdydz(1,7),0.0f,dNdxdydz(2,8),dNdxdydz(1,8),0.0f,dNdxdydz(2,9),dNdxdydz(1,9),0.0f,dNdxdydz(2,10),dNdxdydz(1,10),0.0f,dNdxdydz(2,11),dNdxdydz(1,11),0.0f,dNdxdydz(2,12),dNdxdydz(1,12),0.0f,dNdxdydz(2,13),dNdxdydz(1,13),0.0f,dNdxdydz(2,14),dNdxdydz(1,14),0.0f,dNdxdydz(2,15),dNdxdydz(1,15),0.0f,dNdxdydz(2,16),dNdxdydz(1,16),0.0f,dNdxdydz(2,17),dNdxdydz(1,17),0.0f,dNdxdydz(2,18),dNdxdydz(1,18),0.0f,dNdxdydz(2,19),dNdxdydz(1,19),
                     dNdxdydz(2,0),0.0f,dNdxdydz(0,0),dNdxdydz(2,1),0.0f,dNdxdydz(0,1),dNdxdydz(2,2),0.0f,dNdxdydz(0,2),dNdxdydz(2,3),0.0f,dNdxdydz(0,3),dNdxdydz(2,4),0.0f,dNdxdydz(0,4),dNdxdydz(2,5),0.0f,dNdxdydz(0,5),dNdxdydz(2,6),0.0f,dNdxdydz(0,6),dNdxdydz(2,7),0.0f,dNdxdydz(0,7),dNdxdydz(2,8),0.0f,dNdxdydz(0,8),dNdxdydz(2,9),0.0f,dNdxdydz(0,9),dNdxdydz(2,10),0.0f,dNdxdydz(0,10),dNdxdydz(2,11),0.0f,dNdxdydz(0,11),dNdxdydz(2,12),0.0f,dNdxdydz(0,12),dNdxdydz(2,13),0.0f,dNdxdydz(0,13),dNdxdydz(2,14),0.0f,dNdxdydz(0,14),dNdxdydz(2,15),0.0f,dNdxdydz(0,15),dNdxdydz(2,16),0.0f,dNdxdydz(0,16),dNdxdydz(2,17),0.0f,dNdxdydz(0,17),dNdxdydz(2,18),0.0f,dNdxdydz(0,18),dNdxdydz(2,19),0.0f,dNdxdydz(0,19),
                // compute Gauss point contribution to element Ke matrix
                // § 8.5. The Mass Matrix, chapter 8
                global_N_matrix << N1,0.0f,0.0f,N2,0.0f,0.0f,N3,0.0f,0.0f,N4,0.0f,0.0f,N5,0.0f,0.0f,N6,0.0f,0.0f,N7,0.0f,0.0f,N8,0.0f,0.0f,N9,0.0f,0.0f,N10.0f,0.0f,0.0f,N11,0.0f,0.0f,N12,0.0f,0.0f,N13,0.0f,0.0f,N14,0.0f,0.0f,N15,0.0f,0.0f,N16,0.0f,0.0f,N17,0.0f,0.0f,N18,0.0f,0.0f,N19,0.0f,0.0f,N20.0f,0.0f,0.0f,
                                   0.0f,N1,0.0f,0.0f,N2,0.0f,0.0f,N3,0.0f,0.0f,N4,0.0f,0.0f,N5,0.0f,0.0f,N6,0.0f,0.0f,N7,0.0f,0.0f,N8,0.0f,0.0f,N9,0.0f,0.0f,N10.0f,0.0f,0.0f,N11,0.0f,0.0f,N12,0.0f,0.0f,N13,0.0f,0.0f,N14,0.0f,0.0f,N15,0.0f,0.0f,N16,0.0f,0.0f,N17,0.0f,0.0f,N18,0.0f,0.0f,N19,0.0f,0.0f,N20.0f,0.0f,
                                   0.0f,0.0f,N1,0.0f,0.0f,N2,0.0f,0.0f,N3,0.0f,0.0f,N4,0.0f,0.0f,N5,0.0f,0.0f,N6,0.0f,0.0f,N7,0.0f,0.0f,N8,0.0f,0.0f,N9,0.0f,0.0f,N10.0f,0.0f,0.0f,N11,0.0f,0.0f,N12,0.0f,0.0f,N13,0.0f,0.0f,N14,0.0f,0.0f,N15,0.0f,0.0f,N16,0.0f,0.0f,N17,0.0f,0.0f,N18,0.0f,0.0f,N19,0.0f,0.0f,N20.0f;  
                Ke = Ke + w*(B.transpose()*D*B*J.determinant());
                Me = Me + w*mid->get_density()*global_N_matrix.transpose()*global_N_matrix*J.determinant();
    }
    print_element_info_to_log();
}
