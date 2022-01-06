#include <iostream>
#include <cmath>
#include "../include/CPS4.h"

std::vector<Scalar> CPS4::calculate_stress(dynMatrix D,
                                           dynMatrix u)
{
    // stress is calculated at Gauss points, by interpolating the displacements
    // at the nodes (using shape functions)
    // stress = [force/area] = E*epsilon = Young's modulus * strain = D * strain
    // B is strain-displacement matrix. If we use the Gauss points to evaluate
    //  stresses and strains we save the effort of evaluat the B matrix.
    // Want to evaluate σ at the Gauss integration points used in the element stiffness integration and
    // then extrapolate to the element node points.
    // Felipe "Introduction to the finite element method" chapter 28:
    // for a quadrilateral element with 2x2 Gauss integration:
    // ξ = ξ'/ 3, η = η'/3
    // where ' denotes an imaginary "Gauss element" inside our distorted element.
    // for given quantity w:
    // w(xhi', eta') = [w1' w2' w3' w4'][ N1' N2' N3' N4']
    // where for CPS4: N1' = (1 - ξ')(1 - η')/4
    //                 N2' = (1 + ξ')(1 - η')/4
    //                 N3' = (1 + ξ')(1 + η')/4
    //                 N4' = (1 - ξ')(1 + η')/4
    // c.f CPS4::calculate_Me shape functions.
    // in the case of stress w is σ_xx, σ_yy, τ_xy.
    // size(strain) = 1*Number of Gauss points
    // auto strain = calculate_strain(D, node_displacement);
    std::vector<Scalar> stresses;
    auto strains = calculate_strain(D,u);
    // for(const auto& strain: strains)
    // {
    //     auto stress = D*strain;
    //     // stresses.emplace_back(stress);
    // }

}
std::vector<Scalar> CPS4::calculate_strain(dynMatrix D,
                                           dynMatrix u)
{
    std::vector<Scalar> strains;
    for(const auto& Bi: B)
    {
        // size(B) = 3 x 8, size(u) = 8 x 1 = (nnodes x ndofs) x 1
        // --> size(strain) = 3 x 1 = [ε_x ε_y ε_xy]
        auto strain = Bi*u;
        strains.emplace_back(strain);
    }

}


std::vector<Segment>& CPS4::get_segments(Node* node)
{
    // given a node, what other nodes does it connect to?
    // will be different for each element type.
    // base on node position in connectivity vector
    std::cout << "CPS4::get_segments: node "<< node->id << ": x = " << node->x << ", y = " << node->y << std::endl;
    for(size_t i = 0; i < connectivity.size(); i++)
    {
        if(connectivity.at(i)->id == node->id)
        {
            if (i == 0)
            {
                segments.push_back(std::make_pair(connectivity.at(0), connectivity.at(1)));
                segments.push_back(std::make_pair(connectivity.at(1), connectivity.at(3)));
            }
            else if (i == 1)
            {
                segments.push_back(std::make_pair(connectivity.at(1), connectivity.at(0)));
                segments.push_back(std::make_pair(connectivity.at(1), connectivity.at(2)));
            }
            else if (i == 2)
            {
                segments.push_back(std::make_pair(connectivity.at(2), connectivity.at(1)));
                segments.push_back(std::make_pair(connectivity.at(2), connectivity.at(3)));
            }
            else if (i == 3)
            {
                segments.push_back(std::make_pair(connectivity.at(3), connectivity.at(2)));
                segments.push_back(std::make_pair(connectivity.at(3), connectivity.at(0)));
            }
        }
    }
    return segments;
};

void CPS4::calculate_Ke(){
    setup_coord();
    Mid* mid = pid->get_mid();
    float v = mid->get_v();
    float E = mid->get_E();
    Eigen::Matrix<float,3,3> D;
    D << 1, v, 0,
         v, 1, 0,
         0, 0, 0.5*(1-v);
    D *= E/(1-(v*v));
     // create coord matrix needed to find Jacobian
    // setup_coord();
    Eigen::Matrix<float,2,4> dNdXhidEta;
    Eigen::Matrix<float,2,4> dNdxdy;
    Eigen::Matrix<float,2,2> J;
    Eigen::Matrix<float,3,8> Bi;
    float xhi,eta,w;
    float dN1dXhi,dN2dXhi,dN3dXhi,dN4dXhi;
    float dN1dEta,dN2dEta,dN3dEta,dN4dEta;
    for (unsigned char i = 0; i < gauss_points->size(); i++)
    {            
        xhi = gauss_points->at(i).at(0);
        eta = gauss_points->at(i).at(1);
        w   = gauss_weights->at(i);
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
        if (detJ.back() < 0.1f)
        {
            std::cout << "WARNING: Jacobian determinant less than 0.1 for element #" << id << std::endl;    
        }
        dNdxdy = J.inverse()*dNdXhidEta;
        Bi <<    dNdxdy(0,0),           0, dNdxdy(0,1),           0, dNdxdy(0,2),           0, dNdxdy(0,3),           0,
                          0, dNdxdy(1,0),           0, dNdxdy(1,1),           0, dNdxdy(1,2),           0, dNdxdy(1,3),
                dNdxdy(1,0), dNdxdy(0,0), dNdxdy(1,1), dNdxdy(0,1), dNdxdy(1,2), dNdxdy(0,2), dNdxdy(1,3), dNdxdy(0,3); 
        B.emplace_back(Bi);
        Ke = Ke + w*(Bi.transpose()*D*Bi*detJ.back());
    }
};
void CPS4::calculate_Me(){
    setup_coord();
    // size(N) = dof per node x (nnodes*dof per node)
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
        N <<  N1, 0.f,   N2,  0.f,   N3,  0.f,   N4,  0.f,
             0.f,  N1,  0.f,   N2,  0.f,   N3,  0.f,   N4;
        Me = Me + w*pid->get_mid()->get_density()*N.transpose()*N*detJ.at(i);
    }
};


CPS4::CPS4(unsigned int                        id,
           std::vector<Node*>                  connectivity,
           Pid*                                pid,
           const unsigned short                nnodes,
           const unsigned short                ndofs,
           const unsigned short                vtk_identifier,
           const unsigned short                ngp,
           const unsigned short                dimensions,
           std::string                         element_type):
Element{id,connectivity,pid,nnodes,ndofs,vtk_identifier,ngp,dimensions,element_type}{
    print_element_info_to_log();
}



