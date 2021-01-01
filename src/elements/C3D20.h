#pragma once
#include <vector>
#include <array>
#include <memory>
#include <Eigen/Dense>
#include "../element.h"
#include "../pid.h"
#include "../node.h"


class C3D20 : public Element
{
private:
    unsigned int id;
    static const std::string element_type;
    static const unsigned short nnodes  = 20;
    static const unsigned char ndofs   = 60; // 30*3
    // https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    static const unsigned short vtk_identifier = 25;
    std::vector<unsigned int> dofs_id;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    // dim(coord) = nnodes*ndim
    Eigen::Matrix<float,20,3> coord;
    // dim(K) = (nnodes*ndim) x (nnodes*ndim)
    Eigen::Matrix<float,60,60> Me;
    Eigen::Matrix<float,60,60> Ke;
    Eigen::Matrix<float,60,1> fe;  
    // Eigen::Matrix<float,3,10> shape_functions(Eigen::Matrix<float,1,4> evaluation_points);
    const unsigned short ngp = 27;  
    const unsigned short ngp_per_dim = 3;
    // Eigen::Matrix<float,1,8> xhi;
    // Eigen::Matrix<float,1,8> eta;
    // Eigen::Matrix<float,1,8> my;
    // Gauss points position and weights are from Calculix. Positions are standard but weights are weird, no idea where they found them lol
    std::array<std::array<float, 3>,27> gauss_points = {{
                                            {-0.774596669241483,  -0.774596669241483,  -0.774596669241483},
                                            {0.000000000000000,  -0.774596669241483,  -0.774596669241483},
                                            {0.774596669241483,  -0.774596669241483,  -0.774596669241483},
                                            {-0.774596669241483,   0.000000000000000,  -0.774596669241483},
                                            {0.000000000000000,   0.000000000000000,  -0.774596669241483},
                                            {0.774596669241483,   0.000000000000000,  -0.774596669241483},
                                            {-0.774596669241483,   0.774596669241483,  -0.774596669241483},
                                            {0.000000000000000,   0.774596669241483,  -0.774596669241483},
                                            {0.774596669241483,   0.774596669241483,  -0.774596669241483},
                                            {-0.774596669241483,  -0.774596669241483,   0.000000000000000},
                                            {0.000000000000000,  -0.774596669241483,   0.000000000000000},
                                            {0.774596669241483,  -0.774596669241483,   0.000000000000000},
                                            {-0.774596669241483,   0.000000000000000,   0.000000000000000},
                                            {0.000000000000000,   0.000000000000000,   0.000000000000000},
                                            {0.774596669241483,   0.000000000000000,   0.000000000000000},
                                            {-0.774596669241483,   0.774596669241483,   0.000000000000000},
                                            {0.000000000000000,   0.774596669241483,   0.000000000000000},
                                            {0.774596669241483,   0.774596669241483,   0.000000000000000},
                                            {-0.774596669241483,  -0.774596669241483,   0.774596669241483},
                                            {0.000000000000000,  -0.774596669241483,   0.774596669241483},
                                            {0.774596669241483,  -0.774596669241483,   0.774596669241483},
                                            {-0.774596669241483,   0.000000000000000,   0.774596669241483},
                                            {0.000000000000000,   0.000000000000000,   0.774596669241483},
                                            {0.774596669241483,   0.000000000000000,   0.774596669241483},
                                            {-0.774596669241483,   0.774596669241483,   0.774596669241483},
                                            {0.000000000000000,   0.774596669241483,   0.774596669241483},
                                            {0.774596669241483,   0.774596669241483,   0.774596669241483}   }};

    std::array<float, 27> gauss_weights = { 0.171467764060357,   0.274348422496571,   0.171467764060357,
                                           0.274348422496571,   0.438957475994513,   0.274348422496571,
                                           0.171467764060357,   0.274348422496571,   0.171467764060357,
                                           0.274348422496571,   0.438957475994513,   0.274348422496571,
                                           0.438957475994513,   0.702331961591221,   0.438957475994513,
                                           0.274348422496571,   0.438957475994513,   0.274348422496571,
                                           0.171467764060357,   0.274348422496571,   0.171467764060357,
                                           0.274348422496571,   0.438957475994513,   0.274348422496571,
                                           0.171467764060357,   0.274348422496571,   0.171467764060357};

public:
    std::shared_ptr<Pid> get_pid(){return this->pid;};
    unsigned int get_id(){return id;};
    std::vector<unsigned int> get_element_dof_ids(){return dofs_id;};
    std::vector<std::shared_ptr<Node>> get_connectivity(){return this->connectivity;};
    unsigned short get_element_ndofs(){return ndofs;}
    unsigned short get_element_nnodes(){return nnodes;}
    unsigned short get_vtk_identifier(){return vtk_identifier;}
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Ke(){return Ke;}
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Me(){return Me;}
    std::string get_element_type(){return element_type;}
    C3D20(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~C3D20();
};