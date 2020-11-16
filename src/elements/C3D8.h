#pragma once
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include "../element.h"
#include "../pid.h"
#include "../node.h"


class C3D8 : public Element
{
private:
    unsigned int id;
    static const std::string element_type;
    static const unsigned short nnodes  = 8;
    static const unsigned char ndofs   = 24; // 10*3
    // https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    static const unsigned short vtk_identifier = 12;
    std::vector<unsigned int> dofs_id;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    float A;
    // dim(coord) = nnodes*ndim
    Eigen::Matrix<float,8,3> coord;
    // size(N) = ngp*nnodes
    Eigen::Matrix<float,8,8> N;
    // size(dNdXhi) = ngp*ngp, 4 vectors of (1x4)
    Eigen::Matrix<float,4,4> dNdXhi;
    Eigen::Matrix<float,4,4> dNdEta;

    // dim(dNdxdydz) = ndim x nnodes
    Eigen::Matrix<float,3,8> dNdxdydz;
    Eigen::Matrix<float,3,8> dNdXhidEtadMy;

    // Eigen::Matrix<float,3,10> shape_functions(Eigen::Matrix<float,1,4> evaluation_points);
    bool initialized; // init to be set after we've constructed shape functions for this element, dont need to run for every object instance
    const unsigned short ngp = 8; // 2*2*2 
    const unsigned short ngp_per_dim = 2; // 2*2*2 
    const unsigned short W = 1; 
    Eigen::Matrix<float,1,2> xhi;
    Eigen::Matrix<float,1,2> eta;
    Eigen::Matrix<float,1,2> my;

public:
    Eigen::Matrix<float,24,24> Ke;
    Eigen::Matrix<float,24,1> fe;  
    std::shared_ptr<Pid> get_pid(){return this->pid;};
    unsigned int get_id(){return id;};
    std::vector<unsigned int> get_element_dof_ids(){return dofs_id;};
    std::vector<std::shared_ptr<Node>> get_connectivity(){return this->connectivity;};
    unsigned short get_element_ndofs(){return ndofs;}
    unsigned short get_element_nnodes(){return nnodes;}
    unsigned short get_vtk_identifier(){return vtk_identifier;}
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Ke(){return Ke;}
    std::string get_element_type(){return element_type;}
    C3D8(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~C3D8();
};