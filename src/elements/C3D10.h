#pragma once
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include "../element.h"
#include "../pid.h"
#include "../node.h"


class C3D10 : public Element
{
private:
    unsigned int id;
    static const std::string element_type;
    static const unsigned short nnodes  = 10;
    static const unsigned char ndofs   = 30; // 10*3
    // https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    static const unsigned short vtk_identifier = 24;
    std::vector<unsigned int> dofs_id;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    Eigen::Matrix<float,4,4> J;
    std::vector<float> detJ;
    Eigen::Matrix<float,6,10> B;
    Eigen::Matrix<float,6,6> D;
    float A;
    Eigen::Matrix<float,10,2> coord;
    // size(N) = ngp*nnodes
    Eigen::Matrix<float,4,10> N;
    // size(dNdXhi) = ngp*ngp, 4 vectors of (1x4)
    Eigen::Matrix<float,4,4> dNdXhi;
    Eigen::Matrix<float,4,4> dNdEta;
    Eigen::Matrix<float,3,10> dNdxdydz;
    // Eigen::Matrix<float,3,10> shape_functions(Eigen::Matrix<float,1,4> evaluation_points);
    bool initialized; // init to be set after we've constructed shape functions for this element, dont need to run for every object instance
    const unsigned short ngp = 4; // 
    const unsigned short W = 0.25; 
    Eigen::Matrix<float,4,1> xhi;
    Eigen::Matrix<float,4,1> eta;

public:
    Eigen::Matrix<float,30,30> Ke;
    Eigen::Matrix<float,30,1> fe;  
    std::shared_ptr<Pid> get_pid(){return this->pid;};
    unsigned int get_id(){return id;};
    std::vector<unsigned int> get_element_dof_ids(){return dofs_id;};
    std::vector<std::shared_ptr<Node>> get_connectivity(){return this->connectivity;};
    unsigned short get_element_ndofs(){return ndofs;}
    unsigned short get_element_nnodes(){return nnodes;}
    unsigned short get_vtk_identifier(){return vtk_identifier;}
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Ke(){return Ke;}
    std::string get_element_type(){return element_type;}
    C3D10(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~C3D10();
};