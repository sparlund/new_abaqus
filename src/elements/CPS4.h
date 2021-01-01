#pragma once
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include "../element.h"
#include "../pid.h"
#include "../node.h"

// CPS4 is 4 node quadrilateral element
// *---------* 
// |  x   x  |
// |         |
// |  x   x  |
// *---------*

class CPS4 : public Element
{
private:
    unsigned int id;
    static const std::string element_type;
    static const unsigned short nnodes  = 4;
    static const unsigned char ndofs   = 8; // 4*2
    // https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    static const unsigned short vtk_identifier = 9;
    std::vector<unsigned int> dofs_id;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    Eigen::Matrix<float,2,2> J;
    Eigen::Matrix<float,3,8> B;
    Eigen::Matrix<float,3,3> D;
    float A;
    Eigen::Matrix<float,4,2> coord;
    // size(N) = ngp*nnodes
    Eigen::Matrix<float,4,4> N;
    // size(dNdXhi) = 4 vectors of (1x4)
    Eigen::Matrix<float,4,4> dNdXhi;
    Eigen::Matrix<float,4,4> dNdEta;
    Eigen::Matrix<float,2,4> dNdxdy;
    bool initialized; // init to be set after we've constructed shape functions for this element, dont need to run for every object instance
    const unsigned char ngp = 4; // to use 9 ngp use element CPS8!
    const unsigned char W = 1;
    Eigen::Matrix<float,4,1> xhi;
    Eigen::Matrix<float,4,1> eta;

public:
    Eigen::Matrix<float,8,8> Ke;
    Eigen::Matrix<float,8,8> Me;
    Eigen::Matrix<float,8,1> fe;  
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
    CPS4(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~CPS4();
};







