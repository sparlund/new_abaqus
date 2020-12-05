#pragma once
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <array>
#include "../element.h"
#include "../pid.h"
#include "../node.h"

// CPS3 is 3 node tria shell element
class CPS3 : public Element
{
private:
    unsigned int id;
    static const std::string element_type;
    static const unsigned short nnodes   = 3;
    static const unsigned short ndofs     = 6; // 3*2
    unsigned short vtk_identifier = 5;
    std::vector<unsigned int> dofs_id;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    Eigen::Matrix<float,2,2> J;
    Eigen::Matrix<float,3,6> B;
    float A;
public:
    Eigen::Matrix<float,6,6> Ke;
    Eigen::Matrix<float,6,6> Me;
    Eigen::Matrix<float,6,1> fe;  
    static Eigen::Matrix<float,6,2> N; 
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
    CPS3(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~CPS3();
};







