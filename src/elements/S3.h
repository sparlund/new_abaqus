#pragma once
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include "../element.h"
#include "../pid.h"
#include "../node.h"

// S3 is 3 node tria shell element
class S3 : public Element
{
private:
    unsigned int id;
    static const std::string element_type;
    static const unsigned short nnodes   = 3;
    static const unsigned short ndofs     = 9; // 3*3
    static const unsigned short vtk_identifier = 5;
    std::vector<unsigned int> dofs_id;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    Eigen::Matrix<float,6,6> Ke;
    Eigen::Matrix<float,6,6> Me;
    Eigen::Matrix<float,6,1> fe;
    Eigen::Matrix<float,6,6> C;
    Eigen::Matrix<float,6,6> B;
    float A;
public:
    static Eigen::Matrix<float,6,2> N; 
    std::shared_ptr<Pid> get_pid(){return this->pid;};
    unsigned int get_id(){return id;};
    std::vector<std::shared_ptr<Node>> get_connectivity(){return this->connectivity;};
    unsigned short get_element_ndofs(){return ndofs;}
    unsigned short get_element_nnodes(){return nnodes;}
    std::vector<unsigned int> get_element_dof_ids(){return dofs_id;};
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Ke(){return Ke;}
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Me(){return Me;}
    std::string get_element_type(){return element_type;}
    unsigned short get_vtk_identifier(){return vtk_identifier;}
    S3(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~S3();
};







