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
    static const unsigned char nnodes   = 3;
    static const unsigned char ndofs     = 9; // 3*3
    std::vector<unsigned int> dofs_id;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    Eigen::Matrix<float,6,6> Ke;
    Eigen::Matrix<float,6,1> fe;
    Eigen::Matrix<float,6,6> C;
    Eigen::Matrix<float,6,6> B;
    float A;

public:
    static Eigen::Matrix<float,6,2> N; 
    std::shared_ptr<Pid> get_pid(){return this->pid;};
    unsigned int get_id(){return id;};
    std::vector<std::shared_ptr<Node>> get_connectivity(){return this->connectivity;};
    unsigned char get_element_ndofs(){return ndofs;}
    std::vector<unsigned int> get_element_dof_ids(){return dofs_id;};
    unsigned char get_element_nnodes(){return nnodes;}
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Ke(){return Ke;}
    std::string get_element_type(){return element_type;}
    S3(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~S3();
};







