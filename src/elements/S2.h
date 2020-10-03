#pragma once
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include "../element.h"
#include "../pid.h"
#include "../node.h"

// S2 is 3 node tria shell element
class S2 : public Element
{
private:
    unsigned int id;
    static const unsigned char dofs     = 6; // 3*2
    static const unsigned char nnodes   = 3;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    Eigen::Matrix<float,6,6> Ke;
    Eigen::Matrix<float,6,1> fe;  
    Eigen::Matrix<float,6,6> C;
    Eigen::Matrix<float,3,6> B;
    float A;

public:
    static Eigen::Matrix<float,6,2> N; 
    std::shared_ptr<Pid> get_pid(){return this->pid;};
    std::vector<std::shared_ptr<Node>> get_connectivity(){return this->connectivity;};
    unsigned char get_element_ndofs(){return dofs;}
    unsigned char get_element_nnodes(){return nnodes;}
    S2(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~S2();
};







