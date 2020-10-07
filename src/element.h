#pragma once
#include <vector>
#include <memory>
#include <Eigen/Sparse>
#include "node.h"
#include "pid.h"

// virtual class
class Element{
private:
    
public:
    unsigned int id;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    virtual std::vector<std::shared_ptr<Node>> get_connectivity()=0;
    virtual unsigned char get_element_ndofs()=0;
    virtual unsigned char get_element_nnodes()=0;
    virtual std::shared_ptr<Pid> get_pid()=0;
    virtual Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>get_Ke()=0; 
    virtual ~Element(){};
    Element(){};
};
