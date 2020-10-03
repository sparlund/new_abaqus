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
    std::vector<std::shared_ptr<Node>> get_connectivity(){return this->connectivity;};
    S3(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~S3();
};






