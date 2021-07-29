#pragma once
#include <vector>
#include <memory>
#include "../../external_libs/Eigen/Dense"
#include "../element.h"
#include "../pid.h"
#include "../node.h"

// S3 is 3 node tria shell element
class S3 : public Element
{
public:
    void calculate_Ke();
    void calculate_Me();
    S3(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~S3();
};







