#pragma once
#include "element.h"
#include "node.h"
#include "pid.h"

#include <Eigen/Dense>
#include <memory>
#include <vector>

// S3 is 3 node tria shell element
class S3 : public Element
{
public:
    void calculate_Ke();
    void calculate_Me();
    S3(unsigned int                        id,
         std::vector<Node*>                  connectivity,
         Pid*                                pid);
};







