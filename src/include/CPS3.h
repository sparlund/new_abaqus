#pragma once
#include "element.h"
#include "node.h"
#include "pid.h"

#include <array>
#include <Eigen/Dense>
#include <memory>
#include <vector>

// CPS3 is 3 node tria shell element. this element is shit
class CPS3 : public Element
{
public:
    void calculate_Ke();
    void calculate_Me();
    CPS3(unsigned int                        id,
         std::vector<Node*>                  connectivity,
         Pid*                                pid);
};







