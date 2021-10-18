#pragma once
#include "element.h"
#include "../node.h"
#include "../pid.h"

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
         Pid*                                pid,
         const unsigned short                nnodes,
         const unsigned short                ndofs,
         const unsigned short                vtk_identifier,
         const unsigned short                ngp,
         const unsigned short                dimensions,
         std::string                         element_type);
};







