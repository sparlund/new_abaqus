#pragma once
#include <vector>
#include <memory>
#include "../../external_libs/Eigen/Dense"
#include <array>
#include "../element.h"
#include "../pid.h"
#include "../node.h"

// CPS3 is 3 node tria shell element
class CPS3 : public Element
{  
public:
    void calculate_Ke();
    void calculate_Me();
    CPS3(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~CPS3();
};







