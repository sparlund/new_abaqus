#pragma once
#include <vector>
#include <memory>
#include "../../external_libs/Eigen/Dense"
#include "../element.h"
#include "../pid.h"
#include "../node.h"
#include "../Gauss.h"

// CPS4 is 4 node quadrilateral element
// *---------* 
// |  x   x  |
// |         |
// |  x   x  |
// *---------*

class CPS4 : public Element
{
public: 
    void calculate_Ke();
    void calculate_Me();
    CPS4(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~CPS4();
};







