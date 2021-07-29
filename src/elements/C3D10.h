#pragma once
#include <vector>
#include <memory>
#include "../../external_libs/Eigen/Dense"
#include "../element.h"
#include "../pid.h"
#include "../node.h"
#include "../Gauss.h"

class C3D10 : public Element
{
    // The target rank of Ke is 30 − 6 = 24. Since each Gauss point adds 6 to the rank up to a maximum
    // of 24, the number of Gauss points should be 4 or higher. 
public:
    void calculate_Ke();
    void calculate_Me();
    C3D10(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~C3D10();
};