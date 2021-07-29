#pragma once
#include <vector>
#include <array>
#include <memory>
#include "../../external_libs/Eigen/Dense"
#include "../element.h"
#include "../pid.h"
#include "../node.h"
#include "../Gauss.h"

class C3D20 : public Element
{
public:
    void calculate_Ke();
    void calculate_Me();
    C3D20(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~C3D20();
};