#pragma once
#include <vector>
#include <memory>
#include "../../external_libs/Eigen/Dense"
#include "../element.h"
#include "../pid.h"
#include "../node.h"
#include "../Gauss.h"

class C3D8 : public Element
{
public:
    void calculate_Ke();
    void calculate_Me();
    C3D8(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~C3D8();
};