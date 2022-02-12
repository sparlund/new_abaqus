#pragma once
#include "element.h"
#include "Gauss.h"
#include "node.h"
#include "pid.h"

#include <array>
#include <Eigen/Dense>
#include <memory>
#include <vector>

class C3D20 : public Element
{
private:
    const std::array<std::array<float, 3>,27>*  gauss_points  = &Gauss::_3D::integration_points_3_by_3_by_3;
    const std::array<float, 27>*                gauss_weights = &Gauss::_3D::gauss_weights_3_by_3_by_3;
public:
    void calculate_Ke();
    void calculate_Me();
    C3D20(unsigned int                        id,
          std::vector<Node*>                  connectivity,
          Pid*                                pid);
};