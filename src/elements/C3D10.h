#pragma once
#include "element.h"
#include "../Gauss.h"
#include "../node.h"
#include "../pid.h"

#include <Eigen/Dense>
#include <memory>
#include <vector>

class C3D10 : public Element
{
private:
    const std::array<std::array<float, 3>,4>* gauss_points  = &Gauss::_3D::integration_points_4;
    const std::array<float,4>*                gauss_weights = &Gauss::_3D::gauss_weights_4;
    // The target rank of K e is 30 − 6 = 24. Since each Gauss point adds 6 to the rank up to a maximum
    // of 24, the number of Gauss points should be 4 or higher. 
public:
    void calculate_Ke();
    void calculate_Me();
    C3D10(unsigned int                        id,
          std::vector<Node*>                  connectivity,
          Pid*                                pid,
          const unsigned short                nnodes,
          const unsigned short                ndofs,
          const unsigned short                vtk_identifier,
          const unsigned short                ngp,
          const unsigned short                dimensions,
          std::string                         element_type);
};