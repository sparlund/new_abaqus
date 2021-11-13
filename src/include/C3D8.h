#pragma once
#include "element.h"
#include "Gauss.h"
#include "node.h"
#include "pid.h"

#include <Eigen/Dense>
#include <memory>
#include <vector>

class C3D8 : public Element
{
private:
    const std::array<std::array<float, 3>,8>* gauss_points  = &Gauss::_3D::integration_points_2_by_2_by_2;
    const std::array<float, 8>*               gauss_weights = &Gauss::_3D::gauss_weights_2_by_2_by_2;
public:
    void calculate_Ke();
    void calculate_Me();
    C3D8(unsigned int                        id,
         std::vector<Node*>                  connectivity,
         Pid*                                pid,
         const unsigned short                nnodes,
         const unsigned short                ndofs,
         const unsigned short                vtk_identifier,
         const unsigned short                ngp,
         const unsigned short                dimensions,
         std::string                         element_type);
};