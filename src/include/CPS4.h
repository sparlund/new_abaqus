#pragma once
#include "element.h"
#include "Gauss.h"
#include "node.h"
#include "pid.h"

#include <Eigen/Dense>
#include <memory>
#include <vector>

// CPS4 is 4 node quadrilateral element
// *---------* 
// |  x   x  |
// |         |
// |  x   x  |
// *---------*

class CPS4 : public Element
{
private:
    const std::array<std::array<float, 2>, 4>*  gauss_points  = &Gauss::_2D::integration_points_2_by_2;
    const std::array<float, 4>*                 gauss_weights = &Gauss::_2D::gauss_weights_2_by_2;
public:
    void calculate_Ke();
    void calculate_Me();
    std::vector<Segment>&                    get_segments(Node*) override;
    std::vector<Scalar>                          calculate_stress(dynMatrix,
                                                                 dynMatrix) override;
    std::vector<Scalar>                          calculate_strain(dynMatrix,
                                                                 dynMatrix) override;
    CPS4(unsigned int                        id,
         std::vector<Node*>                  connectivity,
         Pid*                                pid,
         const unsigned short                nnodes,
         const unsigned short                ndofs,
         const unsigned short                vtk_identifier,
         const unsigned short                ngp,
         const unsigned short                dimensions,
         std::string                         element_type);
};







