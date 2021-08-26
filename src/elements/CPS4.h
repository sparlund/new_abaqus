#pragma once
#include <vector>
#include <memory>
#include <Eigen/Dense>
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
private:
    const std::array<std::array<float, 2>, 4>*  gauss_points  = &Gauss::_2D::integration_points_2_by_2;
    const std::array<float, 4>*                 gauss_weights = &Gauss::_2D::gauss_weights_2_by_2;
public:
    void calculate_Ke();
    void calculate_Me();
    CPS4(unsigned int                        id,
         std::vector<std::shared_ptr<Node>>  connectivity,
         std::shared_ptr<Pid>                pid,
         const unsigned short                nnodes,
         const unsigned short                ndofs,
         const unsigned short                vtk_identifier,
         const unsigned short                ngp,
         const unsigned short                dimensions,
         std::string                         element_type);
};







