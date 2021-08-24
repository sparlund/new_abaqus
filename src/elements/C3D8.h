#pragma once
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include "../element.h"
#include "../pid.h"
#include "../node.h"
#include "../Gauss.h"

class C3D8 : public Element
{
protected:
    unsigned int id;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    std::vector<unsigned int> dofs_id;
    std::vector<float> detJ;
    float area, volume, weight;
    static const std::string element_type;
    const std::array<std::array<float, 3>,8>* gauss_points = &Gauss::_3D::integration_points_2_by_2_by_2;
    const std::array<float, 8>* gauss_weights = &Gauss::_3D::gauss_weights_2_by_2_by_2;
public:
    void calculate_Ke();
    void calculate_Me();
    C3D8(unsigned int                        id,
         std::vector<std::shared_ptr<Node>>  connectivity,
         std::shared_ptr<Pid>                pid,
         const unsigned short                nnodes,
         const unsigned short                ndofs,
         const unsigned short                vtk_identifier,
         const unsigned short                ngp,
         const unsigned short                dimensions);
};