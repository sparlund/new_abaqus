#pragma once
#include "dof.h"
#include <array>
#include <memory>
#include <vector>

class Element;
class Node
{
public:
    const unsigned int id;
    double x, y, z;
    const double original_x, original_y, original_z;
    std::vector<std::unique_ptr<Dof>> dofs;
    std::vector<Element*> connected_elements;
    Node(unsigned int global_id, double x_, double y_, double z_ = 0);
};
