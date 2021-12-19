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
    float x, y, z;
    const float original_x, original_y, original_z;
    std::vector<std::unique_ptr<Dof>> dofs;
    std::vector<Element*> connected_elements;
    Node(unsigned int global_id, float x_, float y_, float z_ = 0);
};
