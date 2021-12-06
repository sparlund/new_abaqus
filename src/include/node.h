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
    const float x, y, z;
    std::vector<std::unique_ptr<Dof>> dofs;
    std::vector<Element*> connected_elements;
    Node(unsigned int global_id, float x, float y, float z = 0);
};
