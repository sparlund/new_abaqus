#pragma once
#include "dof.h"
#include <array>
#include <vector>

class Node
{
private:
    static unsigned int node_counter;
public:
    // Node has id and positional data
    const unsigned int id;
    const float x, y, z;
    std::vector<Dof> dofs;  
    unsigned int get_node_counter();
    Node(unsigned int global_id, float x, float y, float z);
};
