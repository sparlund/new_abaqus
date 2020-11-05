#pragma once
#include <array>
#include <vector>
#include "dof.h"

class Node
{
private:
public:
    // Node has id and positional data
    unsigned int id;
    float x, y, z;
    std::vector<Dof> dofs;  
    static unsigned int node_counter;
    Node(unsigned int local_id, float x, float y, float z);
    ~Node();
};
