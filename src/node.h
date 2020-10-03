#pragma once
#include <array>
class Node
{
private:
public:
    // Node has id and positional data
    float x, y, z;
    unsigned int id;
    std::array<float,6> dofs;  
    Node(unsigned int local_id, float x, float y, float z);
    ~Node();
};
