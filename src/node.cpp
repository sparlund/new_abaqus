#include <iostream>
#include "node.h"
#include "dof.h"
Node::Node(unsigned int id, float x, float y, float z):id(id),x(x),y(y),z(z)
{
    // dofs will be added later to the node when we know what type of element will use it!     
    // std::cout << "node " << id << "created!" << std::endl;   
};

Node::~Node(){//std::cout << "node destroyed\n";
};