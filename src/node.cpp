#include <iostream>
#include "node.h"
Node::Node(unsigned int id, float x, float y, float z)
{
    try
    {
        this->x = x;
        this->y = y;
        this->z = z;
        this->id = id;
    }
    catch(const std::exception& e)
    {
        std::cerr << "Error creating node: " << e.what() << '\n';
    }
    // std::cout << "node created\n";
};

Node::~Node(){//std::cout << "node destroyed\n";
};