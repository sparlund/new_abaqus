#include <iostream>
#include "../include/node.h"
#include "../include/dof.h"

Node::Node(unsigned int id, float x_, float y_, float z_):id{id},original_x{x_},original_y{y_},original_z{z_}
{
    x = x_;
    y = y_;
    z = z_;
    // dofs will be added later to the node when we know what type of element will use it!     
    std::cout << "*NODE: id=" << id << ", x=" << x << ", y=" << y << ", z=" << z << std::endl;   
};
