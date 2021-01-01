#include <string>
#include <vector>
#include <iostream>
#include "element.h"
#include "node.h"
#include "dof.h"
#include "misc_string_functions.h"

unsigned int Element::element_counter=0;

float Element::inv_div_by1(float in){
    if (in > std::numeric_limits<float>::min())
    {
        return 1/in;
    }
    else
    {   
        return std::max(1.25e9f,1/std::numeric_limits<float>::min());
    }
}

void Element::print_element_info_to_log(){
    // print info to log file
    std::vector<std::shared_ptr<Node>> nodes = get_connectivity();
    std::cout << "*ELEMENT: type=" << get_element_type() << ", id=" << get_id() << ", nodes=";
    for (unsigned short i = 0; i < get_element_nnodes() ; i++)
    {
        std::cout << nodes.at(i)->id << ",";        
    }
    std::cout << " dofs=";
    for (unsigned short i = 0; i < get_element_nnodes() ; i++)
    {   
        for (unsigned short j = 0; j < get_element_ndofs()/get_element_nnodes() ; j++)
        {
            std::cout << nodes.at(i)->dofs.at(j).id << ",";
        }
    }
    std::cout << std::endl;
}
