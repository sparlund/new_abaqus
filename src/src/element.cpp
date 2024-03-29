#include "../include/element.h"
#include "../include/dof.h"
#include "../include/misc_string_functions.h"
#include "../include/node.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

void Element::calculate_f_internal(dynVector u)
{
    std::cout << "ERROR: element not supported for computing internal element forces." << std::endl;
    exit(0);
}

std::vector<Scalar> Element::calculate_stress(dynMatrix, 
                                             dynMatrix)
{
    std::cout << "ERROR: element not supported for computing stress." << std::endl;
    exit(0);
    return std::vector<Scalar>{};
}                                                                     
std::vector<Scalar> Element::calculate_strain(dynMatrix, 
                                             dynMatrix)
{
    std::cout << "ERROR: element not supported for computing strain" << std::endl;
    exit(0);
    return std::vector<Scalar>{};
}

std::vector<Segment>& Element::get_segments(Node*)
{
    std::cout << "ERROR: element not supported for *CONTACT." << std::endl;
    exit(0);
    return segments;
}

void Element::setup_dofs(){
    // add dofs to each node. can be done first now because now we know how many dofs each node should have    
    for (unsigned int i = 0; i < connectivity.size(); i++)
    {
        // Check if current node already has dofs or if we need to create
        if (connectivity.at(i)->dofs.size() != ndofs/nnodes)
        {
            // create 2 dofs
            auto x  = std::make_unique<Dof>();
            auto y  = std::make_unique<Dof>();
            // put indiviual dof id's in a list for easy access?
            dofs_id.push_back(x->id);
            dofs_id.push_back(y->id);
            // put Dof object itself in list of Dofs for element
            connectivity.at(i)->dofs.emplace_back(std::move(x));
            connectivity.at(i)->dofs.emplace_back(std::move(y));
            if(dimensions == 3){
                auto z  = std::make_unique<Dof>();
                dofs_id.push_back(z->id);
                connectivity.at(i)->dofs.emplace_back(std::move(z));
            }
        }
        else
        {

            // find dofs from node and add to dofs_id vector
            dofs_id.push_back(connectivity.at(i)->dofs.at(0)->id);
            dofs_id.push_back(connectivity.at(i)->dofs.at(1)->id);
            if (dimensions == 3)
            {
                dofs_id.push_back(connectivity.at(i)->dofs.at(2)->id);
            }
            
        }
    }
};

void Element::setup_coord(){
    // coord system is located in the middle of 
    // the element, as weights are from -1 to 1
    for (unsigned short i = 0; i < nnodes; i++)
    {
        coord(i,0) = connectivity.at(i)->x;
        coord(i,1) = connectivity.at(i)->y;
        if (dimensions == 3)
        {
            coord(i,2) = connectivity.at(i)->z;
        }
        
    }
};

Element::Element(unsigned int                        id,
                 std::vector<Node*>                  connectivity,
                 Pid*                                pid,
                 ElementType                         elementType,
                 const unsigned short                nnodes,
                 const unsigned short                ndofs,
                 const unsigned short                vtk_identifier,
                 const unsigned short                ngp,
                 const unsigned short                dimensions):
                id{id},
                connectivity{std::move(connectivity)},
                pid{pid},
                elementType{elementType},
                nnodes{nnodes},
                ndofs{ndofs},
                vtk_identifier{vtk_identifier},
                ngp{ngp},
                dimensions{dimensions}
    {
        Ke.resize(ndofs,ndofs);
        Me.resize(ndofs,ndofs);
        coord.resize(nnodes,dimensions);
        B.resize(ngp);
        for(auto& Bi: B)
        {
            Bi.resize(dimensions, ngp);
        }
        Ke.setZero();
        Me.setZero();
        // dofs can be setup now because an element won't
        // change dofs during analysis but coord cant be setup
        // because elements change shape during analysis.
        // so setup_coord is moved to calculate_{Ke,Me} for each element!
        setup_dofs();
    };

unsigned int Element::element_counter=0;

double Element::inv_div_by1(double in) const {
    if (in > std::numeric_limits<double>::min())
    {
        return 1/in;
    }
    else
    {   
        return std::max(1.25e9d,1/std::numeric_limits<double>::min());
    }
}

void Element::print_element_info_to_log() const {
    // print info to log file
    std::vector<Node*> nodes = get_connectivity();
    std::cout << "*ELEMENT: type=" << elementType << ", PID = " << get_pid()->get_name() << ", id=" << id << ", nodes=";
    for (unsigned short i = 0; i < nnodes ; i++)
    {
        std::cout << nodes.at(i)->id << ",";        
    }
    std::cout << " dofs=";
    for (unsigned short i = 0; i < nnodes ; i++)
    {   
        for (unsigned short j = 0; j < ndofs/nnodes ; j++)
        {
            std::cout << nodes.at(i)->dofs.at(j)->id << ",";
        }
    }
    std::cout << std::endl;
}

std::ostream& operator<<(std::ostream& os, const Element::ElementType t)
{
    switch (t)
    {
        case Element::ElementType::C3D10:
            os << "C3D10";
            break;
        case Element::ElementType::C3D20:
            os << "C3D20";
            break;
        case Element::ElementType::C3D8:
            os << "C3D8";
            break;
        case Element::ElementType::CPS3:
            os << "CPS3";
            break;
        case Element::ElementType::CPS4:
            os << "CPS4";
            break;
        case Element::ElementType::S4:
            os << "S4";
            break;
    }
    return os;
}

