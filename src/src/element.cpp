#include "../include/element.h"
#include "../include/dof.h"
#include "../include/misc_string_functions.h"
#include "../include/node.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

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

Element::Element(unsigned int                       id,
                std::vector<Node*>                  connectivity,
                Pid*                                pid,
                unsigned short                      nnodes,
                unsigned short                      ndofs,
                unsigned short                      vtk_identifier,
                unsigned short                      ngp,
                unsigned short                      dimensions,
                std::string                         element_type):
                id{id},
                connectivity{std::move(connectivity)},
                pid{pid},
                nnodes{nnodes},
                ndofs{ndofs},
                vtk_identifier{vtk_identifier},
                ngp{ngp},
                dimensions{dimensions},
                element_type{std::move(element_type)}{
    Ke.resize(ndofs,ndofs);
    Me.resize(ndofs,ndofs);
    coord.resize(nnodes,dimensions);
    Ke.setZero();
    Me.setZero();
    setup_dofs();
    setup_coord();
    };

std::vector<Node*>                                 Element::get_connectivity() const {return connectivity;}
Pid*                                               Element::get_pid() const {return pid;}
std::vector<unsigned int>                          Element::get_element_dof_ids() const {return dofs_id;}
unsigned short                                     Element::get_element_ndofs() const {return ndofs;}
unsigned short                                     Element::get_element_nnodes() const {return nnodes;}
unsigned int                                       Element::get_id() const {return id;}
unsigned short                                     Element::get_dimensions() const {return dimensions;};
std::string                                        Element::get_element_type() const {return element_type;}
unsigned short                                     Element::get_vtk_identifier() const {return vtk_identifier;}
float                                              Element::get_weight() const {return weight;}
float                                              Element::get_volume() const {return volume;}
unsigned short                                     Element::get_ngp() const {return ngp;};
unsigned int                                       Element::get_element_counter() const {return element_counter;};
Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Element::get_Ke() const {return Ke;};
Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Element::get_Me() const {return Me;}; 

unsigned int Element::element_counter=0;

float Element::inv_div_by1(float in) const {
    if (in > std::numeric_limits<float>::min())
    {
        return 1/in;
    }
    else
    {   
        return std::max(1.25e9f,1/std::numeric_limits<float>::min());
    }
}

void Element::print_element_info_to_log() const {
    // print info to log file
    std::vector<Node*> nodes = get_connectivity();
    std::cout << "*ELEMENT: type=" << get_element_type() << ", PID = " << get_pid()->get_name() << ", id=" << get_id() << ", nodes=";
    for (unsigned short i = 0; i < get_element_nnodes() ; i++)
    {
        std::cout << nodes.at(i)->id << ",";        
    }
    std::cout << " dofs=";
    for (unsigned short i = 0; i < get_element_nnodes() ; i++)
    {   
        for (unsigned short j = 0; j < get_element_ndofs()/get_element_nnodes() ; j++)
        {
            std::cout << nodes.at(i)->dofs.at(j)->id << ",";
        }
    }
    std::cout << std::endl;
}
