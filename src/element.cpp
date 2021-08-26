#include "dof.h"
#include "element.h"
#include "misc_string_functions.h"
#include "node.h"
#include <iostream>
#include <string>
#include <vector>

void Element::setup_dofs(){
    // add dofs to each node. can be done first now because now we know how many dofs each node should have    
    for (unsigned int i = 0; i < connectivity.size(); i++)
    {
        // Check if current node already has dofs or if we need to create
        if (connectivity.at(i)->dofs.size() != ndofs/nnodes)
        {
            // create 3 dofs
            Dof x = Dof();
            Dof y = Dof();
            Dof z = Dof();
            // put Dof object itself in list of Dofs for element
            connectivity.at(i)->dofs.push_back(x);
            connectivity.at(i)->dofs.push_back(y);
            // put indiviual dof id's in a list for easy access?
            dofs_id.push_back(x.id);
            dofs_id.push_back(y.id);
            if(dimensions == 3){
            connectivity.at(i)->dofs.push_back(z);
            dofs_id.push_back(z.id);
            }
        }
        else
        {
            // find dofs from node and add to dofs_id vector
            dofs_id.push_back(connectivity.at(i)->dofs.at(0).id);
            dofs_id.push_back(connectivity.at(i)->dofs.at(1).id);
            if(dimensions == 3){
                dofs_id.push_back(connectivity.at(i)->dofs.at(2).id);
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
        if(dimensions == 3){
            coord(i,2) = connectivity.at(i)->z;
        }
    }
};

Element::Element(unsigned int                       id,
                std::vector<std::shared_ptr<Node>>  connectivity,
                std::shared_ptr<Pid>                pid,
                unsigned short                      nnodes,
                unsigned short                      ndofs,
                unsigned short                      vtk_identifier,
                unsigned short                      ngp,
                unsigned short                      dimensions,
                std::string                         element_type):
                id{id},
                connectivity{connectivity},
                pid{pid},
                nnodes{nnodes},
                ndofs{ndofs},
                vtk_identifier{vtk_identifier},
                ngp{ngp},
                dimensions{dimensions},
                element_type{element_type}{print_element_info_to_log();};

std::vector<std::shared_ptr<Node>>                 Element::get_connectivity(){return connectivity;}
std::shared_ptr<Pid>                               Element::get_pid(){return pid;}
std::vector<unsigned int>                          Element::get_element_dof_ids(){return dofs_id;}
unsigned short                                     Element::get_element_ndofs(){return ndofs;}
unsigned short                                     Element::get_element_nnodes(){return nnodes;}
unsigned int                                       Element::get_id(){return id;}
std::string                                        Element::get_element_type(){return element_type;}
unsigned short                                     Element::get_vtk_identifier(){return vtk_identifier;}
float                                              Element::get_weight(){return weight;}
float                                              Element::get_volume(){return volume;}
unsigned short                                     Element::get_ngp(){return ngp;};
unsigned int                                       Element::get_element_counter(){return element_counter;};
Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Element::get_Ke(){return Ke;};
Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Element::get_Me(){return Me;}; 

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
