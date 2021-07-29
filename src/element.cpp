#include <string>
#include <vector>
#include <iostream>
#include "element.h"
#include "dof.h"
#include "Gauss.h"
#include "misc_string_functions.h"
#include "node.h"


unsigned short                     Element::get_dimensions(){return dimensions;};
std::vector<std::shared_ptr<Node>> Element::get_connectivity(){return connectivity;}
std::shared_ptr<Pid>               Element::get_pid(){return pid;}
std::vector<unsigned int>          Element::get_element_dof_ids(){return dofs_id;}
unsigned short                     Element::get_element_ndofs(){return ndofs;}
unsigned short                     Element::get_element_nnodes(){return nnodes;}
unsigned int                       Element::get_id(){return id;}
std::string                        Element::get_element_type(){return element_type;}
unsigned short                     Element::get_vtk_identifier(){return vtk_identifier;}
float                              Element::get_weight(){return weight;}
float                              Element::get_volume(){return volume;}
unsigned short                     Element::get_ngp(){return ngp;};
Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Element::get_Ke(){return Ke;};
Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Element::get_Me(){return Me;}; 


Element::Element(unsigned short nnodes,
                 unsigned short ndofs,
                 unsigned short vtk_identifier,
                 unsigned short ngp,
                 unsigned short dimensions,
                 std::string element_type,
                 unsigned int id,
                 std::vector<std::shared_ptr<Node>> connectivity,
                 std::shared_ptr<Pid> pid):
                 nnodes{nnodes},
                 ndofs{ndofs},
                 vtk_identifier{vtk_identifier},
                 ngp{ngp},
                 dimensions{dimensions},
                 element_type{element_type},
                 id{id},
                 connectivity{connectivity},
                 pid{pid}
{
    // resize Me, Ke, load vector fe and coordinate vector
    fe.resize(ndofs,1);
    Me.resize(ndofs,ndofs);
    Ke.resize(ndofs,ndofs);
    // resize coord and set values
    coord.resize(nnodes,dimensions);
    // coord system is located in the middle of the element
    for (unsigned short i = 0; i < nnodes; i++)
    {
        coord(i,0) = connectivity.at(i)->x;
        coord(i,1) = connectivity.at(i)->y;
        coord(i,2) = connectivity.at(i)->z;
    }
    // set Gauss integration information
    gauss_points = Gauss::get_gauss_integration_points(dimensions,ngp);
    gauss_weights = Gauss::get_gauss_integration_weights(dimensions,ngp);
    // add dofs to each node
    for (unsigned int i = 0; i < connectivity.size(); i++)
    {
        // Check if current node already has dofs or if we need to create
        if (connectivity.at(i)->dofs.size() != ndofs/nnodes)
        {
            // create 3 dofs
            Dof x = Dof();
            Dof y = Dof();
            // put Dof object itself in list of Dofs for element
            connectivity.at(i)->dofs.push_back(x);
            connectivity.at(i)->dofs.push_back(y);
            // put indiviual dof id's in a list for easy access?
            dofs_id.push_back(x.id);
            dofs_id.push_back(y.id);
            // If 3D element do same thing for z-direction
            if(dimensions == 3){
                Dof z = Dof();
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
    print_element_info_to_log();
}

Element::~Element(){
    --element_counter;
}

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
