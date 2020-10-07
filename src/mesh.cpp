#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <memory>
#include <Eigen/Dense>
#include "mesh.h"
#include "mid.h"
#include "pid.h"
#include "node.h"
#include "dof.h"
#include "element.h"
#include "elements/S3.h"
#include "elements/S2.h"
#include "misc_string_functions.h"

// unsigned int Mesh::ndofs = 0;

void Mesh::assemble(){
    // By the time of assemble we know the number of dofs --> pre-allocate K, f & solution u
    // std::cout << ndofs << "\n";
    K.resize(ndofs,ndofs);
    K.Zero(ndofs,ndofs);    
    f.resize(ndofs,1);
    f.Zero(ndofs,1);
    u.resize(ndofs,1);
    u.Zero(ndofs,1);
    // Loop over all elements!
    std::shared_ptr<Element> current_element;
    std::cout << ":)\n"; 
    for (unsigned int i = 0; i < elements.size(); i++)
    { 
        current_element = elements.at(i);
        std::vector<std::shared_ptr<Node>> current_element_nodes =  current_element->get_connectivity();        
        std::cout << current_element->get_Ke() << "\n";
        // K = K + current_element->get_Ke();
        // f = f + current_element->fe();
    }
    

}



void Mesh::add_mid(std::unordered_map<std::string, std::string> options){
    // Create new MID based on input on *MATERIAL data line.
    // Other functions add mtrl data such as density etc
    // Mid temp_mid(mid_counter,options["NAME"]);
    // mid_counter++;
    
};



Mesh::Mesh(){
    // init object to hold the mesh data.
};


void Mesh::read_file(std::string filename){
    misc::append_newline_to_textfile(filename);
    // Create input stream object
    std::ifstream input_file(filename);
    std::vector<std::string> filename_split = misc::split_on(filename,'/');   
    std::string line;
    unsigned int row_counter = 0;
    while (getline(input_file, line))
    {
        row_counter++;
        if (misc::is_keyword(line) == true){
            std::string keyword = misc::split_on(line, ',').at(0);
            // extract a map of parameters and their values
            std::unordered_map<std::string,std::string> options = misc::options_map(line);

            if (keyword == "*SOLID SECTION" or keyword == "*SHELL SECTION")
            {
                add_pid(options);
            }     
            else if (keyword == "*MATERIAL")
            {
                // a complete material is typically described over several
                // lines so must make a complicated loop here..
                add_mid(options);
                // while (/* condition */)
                // {
                //     /* code */
                // }
                
            }
                   
            else if (keyword == "*NODE" or keyword == "*ELEMENT")
            {
                getline(input_file, line);
                row_counter++;
                bool inner_loop_keyword = true;
                while (inner_loop_keyword == true){
                // Ignore if it's a comment! Still on same keyword.
                    if (misc::is_comment(line) == false){
                        if (keyword == "*NODE")
                        {
                            add_node(line,options);
                        }
                        else if (keyword == "*ELEMENT")
                        {
                            add_element(line,options);
                        }
                    }
                    // Want to peek next line, if it's a keyword or empty line we break
                    // the while loop and start over!
                    unsigned int previous_pos = input_file.tellg();
                    getline(input_file, line);
                    // std::cout << line << ", is  keyword?"<< misc::is_keyword(line) << "\n";
                    if (misc::is_keyword(line) == true or line.empty() == true)
                    {
                        input_file.seekg(previous_pos);
                        inner_loop_keyword = false;
                    }
                }                
            }
            if (input_file.eof() == true)
            {
                break;
            }
            row_counter++;   
        }
    }       
};

void Mesh::add_pid(std::unordered_map<std::string, std::string> options){
    
    // Create random color
    // std::cout << colors.size() << "\n";
    std::string pid_name = options["ELSET"];
    // TODO: find MID pointer instead!!
    std::string mid_name = options["MATERIAL"];

    // Create new pid
    std::shared_ptr<Pid> pid  = std::shared_ptr<Pid>(new Pid(pid_name,mid_name));
    pids.push_back(pid);
    pid_map[pid_name] = pid;

};


void Mesh::add_element(std::string line,std::unordered_map<std::string,std::string> options){
    // these options need to be available to create an element
    std::string type = options["TYPE"];
    std::string pid_name = options["ELSET"];
    std::shared_ptr<Pid> pid = pid_map[pid_name];
    std::vector<std::string> dataline_items = misc::split_on(line,',');
    // std::cout << "added data line: " << line << "\n";
    unsigned int element_id = std::stoi(dataline_items.at(0));
    // Find node pointers for each node
    std::vector<std::shared_ptr<Node>> element_connectivity;
    for (unsigned int i = 1; i < dataline_items.size(); i++)
    {
        unsigned int node_id = std::stoi(dataline_items.at(i));
        element_connectivity.push_back(node_id_2_node_pointer[node_id]);
    }
    std::shared_ptr<Element> element;
    if (type.compare("S3"))
    {
        element = std::shared_ptr<Element>(new S3(element_id,element_connectivity,pid));

    }
    else if (type.compare("S2"))
    {
        element = std::shared_ptr<Element>(new S2(element_id,element_connectivity,pid));
    }
    // Need to add correct amount of dof's to each node, now that we know what
    // type of element have. Before this we don't know how many dofs each node need
    for (unsigned char i = 0; i < element->get_element_nnodes() ; i++)
    {
        for (unsigned char j = 0; j < element->get_element_ndofs(); j++)
        {
            Dof dof;
            // for each node in current element
            // add correct number of dofs to current node
            element_connectivity.at(i)->dofs.push_back(dof);
            ndofs++;
        }
    }
    
    elements.push_back(element);
    // else{
    //     std::cout << "Element of type " << type << " not supported.\n";
    // }
    
    };
void Mesh::add_node(std::string line,std::unordered_map<std::string, std::string> options){
    // aba docs:
    // Data line to define the node
    // Node number.
    // First coordinate of the node.
    // Second coordinate of the node.
    // Third coordinate of the node.
    // First direction cosine of the normal at the node (optional).
    // Second direction cosine of the normal at the node (optional). For nodes entered in a cylindrical or spherical system, this entry is an angle given in degrees.
    // Third direction cosine of the normal at the node (optional). For nodes entered in a spherical system, this entry is an angle given in degrees.
    
    // another node added to the Mesh!
    std::vector<std::string> dataline_items = misc::split_on(line,',');    
    // TODO: add support for more node options like coordinate system and stuff
    unsigned int             id = std::stoi(dataline_items.at(0));
    float                     x = std::stof(dataline_items.at(1));
    float                     y = std::stof(dataline_items.at(2));
    float                     z = std::stof(dataline_items.at(3));
    std::shared_ptr<Node> node  = std::shared_ptr<Node>(new Node(id,x,y,z));
    nodes.push_back(node);
    node_id_2_node_pointer[id] = node;
};

void Mesh::about(){  
    std::cout << "- PID -\n"; 
    for (unsigned int i = 0; i < pids.size(); i++)
    {
        std::cout << pids.at(i)->name << "\n";
    }
    std::cout << "- NODES -\n";   
    for (unsigned int i = 0; i < nodes.size(); i++)
    {
        std::cout << nodes.at(i)->x << ", " << nodes.at(i)->y << ", " << nodes.at(i)->z << "\n";
    }
    std::cout << "- ELEMENTS -\n"; 
    for (unsigned int i = 0; i < elements.size(); i++)
    {
        std::cout << elements.at(i)->id << "\n";
    }
    
};

Mesh::~Mesh(){};


