#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <utility>
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


void Mesh::solve(){
    // Ku=f
    u.resize(ndofs,1);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver;
    solver.compute(K);
    u = solver.solve(f);
    unsigned int node_id = global_2_local_node_id[348];
    std::shared_ptr<Node> node = nodes.at(node_id);
    float x_disp = u(node->dofs.at(0).id);
    std::cout << x_disp << "\n";
}

void Mesh::assemble(){
    // By the time of assemble we know the number of dofs --> pre-allocate K, f & solution u
    ndofs = nodes.at(0)->dofs.at(0).global_dof_id_counter; 
    K.resize(ndofs,ndofs);
    f.resize(ndofs,1);
    std::cout << "ndofs=" << ndofs << "\n";
    // assemble load vector
    for (unsigned int i = 0; i < f_to_be_added.size(); i++)
    {
        f.coeffRef(f_to_be_added.at(i).first) += f_to_be_added.at(i).second;
    }
    // assemble stiffness matrix
    std::shared_ptr<Element> current_element;
    std::cout << "----\n";
    for (unsigned int i = 0; i < elements.size(); i++)
    { 
        current_element = elements.at(i);
        std::vector<unsigned int> dofs = current_element->get_element_dof_ids();
        for(unsigned int j=0;j < dofs.size();j++){
            unsigned int dof_row = dofs.at(j);
            for(unsigned int k=0;k < dofs.size();k++){
                unsigned int dof_column = dofs.at(k);
                Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Ke = current_element->get_Ke();
                K.coeffRef(dof_row,dof_column) += Ke(j,k);
            }
        }
    }
    // boundary conditions:
    for (unsigned int i = 0; i < bc.size(); i++)
    {
        // if current bc==0, make all values in the corresponding row and column in K to zero
        if (bc.at(i).second == 0.0f)
        {
            K.row(bc.at(i).first) *= 0;
            K.col(bc.at(i).first) *= 0;
        }
    }
    std::cout << "K:\n";
    std::cout << K;
    std::cout << "\n";
    std::cout << "f:\n";
    std::cout << f;
    std::cout << "---\n";
}



void Mesh::add_mid(std::unordered_map<std::string, std::string> options){
    // Create new MID based on input on *MATERIAL data line.
    // Other functions add mtrl data such as density etc
    // Mid temp_mid(mid_counter,options["NAME"]);
    // mid_counter++;
    
};



Mesh::Mesh(){
    // init object to hold the mesh data.
}

void Mesh::add_boundary(std::string line,std::unordered_map<std::string, std::string> options){
    std::string type = options["TYPE"];
    std::vector<std::string> data = misc::split_on(line,',');
    if (type=="DISPLACEMENT")
    {
        // node,dof_from,dof_to,magnitude
        unsigned int global_node_id = std::stoi(data.at(0));
        unsigned int dof_from       = std::stoi(data.at(1));
        unsigned int dof_to         = std::stoi(data.at(2));
        float        magnitude      = std::stof(data.at(3));        
        // dof is given above in the local coord system
        std::shared_ptr<Node> node = nodes.at(global_2_local_node_id[global_node_id]);
        for (unsigned int i = dof_from; i <= dof_to; i++)
        {
            // abaqus starts counting dof's at 1, but vectors start at 0
            bc.push_back(std::make_pair(node->dofs.at(i-1).id,magnitude));
        }
        
    }
    
}

void Mesh::add_load(std::string line){
    // for *cload line is: node, dof, magnitude
    std::vector<std::string> data   = misc::split_on(line,',');   
    unsigned int global_node_id     = std::stoi(data.at(0));
    unsigned int local_dof          = std::stoi(data.at(1));
    float magnitude                 = std::stof(data.at(2));
    std::shared_ptr<Node> node = nodes.at(global_2_local_node_id[global_node_id]);
    unsigned int global_dof = node->dofs.at(local_dof-1).id;
    // we don't know complete number of dofs yet,
    // so we can't add directly to global load vector, but have
    // to store it here meanwhile
    f_to_be_added.push_back(std::make_pair(global_dof,magnitude));    
}

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
            else if (keyword=="*BOUNDARY"){
                getline(input_file, line);
                row_counter++;
                bool inner_loop_keyword = true;
                while (inner_loop_keyword == true){
                // Ignore if it's a comment! Still on same keyword.
                    if (misc::is_comment(line) == false){
                        add_boundary(line,options);
                    }
                    // Want to peek next line, if it's a keyword or empty line we break
                    // the while loop and start over!
                    unsigned int previous_pos = input_file.tellg();
                    getline(input_file, line);
                    if (misc::is_keyword(line) == true or line.empty() == true)
                    {
                        input_file.seekg(previous_pos);
                        inner_loop_keyword = false;
                    }
                }
            }
            else if (keyword == "*CLOAD")
            {
                getline(input_file, line);
                row_counter++;
                bool inner_loop_keyword = true;
                while (inner_loop_keyword == true){
                // Ignore if it's a comment! Still on same keyword.
                    if (misc::is_comment(line) == false){
                        add_load(line);
                    }
                    // Want to peek next line, if it's a keyword or empty line we break
                    // the while loop and start over!
                    unsigned int previous_pos = input_file.tellg();
                    getline(input_file, line);
                    if (misc::is_keyword(line) == true or line.empty() == true)
                    {
                        input_file.seekg(previous_pos);
                        inner_loop_keyword = false;
                    }
                }
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
    if (type == "S3")
    {
        element = std::shared_ptr<Element>(new S3(element_id,element_connectivity,pid));

    }
    else if (type == "S2")
    {
        element = std::shared_ptr<Element>(new S2(element_id,element_connectivity,pid));
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
    global_2_local_node_id[id] = nodes.size()-1;  
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
        std::cout << nodes.at(i)->id << ": "<< nodes.at(i)->x << ", " << nodes.at(i)->y << ", " << nodes.at(i)->z << "\n";
    }
    std::cout << "- ELEMENTS -\n"; 
    for (unsigned int i = 0; i < elements.size(); i++)
    {
        std::cout << elements.at(i)->get_id() << "   nodes (";
        std::vector<std::shared_ptr<Node>> connectivity = elements.at(i)->get_connectivity();
        for (unsigned int k = 0; k < connectivity.size()  ; k++)
        {
            std::cout << connectivity.at(k)->id << ", ";
        }
        std::cout << ")\n";
        
    }
    
};

Mesh::~Mesh(){};


