#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <utility>
#include <cmath>
#include <fstream>
#include <ctime>
#include <Eigen/Eigenvalues>
#include <stdlib.h>
#include "mesh.h"
#include "mid.h"
#include "pid.h"
#include "node.h"
#include "dof.h"
#include "element.h"
#include "elements/S3.h"
#include "elements/CPS3.h"
#include "elements/CPS4.h"
#include "elements/C3D10.h"
#include "elements/C3D8.h"
#include "misc_string_functions.h"
// for eigenfrequency solver
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseCholesky.h>

void Mesh::export_2_vtk(){
    // output results to vtk format for the ParaView post processor
    // https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    // http://profs.sci.univr.it/~caliari/aa1314/advanced_numerical_analysis/Pellegrini.pdf
    std::clock_t clock_export;
    float duration_clock_export;
    clock_export = std::clock();
    std::cout << "---    Starting to export VTK results file    ---" << std::endl;
    std::vector<std::string> split1 = misc::split_on(filename,'/');
    std::vector<std::string> split2 = misc::split_on(split1.at(1),'.');
    std::string output_filename = split2.at(0) + ".vtk";
    std::ofstream output(output_filename);
    output << "# vtk DataFile Version 3.0" << std::endl;
    output << "vtk output" << std::endl;
    // ascii or binary format:
    output << "ASCII" << std::endl;
    output << "DATASET UNSTRUCTURED_GRID" << std::endl;
    // no need to use double precision on printing results
    output << "POINTS " << nodes.size() << " float" << std::endl;
    // print nodes to result file
    for (unsigned int i = 0; i < nodes.size() ; i++)
    {
        output << nodes.at(i)->x << " " << nodes.at(i)->y << " " << nodes.at(i)->z << std::endl;
    }
    unsigned int total_amount_of_used_nodes=0;
    for (unsigned int i = 0; i < elements.size(); i++)
    {
        total_amount_of_used_nodes += elements.at(i)->get_element_nnodes();
    }    
    output << "CELLS " << elements.size() << " " << elements.size() + total_amount_of_used_nodes  << std::endl;
    for (unsigned int i = 0; i < elements.size(); i++)
    {
        output << elements.at(i)->get_element_nnodes() << " ";
        std::vector<std::shared_ptr<Node>> connectivity =  elements.at(i)->get_connectivity();
        for (unsigned int  j = 0; j < connectivity.size(); j++)
        {
            output << global_2_local_node_id[connectivity.at(j)->id] << " ";
        }
        output << std::endl;
    }
    
    output << "CELL_TYPES " << elements.size() << std::endl;
    for (unsigned int i = 0; i < elements.size(); i++)
    {
        output << elements.at(i)->get_vtk_identifier() << std::endl;
    }
    // point_data, i.e displacement on nodes
    output << "POINT_DATA " << nodes.size() << std::endl;
    if (this->static_analysis == true)
    {
        output << "FIELD displacement 1" << std::endl;
        // data name, number of values (3 for 3D & 2D, Z is zero for 2D, datatype)
        output << "displacement 3 " << nodes.size() << " float" << std::endl;
        std::shared_ptr<Node> current_node;
        for (unsigned int i = 0; i < nodes.size(); i++)
        {
            current_node = nodes.at(i);        
            output << u(current_node->dofs.at(0).id) << " " << u(current_node->dofs.at(1).id) << " " << u(current_node->dofs.at(2).id) << std::endl;
        }
    }
    if (this->eigenvalue_analysis == true)
    {
        // scale eigenvectors before saving to vtk file
        // normalize against largest value in eigenvector matrix
        eigenvectors /= eigenvectors.maxCoeff();
        for (unsigned int mode = 0; mode < number_of_modes_to_find; mode++)
        {
            output << "FIELD mode" << mode << " 1" << std::endl;
            // data name, number of values (3 for 3D & 2D, Z is zero for 2D, datatype)
            output << "mode" << mode << " 3 " << nodes.size() << " float" << std::endl;
            std::shared_ptr<Node> current_node;
            for (unsigned int i = 0; i < nodes.size(); i++)
            {
                    current_node = nodes.at(i);
                    std::cout << "node id=" << current_node->id << std::endl;
                    std::cout << "dof1=" << current_node->dofs.at(0).id << ", eigenvector value for that dof="<< eigenvectors(current_node->dofs.at(0).id,0) << std::endl;
                    std::cout << "--" << std::endl;
                    output << current_node->x*eigenvectors(current_node->dofs.at(0).id,mode)
                           << " " 
                           << current_node->y*eigenvectors(current_node->dofs.at(1).id,mode)
                           << " " 
                           << current_node->z*eigenvectors(current_node->dofs.at(2).id,mode) << std::endl;
            }
            
        }
        
    }
    
    
    // close file!
    output.close();
    duration_clock_export = ( std::clock() - clock_export ) / (float) CLOCKS_PER_SEC;
    std::cout << "---    Exported to VTK format in " << duration_clock_export << " seconds (wallclock time)   ---" << std::endl;
}

void Mesh::solve(){
    if (eigenvalue_analysis == true)
    {
        solve_eigenfrequency();
    }
    if (static_analysis == true)
    {
        solve_static();
    }
}

void Mesh::solve_eigenfrequency(){
    std::cout << "---    Starting to solve eigenvalue problem    ---" << std::endl;
    std::clock_t clock_solve;
    float duration_clock_solve;
    clock_solve = std::clock();  
    Spectra::SparseSymMatProd<float> opK(K);
    Spectra::SparseCholesky<float> opM(M);
    // ncv must satisfy nev < ncv <= n, n is the size of matrix. But what is a good value?
    unsigned int ncv = (2*number_of_modes_to_find) - 1;

    // ncv that controls the convergence speed of the algorithm.
    // Typically a larger `ncv` means faster convergence, but it may
    // also result in greater memory use and more matrix operations
    // in each iteration. This parameter must satisfy nev < ncv < n,
    // and is advised to take ncv< 2*nev.
    if ( (this->number_of_modes_to_find < ncv && ncv < ndofs) == false)
    {
        std::cout << "ERROR: nev < ncv <= n not satisfied." << std::endl;
    }
    std::cout << ":)" << std::endl;
    Spectra::SymGEigsSolver<float,
                            Spectra::SMALLEST_MAGN,
                            Spectra::SparseSymMatProd<float>,
                            Spectra::SparseCholesky<float>,
                            Spectra::GEIGS_CHOLESKY> geigs(&opK, &opM, number_of_modes_to_find, ncv);
    std::cout << ":))" << std::endl;
    // Initialize and compute
    geigs.init();
    std::cout << ":)))" << std::endl;
    std::cout << geigs.info() << std::endl;
    geigs.compute();
    std::cout << geigs.info() << std::endl;
    std::cout << ":))))" << std::endl;
    // Retrieve results
    if(geigs.info() == Spectra::SUCCESSFUL)
    {
        // omega = sqrt(lambda)
        // frequency = omega/(2*pi) 
        std::cout << ":))))" << std::endl;
        this->eigenvalues = geigs.eigenvalues();
        this->eigenvectors = geigs.eigenvectors();
        std::cout << eigenvalues.size() << std::endl;
        std::cout << eigenvectors.size() << std::endl;
    }
    else
    {
        std::cout << "ERROR: couldn't solve eigenvalue problem" << std::endl;
        return;
    }
    std::cout << "Generalized eigenvalues found:\n" << eigenvalues << std::endl;
    std::cout << "Generalized eigenvectors found:\n" << eigenvectors << std::endl;
    duration_clock_solve = ( std::clock() - clock_solve ) / (float) CLOCKS_PER_SEC;
    std::cout << "---    Solution to eigenvalue problem found in " << duration_clock_solve << " seconds (wallclock time)    ---" << std::endl;
}

void Mesh::solve_static(){
    std::cout << "---    Starting to solve linear problem Ku=f    ---" << std::endl;
    std::clock_t clock_solve;
    float duration_clock_solve;
    clock_solve = std::clock();  
    // Ku=f
    Eigen::SparseLU<Eigen::SparseMatrix<float>> solver;
    solver.compute(K);
    u = solver.solve(f);
    std::cout << "u:" << std::endl;
    std::cout << Eigen::Matrix<float,1,Eigen::Dynamic>(u).transpose() << std::endl;
    duration_clock_solve = ( std::clock() - clock_solve ) / (float) CLOCKS_PER_SEC;
    std::cout << "---    Solution to linear problem found in " << duration_clock_solve << " seconds (wallclock time)    ---" << std::endl;
}


void Mesh::assemble(){
    // By the time of assemble we know the number of dofs --> pre-allocate K, f & solution u
    std::cout << "---    Starting stiffness matrix assembly    ---" << std::endl;
    std::cout << "          Mesh size:" << std::endl;
    std::cout << "               nodes = " << nodes.size() << std::endl;
    std::cout << "            elements = " << elements.size() << std::endl;
    // access arbitrary dof object from arbitrary node object and check the static member to see total ndofs     
    ndofs = nodes.at(0)->dofs.at(0).global_dof_id_counter; 
    std::cout << "  degrees of freedom = " << ndofs << std::endl;
    std::clock_t clock_assemble;
    float duration_assemble;
    clock_assemble = std::clock();
    K.resize(ndofs,ndofs);
    M.resize(ndofs,ndofs);
    f.resize(ndofs,1);
    u.resize(ndofs,1);
    // assemble stiffness- and mass matrix, will alter it later due to boundary conditions
    std::shared_ptr<Element> current_element;

    
    for (unsigned int i = 0; i < elements.size(); i++)
    { 
        current_element = elements.at(i);
        std::vector<unsigned int> dofs = current_element->get_element_dof_ids();
        for(unsigned int j=0;j < dofs.size();j++){
            unsigned int dof_row = dofs.at(j);
            for(unsigned int k=0;k < dofs.size();k++){
                unsigned int dof_column = dofs.at(k);
                Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Ke = current_element->get_Ke();
                Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Me = current_element->get_Me();
                K.coeffRef(dof_row,dof_column) += Ke(j,k);
                M.coeffRef(dof_row,dof_column) += Me(j,k);
            }
        }
    }
    // assemble load vector, will alter due to boundary conditions
    for (unsigned int i = 0; i < f_to_be_added.size(); i++)
    {
        f.coeffRef(f_to_be_added.at(i).first) += f_to_be_added.at(i).second;
    }
    // -----
    // START APPLYING BOUNDARY CONDITIONS!
    // index u=unknown, index d=dircichlet (known!)
    // [K_dd  K_du][u_d]=[f_d] 
    // [K_ud  K_uu][u_u] [f_u]
    // boundary conditions:
    // want to find free indices!
    std::vector<unsigned int> free_dofs;
    std::vector<unsigned int> dirichlet_dofs;
    for (unsigned int i = 0; i < bc.size(); i++)
    {
        dirichlet_dofs.push_back(bc.at(i).first);
    }
    
    for (unsigned int i = 0; i < ndofs; i++)
    {
        // if i not in dirichlet_dofs
        if(std::find(dirichlet_dofs.begin(), dirichlet_dofs.end(), i) != dirichlet_dofs.end()){
        }
        else{
            free_dofs.push_back(i);
        }
    }
    for (unsigned int i = 0; i < bc.size(); i++)
    {
        // if current bc==0, make all values in the corresponding row and column in K to zero
        unsigned int current_global_dof = bc.at(i).first;
        for (unsigned int j = 0; j < ndofs; j++)
        {
            if (j == current_global_dof)
            {
                f.coeffRef(current_global_dof) = bc.at(j).second;
            }
            else
            {
                f.coeffRef(current_global_dof) = 1e8*bc.at(i).second;
            }
        }        
        K.row(current_global_dof) *= 0;
        K.col(current_global_dof) *= 0;
        K.coeffRef(current_global_dof,current_global_dof) = 1e8;
    }
    duration_assemble = ( std::clock() - clock_assemble ) / (float) CLOCKS_PER_SEC;;
    std::cout << "K:" << std::endl;
    std::cout << Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>(K) << std::endl;
    std::cout << "M:" << std::endl;
    std::cout << Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>(M) << std::endl;
    std::cout << "f:" << std::endl;
    std::cout << Eigen::Matrix<float,1,Eigen::Dynamic>(f) << std::endl;
    std::cout << "---    Assembly completed in "<< duration_assemble << "seconds (wallclock time)    ---" << std::endl;
}


void Mesh::add_mid(std::unordered_map<std::string, std::string> options){
    // Create new MID based on input on *MATERIAL data line.
    // Other functions add mtrl data such as density etc
    // MID has to have a name
    std::string mid_name = options["NAME"];
    std::shared_ptr<Mid> mid  = std::shared_ptr<Mid>(new Mid(mid_name));
    mid_name_2_mid_pointer[mid_name] = mid;
    mids.push_back(mid);  
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
        std::cout << "*BOUNDARY: id=" << global_node_id << ", dofs={" << dof_from << "-" << dof_to << "}=";
        for (unsigned int i = 0; i < 3; i++)
        {
            std::cout << nodes.at(global_2_local_node_id[global_node_id])->dofs.at(i).id << ",";
        }
        std::cout <<" magnitude=" << magnitude << std::endl;
    }
    
}

void Mesh::add_load(std::string line){
    // for *cload line is: node, dof, magnitude
    std::vector<std::string> data   = misc::split_on(line,',');   
    unsigned int global_node_id     = std::stoi(data.at(0));
    unsigned int local_dof          = std::stoi(data.at(1));
    float magnitude                 = std::stof(data.at(2));
    // std::cout << global_node_id << std::endl;    
    // std::cout << local_dof  << std::endl;
    // std::cout << magnitude << std::endl;
    try
    {
        std::shared_ptr<Node> node = nodes.at(global_2_local_node_id[global_node_id]);
        unsigned int global_dof = node->dofs.at(local_dof-1).id;
        // we don't know complete number of dofs yet,
        // so we can't add directly to global load vector, but have
        // to store it here meanwhile
        f_to_be_added.push_back(std::make_pair(global_dof,magnitude));   
        std::cout << "*CLOAD: id=" << global_node_id << ", dof={" << local_dof << "}="<< global_dof <<", magnitude=" << magnitude << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cout << "ERROR: can't create *CLOAD because node with id " << global_node_id << " not found, exiting." << std::endl;
        exit(1);
    }
    
}

void Mesh::read_file(std::string filename){
    std::vector<std::string> filename_split = misc::split_on(filename,'/');   
    std::cout << "---    Starting to read input file " << filename_split.back() << "    ---" << std::endl;
    std::clock_t clock_read_file;
    float duration_clock_read_file;
    clock_read_file = std::clock();  
    misc::append_newline_to_textfile(filename);
    this->filename = filename;
    // Create input stream object
    std::ifstream input_file(filename);
    std::string line;
    unsigned int row_counter = 0;
    while (getline(input_file, line))
    {
        row_counter++;
        if (misc::is_keyword(line) == true){
            std::string keyword = misc::split_on(line, ',').at(0);
            // extract a map of parameters and their values
            std::unordered_map<std::string,std::string> options = misc::options_map(line);

            // bool temp = keyword == "*SOLID SECTION";
            // std::cout << keyword << "=" << "*SOLID SECTION" << temp << std::endl;
            if (keyword == "*SOLID SECTION" or keyword == "*SHELL SECTION")
            {
                add_pid(options);
            }     
            else if (keyword == "*MATERIAL")
            {
                // a complete material is typically described over several
                // lines so must make a complicated loop here..
                // Create mid:
                add_mid(options);
                // Add all following options to the last created material card.
                bool material_keywords_loop = true;
                while(material_keywords_loop == true)
                {
                    // Jump to next line 
                    getline(input_file, line);
                    row_counter++;
                    // Check if new line is a keyword, otherwise it's a comment and we skip that line 
                    if (misc::is_keyword(line) == true)
                    {
                        std::string keyword = misc::split_on(line, ',').at(0);
                        if (keyword == "*DENSITY")
                        {
                            // Skip to next line and save the value
                            getline(input_file, line);
                            float density = std::stof(misc::split_on(line, ',').at(0));
                            mids.back()->set_density(density);
                            std::cout << keyword << " = " << density << std::endl;
                        }
                        else if (keyword == "*ELASTIC")
                        {
                            // elastic keyword has to specify an option, such as isotropic, anisotropic..
                            std::unordered_map<std::string,std::string> options = misc::options_map(line);
                            if (options["TYPE"] == "ISOTROPIC")
                            {
                                // Example:
                                // *ELASTIC, TYPE=ISOTROPIC
                                // E,v
                                // Skip to next line and save the value
                                getline(input_file, line);
                                std::vector<std::string> values = misc::split_on(line, ',');
                                float E = std::stof(values.at(0));
                                float v = std::stof(values.at(1));
                                mids.back()->set_E(E);
                                mids.back()->set_v(v);
                                std::cout << keyword << ", " << options["TYPE"] << ", E" << " = " << E << ", v" << " = " << v << std::endl;
                            }
                            
                        }
                        else
                        {
                            // std::cout << keyword << " unsupported material option" << std::endl;
                            // we found keyword either not supported by or it's something
                            // entirely different. set input_file to current line and start reading normal
                            // keywords again
                            // break;

                        }
                    }
                    else
                    {
                        // peep next line, if it's not one of the 
                        // supported material specific keywords we skip adding data to the material
                        unsigned int previous_pos = input_file.tellg();
                        getline(input_file, line);
                        if (misc::is_keyword(line) == true)
                        {
                            std::string keyword = misc::split_on(line, ',').at(0);
                            // if the new keyword is not a supported material keyword we break loop and keep reading file
                            // otherwise keep on truckin'
                            if (misc::is_string_in_string_vector(keyword,mids.back()->supported_material_keywords) == false)
                            {
                            input_file.seekg(previous_pos);
                            material_keywords_loop = false;
                            }                            
                        }
                    }
                }
                
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
            else if (keyword == "*STATIC"){
                // with this flag set a static analis will be ran
                this->static_analysis = true;
            }
            else if (keyword == "*FREQUENCY"){
                // with this flag set a static analis will be ran
                this->eigenvalue_analysis = true;
                // peep next line and save the number of modes to find
                // Jump to next line 
                getline(input_file, line);
                row_counter++;
                std::string number_of_modes_to_find_string = misc::split_on(line, ',').at(0);
                this->number_of_modes_to_find = std::stoi(number_of_modes_to_find_string);
                std::cout << keyword << ": number of eigenvalues to be calculated = " << number_of_modes_to_find << std::endl;
            }
            else if (keyword == "*NODE" or keyword == "*ELEMENT")
            {
                getline(input_file, line);
                row_counter++;
                bool inner_loop_keyword = true;
                while (inner_loop_keyword == true){
                row_counter++;   
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
                    // std::cout << line << ", is  keyword?"<< misc::is_keyword(line) << ", row_counter = " << row_counter <<  "\n";
                    if (misc::is_keyword(line) == true or line.empty() == true)
                    {
                        input_file.seekg(previous_pos);
                        inner_loop_keyword = false;
                    }
                }                
            }
            else
            {
                std::cout << keyword << ": not supported" << std::endl;
            }
            
            if (input_file.eof() == true)
            {
                break;
            }
            row_counter++;   
        }
    }       
    duration_clock_read_file = (std::clock() - clock_read_file ) / (float) CLOCKS_PER_SEC;
    std::cout << "---    Input file read and written information to the log file, in " << duration_clock_read_file << " seconds (wallclock time)    ---" << std::endl;
};

void Mesh::add_pid(std::unordered_map<std::string, std::string> options){    
    std::string pid_name = options["ELSET"];
    std::string mid_name = options["MATERIAL"];
    try
    {
        // Find material
        std::shared_ptr<Mid> mid = mid_name_2_mid_pointer[mid_name];
        // TODO: find MID pointer instead!!
        // Create new pid
        std::shared_ptr<Pid> pid  = std::shared_ptr<Pid>(new Pid(pid_name,mid));
        pids.push_back(pid);
        pid_map[pid_name] = pid;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        std::cout << "ERROR: material with name " << mid_name << "not found, exiting." << std::endl;
        exit(1);
    }

    

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
    else if (type == "CPS3")
    {
        element = std::shared_ptr<Element>(new CPS3(element_id,element_connectivity,pid));
    }
    else if (type == "CPS4")
    {   
        element = std::shared_ptr<Element>(new CPS4(element_id,element_connectivity,pid));
    }
    else if (type == "C3D10")
    {   
        element = std::shared_ptr<Element>(new C3D10(element_id,element_connectivity,pid));
    }
    else if (type == "C3D8")
{   
    element = std::shared_ptr<Element>(new C3D8(element_id,element_connectivity,pid));
}
    else
    {
        std::cout << "*ELEMENT: type=" << type << " not supported" << std::endl;
        return;
    }    
    elements.push_back(element);
    
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
Mesh::~Mesh(){};


