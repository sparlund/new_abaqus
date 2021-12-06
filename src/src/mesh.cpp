#include "../include/mesh.h"
#include "../include/dof.h"
#include "../include/S3.h"
#include "../include/CPS3.h"
#include "../include/CPS4.h"
#include "../include/C3D10.h"
#include "../include/C3D8.h"
#include "../include/C3D20.h"
#include "../include/misc_string_functions.h"
#include "../include/set.h"
#include "../../external_libs/Eigen/Dense"
#include "../../external_libs/Eigen/Sparse"
#include "../../external_libs/Eigen/Cholesky"
#include "../../external_libs/Eigen/Eigenvalues"
// for eigenfrequency solver
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <string>
#include <vector>
#include <array>
#include <complex>
#include <unordered_map>
#include <fstream>
#include <functional>
#include <iostream>
#include <iomanip>
#include <memory>
#include <numeric>
#include <Spectra/SymGEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/Util/GEigsMode.h>
#include <Spectra/Util/SelectionRule.h>
#include <utility>
#include <cmath>
#include <fstream>
#include <ctime>
#include <stdlib.h>

// TODO: move this to misc::
float squirt_and_divide_by_2pi(float in){
    return std::sqrt(in)/(2*3.14159265359);
}

unsigned int Mesh::get_number_of_dofs() const{
    const auto& temp_ndofs = nodes.at(0)->dofs.at(0)->get_global_dof_id_counter();
    return temp_ndofs;
}

void Mesh::print_matrix_to_mtx(const Eigen::SparseMatrix<float>& A,const std::string& output_filename) const
{
    // Print input matrix A in abaqus "COORDINATE" format in a .mtx-file
    std::ofstream mtx;
    mtx.open(output_filename);
    for (int k=0; k<A.outerSize(); ++k){
        for (Eigen::SparseMatrix<float>::InnerIterator it(A,k); it; ++it){
            int row = it.row()+1;   // row index
            int column = it.col()+1;   // col index (here it is equal to k)
            // it.index(); // inner index, here it is equal to it.row()
            // this should only print non-zero but doesn't work correctly? print zeros to, need to add a check
            float v = it.value();
            if (v != 0.0f)
            {
            mtx << row << " " << column << " " << v << std::endl;
            }
            
        }
    }
    mtx.close();
}

void Mesh::export_2_vtk(){
    // output results to vtk format for the ParaView post processor
    // https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    // http://profs.sci.univr.it/~caliari/aa1314/advanced_numerical_analysis/Pellegrini.pdf
    std::clock_t clock_export;
    float duration_clock_export;
    clock_export = std::clock();
    std::cout << "---    Starting to export VTK results file    ---" << std::endl;
    std::string output_filename = analysis_name + ".vtk";
    std::ofstream output(output_filename);
    output << "# vtk DataFile Version 3.0" << std::endl;
    output << "vtk output" << std::endl;
    // ascii or binary format:
    output << "ASCII" << std::endl;
    output << "DATASET UNSTRUCTURED_GRID" << std::endl;
    // no need to use double precision on printing results
    output << "POINTS " << nodes.size() << " float" << std::endl;
    // print nodes to result file
    for(const auto& node: nodes)
    {
        output << node->x << " " << node->y << " " << node->z << std::endl;
    }
    unsigned int total_amount_of_used_nodes=0;
    for(const auto& element: elements)
    {
        total_amount_of_used_nodes += element->get_element_nnodes();
    }    
    output << "CELLS " << elements.size() << " " << elements.size() + total_amount_of_used_nodes  << std::endl;
    for(const auto& element: elements)
    {
        output << element->get_element_nnodes() << " ";
        for (unsigned int  j = 0; j < element->get_connectivity().size(); j++)
        {
            output << global_2_local_node_id[element->get_connectivity().at(j)->id] << " ";
        }
        output << std::endl;
    }
    
    output << "CELL_TYPES " << elements.size() << std::endl;
    for(const auto& element: elements)
    {
        output << element->get_vtk_identifier() << std::endl;
    }
    // point_data, i.e displacement on nodes
    output << "POINT_DATA " << nodes.size() << std::endl;
    if (static_analysis)
    {
        output << "FIELD displacement 1" << std::endl;
        // data name, number of values (3 for 3D & 2D, Z is zero for 2D, datatype)
        output << "displacement 3 " << nodes.size() << " float" << std::endl;
        for(const auto& node: nodes)
        {
            // 2 scenarios available so far: 3 dofs per node & 2 dofs per node
            if (node->dofs.size() == 3)
            {
                output << u(node->dofs.at(0)->id) << " " << u(node->dofs.at(1)->id) << " " << u(node->dofs.at(2)->id) << std::endl;
            }
            else{
                output << u(node->dofs.at(0)->id) << " " << u(node->dofs.at(1)->id) << " 0" << std::endl;
            }
            
        }
    }
    if (eigenvalue_analysis)
    {
        // scale eigenvectors before saving to vtk file
        // normalize against largest value in eigenvector matrix
        eigenvectors /= eigenvectors.maxCoeff();
        for (unsigned int mode = number_of_modes_to_find-1; mode > 1; mode--)
        {
            output << "FIELD mode" << (number_of_modes_to_find - mode) << " 1" << std::endl;
            // data name, number of values (3 for 3D & 2D, Z is zero for 2D, datatype)
            output << "mode" << (number_of_modes_to_find - mode) << " 3 " << nodes.size() << " float" << std::endl;
            // for (unsigned int i = 0; i < nodes.size(); i++)
            for(const auto& node: nodes)
            {
                if (node->dofs.size() == 3)
                {
                    output << eigenvectors(node->dofs.at(0)->id,mode)
                            << " " 
                            << " " 
                            << " " 
                            << " " 
                            << " " 
                            << " " 
                            << " " 
                            << eigenvectors(node->dofs.at(1)->id,mode)
                            << " " 
                            << " " 
                            << " " 
                            << eigenvectors(node->dofs.at(2)->id,mode) << std::endl;
                }
                else
                {
                    output << eigenvectors(node->dofs.at(0)->id,mode)
                        << " " 
                        << " " 
                        << " " 
                        << " " 
                        << " " 
                        << " " 
                        << " " 
                        << eigenvectors(node->dofs.at(1)->id,mode)
                        << " 0" << std::endl;
                }   
            }
        }          
    }    
    output.close();
    duration_clock_export = ( std::clock() - clock_export ) / static_cast<float>(CLOCKS_PER_SEC);
    std::cout << "---    Exported to VTK format in " << duration_clock_export << " seconds (wallclock time)   ---" << std::endl;
}

void Mesh::matrix_generate()
{
    std::cout << "---    Starting to print matrices to mtx format    ---\n" << std::endl;  
    if (mtx_to_print.STIFFNESS_WITH_BC)
    {
        std::string mtx_output_filename = analysis_name + "_STIF2_new_abaqus.mtx";
        print_matrix_to_mtx(K_with_bc,mtx_output_filename);
        std::cout << "Printed stiffness matrix adjusted for boundary conditions to " << mtx_output_filename << std::endl;
    }
    if(mtx_to_print.MASS)
    {
        std::string mtx_output_filename = analysis_name + "_MASS2_new_abaqus.mtx";
        print_matrix_to_mtx(M,mtx_output_filename);
        std::cout << "Printed mass matrix adjusted to " << mtx_output_filename << std::endl;
    }
    if (mtx_to_print.LOAD)
    {
        if(static_analysis)
        {
            std::string mtx_output_filename = analysis_name + "_LOAD2_new_abaqus.mtx";
            print_matrix_to_mtx(f,mtx_output_filename);
            std::cout << "Printed load vector to " << mtx_output_filename << std::endl;
        }
        else
        {
            std::cout << "ERROR: f not used, add *STATIC analysis." << std::endl;
        }
        
    }
    std::cout << "---    Finished printing matrices to mtx format    ---\n" << std::endl;  
}

void Mesh::solve(){
    if (eigenvalue_analysis)
    {
        solve_eigenfrequency();
    }
    if (static_analysis)
    {
        solve_static();
        solve_static_with_contact();
    }
    if(steady_state_dynamics_analysis)
    {
        solve_steady_state_dynamics();
    }
    export_2_vtk();
    matrix_generate();
}

void Mesh::solve_steady_state_dynamics()
{
    std::cout << "---    Starting to solve frequency response function    ---\n" << std::endl;  
    std::clock_t clock_solve;
    auto duration_clock_solve = std::clock();
    // K_eff = -w^2 * M + K;
    auto frequency_step_width = (steady_state_dynamics_upper_limit - steady_state_dynamics_lower_limit) / steady_state_dynamics_number_of_points;
    float w, frequency;
    auto ndofs = get_number_of_dofs();
    // Eigen::SparseVector<float> F{ndofs};
    Eigen::SparseMatrix<std::complex<float>> H{ndofs,ndofs};
    const std::complex<float> i{0.0,1.0};
    Eigen::SparseLU<Eigen::SparseMatrix<float>> solver;
    std::cout << "A" << std::endl;
    for(unsigned int k = 0; k < 1; k++)
    {
        frequency = steady_state_dynamics_lower_limit + k*frequency_step_width;
        w = frequency*2*M_PIl;
        // want to solve (-w^2*M +K)*u = F'
        //               <--------->
        //                    H
        auto F = f*std::exp(i*w);
        std::cout << "E" << std::endl;
        auto t = std::pow(w,2)*M;
        std::cout << "F" << std::endl;
        Eigen::SparseMatrix<float> H = M + K;
        solver.compute(H);
        if(solver.info()!=Eigen::Success)
        {
            std::cout << "decomposition failed!" << std::endl;
            exit(0);
        }
        std::cout << "G" << std::endl;
        auto Ft = F.conjugate().transpose();
        std::cout << "Ft=" << Ft << std::endl;
        std::cout << "I" << std::endl;
        // u = solver.solve(Ft);
        if(solver.info()!=Eigen::Success)
        {
            std::cout << "solver failed!" << std::endl;
            exit(0);
        }
    }
    std::cout << "u=" << u << std::endl;
    duration_clock_solve = ( std::clock() - clock_solve ) / static_cast<float>(CLOCKS_PER_SEC);
    std::cout << "---    Solution to eigenvalue problem found in " << std::setprecision(2) << duration_clock_solve << " seconds (wallclock time)    ---" << std::endl;
}


void Mesh::solve_eigenfrequency(){
    std::cout << "---    Starting to solve eigenvalue problem    ---\n" << std::endl;  
    // Print global stiffness and mass matrices to .mtx format
    
    
    std::clock_t clock_solve;
    clock_solve = std::clock();  
    // Need to modify global stiffness matrix
    // in order to account for boundary conditions
    Eigen::SparseMatrix<float> K_eigen = K;
    for(const auto& i : bc)
    {
        // if current bc==0, make all values in the corresponding row and column in K to zero
        unsigned int current_global_dof = i.first;
        K_eigen.coeffRef(current_global_dof,current_global_dof) = penalty_value;
    }
    Spectra::SymShiftInvert<float, Eigen::Sparse, Eigen::Sparse> opK1(K_eigen,M);
    Spectra::SparseSymMatProd<float> opM1(M);
    unsigned int ncv = (2*number_of_modes_to_find);
    // abaqus defines number of eigenvalues to solve for as 
    // the ones with smallest magntide --> sigma=0
    double sigma = 0;
    Spectra::SymGEigsShiftSolver<float,
                        Spectra::SymShiftInvert<float, Eigen::Sparse, Eigen::Sparse>,
                        Spectra::SparseSymMatProd<float>,
                        Spectra::GEigsMode::ShiftInvert> es(opK1,
                                                            opM1,
                                                            number_of_modes_to_find,
                                                            ncv,
                                                            sigma);
    es.init();
    es.compute(Spectra::SortRule::LargestMagn, 1e4,1e-18);
    // Retrieve results
    this->eigenvalues = es.eigenvalues();
    this->eigenvectors = es.eigenvectors();
    // eigenfrequency = sqrt(eigenvalue)/(2*pi);
    this->eigenfrequencies = eigenvalues.unaryExpr(&squirt_and_divide_by_2pi);

    // Eigen::MatrixXd Mdense = M;
    // Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(A, Bdense);
    auto duration_clock_solve = ( std::clock() - clock_solve ) / static_cast<float>(CLOCKS_PER_SEC);
    std::cout << "                              E I G E N V A L U E    O U T P U T\n" << std::endl;
    std::cout << " MODE NO      EIGENVALUE              FREQUENCY         " << std::endl;//GENERALIZED MASS   COMPOSITE MODAL DAMPING" << std::endl
    std::cout << "                             (RAD/TIME)   (CYCLES/TIME)\n\n" << std::endl;
    int temp_mode_counter_for_text = 1;
    for (unsigned int i = eigenfrequencies.rows()-1; i > 0; --i)
    // for(unsigned int i = 0;i < eigenfrequencies.rows();i++)
    {
        std::cout << std::setw(8) << temp_mode_counter_for_text << "      " << std::setw(10) << std::scientific << eigenvalues[i] << "    " <<  std::setw(19) << std::setprecision(2)  << std::fixed << eigenfrequencies[i] << std::endl;
        temp_mode_counter_for_text++;
    }
    std::cout << "---    Solution to eigenvalue problem found in " << std::setprecision(2) << duration_clock_solve << " seconds (wallclock time)    ---" << std::endl;
}

void Mesh::solve_static_with_contact()
{
    std::cout << "---    Starting to solve non-linear problem with contact mechanics    ---" << std::endl;
    std::clock_t clock_solve;
    float duration_clock_solve;
    clock_solve = std::clock();
    /*
    find which nodes already on contact before load or
    controlled displacement is applied
    for step in number_of_steps:
        Equilibrium iterations below
        for iteration in max_iterations:
            assemble K
            compute internal force vector and assemble in f_int
            if number of nodes in contact != 0:
                f_int = f_int + penatly*displacement at contact nodes
                K = K + penalty*displacement at contact nodes
                        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                    "spring forces"
            Q = nodal forces from prescribed displacement
            R = Q - f_int = residual of nodal forces (free nodes)
            u = solve(K,bc)
            save displacement for this step
            find what nodes are in contact and penetrating too far
            save which nodes are in contact this step
            
            error = norm(R)
            if error < max(abs(R))*relative_tolerance
                convergence reached, exit
            else
                continue
            if iteration == max_iterations:
                max number of iterations reached, did not converge


    
    */
    duration_clock_solve = ( std::clock() - clock_solve ) / (float) CLOCKS_PER_SEC;
    std::cout << "---    Solution to non-linear problem with contact mechanics found in " << duration_clock_solve << " seconds (wallclock time)    ---" << std::endl;
}

void Mesh::solve_static(){
    std::cout << "---    Starting to solve linear problem Ku=f    ---" << std::endl;
    std::clock_t clock_solve;
    float duration_clock_solve;
    clock_solve = std::clock();  
    // Need to alter global stiffness matric and global load vector to account for boundary conditions
    K_with_bc = K;
    for (unsigned int i = 0; i < bc.size(); i++)
    {
        try
        {
            unsigned int current_global_dof = bc.at(i).first;
            for (unsigned int j = 0; j < get_number_of_dofs(); j++)
            {
                if (j == current_global_dof)
                {
                    f.coeffRef(current_global_dof) = bc.at(i).second;
                }
                else
                {
                    f.coeffRef(current_global_dof) = penalty_value*bc.at(i).second;
                }
            }
            // if current bc==0, make all values in the corresponding row and column in K to zero
            K_with_bc.row(current_global_dof) *= 0;
            K_with_bc.col(current_global_dof) *= 0;
            K_with_bc.coeffRef(current_global_dof,current_global_dof) = penalty_value;
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << '\n';
        std::cout << "i=" << i << '\n';
    }
    }
    // std::cout << "K after bc:" << std::endl;
    // std::cout << Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>(K_static) << std::endl;
    // std::cout << "f after bc:" << std::endl;
    // std::cout << Eigen::Matrix<float,1,Eigen::Dynamic>(f) << std::endl;
    // Ku=f
    // Eigen::SparseLU"Eigen::SparseMatrix<float>> solver;
    // solver.compute(K_static);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver(K_with_bc);
    u = solver.solve(f);
    // std::cout << "d" << std::endl;
    // std::cout << "u:" << std::endl;
    // std::cout << Eigen::Matrix<float,1,Eigen::Dynamic>(u).transpose() << std::endl;
    duration_clock_solve = ( std::clock() - clock_solve ) / (float) CLOCKS_PER_SEC;
    std::cout << "---    Solution to linear problem found in " << duration_clock_solve << " seconds (wallclock time)    ---" << std::endl;
}


void Mesh::assemble(){
    // By the time of assemble we know the number of dofs --> pre-allocate K, f & solution u
    std::cout << "---    Starting matrix assembly    ---" << std::endl;
    std::cout << "           Mesh size:" << std::endl;
    std::cout << "               nodes = " << nodes.size() << std::endl;
    std::cout << "            elements = " << elements.size() << std::endl;
    // access arbitrary dof object from arbitrary node object and check the static member to see total ndofs     
    std::cout << "  degrees of freedom = " << get_number_of_dofs() << std::endl << std::endl;
    // Also print weight as a sanity check!
    // float model_total_weight = 0;
    // for (unsigned int i = 0; i < elements.size(); i++)
    // {
    //     model_total_weight = model_total_weight + elements.at(i)->get_weight();
    // }
    // std::cout << "         Model info:" << std::endl;
    // std::cout << "              weight = " << model_total_weight << std::endl;
    
    std::clock_t clock_assemble;
    float duration_assemble;
    clock_assemble = std::clock();
    K.resize(get_number_of_dofs(),get_number_of_dofs());
    M.resize(get_number_of_dofs(),get_number_of_dofs());
    f.resize(get_number_of_dofs(),1);
    u.resize(get_number_of_dofs(),1);
    // assemble stiffness- and mass matrix, will alter it later due to boundary conditions
    std::cout << "Calculating all element stiffness and or mass matrices:" <<std::endl;
    size_t current_element_counter = 1;
    size_t nel = get_number_of_elements();
    for(const auto& element: elements)
    {
        std::cout << "[" << current_element_counter << "/" << nel << "]" << "\r" << std::flush;
        element->calculate_Ke();
        element->calculate_Me();
        current_element_counter++;
    }
    std::cout << std::endl;
    std::cout << "Progress assembly:" <<std::endl;
    current_element_counter = 1;
    std::vector<Eigen::Triplet<float>> K_tripletList;
    std::vector<Eigen::Triplet<float>> M_tripletList;
    K_tripletList.reserve(get_number_of_dofs());
    M_tripletList.reserve(get_number_of_dofs());
    for(const auto& element:elements)
    {
        std::cout << "[" << current_element_counter << "/" << nel << "]\r";
        std::cout.flush();
        auto Ke         = element->get_Ke();
        auto Me         = element->get_Me();
        auto dofs       = element->get_element_dof_ids();
        for(unsigned int j=0;j < dofs.size();j++){
            unsigned int dof_row = dofs.at(j);
            for(unsigned int k=0;k < dofs.size();k++){
                auto dof_column = dofs.at(k);
                K_tripletList.push_back(Eigen::Triplet<float>(dof_row,dof_column,Ke(j,k)));
                M_tripletList.push_back(Eigen::Triplet<float>(dof_row,dof_column,Me(j,k)));
            }
        }
        current_element_counter++;
    }
    K.setFromTriplets(K_tripletList.begin(), K_tripletList.end());
    M.setFromTriplets(M_tripletList.begin(), M_tripletList.end());
    std::cout << std::endl;
    // assemble load vector, will alter due to boundary conditions
    std::cout << "f_to_be_added.size() = " << f_to_be_added.size() << std::endl;
    for (unsigned int i = 0; i < f_to_be_added.size(); i++)
    {
        std::cout << "f_to_be_added.at(i).first = " << f_to_be_added.at(i).first << std::endl;
        std::cout << "f_to_be_added.at(i).second = " << f_to_be_added.at(i).second << std::endl;
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
    for(const auto& i : bc)
    {
        dirichlet_dofs.push_back(i.first);
    }
    for (unsigned int i = 0; i < get_number_of_dofs(); i++)
    {
        // if i not in dirichlet_dofs
        if(std::find(dirichlet_dofs.begin(), dirichlet_dofs.end(), i) != dirichlet_dofs.end()){
        }
        else{
            free_dofs.push_back(i);
        }
    }
    duration_assemble = ( std::clock() - clock_assemble ) / (float) CLOCKS_PER_SEC;;
    // std::cout << "K:" << std::endl;
    // std::cout << Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>(K) << std::endl;
    // std::cout << "M:" << std::endl;
    // std::cout << Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>(M) << std::endl;
    // std::cout << "f:" << std::endl;
    // std::cout << Eigen::Matrix<float,1,Eigen::Dynamic>(f) << std::endl;
    std::cout << "---    Assembly completed in " << std::setprecision(2) << duration_assemble << " seconds (wallclock time)    ---" << std::endl;
}


void Mesh::add_mid(std::unordered_map<std::string, std::string> options){
    // Create new MID based on input on *MATERIAL data line.
    // Other functions add mtrl data such as density etc
    // MID has to have a name
    std::string mid_name = options["NAME"];
    if(misc::is_valid_name(mid_name)){
        auto mid = std::make_unique<Mid>(mid_name);
        mid_map[mid_name] = mid.get();
        mids.push_back(std::move(mid));
    }
    else{
        std::cout << "ERROR: terminate at line " << row_counter << " in file " << current_inputfile << std::endl;
        std::cout << "ERROR: invalid entity name" << std::endl;
        std::terminate();
    }
};

void Mesh::add_boundary(std::string line,std::unordered_map<std::string, std::string> options){
    if(options.find("TYPE") == options.end()){
        std::cout << "ERROR: terminate at line " << row_counter << " in file " << current_inputfile << std::endl;
        std::cout << "ERROR: missing parameter TYPE" << std::endl;
        std::terminate();
    }
    const auto& type = options["TYPE"];
    const auto& data = misc::split_on(line,',');
    if (type=="DISPLACEMENT")
    {
        // Check if first object is a string or int <--> set* or node_id
        // Need to iterate through string. If there are no alpabetical chars it's
        // a node id, otherwise it's a string!
        auto node_id_or_nset_name = data.at(0);
        misc::trim_leading_and_ending_whitespace(node_id_or_nset_name);
        bool node_boundary     = std::find_if(node_id_or_nset_name.begin(),
                                              node_id_or_nset_name.end(),
                                              [](unsigned char c) { return !std::isdigit(c); }) == node_id_or_nset_name.end();
        unsigned int dof_from  = std::stoi(data.at(1));
        unsigned int dof_to    = std::stoi(data.at(2));
        auto         magnitude = std::stof(data.at(3));
        // node,dof_from,dof_to,magnitude
        if (node_boundary)
        {
            unsigned int global_node_id     = std::stoi(data.at(0));
            if(node_map.find(global_node_id) == node_map.end()){
                std::cout << "ERROR: terminate at line " << row_counter << " in file " << current_inputfile << std::endl;
                std::cout << "ERROR: can't find node with id " << global_node_id << std::endl;
                std::terminate();
            }
            const auto& node = node_map[global_node_id];
            std::cout << "*BOUNDARY: id=" << global_node_id << ", dofs=" << dof_from << "-" << dof_to << ", magnitude=" << magnitude << std::endl;
            for(unsigned int j = dof_from; j <= dof_to; j++){
                bc.push_back(std::make_pair(node->dofs.at(j-1)->id,magnitude));
            }
            return;
        }
        else
        {
            const auto& node_set = node_set_from_node_set_name[node_id_or_nset_name];
            auto number_of_nodes = node_set->get_number_of_entities();
            for (unsigned int i = 0; i < number_of_nodes; i++)
            {
                const auto& node           = node_set->get_entity(i);
                const auto& global_node_id = node->id;
                std::cout << "*BOUNDARY: id=" << global_node_id << ", dofs={" << dof_from << "-" << dof_to << "}=" << ", magnitude=" << magnitude << std::endl;
                for (unsigned int j = dof_from; j <= dof_to; j++)
                {
                    // abaqus starts counting dof's at 1, but vectors start at 0
                    const auto& node = node_set->get_entity(i);
                    bc.push_back(std::make_pair(node->dofs.at(j-1)->id,magnitude));
                }
            }
        }
    }
}

void Mesh::add_load(std::string line, std::unordered_map<std::string, std::string> options){
    // for *cload line is: node, dof, magnitude
    auto         data               = misc::split_on(line,',');   
    unsigned int global_node_id     = std::stoi(data.at(0));
    unsigned int local_dof          = std::stoi(data.at(1));
    float        magnitude          = std::stof(data.at(2));
    try
    {
        auto         node       = node_map[global_node_id]; 
        unsigned int global_dof = node->dofs.at(local_dof-1)->id;
        // we don't know complete number of dofs yet,
        // so we can't add directly to global load vector, but have
        // to store it here meanwhile
        f_to_be_added.emplace_back(std::make_pair(global_dof,magnitude));   
        std::cout << "*CLOAD: id=" << global_node_id << ", dof={" << local_dof << "}="<< global_dof <<", magnitude=" << magnitude << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cout << "ERROR: can't create *CLOAD because node with id " << global_node_id << " not found, exiting." << std::endl;
        exit(1);
    }
    
}



void Mesh::read_file_new_method(const std::string& filename){
    std::vector<std::string> filename_split = misc::split_on(filename,'/');   
    std::cout << "---    Starting to read input file " << filename_split.back() << "    ---" << std::endl;
    std::clock_t clock_read_file;
    float duration_clock_read_file;
    clock_read_file = std::clock();  
    misc::append_newline_to_textfile(filename);
    for (auto && keyword : keywords)
    {
        read_file(filename,keyword);
    }
    duration_clock_read_file = (std::clock() - clock_read_file ) / static_cast<float>(CLOCKS_PER_SEC);
    std::cout << "---    Input file "<< filename << " read and written information to the log file, in " << std::setprecision(2) << duration_clock_read_file << " seconds (wallclock time)    ---" << std::endl;

}

void Mesh::read_file(const std::string& filename, const std::string& keyword){
    std::ifstream input_file(filename);
    if (input_file.fail())
    {
        // file could not be opened
        std::cout << "Error: include " << filename << " could not be opened." << std::endl;
        exit(0);
    }
    row_counter = 0;
    current_inputfile = filename;
    std::string line;
    while (getline(input_file, line))
    {
        row_counter++;
        if (misc::is_keyword(line)){
            auto current_line_keyword = misc::split_on(line, ',').at(0);
            if (current_line_keyword == keyword or current_line_keyword == "*INCLUDE"){
                // extract a map of parameters and their values
                auto options = misc::options_map(line);
                if (current_line_keyword == "*INCLUDE")
                {
                    // Call this function recursively if found include.
                    std::string include_filename = options["INPUT"];
                    read_file(include_filename, keyword);
                }
                else if (keyword == "*SOLID SECTION")
                {
                    add_pid(options);
                }
                else if (keyword == "*CONTACT PAIR")
                {
                    // advance to next line, grab sets, then leave.
                    // it could be a comment, new keyword.. then we exit
                    getline(input_file, line);
                    row_counter++;
                    if (!misc::is_comment(line) && !misc::is_keyword(line))
                    {
                        // split on comma and find the two sets.
                        // slave is first
                        auto slave_and_master_set = misc::split_on(line,',');
                        auto slave_set_name = slave_and_master_set.at(0);
                        auto master_set_name = slave_and_master_set.at(1);
                        // find node set entities
                        auto slave_set  = node_set_from_node_set_name[slave_set_name];
                        auto master_set = node_set_from_node_set_name[master_set_name];
                        contact = std::make_unique<Contact>(*master_set,
                                                            *slave_set);
                    }
                    else
                    {
                        std::cout << "ERROR: *CONTACT PAIR not followed by 2 sets." << std::endl;
                    }
                }
                else if(keyword == "*NSET"){
                    add_set(line,options);
                    getline(input_file, line);
                    auto nset = nsets.back().get();
                    bool inner_loop_keyword = true;
                    while(inner_loop_keyword)
                    {
                        row_counter++;
                        if(!misc::is_comment(line) and !misc::is_keyword(line)){
                            // Jump to next line 
                            auto entity_ids = misc::split_on(line, ',');
                            // abaqus documentation says max number of entities per line is 16
                            for (const auto& entity_id : entity_ids)
                            {
                                auto global_node_id = std::stoi(entity_id);
                                auto& node_pointer  = node_map[global_node_id];
                                nset->add_entity(node_pointer);
                            }
                        }
                        // Want to peek next line, if it's a keyword or empty line we break
                        // the while loop and start over!
                        unsigned int previous_pos = input_file.tellg();
                        getline(input_file, line);
                        if (misc::is_keyword(line) or line.empty())
                        {
                            input_file.seekg(previous_pos);
                            inner_loop_keyword = false;
                        }
                    }
                }    
                else if (keyword == "*MATERIAL")
                {
                    // a complete material is typically described over several
                    // lines so must make a complicated loop here..
                    // Create mid:
                    add_mid(options);
                    // Add all following options to the last created material card.
                    bool material_keywords_loop = true;
                    while(material_keywords_loop)
                    {
                        // Jump to next line 
                        getline(input_file, line);
                        row_counter++;
                        // Check if new line is a keyword, otherwise it's a comment and we skip that line 
                        if (misc::is_keyword(line))
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
                                    auto values = misc::split_on(line, ',');
                                    float E = std::stof(values.at(0));
                                    float v = std::stof(values.at(1));
                                    mids.back()->set_E(E);
                                    mids.back()->set_v(v);
                                    std::cout << keyword << ", " << options["TYPE"] << ", E" << " = " << E << ", v" << " = " << v << std::endl;
                                }
                            }
                            else
                            {
                                break;
                            }
                        }
                        else
                        {
                            // peep next line, if it's not one of the 
                            // supported material specific keywords we skip adding data to the material
                            unsigned int previous_pos = input_file.tellg();
                            getline(input_file, line);
                            if (misc::is_keyword(line) or line.empty())
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
                    // TODO:
                    // We've added all the info we want to this specific material, let's set up the constitutive matrices for both 2D and 3D.
                    // The latest createt MID bill be att the back of the material list
                    // mids.back()                
                    // 2D:

                }
                else if (keyword=="*BOUNDARY"){
                    getline(input_file, line);
                    row_counter++;
                    bool inner_loop_keyword = true;
                    while (inner_loop_keyword){
                    // Ignore if it's a comment! Still on same keyword.
                        if (misc::is_comment(line) == false){
                            add_boundary(line,options);
                        }
                        // Want to peek next line, if it's a keyword or empty line we break
                        // the while loop and start over!
                        unsigned int previous_pos = input_file.tellg();
                        getline(input_file, line);
                        if (misc::is_keyword(line) or line.empty())
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
                    while (inner_loop_keyword){
                    // Ignore if it's a comment! Still on same keyword.
                        if (misc::is_comment(line) == false){
                            add_load(line,options);
                        }
                        // Want to peek next line, if it's a keyword or empty line we break
                        // the while loop and start over!
                        unsigned int previous_pos = input_file.tellg();
                        getline(input_file, line);
                        if (misc::is_keyword(line) or line.empty())
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
                else if(keyword == "*MATRIX GENERATE")
                {
                    auto what_matrix_to_generate = misc::split_on(line, ',').at(1);
                    misc::trim_leading_and_ending_whitespace(what_matrix_to_generate);
                    if(what_matrix_to_generate == "STIFFNESS"){
                        mtx_to_print.STIFFNESS_WITH_BC = true;
                    }
                    else if(what_matrix_to_generate == "MASS"){
                        mtx_to_print.MASS = true;
                    }
                    else if(what_matrix_to_generate == "LOAD"){
                        mtx_to_print.LOAD = true;
                    }
                    else
                    {
                        std::cout << "ERROR: terminate at line " << row_counter << " in file " << current_inputfile << std::endl;
                        std::cout << "ERROR: " << what_matrix_to_generate << " not recognized as *MATRIX GENERATE option." << std::endl;
                        std::terminate();
                    }
                }
                else if (keyword == "*STEADY STATE DYNAMICS")
                {
                    this->steady_state_dynamics_analysis = true;
                    getline(input_file, line);
                    row_counter++;
                    auto lowLim_and_uppLim_and_total_noPoints = misc::split_on(line, ',');
                    steady_state_dynamics_lower_limit      = std::stof(lowLim_and_uppLim_and_total_noPoints.at(0));
                    steady_state_dynamics_upper_limit      = std::stof(lowLim_and_uppLim_and_total_noPoints.at(1));
                    steady_state_dynamics_number_of_points = std::stoi(lowLim_and_uppLim_and_total_noPoints.at(2));
                }
                else if (keyword == "*ELEMENT"){
                    getline(input_file, line);
                    row_counter++;
                    bool inner_loop_keyword = true;
                    while (inner_loop_keyword){
                    row_counter++;   
                    // Ignore if it's a comment! Still on same keyword.
                    if (misc::is_comment(line) == false){
                        // If the line ends with a comma (','), the next line will continue to list nodes for that element.
                        // We will ever need 3 lines, so can hard-code for two lines.
                        // Check if last char is a comma:
                        if (line.back() == ',')
                        {
                            // Peek next row in the text file and append it to our data line
                            std::string next_line;
                            getline(input_file, next_line);
                            line += next_line;
                        }
                        add_element(line,options);
                    }
                    // Want to peek next line, if it's a keyword or empty line we break
                    // the while loop and start over!
                    unsigned int previous_pos = input_file.tellg();
                    getline(input_file, line);
                    // std::cout << line << ", is  keyword?"<< misc::is_keyword(line) << ", row_counter = " << row_counter <<  "\n";
                    if (misc::is_keyword(line) or line.empty())
                    {
                        input_file.seekg(previous_pos);
                        inner_loop_keyword = false;
                    }
                }                
                }
                else if (keyword == "*NODE")
                {
                    getline(input_file, line);
                    row_counter++;
                    bool inner_loop_keyword = true;
                    while (inner_loop_keyword){
                    row_counter++;   
                    // Ignore if it's a comment! Still on same keyword.
                        if (misc::is_comment(line) == false){
                            add_node(line,options);
                        }
                        // Want to peek next line, if it's a keyword or empty line we break
                        // the while loop and start over!
                        unsigned int previous_pos = input_file.tellg();
                        getline(input_file, line);
                        if (misc::is_keyword(line) or line.empty())
                        {
                            input_file.seekg(previous_pos);
                            inner_loop_keyword = false;
                        }
                    }                
                }
            }
            if (input_file.eof())
            {
                break;
            }
            row_counter++;   
        }
    }       
};

void Mesh::add_set(std::string line,std::unordered_map<std::string, std::string> options){
    // create new node set, add it to member variable map node_sets
    const auto set_name  = options["NSET"];
    auto       node_set  = std::make_unique<Set<Node*>>(set_name);
    // Add node set to vector of node sets
    nsets.push_back(std::move(node_set));
    // Add entry to map to find node set by name
    node_set_from_node_set_name[set_name] = node_set.get();
}

void Mesh::add_pid(std::unordered_map<std::string, std::string> options){    
    std::string pid_name = options["ELSET"];
    if(misc::is_valid_name(pid_name)){
        std::string mid_name = options["MATERIAL"];
        // Find material
        if(mid_map.find(mid_name) == mid_map.end()){
            std::cout << "ERROR: terminate at line " << row_counter << " in file " << current_inputfile << std::endl;
            std::cout << "ERROR: invalid material name" << std::endl;
            std::terminate();
        }
        auto mid = mid_map[mid_name];
        // Create new pid
        auto pid  = std::make_unique<Pid>(pid_name,mid);
        pid_map[pid_name] = pid.get();
        pids.push_back(std::move(pid));
    }
    else{
        std::cout << "ERROR: terminate at line " << row_counter << " in file " << current_inputfile << std::endl;
        std::cout << "ERROR: invalid entity name" << std::endl;
        std::terminate();
    }
};
void Mesh::add_element(std::string line,std::unordered_map<std::string,std::string> options){
    // these options need to be available to create an element
    auto type = options["TYPE"];
    auto pid_name = options["ELSET"];
    auto pid = pid_map[pid_name];
    auto dataline_items = misc::split_on(line,',');
    unsigned int element_id = std::stoi(dataline_items.at(0));
    // Find node pointers for each node
    std::vector<Node*> element_connectivity;
    // for(const auto& item: dataline_items)
    for(size_t i = 1; i < dataline_items.size();i++)
    {
        unsigned int global_node_id = std::stoi(dataline_items.at(i));
        element_connectivity.push_back(node_map[global_node_id]);
    }
    std::unique_ptr<Element> element;
    if (type == "S3")
    {
        element = std::make_unique<S3>(element_id,element_connectivity,pid,3,3*3,5,1,2,"S3");
    }
    else if (type == "CPS3")
    {
        element = std::make_unique<CPS3>(element_id,element_connectivity,pid,3,3*2,6,1,2,"CPS3");
    }
    else if (type == "CPS4")
    {   
        element = std::make_unique<CPS4>(element_id,element_connectivity,pid,4,4*2,9,4,2,"CPS4");
    }
    else if (type == "C3D10")
    {   
        element = std::make_unique<C3D10>(element_id,element_connectivity,pid,10,30,24,4,3,"C3D10");
    }
    else if (type == "C3D8")
    {      
        element = std::make_unique<C3D8>(element_id,element_connectivity,pid,8,24,12,8,3,"C3D8");
    }
    else if (type == "C3D20")
    {      
        element = std::make_unique<C3D20>(element_id,element_connectivity,pid,20,60,25,27,3,"C3D20");
    }
    else
    {
        std::cout << "*ELEMENT: type=" << type << " not supported" << std::endl;
        return;
    }
    // save what nodes connect to this element. needed for 
    // generating segments for contact
    for (const auto& node : element_connectivity)
    {
        node->connected_elements.push_back(element.get());
    }
    elements.push_back(std::move(element));
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
    
    auto               dataline_items = misc::split_on(line,',');
    // TODO: add support for more node options like coordinate system and stuff
    unsigned int            global_id = std::stoi(dataline_items.at(0));
    auto                            x = std::stof(dataline_items.at(1));
    auto                            y = std::stof(dataline_items.at(2));
    auto                            z = std::stof(dataline_items.at(3));
    auto                         node = std::make_unique<Node>(global_id,x,y,z);
    node_map[global_id]               = node.get();
    nodes.push_back(std::move(node));
    global_2_local_node_id[global_id] = nodes.size()-1;
};