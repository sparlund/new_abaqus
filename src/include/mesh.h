#pragma once
#include "element.h"
#include "mid.h"
#include "misc_string_functions.h"
#include "node.h"
#include "pid.h"
#include "set.h"
#include "contact.h"
#include "../../external_libs/Eigen/Dense"
#include "../../external_libs/Eigen/Sparse"

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <utility>

using Segment = std::pair<Node*, Node*>;
class Mesh
{
private:
    unsigned int row_counter = 0;
    std::string current_inputfile = "";
    // the penalty value is used when solving the ODE. It's a relatively large arbitrary number
    const float penalty_value = 1e36;
    // mesh is the owner of all nodes, elements, pids and mids
    std::vector<std::unique_ptr<Node>>                  nodes;
    std::vector<std::unique_ptr<Element>>               elements;
    std::vector<std::unique_ptr<Pid>>                   pids;
    std::vector<std::unique_ptr<Mid>>                   mids;
    std::vector<std::unique_ptr<Set<Node*>>>            nsets;
    std::vector<std::unique_ptr<Set<Element*>>>         esets;
    std::vector<std::pair<unsigned int,float>>          bc;
    std::unordered_map<unsigned int,unsigned int>       global_2_local_node_id;
    // global node id to pointer
    std::unordered_map<unsigned int, Node*>             node_map;
    std::unordered_map<std::string, Mid*>               mid_map;
    std::unordered_map<std::string, Pid*>               pid_map;
    std::unordered_map<std::string, Set<Node*>*>        node_set_from_node_set_name; 
    std::vector<std::pair<unsigned int,float>>          f_to_be_added;
    // The following method are used to add entities from the input file to create the FE model
    void add_node(std::string line,std::unordered_map<std::string, std::string> options);
    void add_element(std::string line,std::unordered_map<std::string, std::string> options);
    void add_load(std::string line, std::unordered_map<std::string, std::string> options);
    void add_boundary(std::string line,std::unordered_map<std::string, std::string> options);
    void add_pid(std::unordered_map<std::string, std::string> options);
    void add_mid(std::unordered_map<std::string, std::string> options);
    void add_set(std::string line,std::unordered_map<std::string, std::string> options);
    bool static_analysis=false;
    bool eigenvalue_analysis=false;
    std::unique_ptr<Contact> contact;
    bool contact_mechanics_enabled=false;
    bool steady_state_dynamics_analysis=false;
    float steady_state_dynamics_lower_limit = 0;
    float steady_state_dynamics_upper_limit = 20e3;
    int   steady_state_dynamics_number_of_points = 1;
    std::string eigenvalue_solution_method;
    std::string analysis_name;
    // The supported keywords needs to be added in this order!
    const std::vector<std::string> keywords = {"*NODE",
                                               "*NSET",
                                               "*MATERIAL",
                                               "*SOLID SECTION",
                                               "*ELEMENT",
                                               "*BOUNDARY",
                                               "*CLOAD",
                                               "*STATIC",
                                               "*FREQUENCY",
                                               "*STEADY STATE DYNAMICS",
                                               "*MATRIX GENERATE"};
    // global stiffness matrix
    Eigen::SparseMatrix<float> K;
    Eigen::SparseMatrix<float> K_with_bc;
    // global mass matrix
    Eigen::SparseMatrix<float> M;
    // global load vector
    Eigen::SparseVector<float> f;
public:    
    auto          get_K() const {return &K;};
    auto          get_M() const {return &M;};
    auto          get_f() const {return &f;};
    void          set_analysis_name(std::string string_analysis_name){this->analysis_name = string_analysis_name;return;}
    Pid*          get_pid_by_name(std::string pid_name) {return pid_map[pid_name];}
    Mid*          get_mid_by_name(std::string mid_name) {return mid_map[mid_name];}
    auto          get_bc()    const {return &bc;}
    auto          get_nsets() const {return &nsets;}
    auto          get_esets() const {return &esets;}
    auto          get_f_to_be_added() const {return f_to_be_added;}
    auto          get_number_of_pids() const {return pids.size();};
    auto          get_number_of_mids() const {return mids.size();};
    auto          get_number_of_nodes() const {return nodes.size();}
    auto          get_number_of_elements() const {return elements.size();}
    unsigned int  get_number_of_dofs() const;
    unsigned int  number_of_modes_to_find;
    // size(eigenvalues) = number_of_modes*1
    Eigen::Matrix<float,Eigen::Dynamic,1> eigenvalues;
    Eigen::Matrix<float,Eigen::Dynamic,1> eigenfrequencies;
    // size(eigenvectors) = ndofs*number_of_modes
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> eigenvectors;
    // solution to Ku=f
    Eigen::Matrix<float,Eigen::Dynamic,1> u;
    struct Mtx_to_print
    {
        bool STIFFNESS = false;
        bool STIFFNESS_WITH_BC = false;
        bool MASS = false;
        bool LOAD = false;
    };
    Mtx_to_print mtx_to_print;
    void matrix_generate();
    void print_matrix_to_mtx(const Eigen::SparseMatrix<float>&,const std::string&) const;
    void assemble();
    void solve();
    void solve_static();
    void solve_static_with_contact();
    void solve_eigenfrequency();
    void solve_steady_state_dynamics();
    void export_2_vtk();
    void read_file(const std::string& filename, const std::string& keyword);
    void read_file_new_method(const std::string& filename);
    Mesh() = default;
    ~Mesh() = default;
};