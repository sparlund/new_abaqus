#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include "../external_libs/Eigen/Dense"
#include "../external_libs/Eigen/Sparse"
#include <utility>
#include "mid.h"
#include "pid.h"
#include "node.h"
#include "element.h"
#include "misc_string_functions.h"
#include "set.h"

class Mesh
{
private:
    // static counters, if we should decide to
    // add another Mesh object they will still
    // be available
    static unsigned int node_counter;
    static unsigned int element_counter;
    static unsigned int pid_counter;
    static unsigned int mid_counter;
    const float penalty_value = 1e36;
    // arrays of pointers to nodes and elements
    std::vector<std::shared_ptr<Node>> nodes;
    std::unordered_map<unsigned int,unsigned int> global_2_local_node_id;
    std::unordered_map<unsigned int,std::shared_ptr<Node>> node_id_2_node_pointer;
    std::unordered_map<std::string, std::shared_ptr<Mid>> mid_name_2_mid_pointer;
    std::vector<std::shared_ptr<Element>> elements;
    // array of PID's belonging to the Mesh
    std::vector<std::shared_ptr<Pid>> pids;
    // dict of pid name to pid object
    std::unordered_map<std::string,std::shared_ptr<Pid>> pid_map;

    // If elements are created before the PID we want to create a boilerplate PID with the
    // specified name
    std::vector<std::string> pid_2_create; 
    // array of MID's belonging to the Mesh
    std::vector<std::shared_ptr<Mid>> mids;
    // array of elements and their connectivity
    // Method to add Node
    void add_node(std::string line,std::unordered_map<std::string, std::string> options);
    void add_element(std::string line,std::unordered_map<std::string, std::string> options);
    void add_load(std::string line);
    void add_boundary(std::string line,std::unordered_map<std::string, std::string> options);
    void add_pid(std::unordered_map<std::string, std::string> options);
    void add_mid(std::unordered_map<std::string, std::string> options);
    void add_set(std::string line,std::unordered_map<std::string, std::string> options);
    std::vector<std::shared_ptr<Set<Node>>> node_sets;
    std::unordered_map<std::string,std::shared_ptr<Set<Node>>> node_set_from_node_set_name; 
    std::vector<std::shared_ptr<Set<Element>>> element_sets;
    bool static_analysis=false;
    bool eigenvalue_analysis=false;
    std::string eigenvalue_solution_method;
    std::string analysis_name;
    void print_matrix_to_mtx(Eigen::SparseMatrix<float>,std::string);
public:    
    void set_analysis_name(std::string string_analysis_name){this->analysis_name = string_analysis_name;return;}
    unsigned int get_pid_counter(){return pid_counter;};
    unsigned long get_number_of_elements(){return this->elements.size();}
    // ndofs is counter for total number of degrees of freedom for the mesh
    unsigned int ndofs = 0;
    // global stiffness matrix
    Eigen::SparseMatrix<float> K;
    // global mass matrix
    Eigen::SparseMatrix<float> M;
    // global load vector
    Eigen::SparseVector<float> f;
    unsigned int number_of_modes_to_find;
    // size(eigenvalues) = number_of_modes*1
    Eigen::Matrix<float,Eigen::Dynamic,1> eigenvalues;
    Eigen::Matrix<float,Eigen::Dynamic,1> eigenfrequencies;
    // size(eigenvectors) = ndofs*number_of_modes
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> eigenvectors;
    // solution to Ku=f
    Eigen::Matrix<float,Eigen::Dynamic,1> u;
    std::vector<std::pair<unsigned int,float>> f_to_be_added;
    // bc: global dof, value
    std::vector<std::pair<unsigned int,float>> bc;
    void assemble();
    void solve();
    void solve_static();
    void solve_eigenfrequency();
    void export_2_vtk();
    Mesh();
    void read_file(std::string filename);
    std::vector<unsigned int> Mesh_connectivity;
    ~Mesh();
};