#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <utility>
#include "mid.h"
#include "pid.h"
#include "node.h"
#include "element.h"
#include "misc_string_functions.h"

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
    std::string filename;
public:    
    unsigned int get_pid_counter(){return pid_counter;};
    // getters and setters
    unsigned long get_number_of_elements(){return this->elements.size();}
    // Element* get_element(unsigned long id){return this->elements.at(id);}
    // unsigned long get_number_of_nodes(){return this->nodes.size();}
    // Node* get_node(unsigned long id){return this->nodes.at(id);}
    // Pid* get_pid(unsigned int pid_id) {return this->pids.at(pid_id);}
    // Mid* get_mid(unsigned int mid_id) {return this->mids.at(mid_id);}

    // ndofs is counter for total number of degrees of freedom for the mesh
    unsigned int ndofs = 0;
    Eigen::SparseMatrix<float> K;
    Eigen::SparseVector<float> f;
    Eigen::Matrix<float,Eigen::Dynamic,1> u;
    std::vector<std::pair<unsigned int,float>> f_to_be_added;
    // bc: global dof, value
    std::vector<std::pair<unsigned int,float>> bc;
    void assemble();
    void solve();
    void export_2_vtk();


    Mesh();
    void read_file(std::string filename);
    std::vector<unsigned int> Mesh_connectivity;
    ~Mesh();
};