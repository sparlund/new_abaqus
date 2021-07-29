#pragma once
#include <vector>
#include <memory>
#include <string>
#include "../external_libs/Eigen/Sparse"
#include "node.h"
#include "pid.h"

class Element{
protected:
    unsigned int id;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    std::vector<unsigned int> dofs_id;
    std::vector<float> detJ;
    float area, volume, weight;
    const unsigned short nnodes,ndofs,vtk_identifier,ngp,dimensions;
    std::string element_type;
    const std::vector<std::vector<float>>* gauss_points;
    const std::vector<float>* gauss_weights;
    static unsigned int element_counter;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Ke;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Me;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> coord; 
    Eigen::Matrix<float,Eigen::Dynamic,1> fe;
public:
    virtual void calculate_Ke()=0;
    virtual void calculate_Me()=0;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Ke();
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Me();
    std::vector<std::shared_ptr<Node>> get_connectivity();
    std::shared_ptr<Pid>               get_pid();
    std::vector<unsigned int>          get_element_dof_ids();
    unsigned short                     get_dimensions();
    unsigned short                     get_element_ndofs();
    unsigned short                     get_element_nnodes();
    unsigned int                       get_id();
    std::string                        get_element_type();
    unsigned short                     get_vtk_identifier();
    float                              get_weight();
    float                              get_volume();
    unsigned short                     get_ngp();
    void print_element_info_to_log();
    // small function used when going from eigenmode to eigenfrequency, TODO move to misc::
    static float inv_div_by1(float in);
    virtual ~Element();
    Element(unsigned short nnodes,
            unsigned short ndofs,
            unsigned short vtk_identifier,
            unsigned short ngp,
            unsigned short dimensions,
            std::string element_type,
            unsigned int id,
            std::vector<std::shared_ptr<Node>> connectivity,
            std::shared_ptr<Pid> pid);
};
