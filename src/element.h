#pragma once
#include <vector>
#include <memory>
#include <string>
#include <Eigen/Sparse>
#include "node.h"
#include "pid.h"

// virtual class
class Element{
protected:
    // unsigned int id;
    // std::vector<std::shared_ptr<Node>> connectivity;
    // std::shared_ptr<Pid> pid;
    // std::vector<unsigned int> dofs_id;
    // std::vector<float> detJ;
    // std::string element_type;
    // float area, volume, weight;
    // unsigned short nnodes,ndofs,vtk_identifier;
    // the following struct should be a static member of all subclasses
    struct Element_type_description
    {

    };
public:
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Ke;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Me;
    static unsigned int element_counter;
    virtual void calculate_Ke()=0;
    virtual void calculate_Me()=0;
    virtual Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Ke()=0;
    virtual Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Me()=0; 
    virtual std::vector<std::shared_ptr<Node>> get_connectivity()=0;
    virtual std::shared_ptr<Pid>               get_pid()=0;
    virtual std::vector<unsigned int>          get_element_dof_ids()=0;
    virtual unsigned short                     get_element_ndofs()=0;
    virtual unsigned short                     get_element_nnodes()=0;
    virtual unsigned int                       get_id()=0;
    virtual std::string                        get_element_type()=0;
    virtual unsigned short                     get_vtk_identifier()=0;
    virtual float                              get_weight()=0;
    virtual float                              get_volume()=0;
    virtual unsigned short                     get_ngp()=0;
    void print_element_info_to_log();
    // small function used when going from eigenmode to eigenfrequency
    static float inv_div_by1(float in);
    virtual ~Element(){};
    Element(){};
};
