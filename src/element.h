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
    const unsigned int                                  id;
    std::vector<std::shared_ptr<Node>>                  connectivity;
    std::shared_ptr<Pid>                                pid;
    std::vector<unsigned int>                           dofs_id;
    std::vector<float>                                  detJ;
    std::string                                         element_type;
    float                                               area, volume, weight;
    static unsigned int                                 element_counter;
    const unsigned short                                nnodes;
    const unsigned short                                ndofs;
    const unsigned short                                vtk_identifier;
    const unsigned short                                ngp;
    const unsigned short                                dimensions;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>  Ke;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>  Me;
    Eigen::Matrix<float,Eigen::Dynamic,1>               fe;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>  coord; 
public:
    virtual void                                       calculate_Ke()=0;
    virtual void                                       calculate_Me()=0;
    std::vector<std::shared_ptr<Node>>                 get_connectivity();
    std::shared_ptr<Pid>                               get_pid();
    std::vector<unsigned int>                          get_element_dof_ids();
    unsigned short                                     get_element_ndofs();
    unsigned short                                     get_element_nnodes();
    unsigned int                                       get_id();
    std::string                                        get_element_type();
    unsigned short                                     get_vtk_identifier();
    float                                              get_weight();
    float                                              get_volume();
    unsigned short                                     get_ngp();
    unsigned int                                       get_element_counter();
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Ke();
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Me(); 
    void print_element_info_to_log();
    // small function used when going from eigenmode to eigenfrequency
    static float inv_div_by1(float in);
    virtual ~Element(){};
    Element(unsigned int                        id,
            std::vector<std::shared_ptr<Node>>  connectivity,
            std::shared_ptr<Pid>                pid,
            const unsigned short                nnodes,
            const unsigned short                ndofs,
            const unsigned short                vtk_identifier,
            const unsigned short                ngp,
            const unsigned short                dimensions);
};
