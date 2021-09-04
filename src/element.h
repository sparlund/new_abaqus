#pragma once
#include "node.h"
#include "pid.h"
#include <Eigen/Sparse>
#include <memory>
#include <string>
#include <vector>

// virtual class
class Element{
protected:
    const unsigned int                                  id;
    std::vector<Node*>                                  connectivity;
    Pid*                                                pid;
    std::vector<unsigned int>                           dofs_id;
    std::vector<float>                                  detJ;
    float                                               area, volume, weight;
    static unsigned int                                 element_counter;
    const unsigned short                                nnodes;
    const unsigned short                                ndofs;
    const unsigned short                                vtk_identifier;
    const unsigned short                                ngp;
    const unsigned short                                dimensions;
    std::string                                         element_type;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>  Ke;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>  Me;
    Eigen::Matrix<float,Eigen::Dynamic,1>               fe;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>  coord; 
    void                                                setup_coord();
    void                                                setup_dofs();
public:
    virtual void                                        calculate_Ke()=0;
    virtual void                                        calculate_Me()=0;
    std::vector<Node*>                                  get_connectivity() const ;
    Pid*                                                get_pid() const ;
    std::vector<unsigned int>                           get_element_dof_ids() const ;
    unsigned short                                      get_element_ndofs() const ;
    unsigned short                                      get_element_nnodes() const ;
    unsigned int                                        get_id() const ;
    std::string                                         get_element_type() const ;
    unsigned short                                      get_vtk_identifier() const ;
    float                                               get_weight() const ;
    float                                               get_volume() const ;
    unsigned short                                      get_ngp() const ;
    unsigned int                                        get_element_counter() const ;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>  get_Ke() const;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>  get_Me() const;  
    void                                                print_element_info_to_log() const ;
    // small function used when going from eigenmode to eigenfrequency
    float inv_div_by1(float in) const;
    virtual ~Element(){};
    Element(unsigned int                        id,
            std::vector<Node*>                  connectivity,
            Pid*                                pid,
            const unsigned short                nnodes,
            const unsigned short                ndofs,
            const unsigned short                vtk_identifier,
            const unsigned short                ngp,
            const unsigned short                dimensions,
            std::string                         element_type);
};
