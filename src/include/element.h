#pragma once
#include "node.h"
#include "pid.h"

#include <Eigen/Sparse>
#include <memory>
#include <string>
#include <vector>
#include<utility>

// ABC
using Segment    = std::pair<Node*, Node*>;
using dynMatrix  = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>;
using Scalar     = Eigen::Product<Eigen::MatrixXf, Eigen::MatrixXf, 0>;
class Element{
protected:
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
    dynMatrix  Ke;
    dynMatrix  Me;
    Eigen::Matrix<float,Eigen::Dynamic,1>               fe;
    dynMatrix  coord; 
    std::vector<Segment>                                segments;
    std::vector<dynMatrix>  B;
    void                                                setup_coord();
    void                                                setup_dofs();
public:
    const unsigned int                                  id;
    virtual void                                        calculate_Ke()=0;
    virtual void                                        calculate_Me()=0;
    virtual std::vector<Scalar>                         calculate_stress(dynMatrix,dynMatrix);
    virtual std::vector<Scalar>                         calculate_strain(dynMatrix,dynMatrix);
    virtual std::vector<Segment>&                       get_segments(Node*);
    std::vector<Node*>                                  get_connectivity() const ;
    Pid*                                                get_pid() const ;
    std::vector<unsigned int>                           get_element_dof_ids() const ;
    unsigned short                                      get_element_ndofs() const ;
    unsigned short                                      get_element_nnodes() const ;
    unsigned short                                      get_dimensions() const;
    std::string                                         get_element_type() const ;
    unsigned short                                      get_vtk_identifier() const ;
    float                                               get_weight() const ;
    float                                               get_volume() const ;
    unsigned short                                      get_ngp() const ;
    unsigned int                                        get_element_counter() const ;
    dynMatrix                                           get_Ke() const;
    dynMatrix                                           get_Me() const;
    void                                                print_element_info_to_log() const ;
    // small function used when going from eigenmode to eigenfrequency
    float inv_div_by1(float in) const;
    virtual ~Element() = default;
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
