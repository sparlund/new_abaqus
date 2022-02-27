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
using dynMatrix  = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using dynVector  = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Scalar     = Eigen::Product<Eigen::MatrixXd, Eigen::MatrixXd, 0>;

class Element
{
protected:
    enum class ElementType
    {
        C3D10,
        C3D20,
        C3D8,
        CPS3,
        CPS4,
        S4,
    };
    friend std::ostream& operator<<(std::ostream& os, const Element::ElementType t);
    std::vector<Node*>                                  connectivity;
    Pid*                                                pid;
    std::vector<unsigned int>                           dofs_id;
    std::vector<double>                                  detJ;
    static unsigned int                                 element_counter;
    dynMatrix Ke;
    dynMatrix Me;
    dynVector fe;
    dynVector f_int;
    dynMatrix  coord; 
    std::vector<Segment>                                segments;
    std::vector<dynMatrix>  B;
    void                                                setup_coord();
    void                                                setup_dofs();
public:
    // TODO: these are not calculated at element conception for every element!
    const double                                         area = 0, volume = 0, weight = 0;
    const unsigned short                                nnodes;
    const unsigned short                                ndofs;
    const unsigned short                                vtk_identifier;
    const unsigned short                                ngp;
    const unsigned short                                dimensions;
    const ElementType                                   elementType;
    const unsigned int                                  id;
    // Let's not care about thickness too much yet
    const double                                         t = 1.0;
    virtual void                                        calculate_Ke()=0;
    virtual void                                        calculate_Me()=0;
    virtual void                                        calculate_f_internal(dynVector u);
    virtual std::vector<Scalar>                         calculate_stress(dynMatrix,dynMatrix);
    virtual std::vector<Scalar>                         calculate_strain(dynMatrix,dynMatrix);
    virtual std::vector<Segment>&                       get_segments(Node*);
    unsigned int                                        get_element_counter() const {return element_counter;};
    std::vector<Node*>                                  get_connectivity() const {return connectivity;};
    Pid*                                                get_pid() const {return pid;};
    std::vector<unsigned int>                           get_element_dof_ids() const { return dofs_id;};
    dynMatrix                                           get_Ke() const {return Ke;};
    dynMatrix                                           get_Me() const {return Me;};
    void                                                print_element_info_to_log() const ;
    // small function used when going from eigenmode to eigenfrequency
    double inv_div_by1(double in) const;
    virtual ~Element() = default;
    Element(unsigned int                        id,
            std::vector<Node*>                  connectivity,
            Pid*                                pid,
            ElementType                         elementType,
            const unsigned short                nnodes,
            const unsigned short                ndofs,
            const unsigned short                vtk_identifier,
            const unsigned short                ngp,
            const unsigned short                dimensions);
};
