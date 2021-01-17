#pragma once
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include "../element.h"
#include "../pid.h"
#include "../node.h"
#include "../Gauss.h"

class C3D8 : public Element
{
protected:
    unsigned int id;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    std::vector<unsigned int> dofs_id;
    std::vector<float> detJ;
    float area, volume, weight;
    static const std::string element_type;
    static const unsigned short nnodes          = 8;
    static const unsigned short ndofs           = 24; 
    static const unsigned short vtk_identifier  = 12;
    static const unsigned short ngp             = 8;
    static const unsigned short dimensions      = 3;
    Eigen::Matrix<float,ndofs,1> fe;
    Eigen::Matrix<float,ndofs,ndofs> Me;
    Eigen::Matrix<float,ndofs,ndofs> Ke;
    Eigen::Matrix<float,nnodes,dimensions> coord; 
    // Eigen::Matrix<float,3,10> shape_functions(Eigen::Matrix<float,1,4> evaluation_points);
    const std::array<std::array<float, dimensions>,ngp>* gauss_points = &Gauss::_3D::integration_points_2_by_2_by_2;
    const std::array<float, ngp>* gauss_weights = &Gauss::_3D::gauss_weights_2_by_2_by_2;
public:
    std::vector<std::shared_ptr<Node>> get_connectivity(){return this->connectivity;}
    std::shared_ptr<Pid>               get_pid(){return this->pid;}
    std::vector<unsigned int>          get_element_dof_ids(){return dofs_id;}
    unsigned short                     get_element_ndofs(){return ndofs;}
    unsigned short                     get_element_nnodes(){return nnodes;}
    unsigned int                       get_id(){return id;}
    std::string                        get_element_type(){return element_type;}
    unsigned short                     get_vtk_identifier(){return vtk_identifier;}
    float                              get_weight(){return weight;}
    float                              get_volume(){return volume;}
    unsigned short                     get_ngp(){return ngp;};
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Ke(){return Ke;};
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Me(){return Me;}; 
    void calculate_Ke();
    void calculate_Me();
    C3D8(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~C3D8();
};