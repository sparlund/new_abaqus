#pragma once
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include "../element.h"
#include "../pid.h"
#include "../node.h"

// S3 is 3 node tria shell element
class S3 : public Element
{
protected:
    unsigned int id;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    std::vector<unsigned int> dofs_id;
    std::vector<float> detJ;
    float area, volume, weight;
    static const std::string element_type;
    static const unsigned short nnodes          = 3;
    static const unsigned short ndofs           = 9; // 3*3
    static const unsigned short vtk_identifier  = 5;
    static const unsigned short ngp             = 1;
    Eigen::Matrix<float,ndofs,1> fe;
    Eigen::Matrix<float,ndofs,ndofs> Me;
    Eigen::Matrix<float,ndofs,ndofs> Ke;
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
    S3(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~S3();
};







