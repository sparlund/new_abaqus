#pragma once
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include "../element.h"
#include "../pid.h"
#include "../node.h"


class C3D10 : public Element
{
private:
    unsigned int id;
    static const std::string element_type;
    static const unsigned short nnodes  = 10;
    static const unsigned char ndofs   = 30; // 10*3
    // https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    static const unsigned short vtk_identifier = 24;
    std::vector<unsigned int> dofs_id;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    const unsigned short ngp = 4;   // The target rank of K e is 30 − 6 = 24. Since each Gauss point adds 6 to the rank up to a maximum
                                    // of 24, the number of Gauss points should be 4 or higher. 
    // dim(coord) = nnodes*ndim
    Eigen::Matrix<float,10,3> coord;
public:
    Eigen::Matrix<float,30,30> Ke;
    Eigen::Matrix<float,30,30> Me;
    Eigen::Matrix<float,30,1> fe;  
    std::shared_ptr<Pid> get_pid(){return this->pid;};
    unsigned int get_id(){return id;};
    std::vector<unsigned int> get_element_dof_ids(){return dofs_id;};
    std::vector<std::shared_ptr<Node>> get_connectivity(){return this->connectivity;};
    unsigned short get_element_ndofs(){return ndofs;}
    unsigned short get_element_nnodes(){return nnodes;}
    unsigned short get_vtk_identifier(){return vtk_identifier;}
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Ke(){return Ke;}
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> get_Me(){return Me;}
    std::string get_element_type(){return element_type;}
    C3D10(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid);
    ~C3D10();
};