#pragma once
#include <vector>
#include <memory>
#include <string>
#include <Eigen/Sparse>
#include "node.h"
#include "pid.h"

// virtual class
class Element{
private:
    
public:
    static unsigned int element_counter;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    virtual std::vector<std::shared_ptr<Node>> get_connectivity()=0;
    virtual Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>get_Ke()=0;
    virtual Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>get_Me()=0; 
    virtual std::shared_ptr<Pid> get_pid()=0;
    virtual std::vector<unsigned int> get_element_dof_ids()=0;
    virtual unsigned short get_element_ndofs()=0;
    virtual unsigned short get_element_nnodes()=0;
    virtual unsigned int   get_id()=0;
    virtual std::string    get_element_type()=0;
    virtual unsigned short get_vtk_identifier()=0;
    void print_element_info_to_log();
    static float inv_div_by1(float in);
    virtual ~Element(){};
    Element(){};
};
