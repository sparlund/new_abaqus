#include "S3.h"

S3::S3(unsigned int                        id,
       std::vector<Node*>                  connectivity,
       Pid*                                pid,
       const unsigned short                nnodes,
       const unsigned short                ndofs,
       const unsigned short                vtk_identifier,
       const unsigned short                ngp,
       const unsigned short                dimensions,
       std::string                         element_type):
Element{id,connectivity,pid,nnodes,ndofs,vtk_identifier,ngp,dimensions,element_type}{}

// TODO: calculate stiffness- and mass matrices
void S3::calculate_Ke(){};
void S3::calculate_Me(){};
