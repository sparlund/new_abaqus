#include "../include/S3.h"

S3::S3(unsigned int                        id,
       std::vector<Node*>                  connectivity,
       Pid*                                pid):
Element{id,connectivity,pid,ElementType::S4,3,3*3,5,1,2}{}

// TODO: calculate stiffness- and mass matrices
void S3::calculate_Ke(){};
void S3::calculate_Me(){};
