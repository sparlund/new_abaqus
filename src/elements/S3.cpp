#include "S3.h"

S3::~S3(){}
S3::S3(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid):
    Element(3, // nnodes
            9, // ngp
            5, // vtk identifier
            1, // ngp
            2, // dimensions
            "S3",
            id,
            connectivity,
            pid){
    
}

void S3::calculate_Ke(){};
void S3::calculate_Me(){};
