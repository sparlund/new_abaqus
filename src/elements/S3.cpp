#include "S3.h"

const std::string S3::element_type = "S3";
S3::~S3(){}
S3::S3(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid):
id(id),connectivity(connectivity),pid(pid){
    
}

void S3::calculate_Ke(){};
void S3::calculate_Me(){};
