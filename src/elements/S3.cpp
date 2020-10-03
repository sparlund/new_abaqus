#include "S3.h"



S3::~S3(){}
S3::S3(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid):
    id(id),connectivity(connectivity),pid(pid){      
    }
