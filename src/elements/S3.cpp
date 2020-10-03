#include "S3.h"


S3::~S3(){}
S3::S3(unsigned int id, std::vector<std::shared_ptr<Node>> connectivity,std::shared_ptr<Pid> pid):
    id(id),connectivity(connectivity),pid(pid){
        // Want to find elements addition to the stiffness- & load matrix
        C <<  1,  ex(1), ey(1), 0,    0,       0,  
              0,      0,     0, 1,ex(1),   ey(1),
              1,  ex(2), ey(2), 0,    0,       0,  
              0,      0,     0, 1,ex(2),   ey(2),
              1,  ex(3), ey(3), 0,    0,       0,  
              0,      0,     0, 1,ex(3),   ey(3);
    }
