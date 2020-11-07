#include <iostream>
#include "pid.h"
#include "mid.h"


Pid::Pid(std::string name,std::shared_ptr<Mid> mid):name(name),mid(mid){
    std::cout << "*SECTION: ELSET = " << this->name << ", MATERIAL = " << this->mid->get_name() << std::endl;
}

Pid::~Pid(){}