#include "pid.h"
#include "mid.h"
#include <iostream>

std::shared_ptr<Mid> Pid::get_mid(){return this->mid;};

std::unordered_map<std::string,std::string>* Pid::get_options(){return &options;};

Pid::Pid(std::string name,std::shared_ptr<Mid> mid):name(name),mid(mid){
    std::cout << "*SECTION: ELSET = " << this->name << ", MATERIAL = " << this->mid->get_name() << std::endl;
}