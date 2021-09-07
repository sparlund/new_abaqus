#include "pid.h"
#include "mid.h"
#include <iostream>

std::string Pid::get_name() const {
    return name;
};

Mid* Pid::get_mid(){return this->mid;};

std::unordered_map<std::string,std::string>* Pid::get_options(){return &options;};

Pid::Pid(std::string name, Mid* mid):name{std::move(name)},mid{mid}{
    std::cout << "*SECTION: ELSET = " << this->name << ", MATERIAL = " << this->mid->get_name() << std::endl;
}