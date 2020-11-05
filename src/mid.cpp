#include <iostream>
#include "mid.h"
Mid::Mid(int id, std::string name):id(id),name(name){
    std::cout << "*MATERIAL: name=" << name << std::endl;
};

Mid::~Mid(){};
