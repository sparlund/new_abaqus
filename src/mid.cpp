#include <iostream>
#include "mid.h"

std::vector<std::string> Mid::supported_material_keywords = { "xyzzy", "plugh", "abracadabra" };
Mid::Mid(std::string name):name(name){
    std::cout << "*MATERIAL: name=" << name << std::endl;
};

Mid::~Mid(){};
