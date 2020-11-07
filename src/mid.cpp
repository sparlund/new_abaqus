#include <iostream>
#include "mid.h"

const std::vector<std::string> Mid::supported_material_keywords = { "*ELASTIC",
                                                              "*DENSITY"};
Mid::Mid(std::string name):name(name){
    std::cout << "*MATERIAL: name=" << name << std::endl;
};

Mid::~Mid(){};
