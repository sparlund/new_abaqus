#include "../include/mid.h"
#include <iostream>

const std::vector<std::string> Mid::supported_material_keywords = { "*ELASTIC",
                                                                    "*DENSITY"};
Mid::Mid(const std::string& name):linear{true},density{0.d},E{0.d},v{0.d},name{name}
{
    std::cout << "*MATERIAL: name=" << name << std::endl;
};