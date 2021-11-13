#include "../include/mid.h"
#include <iostream>

const std::vector<std::string> Mid::supported_material_keywords = { "*ELASTIC",
                                                                    "*DENSITY"};
Mid::Mid(const std::string& name):linear{true},density{0.0f},E{0.0f},v{0.0f},name{name}
{
    std::cout << "*MATERIAL: name=" << name << std::endl;
};