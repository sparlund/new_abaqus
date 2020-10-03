#pragma once
#include <string>
#include <vector>
#include <unordered_map> 
#include <array> 
#include "mid.h"
class Pid
{
private:
    Mid* mid;
    // 
public:
    Mid* get_mid(){return this->mid;};
    std::string name;
    std::string material_name;
    Pid(std::string name,std::string material_name);
    std::unordered_map<std::string,std::string> options;
    ~Pid();
};
