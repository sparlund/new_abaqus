#pragma once
#include <string>
#include <vector>
#include <unordered_map> 
#include <array> 
#include <memory>
#include "mid.h"
class Pid
{
private:
    std::string name;
    Mid* mid;
    std::unordered_map<std::string,std::string> options;
public:
    Pid(std::string name, Mid* mid);
    std::unordered_map<std::string,std::string>* get_options();
    Mid* get_mid();
    std::string get_name() const;
};
