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
    std::shared_ptr<Mid> mid;
    std::unordered_map<std::string,std::string> options;
public:
    Pid(std::string name,std::shared_ptr<Mid> mid);
    std::unordered_map<std::string,std::string>* get_options();
    std::shared_ptr<Mid> get_mid();
};
