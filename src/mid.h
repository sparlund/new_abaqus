#pragma once
#include <string>
#include <vector>
#include <unordered_map>
		
class Mid
{
private:
    /* data */
public:
    unsigned int id;
    std::string name;
    float density;
    float E;
    float v;
    static std::vector<std::string> supported_material_keywords;
    Mid(std::string name);
    ~Mid();
};
