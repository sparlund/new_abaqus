#pragma once
#include <string>
#include <vector>
#include <unordered_map>
		
class Mid
{
private:
    float density;
    float E;
    float v;
    std::string name;
public:
    void set_density(float density){this->density = density;};
    void set_E(float E){this->E = E; };
    void set_v(float v){this->v = v;};
    // 
    std::string get_name(){return name;};
    float get_density(){return density;};
    float get_E(){return E;};
    float get_v(){return v;};
    static const std::vector<std::string> supported_material_keywords;
    Mid(std::string name);
    ~Mid();
};
