#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>

class Mid
{
private:
    // If the input file contains anything that indicates a non linear material model, set this boolean false
    bool linear=true;
    float density;
    float E;
    float v;
    std::string name;
    // Constitutive matrix (linear continuum mechanics)
    Eigen::Matrix<float,6,6> D_3D_linear_continuum_mechanics;
    Eigen::Matrix<float,3,3> D_2D_linear_continuum_mechanics;
public:
    void set_density(float density){this->density = density;};
    void set_E(float E){this->E = E; };
    void set_v(float v){this->v = v;};
    // 
    std::string get_name(){return name;};
    float get_density(){return density;};
    float get_E(){return E;};
    float get_v(){return v;};
    // fill this vector with the supported material keywords in mid.cpp
    static const std::vector<std::string> supported_material_keywords;
    Mid(std::string name);
    ~Mid();
};
