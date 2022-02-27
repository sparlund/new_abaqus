#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include "../../external_libs/Eigen/Dense"

class Mid
{
private:
    // If the input file contains anything that indicates a non linear material model, set this boolean false
    bool        linear;
    double       density;
    double       E;
    double       v;
    std::string name;
public:
    // Constitutive matrix (linear continuum mechanics)
    Eigen::Matrix<double,6,6> D_3D_linear_continuum_mechanics;
    Eigen::Matrix<double,3,3> D_2D_linear_continuum_mechanics;
    void        set_density(double density){this->density = density;};
    void        set_E(double E){this->E = E; };
    void        set_v(double v){this->v = v;};
    std::string get_name(){return name;};
    double       get_density(){return density;};
    double       get_E(){return E;};
    double       get_v(){return v;};
    static const std::vector<std::string> supported_material_keywords;
    Mid(const std::string& name);
};
