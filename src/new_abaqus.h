#pragma once
#include "mesh.h"

// this class should hold a mesh and results
class New_abaqus
{
private:
public:
    void solve();
    Mesh mesh;
    New_abaqus(/* args */);
    ~New_abaqus();
};


