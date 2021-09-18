#pragma once
#include <iostream>
class Dof
{
private:
    float value;
    static unsigned int global_dof_id_counter;
public:
    const unsigned int id;
    unsigned int get_global_dof_id_counter();
    float get_value();
    Dof();
    ~Dof();
};

