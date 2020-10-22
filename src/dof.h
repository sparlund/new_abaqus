#pragma once

class Dof
{
private:
    /* data */
public:
    static unsigned int global_dof_id_counter;
    unsigned int id;
    float value;
    Dof();
    ~Dof();
};

