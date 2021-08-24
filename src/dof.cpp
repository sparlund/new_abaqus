#include "dof.h"
#include <iostream>

unsigned int Dof::global_dof_id_counter = 0;
unsigned int Dof::get_global_dof_id_counter(){return global_dof_id_counter;};
float Dof::get_value(){return value;};

Dof::Dof():value{0},id{Dof::global_dof_id_counter}
{
    // std::cout << "dof created, id=" << id << "\n";
    global_dof_id_counter++;
}
