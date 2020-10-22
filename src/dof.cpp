#include "dof.h"
#include <iostream>
unsigned int Dof::global_dof_id_counter = 0;
Dof::Dof()
{
    
    id = global_dof_id_counter;
    std::cout << "dof created, id=" << id << "\n";
    global_dof_id_counter++;
}

Dof::~Dof()
{
}