#include "dof.h"

unsigned int Dof::dof_id_counter = 0;
Dof::Dof()
{
    
    id = dof_id_counter + 1;
    dof_id_counter++;
}

Dof::~Dof()
{
}