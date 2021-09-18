#include "dof.h"
#include <iostream>

unsigned int Dof::global_dof_id_counter = 0;
unsigned int Dof::get_global_dof_id_counter(){return global_dof_id_counter;};
float Dof::get_value(){return value;};

Dof::~Dof(){
    --global_dof_id_counter;
}

Dof::Dof():value{0.f},id{Dof::global_dof_id_counter}{
    ++global_dof_id_counter;
}
