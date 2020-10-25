#include <string>
#include <iostream>
#include "src/new_abaqus.h"
#include "src/mesh.h"



int main(int argc, char const *argv[])
{
    // init new_abaqus model object
    New_abaqus na;
    if (argc > 1)
    {
        for (int i = 1; i < argc; i++)
        {
            na.mesh.read_file(std::string(argv[1]));
        }
    }
    na.mesh.assemble();
    na.mesh.solve();
    na.mesh.export_2_vtk();
    // na.mesh.about();


    return 0;
}




