#include <string>
#include <iostream>
#include <cstdio>
#include <vector>
#include "src/new_abaqus.h"
#include "src/mesh.h"
#include "src/misc_string_functions.h"

int main(int argc, char const *argv[])
{
    // init new_abaqus model object
    New_abaqus na;
    // split input first on '/' then '.', and extract filename
    std::string input_filename = std::string(argv[1]);
    std::string logfile_filename = misc::split_on(misc::split_on(input_filename,'/').back(),'.').at(0) + ".log";
    // re-direct stdout & stderr  from terminal to log-file
    std::freopen(logfile_filename.c_str(), "w", stdout );
    std::freopen(logfile_filename.c_str(), "w", stderr );
    na.mesh.read_file(input_filename);
    na.mesh.assemble();
    na.mesh.solve();
    na.mesh.export_2_vtk();
    return 0;
}




