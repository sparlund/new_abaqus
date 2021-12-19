#include <string>
#include <iostream>
#include <cstdio>
#include <vector>
#include "../include/new_abaqus.h"
#include "../include/mesh.h"
#include "../include/misc_string_functions.h"

int main(int argc, char const *argv[])
{
    
    // TODO: move some of the below to inside the new abaqus object
    // split input first on '/' then '.', and extract filename
    std::string input_complete_filename = std::string(argv[1]);
    std::string analysis_name = misc::split_on(misc::split_on(input_complete_filename,'/').back(),'.').at(0);
    std::string logfile_filename = analysis_name + ".log";
    // re-direct stdout & stderr  from terminal to log-file
    std::freopen(logfile_filename.c_str(), "wb", stdout );
    std::freopen(logfile_filename.c_str(), "wb", stderr );
    std::cout << std::endl;
    // init new_abaqus model object
    New_abaqus na;
    na.mesh.set_analysis_name(analysis_name);
    na.mesh.read_file_new_method(input_complete_filename);
    na.mesh.assemble();
    na.mesh.solve();
    return 0;
}




