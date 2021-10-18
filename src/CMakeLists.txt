add_executable(new_abaqus 
               dof.cpp
               new_abaqus.cpp
               set.cpp
               mesh.cpp
               node.cpp
               pid.cpp
               Gauss.cpp
               misc_string_functions.cpp
               staring_point.cpp
               mid.cpp
)
target_link_libraries(new_abaqus elements)
target_compile_options(new_abaqus PRIVATE -O2)

