add_executable(new_abaqus 
               new_abaqus.cpp
               staring_point.cpp
)
target_compile_options(new_abaqus PRIVATE -O2)
target_link_libraries(new_abaqus
                      elements
                      properties
                      mesh
)
add_library(elements SHARED
            element.cpp
            C3D10.cpp
            C3D20.cpp
            C3D8.cpp
            CPS3.cpp
            CPS4.cpp
            S3.cpp
)
add_library(properties SHARED
            pid.cpp
            mid.cpp
            set.cpp
)
add_library(mesh SHARED
            contact.cpp
            Gauss.cpp
            mesh.cpp
            dof.cpp
            node.cpp
            misc_string_functions.cpp
)
target_link_libraries(mesh
                      elements
                      properties)
target_compile_options(elements PRIVATE -O2)
target_compile_options(properties PRIVATE -O2)
target_compile_options(mesh PRIVATE -O2)