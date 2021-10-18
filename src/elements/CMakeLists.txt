add_library(elements STATIC 
            element.cpp
            C3D10.cpp
            C3D20.cpp
            C3D8.cpp
            CPS3.cpp
            CPS4.cpp
            S3.cpp
)
target_compile_options(elements PRIVATE -O2)
