include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
  )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
  
enable_testing()
include(CodeCoverage.cmake)
append_coverage_compiler_flags()

add_executable(
  UT
  tests.cpp
  "../src/src/C3D10.cpp"
  "../src/src/C3D20.cpp"
  "../src/src/C3D8.cpp"
  "../src/src/contact.cpp"
  "../src/src/CPS3.cpp"
  "../src/src/CPS4.cpp"
  "../src/src/dof.cpp"
  "../src/src/element.cpp"
  "../src/src/Gauss.cpp"
  "../src/src/mesh.cpp"
  "../src/src/mid.cpp"
  "../src/src/misc_string_functions.cpp"
  "../src/src/node.cpp"
  "../src/src/pid.cpp"
  "../src/src/S3.cpp"
  "../src/src/set.cpp"
)
target_link_libraries(UT
  gtest_main
)
target_compile_options(UT PRIVATE -O2)
include(GoogleTest)
SETUP_TARGET_FOR_COVERAGE_LCOV(
        NAME UT_coverage
        EXECUTABLE UT
        EXCLUDE "external_libs/**" "_deps/**" "tests/**" "/usr/include/**" "/usr/local/include/**" "usr/lib/**"
)