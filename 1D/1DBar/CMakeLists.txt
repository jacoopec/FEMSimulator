cmake_minimum_required(VERSION 3.10)
project(FEM_1D)
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add executable and source files
add_executable(fem_program
    main.cpp
    fem_utils.cpp
)

# Optional: show compile commands (for debugging)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
