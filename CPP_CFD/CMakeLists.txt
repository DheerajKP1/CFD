cmake_minimum_required(VERSION 3.0)
project(CFD VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(Python3_FIND_STRATEGY LOCATION)
# set(Python3_ROOT_DIR "C:/Users/dheeraj/miniconda3")
find_package(MPI REQUIRED COMPONENTS CXX)
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
# Find OpenMP
find_package(OpenMP REQUIRED)

add_executable(CFD CFD.cpp)
target_link_directories(CFD PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}$)
# Link OpenMP
target_link_libraries(CFD PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(CFD PRIVATE MPI::MPI_CXX)
target_link_libraries(CFD PRIVATE Python3::Python Python3::NumPy)