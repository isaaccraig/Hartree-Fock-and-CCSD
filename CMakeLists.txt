cmake_minimum_required(VERSION 3.9)
project(HartreeFockMP2)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
set(CMAKE_CXX_STANDARD 11)

add_executable(HartreeFockMP2
        HartreeFock.cpp
        HartreeFock.hpp
        main.cpp
        QuantumUtils.cpp
        QuantumUtils.hpp
        Read.cpp
        Read.hpp)

target_link_libraries (HartreeFockMP2 Eigen3::Eigen)