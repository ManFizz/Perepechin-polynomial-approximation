cmake_minimum_required(VERSION 3.27)
project(VKR_1)

set(CMAKE_CXX_STANDARD 23)

add_executable(VKR_1
        src/main.cpp
        include/Chebyshev.hpp
        include/Chebyshev.cpp
        include/Legendre.cpp
        include/Legendre.hpp
        include/DataResult.h
        include/Taylor.cpp
        include/Taylor.hpp
)