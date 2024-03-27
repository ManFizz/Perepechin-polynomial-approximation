cmake_minimum_required(VERSION 3.27)
project(VKR_1)

set(CMAKE_CXX_STANDARD 23)

add_executable(VKR_1
        src/main.cpp
        src/LegendrePolynomial.cpp
        include/LegendrePolynomials.h
        src/GramSchmidtProcess.cpp
        include/GramSchmidtProcess.h
)
