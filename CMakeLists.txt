cmake_minimum_required(VERSION 3.10)
project(VKR)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${CMAKE_SOURCE_DIR}/lib)
link_directories(${CMAKE_SOURCE_DIR}/lib)

file(GLOB SOURCES "src/main.cpp" "lib/bignum.cpp")

add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(VKR chsvlib)