cmake_minimum_required(VERSION 3.10)

Project(ISING)

add_executable(${PROJECT_NAME} p3_auan.cxx)
SET(GCC_COVERAGE_COMPILE_FLAGS " -march=native -O3")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
