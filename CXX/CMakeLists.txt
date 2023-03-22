cmake_minimum_required(VERSION 3.10)

Project(ISING)

add_executable(${PROJECT_NAME} p3_auan.cxx)
add_compile_options(-march=native -Ofast)
