cmake_minimum_required(VERSION 4.1)
project(LINAL1)

set(CMAKE_CXX_STANDARD 20)

add_executable(LINAL1 main.cpp
        methods/gauss.cpp
        methods/gauss.h
        matrix/matrix.cpp
        matrix/matrix.h
        methods/lu.cpp
        methods/lu.h
        utilits/utilits.cpp
        utilits/utilits.h)
