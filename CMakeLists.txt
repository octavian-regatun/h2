cmake_minimum_required(VERSION 3.27)
project(h2)

set(CMAKE_CXX_STANDARD 23)

add_executable(h2 main.cpp
        Chromosome.cpp
        Chromosome.h
        Population.cpp
        Population.h
        MathFunction.cpp
        MathFunction.h
)
