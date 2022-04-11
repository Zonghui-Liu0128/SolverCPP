cmake_minimum_required(VERSION 3.21)
project(LinearSolver)

set(CMAKE_CXX_STANDARD 23)

add_executable(LinearSolver main.cpp Matrix.h Matrix.cpp
        DenseSolver.h DenseSolver.cpp main.cpp
        test.h TestRunner.h TestRunner.cpp utilities.h
        compare.h CSRMatrix.cpp CSRMatrix.h
        sparseSolver.cpp sparseSolver.h)
