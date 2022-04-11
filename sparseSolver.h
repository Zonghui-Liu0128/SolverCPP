//
//  sparseSolver.hpp
//  advanced_programming
//
//  Created by 张致昱 on 2022/1/27.
//


//#ifndef _SOLVER_h_
//#define _SOLVER_h_
#pragma once


//#include "Matrix.h"
//#include "Matrix.cpp"
#include "CSRMatrix.h"
//#include "CSRMatrix.cpp"
#include <memory>
#include <set>
#include <vector>

using namespace std;




template <class T>
class sparseSolver
{
public:
    sparseSolver();
    sparseSolver(CSRMatrix<T> &A, unique_ptr<T[]> &b, unique_ptr<T[]> &x,int size);

    // destructor
    ~sparseSolver();
    
    void cholesky_decomp();
    float solve_cholesky();
    float solve_CG(float & tolerance);
    CSRMatrix<T> A;
    unique_ptr<T[]> b = nullptr;
    unique_ptr<T[]> x = nullptr;
    int size;

};

//#endif
