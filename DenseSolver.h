#ifndef LINEARSOLVER_DENSESOLVER_H
#define LINEARSOLVER_DENSESOLVER_H

#include <iostream>
#include "Matrix.h"
#include <vector>

using namespace std;

template<class T>
class DenseSolver {
public:
    Matrix<T> A;
    vector<T> b{};
    int size = -1;

    // Constructor
    DenseSolver(Matrix<T> &A, vector<T> &b);

    // Constructor - creats a random matrix of dimensions sizexsize
    DenseSolver(int size);

    // Copy constructor
    DenseSolver(const DenseSolver<T> &DS2);

    // Destructor
    virtual ~DenseSolver();

    /** Gauss Elimination and backward substitution(tested) **/
    void GaussSolver();

    /** LU Solver with partial pivoting **/
    vector<int> LUDecomposition(Matrix<T> &LU);

    void LUSolver(Matrix<T> &LU, vector<int> &max_row_index, vector<T> &x, vector<T> &b);

    /** Jacobi(tested) **/
    void Jacobi(vector<T> &x, double &tol, int &it_max);

    /** Gauss Seidel(tested) **/
    void GaussSeidel(vector<T> &x, double &tol, int &it_max);


    /** Conjugate gradient method **/
    void update_CG(int size, T scalar, vector<T> &x, vector<T> &y);
    void CG(vector<T> &x, double tol);



};

#endif //LINEARSOLVER_DENSESOLVER_H