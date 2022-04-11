#include <iostream>
#include <vector>
#include <memory>
#include "Matrix.h"
#include "Matrix.cpp"
#include "DenseSolver.h"
#include "DenseSolver.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include "sparseSolver.h"
#include "sparseSolver.cpp"
#include "test.h"
#include "compare.h"

using namespace std;

int main() {
/** -----------------------------Menu----------------------------- **/
    int dense_sparse = 1;

    while (dense_sparse != 0) {
        cout << "\nPlease choose the type of matrix you want to solve" << endl;
        cout << "1: Dense" << endl << "2: Sparse" << endl << "3: Test all solvers(Dense and Sparse)" << endl
             << "0: exit" << endl;
        cin >> dense_sparse;
        if (dense_sparse == 1) {
            int dense_choose;
            cout << "Please choose the dense method you want to use" << endl;
            cout << "1: Gauss Solver" << endl << "2: LU Solver with partial pivoting" << endl
                 << "3: Jacobi method" << endl << "4: Gauss-Seidel method" << endl << "5: Test all the dense solver"
                 << endl
                 << "6: Compare the performance of dense solvers" << endl << "0: exit" << endl;
            cin >> dense_choose;
            if (dense_choose == 1) {
                run_GaussSolver_test();
            } else if (dense_choose == 2) {
                run_LUSolver_test();
            } else if (dense_choose == 3) {
                run_Jacobi_test();
            } else if (dense_choose == 4) {
                run_GaussSeidel_test();
            } else if (dense_choose == 5) {
                run_GaussSolver_test();
                run_LUSolver_test();
                run_Jacobi_test();
                run_GaussSeidel_test();
            } else if (dense_choose == 6) {
                run_performance();
            }
        } else if (dense_sparse == 2) {
            int sparse_choose;
            cout << "Please choose the sparse method you want to use" << endl;
            cout << "1: Cholesky method" << endl << "2: Conjugate Gradient method" << endl << "3: Test all the methods"
                 << endl
                 << "0: exit" << endl;
            cin >> sparse_choose;
            if (sparse_choose == 1) {
                run_cholesky_tests();
            } else if (sparse_choose == 2) {
                run_CG_tests();
            } else if (sparse_choose == 3) {
                run_cholesky_tests();
                run_CG_tests();
            }
        } else if (dense_sparse == 3) {
            run_tests();
        }
    }
    return 0;
}