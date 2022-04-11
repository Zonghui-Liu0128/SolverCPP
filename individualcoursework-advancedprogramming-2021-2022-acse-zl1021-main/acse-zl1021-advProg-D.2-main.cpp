/**
 * GitHub ID: acse-zl1021
 * Operating System (including version): macOS Monterey(version 12.1)
 * Compiler (including version):
    * g++ 4.2.1
    * clang 13.0.0(clang-1300.0.29.30)
    * x86_64-apple-darwin21.2.0
 * IDE (including version): CLion 2021.3.2
 */
#include<iostream>
#include "CSRMatrix.h"

int main() {
    //initialize the sparse matrix
    double values[] = {2, 1, 5, 1, 12, 1};
    int row_pos[] = {0, 2, 3, 6};
    int col_idx[] = {0, 2, 1, 0, 2, 3};
    CSRMatrix<double> A(3, 4, 6, values, row_pos, col_idx);
    auto res = A.getInverseDiagonal();
    for (int i = 0; i < 3; ++i) {
        std::cout << res[i] << ' ';
    }
    return 0;
}

