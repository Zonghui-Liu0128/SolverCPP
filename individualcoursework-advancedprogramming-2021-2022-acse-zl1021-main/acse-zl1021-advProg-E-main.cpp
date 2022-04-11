/**
 * GitHub ID: acse-zl1021
 * Operating System (including version): macOS Monterey(version 12.1)
 * Compiler (including version):
    * g++ 4.2.1
    * clang 13.0.0(clang-1300.0.29.30)
    * x86_64-apple-darwin21.2.0
 * IDE (including version): CLion 2021.3.2
 */
#include <iostream>
#include "acse-zl1021-advProg-E-Matrix.h"


int main() {
    //the number of rows and cols
    int rows = 10;
    int cols = 10;
    int size_values = rows * cols;

    double A_data[size_values], x_data[rows];

    //initialize the data in matrix A and vector x
    for (int i = 0; i < size_values; ++i) {
        A_data[i] = i * 2;
    }
    for (int i = 0; i < rows; ++i) {
        x_data[i] = i;
    }

    //instantiate the matrix A and vector x
    acse<double> A(rows, cols, A_data);
    acse<double> x(rows, 1, x_data);

    //print matrix A
    A.printMatrix();

    //deep copy the matrix A
    acse<double> A_copy(A);

    //initialize the vector of result
    acse<double> res(rows, 1, true);

    //compute output = AT * input
    A_copy.matVecMultTranspose(x.values, res.values);

    //print the result vector
    std::cout<<std::endl;
    res.printValues();
    return 0;
}