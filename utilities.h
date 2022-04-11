#pragma once

#include <iostream>
#include <vector>
#include "Matrix.h"

template<typename T>
void printVector(vector<T> vec) {
    cout << "Vector is: \n";
    for (int i = 0; i < vec.size(); i++) {
        cout << vec[i] << "\t";
    }
    cout << "\n";
}

template<typename T>
void checkDimensions(Matrix<T> &M1, std::vector<T> &vec) {

    // Check if square matrix
    if (M1.cols != M1.rows) {
        throw std::invalid_argument("Only implemented for square matrix");
    }

    if (M1.cols != vec.size()) {
        throw std::invalid_argument("Dimensions don't match");
    }
}

/** Vector-Vector Multiplication **/
template<typename T>
double vecVecMult(vector<T> v1, vector<T> v2) {
    double result = 0;
    if (v1.size() == v2.size()) {
        for (int i = 0; i < v1.size(); i++) {
            result += v1[i] * v2[i];
        }
    }
    return result;
}

/** Calculate the residual(return value of residual) **/
template<class T>
T calResidual_vector(Matrix<T> &A, vector<T> &x, vector<T> &b) {
    T residual = 0;
    vector<T> output_b(A.cols, 0);

    // output_b = A * x
    A.matVecMult_reference(x, output_b);
    for (int i = 0; i < A.rows; i++) {
        residual += pow(output_b[i] - b[i], 2.0);
    }
    return sqrt(residual);
}

/** Calculate the residual(return a residual vector) **/
template<class T>
vector<T> calResidual_residual(Matrix<T> &A, vector<T> &x, vector<T> &b) {
    vector<T> residual(A.cols, 0);
    vector<T> output_b(A.cols, 0);
    A.matVecMult_reference(x, output_b);
    for (int i = 0; i < A.rows; i++) {
        residual[i] += pow(output_b[i] - b[i], 2.0);
    }
    return residual;
}

template<class T>
T dotProduce(vector<T> &a, vector<T> &b, int size) {
    T result = 0;
    for (int i = 0; i < size; ++i) {
        result += a[i] * b[i];
    }
    return result;
}
