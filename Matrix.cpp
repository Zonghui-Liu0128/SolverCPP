#include <iostream>
#include <vector>
#include <memory>
#include "Matrix.h"

using namespace std;

// Default constructor
template<class T>
Matrix<T>::Matrix() {}

// Constructor - preallocate ourselves
template<class T>
Matrix<T>::Matrix(int rows, int cols, bool preallocate) : rows(rows), cols(cols), size_of_values(rows * cols),
                                                          preallocated(preallocate) {
    // If we want to handle memory ourselves in the class
    if (this->preallocated) {
        shared_ptr<T[]> data_ptr(new T[this->size_of_values]);
        this->values = data_ptr;
    }
}

// Constructor - setting the value of our pointer
template<class T>
Matrix<T>::Matrix(int rows, int cols, shared_ptr<T[]> values_ptr) : rows(rows), cols(cols),
                                                                    size_of_values(rows * cols),
                                                                    values(values_ptr) {}


// Destructor
template<class T>
Matrix<T>::~Matrix() {}

// Print the values in the matrix
template<class T>
void Matrix<T>::printValues() {
    cout << "Printing values:" << endl;
    for (int i = 0; i < this->size_of_values; i++) {
        cout << this->values[i] << "\t";
    }
    cout << endl << endl;
}

// Print the whole matrix
template<class T>
void Matrix<T>::printMatrix() {
    cout << "Printing matrix:" << endl;
    for (int j = 0; j < this->cols; j++) {
        for (int i = 0; i < this->rows; i++) {
            // We have explicitly used a row-major ordering here
            cout << this->values[i + j * this->rows] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

// Function to perform a matrix-vector multiplication
template<class T>
T Matrix<T>::matVecMult_pointer(Matrix<T> &A, T *x, T *out) {
    int n = A.rows;

    for (int i = 0; i < A.rows; i++) {
        // we make sure that our output vector is always starting from 0
        out[i] = 0;
        for (int j = 0; j < A.cols; j++) {
            out[i] += A.values[i * A.cols + j] * x[j];
        }
    }
    return *out;
}

// Do matrix vector multiplication for rowmajor
// output =  this * vec
template<class T>
void Matrix<T>::matVecMult_reference(vector<T> &vec, vector<T> &output) {
//    this->printMatrix();
    for (int i = 0; i < this->rows; i++) {
        output[i] = 0;
        for (int j = 0; j < this->cols; j++) {
            output[i] += this->values[i * this->cols + j] * vec[j];
        }
    }
}

// Multiplication: matrix-matrix, output =  this * mat_right
template<class T>
void Matrix<T>::matMatMult(Matrix &mat_right, Matrix &output) {
    // Check our dimensions match
    if (this->cols != mat_right.rows) {
        cerr << "Input dimensions for matrices don't match" << endl;
        return;
    }

    // Check if our output matrix has had space allocated to it
    if (output.values != nullptr) {
        // Check our dimensions match
        if (this->rows != output.rows || mat_right.cols != output.cols) {
            cerr << "Input dimensions for matrices don't match" << endl;
            return;
        }
    }
        // The output hasn't been preallocated, so we are going to do that
    else {
        shared_ptr<T[]> vals(new T[this->rows * mat_right.cols]);
        output.values = vals;
        output.preallocated = true;
    }

    // Set values to zero before hand
    for (int i = 0; i < output.size_of_values; i++) {
        output.values[i] = 0;
    }

    // Now we can do our matrix-matrix multiplication
    for (int i = 0; i < this->rows; i++) {
        for (int k = 0; k < this->cols; k++) {
            for (int j = 0; j < mat_right.cols; j++) {
                output.values[i * output.cols + j] +=
                        this->values[i * this->cols + k] * mat_right.values[k * mat_right.cols + j];
            }
        }
    }
}
