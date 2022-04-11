#include <iostream>
#include "acse-zl1021-advProg-E-Matrix.h"

// Constructor - using an initialisation list here
template<class T>
acse<T>::acse(int rows, int cols, bool preallocate): rows(rows), cols(cols), size_of_values(rows * cols),
                                                     preallocated(preallocate) {
    // If we want to handle memory ourselves
    if (this->preallocated) {
        // Must remember to delete this in the destructor
        this->values = new T[size_of_values];
    }
}

// Constructor - now just setting the value of our T pointer
template<class T>
acse<T>::acse(int rows, int cols, T *values_ptr): rows(rows), cols(cols), size_of_values(rows * cols),
                                                  values(values_ptr) {}

// Part1: Constructor - Copy(deep copy)
template<class T>
acse<T>::acse(const acse<T> &other): acse(other.rows, other.cols, true) {
    /**Part4: cblas_dcopy(other.size_of_values, other.values, 1, values, 1);**/
    std::copy(other.values, other.values + size_of_values, values);
}

// destructor
template<class T>
acse<T>::~acse() {
    // Delete the values array
    if (this->preallocated) {
        delete[] this->values;
    }
}

// Just print out the values in our values array
template<class T>
void acse<T>::printValues() {
    std::cout << "Printing values" << std::endl;
    for (int i = 0; i < this->size_of_values; i++) {
        std::cout << this->values[i] << " ";
    }
    std::cout << std::endl;
}

// Explicitly print out the values in values array as if they are a matrix
template<class T>
void acse<T>::printMatrix() {
    std::cout << "Printing matrix" << std::endl;
    for (int j = 0; j < this->rows; j++) {
        std::cout << std::endl;
        for (int i = 0; i < this->cols; i++) {
            // There is an explicit major ordering assumption here
            std::cout << this->values[i + j * this->cols] << " ";
        }
    }
    std::cout << std::endl;
}

// Do matrix matrix multiplication
// There is an explicit major ordering assumption in this function
// output = mat_left * this
template<class T>
void acse<T>::matMatMult(acse &mat_left, acse &output) {

    // Check our dimensions match
    if (this->cols != output.cols) {
        std::cerr << "Input dimensions for matrices don't match" << std::endl;
        return;
    }

    // Check if our output matrix has had space allocated to it
    if (output.values != nullptr) {
        // Check our dimensions match
        if (this->rows != mat_left.cols || mat_left.rows != output.rows) {
            std::cerr << "Input dimensions for matrices don't match" << std::endl;
            return;
        }
    }
        // The output hasn't been preallocated, so we are going to do that
    else {
        output.values = new T[this->rows * mat_left.cols];
        // Don't forget to set preallocate to true now it is protected
        output.preallocated = true;
    }

    // Set values to zero before hand
    for (int i = 0; i < output.size_of_values; i++) {
        output.values[i] = 0;
    }

    // Now we can do our matrix-matrix multiplication
    // CHANGE THIS FOR LOOP ORDERING AROUND
    // AND CHECK THE TIME SPENT
    // Does the ordering matter for performance. Why??
    for (int i = 0; i < this->rows; i++) {
        for (int k = 0; k < this->cols; k++) {
            for (int j = 0; j < mat_left.cols; j++) {
                output.values[i * output.cols + j] +=
                        this->values[i * this->cols + k] * mat_left.values[k * mat_left.cols + j];
            }
        }
    }
}

//Part2: compute output = AT * input
template<class T>
void acse<T>::matVecMultTranspose(const T *input, T *output) {
    // fill the res matrix with 0
    std::fill(output, output + rows, 0.);

    // every element do
    for (int i = 0; i < cols; ++i) {
        /**Part4: cblas_daxpy(rows, input[i], values + i*rows, 1, output, 1);**/
        for (int j = 0; j < rows; ++j) {
            //using row-majoring
            output[j] += values[i * rows + j] * input[i];
        }
    }
}

template
class acse<double>;