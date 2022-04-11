#pragma once

#include <memory>
#include <vector>

using namespace std;

template<class T>
class Matrix {
public:
    // Default constructor
    Matrix();

    // Constructor where we want to preallocate ourselves
    Matrix(int rows, int cols, bool preallocate);

    // Constructor where we already have allocated memory outside
    // We use a shared pointer to the values
    Matrix(int rows, int cols, shared_ptr<T[]> values_ptr);


    // Destructor
    virtual ~Matrix();

    // Print out the values in our matrix
    void printValues();

    virtual void printMatrix();

    // Multiplication: matrix-matrix
    void matMatMult(Matrix<T> &mat_left, Matrix<T> &output);

    // Multiplication: matrix-vector
    T matVecMult_pointer(Matrix<T> &A, T *x, T *out);

    void matVecMult_reference(vector<T> &vec, vector<T> &output);

    // Declare the shared pointer to the values(data)
    shared_ptr<T[]> values;
    int rows = -1;
    int cols = -1;
    int size_of_values = -1;

// We want our subclass to know about this
protected:
    bool preallocated = false;

};
