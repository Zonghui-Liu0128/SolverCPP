#pragma once

template<class T>
class acse {
public:

    // constructor where we want to preallocate ourselves
    acse(int rows, int cols, bool preallocate);

    // constructor where we already have allocated memory outside
    acse(int rows, int cols, T *values_ptr);

    // Part1: copy constructor(deep copy)
    acse(const acse<T> &other);

    // destructor
    virtual ~acse();

    // Print out the values in our matrix
    void printValues();

    virtual void printMatrix();

    // Perform some operations with our matrix
    void matMatMult(acse<T> &mat_left, acse<T> &output);

    //Part2: compute output = AT * input
    void matVecMultTranspose(const T *input, T *output);

    // Explicitly using the C++11 nullptr here
    T *values = nullptr;
    int rows = -1;
    int cols = -1;

// We want our subclass to know about this
protected:
    bool preallocated = false;

// Private variables - there is no need for other classes 
// to know about these variables
private:

    int size_of_values = -1;
};