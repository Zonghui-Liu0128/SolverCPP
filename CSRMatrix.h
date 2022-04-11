#ifndef _CSRMATRIX_H_
#define _CSRMATRIX_H_

#include "Matrix.h"
#include <memory>


using namespace std;



template <class T>
class CSRMatrix: public Matrix<T>
{
public:

   // constructor where we want to preallocate ourselves
   CSRMatrix(int rows, int cols, int nnzs, bool preallocate);
   // constructor where we already have allocated memory outside
   CSRMatrix(int rows, int cols, int nnzs, unique_ptr<T[]> &values_ptr, unique_ptr<int[]> &row_position, unique_ptr<int[]> &col_index);
   CSRMatrix(CSRMatrix & m1);
    CSRMatrix();



   // destructor
   ~CSRMatrix();

   // Print out the values in our matrix
	virtual void printMatrix();

   // Perform some operations with our matrix
   void matMatMult(CSRMatrix<T>& mat_left, CSRMatrix<T>& output);
   // Perform some operations with our matrix
   void matVecMult(T *input, T *output);

    void dense2sparse(Matrix<T> m1);
    void m_transpose();

    

   // Explicitly using the C++11 nullptr here
    unique_ptr<int[]> row_position = nullptr;
    unique_ptr<int[]> col_index = nullptr;
    unique_ptr<T[]> values = nullptr;

   // How many non-zero entries we have in the matrix
   int nnzs=-1;

// Private variables - there is no need for other classes 
// to know about these variables
private:
   
};

#endif
