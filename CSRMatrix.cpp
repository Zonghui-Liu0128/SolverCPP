#include <iostream>
#include "CSRMatrix.h"
#include <memory>
#include <vector>
using namespace std;

// Constructor - using an initialisation list here
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, bool preallocate): Matrix<T>(rows, cols, false), nnzs(nnzs)
{
   // If we don't pass false in the initialisation list base constructor, it would allocate values to be of size
   // rows * cols in our base matrix class
   // So then we need to set it to the real value we had passed in
   this->preallocated = preallocate;

   // If we want to handle memory ourselves
   if (this->preallocated)
   {
      // Must remember to delete this in the destructor
       unique_ptr<T[]> values(new T[this->nnzs]);
       unique_ptr<int[]> row_position(new int[this->rows + 1]);
       unique_ptr<int[]> col_index(new int[this->nnzs]);
      this->values.reset(values.release());
      this->row_position.reset(row_position.release());
      this->col_index.reset(col_index.release());
   }
}
//
// Constructor - now just setting the value of our T pointer
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, unique_ptr<T[]> &values_ptr, unique_ptr<int[]> &row_position, unique_ptr<int[]> &col_index)
{
    this->rows = rows;
    this->cols = cols;
    this->nnzs = nnzs;
    unique_ptr<T[]> values_(new T[this->nnzs]);
    unique_ptr<int[]> row_position_(new int[this->rows + 1]);
    unique_ptr<int[]> col_index_(new int[this->nnzs]);
    for(int i=0;i<this->rows+1;i++){
        row_position_[i] = row_position[i];
    }
    for(int i=0;i<this->nnzs;i++){
        values_[i] = values_ptr[i];
        col_index_[i] = col_index[i];
    }

    this->values.reset(values_.release());
    this->row_position.reset(row_position_.release());
    this->col_index.reset(col_index_.release());
}



//
//copy constructor, deep copy
template <class T>
CSRMatrix<T>::CSRMatrix(CSRMatrix & m1)
{
    
    this->rows = m1.rows;
    this->cols = m1.cols;
    this->nnzs = m1.nnzs;
    unique_ptr<T[]> values_(new T[m1.nnzs]);
    unique_ptr<int[]> row_position_(new int[m1.rows + 1]);
    unique_ptr<int[]> col_index_(new int[m1.nnzs]);
    for(int i=0;i<this->rows+1;i++){
        row_position_[i] = m1.row_position[i];
    }
    for(int i=0;i<this->nnzs;i++){
        values_[i] = m1.values[i];
        col_index_[i] = m1.col_index[i];
    }

    this->values.reset(values_.release());
    this->row_position.reset(row_position_.release());
    this->col_index.reset(col_index_.release());
}



template <class T>
CSRMatrix<T>::CSRMatrix()
{}



// destructor
template <class T>
CSRMatrix<T>::~CSRMatrix()
{
   // Delete the values array
   // The super destructor is called after we finish here
   // This will delete this->values if preallocated is true
}

template <class T>
void CSRMatrix<T>::m_transpose()
{
    unique_ptr<T[]> values_temp(new T[this->nnzs]);
    unique_ptr<int[]> row_position_temp(new int[this->rows + 1]);
    unique_ptr<int[]> col_index_temp(new int[this->nnzs]);
    vector<int> v_row_position;
    vector<int> v_col_index;
    vector<T> v_values;
    
    v_row_position.push_back(0);
    int num = 0;
    for(int i=0;i<this->cols;i++){
        
        for(int j=0;j<this->nnzs;j++){
            if(this->col_index[j] == i){
                num+=1;
                v_values.push_back(this->values[j]);
                for(int k=0;k<this->rows;k++){
                    if(j>=this->row_position[k] and j<this->row_position[k+1]){
                        v_col_index.push_back(k);
                    }
                }
            }

        }
        v_row_position.push_back(num);
    }

    for(int i=0;i<this->rows+1;i++){
        row_position_temp[i] = v_row_position[i];
    }
    for(int i=0;i<this->nnzs;i++){
        values_temp[i] = v_values[i];
        col_index_temp[i] = v_col_index[i];
    }
    
    this->values.reset(values_temp.release());
    this->row_position.reset(row_position_temp.release());
    this->col_index.reset(col_index_temp.release());
}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void CSRMatrix<T>::printMatrix() 
{ 
   std::cout << "Printing matrix" << std::endl;
   std::cout << "Values: ";
   for (int j = 0; j< this->nnzs; j++)
   {  
      std::cout << this->values[j] << " ";      
   }
   std::cout << std::endl;
   std::cout << "row_position: ";
   for (int j = 0; j< this->rows+1; j++)
   {  
      std::cout << this->row_position[j] << " ";      
   }
   std::cout << std::endl;   
   std::cout << "col_index: ";
   for (int j = 0; j< this->nnzs; j++)
   {  
      std::cout << this->col_index[j] << " ";      
   }
   std::cout << std::endl;   
}

// Do a matrix-vector product
// output = this * input
template<class T>
void CSRMatrix<T>::matVecMult(T *input, T *output)
{
   if (input == nullptr || output == nullptr)
   {
      std::cerr << "Input or output haven't been created" << std::endl;
      return;
   }

   // Set the output to zero
   for (int i = 0; i < this->rows; i++)
   {
      output[i] = 0.0;
   }

   // Loop over each row
   for (int i = 0; i < this->rows; i++)
   {
      // Loop over all the entries in this col
      for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
      {
         // This is an example of indirect addressing
         // Can make it harder for the compiler to vectorise!
         output[i] += this->values[val_index] * input[this->col_index[val_index]];

      }
   }
}


// Do matrix matrix multiplication
// output = mat_left * this
template <class T>
void CSRMatrix<T>::matMatMult(CSRMatrix<T>& mat_left, CSRMatrix<T>& output)
{

   // Check our dimensions match
   if (this->cols != output.cols)
   {
      std::cerr << "Input dimensions for matrices don't match" << std::endl;
      return;
   }

   // Check if our output matrix has had space allocated to it
   if (output.values != nullptr) 
   {
      // Check our dimensions match
      if (this->rows != mat_left.cols || mat_left.rows != output.rows)
      {
         std::cerr << "Input dimensions for matrices don't match" << std::endl;
         return;
      }      
   }
   // The output hasn't been preallocated, so we are going to do that
   else
   {
      std::cerr << "OUTPUT HASN'T BEEN ALLOCATED" << std::endl;
   }

   // HOW DO WE SET THE SPARSITY OF OUR OUTPUT MATRIX HERE??
}

template<class T>
void CSRMatrix<T>::dense2sparse(Matrix<T> m1)
{
    int cur = 0;
    int j = 0;

    this->row_position[0] = 0;
    for (int i=0;i<m1.cols*m1.rows;i++){
        if (m1.values[i]!=0){
            this->values[cur] = m1.values[i];
            
            this->col_index[cur] = i % m1.rows;
            cur++;
        }
        if(i%m1.cols==m1.cols-1){
            j++;
            this->row_position[j] = cur;
        }
    }
    return ;
}

