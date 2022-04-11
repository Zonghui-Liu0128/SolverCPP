#include <memory>
#include <set>
#include <vector>
#include "sparseSolver.h"
#include <math.h>
#include "iostream"

// Default constructor
template <class T>
sparseSolver<T>::sparseSolver()
{}

template <class T>
sparseSolver<T>::sparseSolver(CSRMatrix<T> &A, unique_ptr<T[]> &b, unique_ptr<T[]> &x, int size) : A(A), size(size)
{
    unique_ptr<T[]> b_(new T[size]);
    unique_ptr<T[]> x_(new T[size]);
    for(int i=0;i<size;i++){
        b_[i] = b[i];
        x_[i] = x[i];
    }
    this->b.reset(b_.release());
    this->x.reset(x_.release());
}

// destructor
template <class T>
sparseSolver<T>::~sparseSolver()
{}

//this method do the decomposition, change the matrix to lower triangular matrix
template <class T>
void sparseSolver<T>::cholesky_decomp(){
    /*here I use vector, because some of the 0 values may change to non-0, so I do not
     know how many elements are there in the matrix. (I need an array with dynamic memory)
     */
    vector<float> values_;
    vector<int> col_index_;
    set<int> col_index_temp;
    unique_ptr<int[]> row_position(new int[this->A.rows+1]);
    
    //this circulation loops through all the lines
    for(int i = 0;i < this->A.rows;i++){
        int row_start = this->A.row_position[i];
        int row_end = this->A.row_position[i+1];
        int col_num = 0;
        //judge whether there are elements before the specific element
        bool row_num_before = false;
        //the column index of each element
        col_num = this->A.col_index[row_start];
        int j = row_start;
        
        //this loop finds all the position appears at lower triangular matrix
        while(j<row_end and col_num<=i){
            for(int k = 0;k <= col_num;k++){
                if(j!=row_start and this->A.col_index[j]==k){
                    row_num_before = true;
                }
                
                for(int n=0;n<k;n++){
                    for(int m=this->A.col_index[j-1];m<col_num;m++){
                        if(row_num_before and this->A.row_position[m]==k){
                            col_index_temp.emplace(k);
                        }
                    }
                }
                
                col_index_temp.emplace(col_num);
            }
            j++;
            col_num = this->A.col_index[j];
        }
        
        if(i==1){
            row_position[i] = 1;
        }
        for (auto iter = col_index_temp.begin(); iter != col_index_temp.end(); ++iter) {
            col_index_.push_back(*iter);
            }

        row_position[i+1] = row_position[i] + col_index_temp.size();
        col_index_temp.clear();
        int row_start_2 = row_position[i];
        int row_end_2 = row_position[i+1];
        
        //this loop gives value to each element in different position of lower triangular matrix
        //find the first element in the matrix
        if(i==0){
            values_.push_back(sqrt(this->A.values[0]));
        }
        else{
            for(int j = row_start_2;j < row_end_2;j++){
                col_num = col_index_[j];
                //give the value in the first column
                if(col_num==0){
                    for(int b = row_start;b < row_end;b++){
                        if(this->A.col_index[b] == col_index_[j]){
                            values_.push_back(this->A.values[b]/values_[0]);
                        }
                    }
                }
                
                //gives the value except elements in the first column and the diagonal
                else if(col_num != i){
                    T num = 0;
                    int row_col_start = row_position[col_index_[j]];
                    int row_col_end = row_position[col_index_[j]+1];;
                    for(int k = row_start_2;k < j;k++){
                        for(int n = row_col_start;n < row_col_end;n++){
                            if(col_index_[k] == col_index_[n]){
                                num -= values_[k] * values_[n];
                            }
                        }
                    }
                    for(int b = row_start;b < row_end;b++){
                        if(this->A.col_index[b] == col_index_[j]){
                            num += this->A.values[b];
                        }
                    }
                    for(int m = row_col_start;m < row_col_end;m++){
                        if(i-1==col_index_[m]){
                            num = num / values_[m];
                        }
                    }
                    values_.push_back(num);
                }
                
                //gives the values to the elements on the diagonal
                if(col_num==i){
                    T num = 0;
                    for(int k = row_start_2;k < row_end_2 and col_index_[k]!=i;k++){
                        num-=pow(values_[k], 2);
                    }
                    for(int m = row_start;m< row_end;m++){
                        if(this->A.col_index[m]==i){
                            num += this->A.values[m];
                        }
                    }
                    values_.push_back(sqrt(num));
                }
                
            }
        }
    }
    unique_ptr<int[]> col_index_final(new int[values_.size()]);
    unique_ptr<T[]> values_final(new T[values_.size()]);
    
    //gives values to arrays
    for(int i=0;i<values_.size();i++){
        col_index_final[i] = col_index_[i];
        values_final[i] = values_[i];
    }
    
    //gives the values of lower triangular matrix to A
    this->A.nnzs = values_.size();
    this->A.values.reset(values_final.release());
    this->A.row_position.reset(row_position.release());
    this->A.col_index.reset(col_index_final.release());
}


//this method calculate the value with lower triangular matrix and upper triangular matrix
template <class T>
float sparseSolver<T>::solve_cholesky(){
    this->cholesky_decomp();
    
    //forward substitution
    float x[this->A.rows];
    x[0] = this->b[0]/this->A.values[0];
    for(int i=1;i<this->A.rows;i++){
        x[i] = this->b[i];
        int row_start = this->A.row_position[i];
        int row_end = this->A.row_position[i+1];
        for(int j=row_start;j<row_end-1;j++){
            x[i] = x[i] - this->A.values[j] * x[this->A.col_index[j]];
            
        }
        x[i]/=this->A.values[row_end-1];
    }
    
    //backward substitution
    this->A.m_transpose();
    float x_new[this->A.rows];
    x_new[this->A.rows-1] = x[this->A.rows-1]/this->A.values[this->A.nnzs-1];
    for(int i=this->A.rows-2;i>=0;i--){
        x_new[i] = x[i];
        int row_start = this->A.row_position[i];
        int row_end = this->A.row_position[i+1];
        for(int j=row_start+1;j<row_end;j++){
            x_new[i] = x_new[i] - this->A.values[j] * x_new[this->A.col_index[j]];
        }
        x_new[i] /= this->A.values[row_start];
    
    }
    float residual = 0;
    //calculate the residual and return the residual
//    cout << "The calculated value is:" << endl;
//    for(int i=0;i<this->A.rows;i++){
//        cout << x_new[i] << " ";
//        residual += pow((x_new[i] - this->x[i]), 2);
//    }
//    cout << endl;
    residual = sqrt(residual);
    return residual;
}

template <typename T>
float sparseSolver<T>::solve_CG(float & tol){
    //the maximum number of iterations
    int max_cycle=100;
    int size = this->A.rows;
    //the value of r dot r
    float r_dotr = 0;
    //the value of r dot r after the update of r value
    float R_dotR = 0;
    //the value of residue, and it is used to decide whether exit the cycle or not
    float r_residual = 0;
    //the value of A dot A_P
    float A_dotP = 0;
    float alpha, beta = 0;
    //the array to store the value of A * P
    float A_P[size];
    float R[size];
    float P[size];
    //the initial value of X
    float X[6] = {9, 0, -2, 3, -2, 5};
    
    //calculate the value of A * P
    this->A.matVecMult(X, A_P);
    //renew the value of R
    for (int i = 0;i < size;i++){
        R[i] = A_P[i] - this->b[i];
        P[i] = -R[i];
    }
    //
    for (int j=0;j<size;j++){
        P[j]=beta * P[j] - R[j];
        r_residual += pow(R[j], 2);
    }
    //calculate the value of residual
    r_residual = sqrt(r_residual);
    
    
    for(int i=0;i<max_cycle;i++){
        r_dotr = 0;
        A_dotP = 0;
        R_dotR = 0;
        r_residual = 0;
        //renew the value of A * P
        this->A.matVecMult(P, A_P);
        //renew the value of R dot R
        for(int j=0;j<size;j++){
            r_dotr += R[j] * R[j];
        }
        //renew the value of P dot A_P
        for(int j=0;j<size;j++){
            A_dotP += P[j] * A_P[j];
        }
        //calculate alpha
        alpha = r_dotr / A_dotP;
        for(int j=0;j<size;j++){
            X[j] += alpha * P[j];
            R[j] += alpha * A_P[j];
            R_dotR += R[j] * R[j];
        }
        //calculate beta
        beta = R_dotR/r_dotr;
        //renew the value of P and residual
        for (int j=0;j<size;j++){
            P[j] = beta * P[j] - R[j];
            r_residual += pow(R[i], 2);
        }
        r_residual = sqrt(r_residual);
        //if the residual<tolerance, then we could exit the cycle
        if(r_residual < tol){
            break;
        }
    }
    return r_residual;
}

