#include <iostream>
#include <cmath>
#include "DenseSolver.h"
#include "Matrix.h"
#include "utilities.h"
#include <stdexcept>
#include <vector>

using namespace std;

/** Constructor
 * @param A: the coefficient matrix of the linear system
 * @param b: the constant vector on the right of the equation
 * ToDo: set the linear system to the dense solver and check the dimensions match or not
 */
template<class T>
DenseSolver<T>::DenseSolver(Matrix<T> &A, vector<T> &b) : A(A), b(b) {
    // Check dimensions of A and b
    if (A.cols != b.size()) {
        cerr << "Input dimensions don't match" << endl;
        return;
    }
}


/** Constructor(creates a random matrix)
 * @param size: the rows and cols of the coefficient matrix which should be a square matrix
 * ToDo: generate a random diagonally dominant square matrix A using random seed and generate a constant vector b
 */
template<class T>
DenseSolver<T>::DenseSolver(int size) : size(size) {
    // Set random seed
    srand(time(NULL));

    // Create random diagonally dominant square matrices
    A = Matrix<T>(size, size, true);
    b.reserve(size);
    for (int i = 0; i < size; i++) {
        b.push_back(T(rand() % 10 + 1));
        for (int j = 0; j < size; j++) {
            // Set the elements on the diagonal
            if (i == j) {
                A.values[i * size + j] = T(rand() % 100000 + 10);
            } else {
                // Set the other elements
                A.values[i * size + j] = T(rand() % 10);
            }
        }
    }
}


/** Constructor(copy from another dense solver object)
 * @param DS2: another dense solver object
 * ToDo: set the linear system using another dense solver object(A and b)
 */
template<class T>
DenseSolver<T>::DenseSolver(const DenseSolver<T> &DS2) {
    A = DS2.A;
    b = DS2.b;
}


/** Destructor **/
template<class T>
DenseSolver<T>::~DenseSolver() {}


/** Gaussian elimination and back-substitution(basic but important method)
 * ToDo:
 * 1) Gauss Elimination: we will get an upper-triangular matrix
 * 2) Back-substitution: from the last row to the first row, we will
 *    get the solution vector and it will be stored in the vector of b
 * solution location: stored in the vector of b
 * **/
template<class T>
void DenseSolver<T>::GaussSolver() {
    /** Gauss Elimination **/
    int n = A.cols;
    for (int i = 0; i < n - 1; i++) {
        for (int k = i + 1; k < n; k++) {
            // Get an nonzero scalar alpha
            double alpha = A.values[k * A.cols + i] / A.values[i * A.cols + i];
            for (int j = i; j < n; j++) {
                // Add a multiple of one row to another row for the coefficient matrix A
                A.values[k * A.cols + j] = A.values[k * A.cols + j] - alpha * A.values[i * A.cols + j];
            }
            // Add a multiple of one row to another row for the constant vector b
            b[k] = b[k] - alpha * b[i];
        }
    }

    /** Back-substitution **/
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            if (i != j) {
                b[i] = b[i] - A.values[i * A.rows + j] * b[j];
            }
        }
        b[i] = b[i] / A.values[i * A.cols + i];
    }
}


/** Jacobi iteration method
 * ToDo:
 * 1) Initialize the guess vector
 * 2) Update the guess vector statically based on the residual until reach the ending conditions
 * @param x: the solution vector(store the guess values of each iteration temporarily
 * @param tol: the tolerance can be allowed(one of the iterations ending condition)
 * @param it_max: the biggest iteration times can be allowed(one of the iterations ending condition)
 */
template<class T>
void DenseSolver<T>::Jacobi(vector<T> &x, double &tol, int &it_max) {
    T residual = 0;
    T sum = 0;
    vector<T> output_b(x.size(), 0);
    vector<T> x_old(x.size(), 0);

    /** Initialize the guess vector **/
    for (int i = 0; i < x.size(); i++) {
        x[i] = 0;
    }

    int cnt;
    for (cnt = 0; cnt < it_max; ++cnt) {
        for (int i = 0; i < A.rows; ++i) {
            sum = 0;
            for (int j = 0; j < A.cols; ++j) {
                if (j != i) {
                    /** Update the guess vector statically(using the guess vector of last iteration) **/
                    sum += A.values[i * A.cols + j] * x_old[j];
                }
            }
            // Update the guess vector
            x[i] = (1.0 / A.values[i + i * A.rows]) * (b[i] - sum);;
        }
        residual = calResidual_vector(A, x, b);
        if (residual < tol) {
            break;
        }
        // After this time iteration, update the guess vector of last iteration
        x_old = x;
    }
}


/** Gauss-Seidel iteration method
 * ToDo: the implement of the solver is similar to the Jacobi solver
 * 1) Initialize the guess vector
 * 2) Update the guess vector dynamically based on the residual until reach the ending conditions
 * @param x: the solution vector(store the guess values of each iteration temporarily
 * @param tol: the tolerance can be allowed(one of the iterations ending condition)
 * @param it_max: the biggest iteration times can be allowed(one of the iterations ending condition)
 */
template<class T>
void DenseSolver<T>::GaussSeidel(vector<T> &x, double &tol, int &it_max) {
    T residual = 0;
    T sum = 0;
    vector<T> output_b(x.size(), 0);
    for (int i = 0; i < x.size(); i++) {
        x[i] = 0;
    }
    int cnt;
    for (cnt = 0; cnt < it_max; ++cnt) {
        for (int i = 0; i < A.rows; ++i) {
            sum = 0;
            for (int j = 0; j < A.cols; ++j) {
                if (j != i) {
                    // Note that we use the latest guess values to update the guess vector
                    sum += A.values[i * A.cols + j] * x[j];
                }
            }
            x[i] = (1.0 / A.values[i + i * A.rows]) * (b[i] - sum);;
        }
        residual = calResidual_vector(A, x, b);
        if (residual < tol) {
            break;
        }
    }
}


/** LU Decomposition
 * @param LU: the coefficient matrix and then the Lower and Upper one after decomposation
 * @return max_row_index: the index of the biggest value of each row
 * ToDo:
 * 1) Find the biggest value's index of each row in the coefficient matrx
 * 2) Implicit pivoting which can be used to make it independent of scaling of equations
 * 3) Partial pivoting which can be used to ensure the stability of the LUSolver
 */
template<class T>
vector<int> DenseSolver<T>::LUDecomposition(Matrix<T> &LU) {
    int max_index;

    vector<int> max_row_index(A.rows);
    vector<T> scaling(A.rows);

    // Copy values into LU from A
    for (int i = 0; i < A.size_of_values; i++) {
        LU.values[i] = A.values[i];
    }

    /** Find the implicit scaling(scaling = 1 / max abs value in each row)
     *  and store scaling factor into a vector **/
    for (int i = 0; i < A.rows; i++) {
        double max = 0.0;
        for (int j = 0; j < A.cols; j++) {
            double temp = abs(LU.values[i * LU.cols + j]);
            if (temp > max)
                max = temp;
        }

        /** Check the matrix is sigular or not **/
        // Note: max == 0 means the pivot is zero which means the matrix is singular,
        // LU Decomposation can not solve it
        if (max == 0)
            throw std::invalid_argument("Matrix is singular");
        scaling[i] = 1.0 / max;
    }

    for (int k = 0; k < A.rows; k++) {
        double max = 0.0;
        for (int i = k; i < A.cols; i++) {
            double temp = scaling[i] * abs(LU.values[i * LU.cols + k]);
            /** Store best pivot row so far **/
            if (temp > max) {
                max = temp;
                max_index = i;
            }
        }

        /** Swap rows(if the row is not the best pivot row) **/
        if (k != max_index) {
            for (int j = 0; j < A.rows; j++) {
                // Swap rows using temp
                double temp = LU.values[max_index * LU.cols + j];
                LU.values[max_index * LU.cols + j] = LU.values[k * LU.cols + j];
                LU.values[k * LU.cols + j] = temp;
            }
            // Update scaling
            scaling[max_index] = scaling[k];
        }
        // Update the col number/index of each row which has the biggest abs value
        max_row_index[k] = max_index;


        /** Inner loop of LU decomposition **/
        for (int i = k + 1; i < A.rows; i++) {
            // Divide by pivot element
            LU.values[i * LU.cols + k] /= LU.values[k * LU.cols + k];
            double temp = LU.values[i * LU.cols + k];

            for (int j = k + 1; j < A.rows; j++) {
                LU.values[i * LU.cols + j] -= temp * LU.values[k * LU.cols + j];
            }
        }
    }
    return max_row_index;
}


/** LUSolver: solve the linear equations L*y=b and the U*x=y
 * @param LU: matrix of A
 * @param max_row_index: the index of the biggest value of each row(get from LUDecomposation function)
 * @param x: the solution vector, we store the solution of L*y=b
 *           and then each value of this vector will be overwritten by the solution of U*x=y
 * @param b: the constant vector of the linear system
 */
template<class T>
void DenseSolver<T>::LUSolver(Matrix<T> &LU, vector<int> &max_row_index, vector<T> &x, vector<T> &b) {
    for (int i = 0; i < LU.rows; i++) {
        x[i] = b[i];
    }

    /** Forward substitution to solve the linear equation L*y=b,
     * note that the solution of y will be stored in the vector of x temporarily **/
    for (int i = 0; i < LU.rows; i++) {
        int max_index = max_row_index[i];
        T tmp = x[max_index];
        x[max_index] = x[i];
        for (int j = 0; j < i; j++) {
            tmp -= LU.values[i * LU.cols + j] * x[j];
        }
        x[i] = tmp;
    }

    /** Backward substitution to solve the linear equation U*x = y,
     *  note that the values of y are stored in the vector of x **/
    for (int i = LU.rows - 1; i >= 0; i--) {
        T tmp = x[i];
        for (int j = i + 1; j < LU.rows; j++) {
            tmp -= LU.values[i * LU.cols + j] * x[j];
        }
        x[i] = tmp / LU.values[i * LU.cols + i];
    }
}

/** Conjugate gradient method **/
template<class T>
void DenseSolver<T>::update_CG(int size, T scalar, vector<T> &x, vector<T> &y) {
    for (int i = 0; i < size; i++) {
        y[i] += x[i] * scalar;
    }
}

template<class T>
void DenseSolver<T>::CG(vector<T> &x, double tol) {
    int size = A.cols;
    vector<T> residual(b.size(), 0);
    vector<T> direction(b.size(), 0);
    vector<T> guess(b.size(), 0);
    vector<T> Ax(b.size(), 0);;
    T alpha, beta;
    int cnt = 0;

    /** Initialization **/
    // Initialize x(0) = 0
    for (int i = 0; i < size; ++i) {
        x[i] = 0;
    }
    // residual(0) = Ax(0) - b
    residual = calResidual_residual(A, x, b);
    // guess = residual_0(accelerate)
    direction = residual;

    for (int i = 0; i < size; ++i) {
        // initialise the guess vector
        A.matVecMult_reference(direction, guess);
        alpha = dotProduce(direction, residual, size) /
                dotProduce(direction, guess, size);
        // update x_k+1
        update_CG(size, alpha, direction, x);
        // update residual_k+1
        residual = calResidual_residual(A, x, b);

        // end the iteration
        if (sqrt(dotProduce(residual, residual, size)) < tol) {
            cout << "The algorithm converged in: " << cnt << " iteration(s)." << endl;
            break;
        } else {
            beta = -dotProduce(residual, guess, size) /
                   dotProduce(direction, guess, size);
            for (int j = 0; j < size; ++j) {
                guess[i] = residual[i] + direction[i] * beta;
            }
            cnt++;
        }
    }
}