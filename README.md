# ACSE-5 Group Project

## Prerequisites

This library requires a compiler with C++23 feature support.

### The Environment We Used
- Operating System (including version): macOS Monterey(version 12.1)
- Compiler (including version):
    * g++ 4.2.1
    * clang 13.0.0(clang-1300.0.29.30)
    * x86_64-apple-darwin21.2.0
 * IDE (including version): CLion 2021.3.2

### Quickstart

Run main.cpp, users can choose the type of matrix to solve (1.Dense. 2.Sparse. 3.Test all solvers(Dense and Sparse). 0. exit). 

When the usere choose 1(Dense Matrix), they will get a sub-menu:
- Please choose the dense method you want to use:
1. Gauss Solver
2. LU Solver with partial pivoting
3. Jacobi method
4. Gauss-Seidel method
5. Test all the dense solver
6. Compare the performance of dense solvers
0. exit

When the usere choose 2(Sparse Matrix), they will get another sub-menu:
- Please choose the sparse method you want to use:
1. Cholesky method
2. Conjugate Gradient method
3. Test all the methods
0. exit

## Matrix

Matrix is a template class, so the values can be of any type T. The values of the matrix are stored in form of dynamically allocated array. The memory is managed by a smart point(shared pointer). 

There are three types of constructors: the constructor allocates memory inside the class, the constructor allocates memory outside and using a shared pointer to the values array and the default constructor.

### Properties

- `rows`(`int`): the number of rows
- `cols`(`int`): the number of columns
- `size_of_values`(`int`): the size of values in the matrix
- `values`(`std::shared_ptr<T[]>`): the shared pointer to the array of values
- `preallocated`(`bool`): the flag to make sure the memory will be allocated inside or outside the class

### Methods

- `void printValues()`
- `void printMatrix()`
- `void matMatMult(Matrix<T> &mat_left, Matrix<T> &output)`
- `T matVecMult_pointer(Matrix<T> &A, T *x, T *out)`
- `void matVecMult_reference(vector<T> &vec, vector<T> &output)`


## CSRMatrix

This class inherits from matrix. The real Symmetric positive definite matrix is required when using the Cholesky method and the Conjugate Gradient method.

### Additional Properties

`values(unique_ptr<T[]>)`: The number of non-zero values in the matrix

`row_position(unique_ptr<int[]>)`: used to indicate how many elements there are in each row

`col_index(unique_ptr<int[]>)`: Each element is the column index of the corresponding value.

### Methods

`CSRMatrix(CSRMatrix & m1)`

`void dense2sparse(Matrix<T> m1)`

`void m_transpose()`

## DenseSolver

It implements algorithms to solve the linear equation `A`**`x`**`=`**`b`** for a _dense matrix_ `A` of type Matrix<T>.

### Properties
- `A`(`Matrix<T>`): A matrix object.
- **`b`**(`std::vector<T>`): a constant vector representing the right hand side of the equation.
- `size`(`int`): the rows/cols of the coefficient matrix.

### Usage
- To solve a linear system using the **Gauss Solver**:

```cpp
// Define A and b beforehand
DenseSolver<double> denseSolver = DenseSolver(A, b);
denseSolver.GaussSolver();
```

- To solve a linear system using the **LU Solver with partial pivoting**:

```cpp
// Define A and b beforehand
DenseSolver<double> denseSolver = DenseSolver(A, b);

// Define matrix object LU
Matrix<double> LU(size, size, true);
vector<int> piv = denseSolver.LUDecomposition(LU);
denseSolver.LUSolver(LU, piv, x, b);
```

- To solve a linear system using the **Jacobi algorithm**:

```cpp
// Set the ending conditions
double tol = 1e-6; 
int it_max = 1000; 

// define A and b beforehand
DenseSolver<double> denseSolver = DenseSolver(A, b);
vector<double> x(b.size(), 0);

denseSolver.Jacobi(x, tol, it_max);
```

- To solve a linear system using the **Gauss-Seidel algorithm**:

```cpp
// Set the ending conditions
double tol = 1e-6; 
int it_max = 1000; 

// define A and b beforehand
DenseSolver<double> denseSolver = DenseSolver(A, b);
vector<double> x(b.size(), 0);

denseSolver.GaussSeidel(x, tol, it_max);
```

### Methods

- `void GaussSolver()`
- `vector<int> LUDecomposition(Matrix<T> &LU)`
- `void LUSolver(Matrix<T> &LU, vector<int> &max_row_index, vector<T> &x, vector<T> &b)`
- `void Jacobi(vector<T> &x, double &tol, int &it_max)`
- `void GaussSeidel(vector<T> &x, double &tol, int &it_max)`

**Constructors**:

- `DenseSolver(Matrix<T> &A, vector<T> &b)`
- `DenseSolver(int size)`
- `DenseSolver(const DenseSolver<T> &DS2)`

The linear system(`A` and **`b`**) cen be set outside, randomly and from another linear system.

## SparseSolver

It implements algorithms to solve the equation `A`**`x`**`=`**`b`** for a _sparse matrix_ `A` of type CSRMatrix<T>.

### Usage

- To solve a linear system using the **Cholesky Solver**:

  ```cpp
  //Define the matrix
  sparseSolver<float> solver1(CSRm2, value_b, result, 10);
  residual = solver1.solve_cholesky();
  ```

  

- To solve a linear system using the **CG Solver**:

  ```cpp
  //Define the matrix
  sparseSolver<float> solver1(CSRm2, value_b, result, 6);
  residual = solver1.solve_CG(input_value);
  ```

  

### Methods

`void cholesky_decomp()`

`void solve_cholesky()`

`void solve_CG(float & tolerance)`


## Test Framework

Each function has a test, and run_test() can run all the test.

### TestRunner

The class will give a summary about how many tests have passed and how many have failed. Also the class will print the testing information(the status of each test) to users in the terminal.


### test

We test all the solvers(dense and sparse included) in this file. 

e.g. the unit test of the Jacobi algorithm in `test.h`:

```cpp
bool test_Jacobi() {
    /* the detail of test */
    bool flag = (result == expected);
    return flag;
}
```

In `run_tests()` function, call the `TestRunner::test()` method with the memory address of the unit test function and a title which will be displayed to the user in terminal. In that function we also test other important function(`matVecMult`, `calResidual`,.etc)

```cpp
test_runner_matrix.test(&test_matVecMult, "matrix-vector multiplication.");
test_runner_dense_solver.test(&test_Jacobi, "Jacobi method(dense).");
test_runner_sparse_solver.test(&test_Cholesky_solver, "Cholesky Solver method(sparse).");
test_runner_utils.test(&test_residual_calculation, "Calculate the residual.");
```

## Performance Comparasion

Since we have implemented multiple solutions for both dense and sparse matrices respectively, it is a good idea to compare the performance(running time) of different solvers.

Unfortunately, due to time limitation we have only completed a performance comparison and analysis of solvers for the dense matrix. 

### Usage

e.g. get the running time in different size of matrices using the Jacobi algorithm:

```cpp
// Randomly generate 10 diagonally dominant matrices(size=100,200...1000)
auto *denseSolver = new DenseSolver<double>(size);

// Get running time of the solver
auto t_begin = std::chrono::high_resolution_clock::now();
denseSolver->Jacobi(x, tol, it_max);
auto t_end = std::chrono::high_resolution_clock::now();

// Write to the txt file
if (myfile.is_open()) {
    myfile << size << "," << time << endl;
}
```

Then we will get the running time of different solvers in different size, the data is stored in txt file.

Next, the data in txt files will be used to generate a comparison graph which will be analyzed in our report.