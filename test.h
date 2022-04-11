#include <iostream>
#include <memory>
#include <fstream>
#include <cmath>
#include <vector>
#include "Matrix.h"
#include "CSRMatrix.h"
#include "DenseSolver.h"
#include "sparseSolver.h"
#include "TestRunner.h"
#include "utilities.h"

/** Test all important function we will use in solvers **/
bool test_matVecMult() {
    int size = 4;

    shared_ptr<double[]> values(new double[size * size]{1, 2, 3, 4, -4, -3, -2, -1, 0, 1, 2, 0, 1, 2, 3, 4});

    Matrix<double> m = Matrix<double>(size, size, values);

    double *v = new double[size]{1, 2, -3, 2};
    double *result = new double[size]{0};

    m.matVecMult_pointer(m, v, result);

    double *expected = new double[size]{4, -6, -4, 4};

    for (int i = 0; i < size; i++) {
        if (result[i] != expected[i]) {
            cerr << "Result doesn't match expected values!" << endl;
            return false;
        }
    }
    return true;
}

bool test_check_dimensions_matching() {
    Matrix<int> m = Matrix<int>(3, 3, true);
    vector<int> v(3, 0);
    try {
        checkDimensions(m, v);
        return true;
    }
    catch (const std::exception &e) {
        return false;
    }
}

bool test_residual_calculation() {
    int size = 4;

    shared_ptr<double[]> values(
            new double[size * size]{1., 2., 3., 4., -4., -3., -2., -1., 0., 1., 2., 0., 1., 2., 3., 4.});

    Matrix<double> A = Matrix<double>(size, size, values);

    vector<double> b{4., -6., -4., 4.};
    vector<double> x{1., 1., -2., 2.};

    double residual = calResidual_vector(A, x, b);

    double expected = 2.;

    return residual == expected;
}

/** Test all the dense solver **/
bool test_Gauss_solver() {
    int size = 3;
    double tol = 1.0e-6;
    shared_ptr<double[]> dense_values(
            new double[size * size]{2., 1., -1., -3., -2., 4., 5., -3., -1.});
    Matrix<double> A = Matrix<double>(size, size, dense_values);
    vector<double> b = {1., 5., -4.};

    DenseSolver<double> denseSolver = DenseSolver<double>(A, b);

    auto t_begin = std::chrono::high_resolution_clock::now();
    denseSolver.GaussSolver();
    auto t_end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(t_end - t_begin).count();
    cerr << "Time taken for Gauss Solver: " << time << " s" << endl;

    shared_ptr<double[]> dense_values_intial(
            new double[size * size]{2., 1., -1., -3., -2., 4., 5., -3., -1.});
    Matrix<double> A_initial = Matrix<double>(size, size, dense_values_intial);
    double residual = 0;
    residual = calResidual_vector(A_initial, denseSolver.b, b);
    cerr << "The residual of Gauss Solver: " << residual << endl << endl;
    return TestRunner::assertBelowTolerance(residual, tol);
}

bool test_Jacobi() {
    int size = 3;
    double tol = 1e-6;
    int it_max = 1000;
    vector<double> b = {5, 15, 8};
    shared_ptr<double[]> dense_values(
            new double[size * size]{2., 1., 1., 3., 5., 2., 2., 1., 4.});
    Matrix<double> A = Matrix<double>(size, size, dense_values);
    DenseSolver<double> denseSolver = DenseSolver(A, b);
    vector<double> x(b.size(), 0);
    auto t_begin = std::chrono::high_resolution_clock::now();
    denseSolver.Jacobi(x, tol, it_max);
    auto t_end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(t_end - t_begin).count();
    cerr << "Time taken for Jacobi: " << time << " s" << endl;
    double residual = 0;
    residual = calResidual_vector(A, x, b);
    cerr << "The residual of Jacobi: " << residual << endl << endl;
    return TestRunner::assertBelowTolerance(residual, tol);
}

bool test_GaussSeidel() {
    int size = 3;
    double tol = 1e-6;
    int it_max = 1000;
    vector<double> b = {5, 15, 8};
    shared_ptr<double[]> dense_values(
            new double[size * size]{2., 1., 1., 3., 5., 2., 2., 1., 4.});
    Matrix<double> A = Matrix<double>(size, size, dense_values);
    DenseSolver<double> denseSolver = DenseSolver(A, b);
    vector<double> x(b.size(), 0);
    auto t_begin = std::chrono::high_resolution_clock::now();
    denseSolver.GaussSeidel(x, tol, it_max);
    auto t_end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(t_end - t_begin).count();
    cerr << "Time taken for Gauss-Seidel: " << time << " s" << endl;
    double residual = 0;
    residual = calResidual_vector(A, x, b);
    cerr << "The residual of Gauss-Seidel: " << residual << endl << endl;
    return TestRunner::assertBelowTolerance(residual, tol);
}

bool test_LUSolver() {
    int size = 3;
    double tol = 1.0e-6;
    vector<double> x(size, 0);

    shared_ptr<double[]> init_dense_values(
            new double[size * size]{2., 1., 1., 3., 5., 2., 2., 1., 4.});

    Matrix<double> dense_mat = Matrix<double>(size, size, init_dense_values);
    vector<double> b = {5, 15, 8};

    DenseSolver<double> dense_solver = DenseSolver<double>(dense_mat, b);

    Matrix<double> LU(size, size, true);

    auto t_begin = std::chrono::high_resolution_clock::now();
    vector<int> piv = dense_solver.LUDecomposition(LU);
    dense_solver.LUSolver(LU, piv, x, b);
    auto t_end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration<double>(t_end - t_begin).count();
    cerr << "Time taken for LU Solver with partial pivoting: " << duration << " s " << endl;

    double residual = calResidual_vector(dense_mat, x, b);
    cerr << "The residual of LU Solver: " << residual << endl << endl;
    return TestRunner::assertBelowTolerance(residual, tol);
}

/** Test all the sparse solver **/
bool test_Cholesky_solver() {
    ifstream fin;
    int buf_int = 0;
    float buf_float = 0;
    //get the value from the file
    fin.open("../data/Cholesky_data.txt", ios::in);
    if (!fin) {
        cout << "wrong" << endl;
    }
    unique_ptr<int[]> init_row_position1(new int[11]{0});
    for (int i = 0; i < 11; i++) {
        fin >> buf_int;
        init_row_position1[i] = buf_int;
    }
    unique_ptr<int[]> init_col_index1(new int[20]{0});
    for (int i = 0; i < 20; i++) {
        fin >> buf_int;
        init_col_index1[i] = buf_int;
    }
    unique_ptr<float[]> init_sparse_values1(new float[20]{0});
    for (int i = 0; i < 20; i++) {
        fin >> buf_float;
        init_sparse_values1[i] = buf_float;
    }
    unique_ptr<float[]> value_b(new float[10]{0});
    for (int i = 0; i < 10; i++) {
        fin >> buf_float;
        value_b[i] = buf_float;
    }
    unique_ptr<float[]> result(new float[10]{0});
    for (int i = 0; i < 10; i++) {
        fin >> buf_float;
        result[i] = buf_float;
    }

    //initialize the value of matrix and solver
    CSRMatrix<float> CSRm2(10, 10, 20, init_sparse_values1, init_row_position1, init_col_index1);
    sparseSolver<float> solver1(CSRm2, value_b, result, 10);

//    cout << "The matrix is:" << endl;
//    CSRm2.printMatrix();
//    cout << "The result is:" << endl;
//
//    for(int i=0;i<10;i++){
//        cout << result[i] << " ";
//    }
//    cout << endl;
    float residual = 0;
    auto t_begin = std::chrono::high_resolution_clock::now();
    residual = solver1.solve_cholesky();
    auto t_end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(t_end - t_begin).count();
    cerr << "Time taken for 10*10 sparse matrix(Cholesky): " << time << " s" << endl;
    cerr << "The residual is " << residual << endl << endl;
    return TestRunner::assertBelowTolerance(residual, 0.1);
}

bool test_CG_solver() {
    ifstream fin;
    int buf_int = 0;
    float buf_float = 0;
    //get the value from the file
    fin.open("../data/CG_data.txt", ios::in);
    if (!fin) {
        cout << "wrong" << endl;
    }
    unique_ptr<int[]> init_row_position1(new int[7]{0});
    for (int i = 0; i < 7; i++) {
        fin >> buf_int;
        init_row_position1[i] = buf_int;
    }
    unique_ptr<int[]> init_col_index1(new int[36]{0});
    for (int i = 0; i < 36; i++) {
        fin >> buf_int;
        init_col_index1[i] = buf_int;
    }
    unique_ptr<float[]> init_sparse_values1(new float[36]{0});
    for (int i = 0; i < 36; i++) {
        fin >> buf_float;
        init_sparse_values1[i] = buf_float;
    }
    unique_ptr<float[]> value_b(new float[6]{0});
    for (int i = 0; i < 6; i++) {
        fin >> buf_float;
        value_b[i] = buf_float;
    }
    unique_ptr<float[]> result(new float[6]{0});
    for (int i = 0; i < 6; i++) {
        fin >> buf_float;
        result[i] = buf_float;
    }
    fin.close();

    //initialize the value of matrix and solver
    CSRMatrix<float> CSRm2(6, 6, 36, init_sparse_values1, init_row_position1, init_col_index1);
    sparseSolver<float> solver1(CSRm2, value_b, result, 6);
    float input_value = 0.01;

//    cout << "The matrix is:" << endl;
//    CSRm2.printMatrix();


    float residual = 0;
    auto t_begin = std::chrono::high_resolution_clock::now();
    residual = solver1.solve_CG(input_value);
    auto t_end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(t_end - t_begin).count();
    cerr << "Time taken for 10*10 sparse matrix(CG-sparse): " << time << " s" << endl;
    cerr << "The residual is " << residual << endl;
    return TestRunner::assertBelowTolerance(residual, 1);
}

/** Test script for each solver **/
void run_GaussSolver_test() {
    TestRunner test_runner_dense_solver = TestRunner("DenseSolver");
    test_runner_dense_solver.test(&test_Gauss_solver, "Gauss Solver method(dense).");
}

void run_LUSolver_test() {
    TestRunner test_runner_dense_solver = TestRunner("DenseSolver");
    test_runner_dense_solver.test(&test_LUSolver, "LU Solver with partial pivoting method(dense).");
}

void run_Jacobi_test() {
    TestRunner test_runner_dense_solver = TestRunner("DenseSolver");
    test_runner_dense_solver.test(&test_Jacobi, "Jacobi method(dense).");
}

void run_GaussSeidel_test() {
    TestRunner test_runner_dense_solver = TestRunner("DenseSolver");
    test_runner_dense_solver.test(&test_GaussSeidel, "Gauss-Seidel method(dense).");
}

void run_cholesky_tests() {
    TestRunner test_runner_sparse_solver = TestRunner("CSRMatrix");
    test_runner_sparse_solver.test(&test_Cholesky_solver, "Cholesky Solver method.");
}

void run_CG_tests() {
    TestRunner test_runner_sparse_solver = TestRunner("CSRMatrix");
    test_runner_sparse_solver.test(&test_CG_solver, "Conjugate Gradient Solver method.");
}

/** Test script for all functions **/
void run_tests() {
    // MATRIX
    TestRunner test_runner_matrix = TestRunner("Matrix");
    test_runner_matrix.test(&test_matVecMult, "matrix-vector multiplication.");

    // DENSE SOLVER
    TestRunner test_runner_dense_solver = TestRunner("DenseSolver");
    test_runner_dense_solver.test(&test_Gauss_solver, "Gauss Solver method(dense).");
    test_runner_dense_solver.test(&test_LUSolver, "LU Solver with partial pivoting method(dense).");
    test_runner_dense_solver.test(&test_Jacobi, "Jacobi method(dense).");
    test_runner_dense_solver.test(&test_GaussSeidel, "Gauss-Seidel method(dense).");

    // SPARSE SOLVER
    TestRunner test_runner_sparse_solver = TestRunner("SparseSolver");
    test_runner_sparse_solver.test(&test_Cholesky_solver, "Cholesky Solver method(sparse).");
    test_runner_sparse_solver.test(&test_CG_solver, "Conjugate Gradient Solver method(sparse).");

    // UTILITIES
    TestRunner test_runner_utils = TestRunner("Utilities");
    test_runner_utils.test(&test_check_dimensions_matching, "checkDimensions for matching matrices.");
    test_runner_utils.test(&test_residual_calculation, "Calculate the residual.");
}