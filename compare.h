#pragma once

#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>
#include <fstream>
#include <string>
#include "Matrix.h"
#include "DenseSolver.h"
#include "TestRunner.h"
#include "utilities.h"

void Gauss_performance(int minsize, int maxsize) {
    string filename;
    double tol = 1.0e-6;
    filename = "../data/GaussSolver_dense_100-1000.txt";
    ofstream myFile;
    myFile.open(filename);

    int size = minsize;
    while (size <= maxsize) {
        vector<double> x(size, 0);
        auto *denseSolver = new DenseSolver<double>(size);
        auto t_begin = std::chrono::high_resolution_clock::now();
        denseSolver->GaussSolver();
        auto t_end = std::chrono::high_resolution_clock::now();

        auto time = std::chrono::duration<double>(t_end - t_begin).count();
        cout << "Gauss Solver for size " << size << ", time = " << time << " s " << endl;

        // Write to file
        if (myFile.is_open()) {
            myFile << size << "," << time << endl;
        } else
            cout << "Unable to open file";
        size += 100;
    }
    myFile.close();
    cout << endl;
}

void LUSolver_performance(int minsize, int maxsize) {
    string filename;
    double tol = 1.0e-6;
    filename = "../data/LUSolver_dense_100-1000.txt";
    ofstream myfile;
    myfile.open(filename);
    int size = minsize;
    while (size <= maxsize) {
        vector<double> x(size, 0);
        vector<double> b(size, 0);
        auto *denseSolver = new DenseSolver<double>(size);

        auto *LU = new Matrix<double>(size, size, true);

        auto t_begin = std::chrono::high_resolution_clock::now();
        vector<int> piv = denseSolver->LUDecomposition(*LU);
        denseSolver->LUSolver(*LU, piv, x, b);
        auto t_end = std::chrono::high_resolution_clock::now();

        auto time = std::chrono::duration<double>(t_end - t_begin).count();
        cout << "LUSolver for size " << size << ", time = " << time << " s " << endl;

        // Write to file
        if (myfile.is_open()) {
            myfile << size << "," << time << endl;
        } else
            std::cout << "Unable to open file";

        if (calResidual_vector(denseSolver->A, x, b) > tol) {
            throw "LU residual is above 1e-6";
        }
        // Delete objects to save memory usage
        delete denseSolver;
        delete LU;
        size += 100;
    }
    myfile.close();
    cout << endl;
}

void Jacobi_performance(int minsize, int maxsize) {
    string filename;
    double tol = 1.0e-6;
    int it_max = 1000;
    filename = "../data/Jacobi_dense_100-1000.txt";
    ofstream myfile;
    myfile.open(filename);
    int size = minsize;
    while (size <= maxsize) {
        vector<double> x(size, 0);
        vector<double> b(size, 0);
        auto *denseSolver = new DenseSolver<double>(size);

        auto t_begin = std::chrono::high_resolution_clock::now();
        denseSolver->Jacobi(x, tol, it_max);
        auto t_end = std::chrono::high_resolution_clock::now();

        auto time = std::chrono::duration<double>(t_end - t_begin).count();
        cout << "Jacobi for size " << size << ", time = " << time << " s " << endl;

        // Write to file
        if (myfile.is_open()) {
            myfile << size << "," << time << endl;
        } else
            cout << "Unable to open file";
        // Delete objects to save memory usage
        delete denseSolver;
        size += 100;
    }
    myfile.close();
    cout << endl;
}

void GS_performance(int minsize, int maxsize) {
    string filename;
    double tol = 1.0e-6;
    int it_max = 1000;
    filename = "../data/GaussSeidel_dense_100-1000.txt";
    ofstream myfile;
    myfile.open(filename);
    int size = minsize;
    while (size <= maxsize) {
        vector<double> x(size, 0);
        vector<double> b(size, 0);
        auto *denseSolver = new DenseSolver<double>(size);

        auto t_begin = std::chrono::high_resolution_clock::now();
        denseSolver->Jacobi(x, tol, it_max);
        auto t_end = std::chrono::high_resolution_clock::now();

        auto time = std::chrono::duration<double>(t_end - t_begin).count();
        cout << "Gauss-Seidel for size " << size << ", time = " << time << " s " << endl;

        // Write to file
        if (myfile.is_open()) {
            myfile << size << "," << time << endl;
        } else
            cout << "Unable to open file";
        // Delete objects to save memory usage
        delete denseSolver;
        size += 100;
    }
    myfile.close();
    cout << endl;
}

void run_performance() {
    int minsize = 100;
    int maxsize = 1000;

    // the performance of Gauss Solver(dense)
    Gauss_performance(minsize, maxsize);

    // the performance of LUSolver(dense)
    LUSolver_performance(minsize, maxsize);

    // the performance of Jacobi(dense)
    Jacobi_performance(minsize, maxsize);

    // the performance of Gauss-Seidel(dense)
    GS_performance(minsize, maxsize);
}
