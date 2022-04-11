#include "TestRunner.h"
#include <iostream>
#include <string>
#include <cmath>

using namespace std;
// Colourful definition
#define RED "\033[31m"
#define GREEN "\033[32m"
#define RESET "\033[0m"
#define YELLOW "\033[1;33m"

TestRunner::TestRunner(string test_name) : title(test_name) {}

TestRunner::~TestRunner() {
    // Print the test summary
    this->completeRun();
}

void TestRunner::test(bool (*test_ptr)(), string title) {
    int total = this->testsFailed + this->testsSucceeded;

    // title displayed before any terminal outputs
    cout << endl << YELLOW << "Test " << total + 1 << ": " << title << RESET << endl;


    bool test_result = test_ptr();

    if (test_result) {
        this->testsSucceeded += 1;
        cout << GREEN << "Passed" << RESET << endl;
    } else {
        this->testsFailed += 1;
        cout << RED << "Failed" << RESET << endl;
    }
}

void TestRunner::completeRun() {
    int total = this->testsFailed + this->testsSucceeded;
    cout << endl;
    // test successfully
    if (this->testsFailed == 0) {
        cout << GREEN << this->title << ": " << total << "/" << total << " tests passed." << RESET << endl;
    } else {
        cout << RED << this->title << ": "
             << " " << this->testsFailed << "/" << total << " tests failed." << RESET << endl
             << endl;
    }
}

void TestRunner::testError(std::string message) {
    cerr << RED << message << RESET << endl;
}


bool TestRunner::assertBelowTolerance(double val, double tol) {

    if (isnan(val)) {
        TestRunner::testError("Something went wrong");
        return false;
    }

    if (val > tol) {
        TestRunner::testError("Something went wrong(value is above tolerance)");
        return false;
    }
    return true;
}
