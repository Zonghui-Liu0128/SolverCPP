#pragma once

#include <iostream>
#include <string>

using namespace std;

class TestRunner {
private:
    // keeps track of every time test method is run & the outcome
    int testsFailed = 0;
    int testsSucceeded = 0;
    string title = "";

    // Test ended, print the summary
    void completeRun();

public:
    TestRunner(string test_name);

    ~TestRunner();

    // runs the function pointed at by test_ptr
    void test(bool (*test_ptr)(), string title);

    static bool assertBelowTolerance(double val, double tol);

    static void testError(std::string message);
};
