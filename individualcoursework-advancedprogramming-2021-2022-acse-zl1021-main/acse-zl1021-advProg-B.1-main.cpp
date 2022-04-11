/**
 * GitHub ID: acse-zl1021
 * Operating System (including version): macOS Monterey(version 12.1)
 * Compiler (including version):
    * g++ 4.2.1
    * clang 13.0.0(clang-1300.0.29.30)
    * x86_64-apple-darwin21.2.0
 * IDE (including version): CLion 2021.3.2
 */
#include <iostream>
#include <vector>

using namespace std;

//over-write the operator using template
template<typename T>
vector<T> operator+(const vector<T> &a, const vector<T> &b) {
    vector<T> res;
    for (uint i = 0; i < a.size(); i++) {
        res.push_back(a[i] + b[i]);
    }
    return res;
}

int main() {
    int rows = 3;
    int cols = 3;

    vector<vector<double>> num1;
    vector<double> v1;
    vector<vector<double>> num2;
    vector<double> v2;

    num1.clear();
    num2.clear();

    double tmp = 0;

    cout<<"Please input the 1st vector: "<<endl;
    for (int i = 0; i < rows; ++i) {
        v1.clear();
        for (int j = 0; j < cols; ++j) {
            cin>>tmp;
            v1.push_back(tmp);
        }
        num1.push_back(v1);
    }

    cout<<"Please input the 2nd vector: "<<endl;
    for (int i = 0; i < rows; ++i) {
        v2.clear();
        for (int j = 0; j < cols; ++j) {
            cin>>tmp;
            v2.push_back(tmp);
        }
        num2.push_back(v2);
    }

//    vector<vector<double>> nums1 = {
//            {1,    2,  3},
//            {4.75, 5,  6},
//            {1,    -2, 3}};
//    vector<vector<double>> nums2 = {
//            {1, 0, 0},
//            {0, 1, 0},
//            {0, 0, 1}};
//    if (num1.size() != num2.size()) {
//        cerr << "Size Error!" << endl;
//    }
    vector<vector<double>> nums3 = num1 + num2;

    cout<<"VECTOR1 + VECTOR2 = "<<endl;
    for (uint i = 0; i < nums3.size(); i++) {
        for (uint j = 0; j < nums3[0].size(); j++) {
            cout << nums3[i][j] << "\t";
        }
    }
    return 0;
}
