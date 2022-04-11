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


int main() {
    vector<string> vec{"Tokyo", "Berlin", "Rio", "Denver", "Helsinki", "Nairobi"};

    sort(vec.begin(), vec.end());

    cout << "Sorting...."<<endl;
    for (vector<string>::iterator iter = vec.begin(); iter != vec.end() ; ++iter) {
        cout<<*iter<<"\t";
    }

    return 0;
}