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


bool isPrime(int num) {
    bool flag = true;

    //"0" and "1" situation
    if (num == 0 || num == 1) {
        flag = false;
    } else {
        for (int i = 2; i <= num / 2; ++i) {
            if (num % i == 0) {
                flag = false;
                break;
            }
        }
    }
    if (flag)
        return true;
    else
        return false;
}

int main() {
    vector<int> vec;
    //counter for the number of the prime number
    int cnt = 0;
    //the number need to check
    int num = 0;
    while (cnt < 1000) {
        // the number is prime
        if (isPrime(num)) {
            vec.push_back(num);
            cnt++;
            num++;
        } else {
            num++;
        }
    }
    //print the elements out
    int cnt_l = 0;
    for (vector<int>::iterator iter = vec.begin(); iter < vec.end(); ++iter) {
        cnt_l ++;
        cout << *iter << ",";
        if (cnt_l == 25){
            cout<<endl;
            cnt_l = 0;
        }
    }
}