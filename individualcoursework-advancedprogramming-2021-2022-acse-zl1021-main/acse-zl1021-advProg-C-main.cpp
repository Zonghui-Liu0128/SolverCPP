/**
 * GitHub ID: acse-zl1021
 * Operating System (including version): macOS Monterey(version 12.1)
 * Compiler (including version):
    * g++ 4.2.1
    * clang 13.0.0(clang-1300.0.29.30)
    * x86_64-apple-darwin21.2.0
 * IDE (including version): CLion 2021.3.2
 */

#include <algorithm>
#include <iostream>
#include <cstdlib>

using namespace std;


double VolumeOfTetra(double a[], double b[], double c[], double d[]) {
    double temp1[3], temp2[3], temp3[3], temp4[3];

    //t1 = a - d
    for (int i = 0; i < 3; i++) {
        temp1[i] = a[i] - d[i];
    }

    //t2 = b - d
    for (int i = 0; i < 3; i++) {
        temp2[i] = b[i] - d[i];
    }

    //t3 = c - d
    for (int i = 0; i < 3; i++) {
        temp3[i] = c[i] - d[i];
    }

    //compute (b - d) * (c - d)
    temp4[0] = temp2[1] * temp3[2] - temp3[1] * temp2[2];
    temp4[1] = temp3[0] * temp2[2] - temp2[0] * temp3[2];
    temp4[2] = temp2[0] * temp3[1] - temp3[0] * temp2[1];

    double tmp = 0;
    for (int i = 0; i < 3; i++) {
        tmp += temp1[i] * temp4[i];
    }
    return abs(tmp) / 6.0;
}



double VolumeOfHexa(double h0[], double h1[], double h2[], double h3[], double h4[], double h5[], double h6[], double h7[]) {
    double sub1 = VolumeOfTetra(h0, h1, h3, h4);
    double sub2 = VolumeOfTetra(h4, h1, h3, h5);
    double sub3 = VolumeOfTetra(h4, h5, h3, h7);
    double sub4 = VolumeOfTetra(h1, h2, h3, h6);
    double sub5 = VolumeOfTetra(h3, h1, h6, h5);
    double sub6 = VolumeOfTetra(h5, h6, h3, h7);

    //adding the volume of sub-tetrahedra
    return sub1 + sub2 + sub3 + sub4 + sub5 + sub6;
}


int main() {
    cout<<"Compute the volume of a tetrahedron:"<<endl;
    double p1[3] = {1, 0, 1}, p2[3] = {2, 2, 3}, p3[3] = {3, 9, 7}, p4[3] = {16, 13, 11};
    cout << VolumeOfTetra(p1, p2, p3, p4) << endl;

    cout<<endl<<"Creating 1,000 random hexahedra and compute their volume..."<<endl;
    double h0[3], h1[3], h2[3], h3[3], h4[3], h5[3], h6[3], h7[3];
    double sum = 0;
    for (int i = 0; i < 10000; i++) {
        //initialize the 8 points randomly
        srand((int) time(0));
        for (int i = 0; i < 3; i++) {
            h0[i] = rand() % 100;
            h1[i] = rand() % 100;
            h2[i] = rand() % 100;
            h3[i] = rand() % 100;
            h4[i] = rand() % 100;
            h5[i] = rand() % 100;
            h6[i] = rand() % 100;
            h7[i] = rand() % 100;
        }
        double v_hexa = VolumeOfHexa(h0, h1, h2, h3, h4, h5, h6, h7);
        sum += v_hexa;
    }
    cout<<endl<<"The sum volume of the 1000 hexahedra"<<endl;
    cout << sum << endl;
}
