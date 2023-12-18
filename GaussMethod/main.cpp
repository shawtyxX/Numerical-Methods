#include <iostream>
#include <fstream>
#include "gauss.h"
#include "errors.h"

using namespace std;

int main(){
    ifstream eqfile = ifstream("equations.txt");
    int size;
    eqfile >> size;

    double** a = new double*[size];
    double*  b = new double [size];

    cout << "[INPUT]" << endl;
    for(int i = 0; i < size; i++){
        a[i] = new double[size];
        for(int j = 0; j < size; j++){
            eqfile >> a[i][j];
            cout << a[i][j] << "\t";
        }
        eqfile >> b[i];
        cout << "|\t" << b[i] << endl;
    }
    cout << endl;

    double* solutions = new double[3]{};
    gauss(a, b, solutions, 3);
    cout << "[SOLUTIONS]" << endl;
    for(int i = 0; i < size; i++){
        cout << "x" << i << " = " << solutions[i] << endl;
    }
    cout << endl;

    //Norm
    cout << "[DIFFERENT ERRORS STUFF]" << endl;
    double* f_vec = F_vector(a, b, solutions, size);
    double _norm = norm(f_vec, size);
    cout << "F_vector = (";
    for(int i = 0; i < size; i++){
        cout << f_vec[i] << (i == size-1 ? "": ", "); // prettify the output 
    }
    cout << ")" << endl;
    cout << "Norm = " << _norm << endl;

    for(int i = 0; i < size; i++){
        delete[] a[i];
    }
    delete[] a;
    delete[] b;
    delete[] solutions;
    delete[] f_vec;
    return 0;
}