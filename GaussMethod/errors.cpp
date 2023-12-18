#include "errors.h"
#include <cmath>

double* F_vector(double** a, double* b, double* solutions, int size) {
    double* f = new double[size];
    for (int i = 0; i < size; i++) {
        f[i] = -b[i];
        for (int j = 0; j < size; j++) {
            f[i] += a[i][j] * solutions[j];
        }
    }
    return f;
}
double norm(double* f, int n){
    double norma = std::abs(f[0]);
    for (int i = 0; i < n; i++) {
        norma = std::max(f[i], std::abs(norma));
    }
    return norma;
}