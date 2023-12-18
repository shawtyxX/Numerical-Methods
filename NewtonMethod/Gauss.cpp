#include "Gauss.h"
#include <cmath>
#include <algorithm>
#include <iostream>
using namespace std;
void Gauss(double** A, double* B, double* X, int N)
{

    for (int k = 0; k < N; k++)
    {
        for (int i = k + 1; i < N; i++)
        {
            if (abs(A[i][k]) > abs(A[k][k]))
            {
                swap(A[i], A[k]);
                swap(B[i], B[k]);   
            }
        }
        float A_Main = A[k][k];
        if (A_Main == 0)
        {
            cout << "error\n";
            exit(0);
        }
        for (int i = k; i < N; i++)
        {
            A[k][i] /= A_Main;
        }
        B[k] /= A_Main;
        for (int i = k + 1; i < N; i++)
        {
            float s = A[i][k];
            for (int j = k; j < N; j++)
            {
                A[i][j] -= s * A[k][j];
            }
            B[i] -= s * B[k];
        }
    }
    for (int k = N - 1; k >= 0; k--) 
    {
        X[k] = B[k];
        for (int i = N - 1; i > k; i--)
        {
            X[k] -= A[k][i] * X[i];
        }
    }

}
