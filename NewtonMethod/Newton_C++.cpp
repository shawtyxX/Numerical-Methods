#include <iostream>
#include <iomanip>
#include <cmath>
#include "Gauss.h"

using namespace std;

double E1 = 1e-9;
double E2 = 1e-9;

double f1(double x1, double x2)
{
    return pow(x1, 3) + pow(x2, 3) - 6 * x1 + 3;
}

double f2(double x1, double x2)
{
    return pow(x1, 3) - pow(x2, 3) + 6 * x2 + 2;
}
double df1dx1(double x1, double x2)
{
    return 3 * pow(x1, 2) - 6;
}

double df1dx2(double x1, double x2)
{
    return 3 * pow(x2, 2);
}

double df2dx1(double x1, double x2)
{
    return 3 * pow(x1, 2);
}

double df2dx2(double x1, double x2)
{
    return -3 * pow(x2, 2) + 6;
}

void Jcobian(double x1, double x2, double **J)
{
	J[0][0] = df1dx1(x1, x2);
	J[0][1] = df1dx2(x1, x2);
	J[1][0] = df2dx1(x1, x2);
	J[1][1] = df2dx2(x1, x2);
}
void Jcobian_M(double x1, double x2, double **J, double M)
{
	J[0][0] = (f1(x1 + M * x1, x2) - f1(x1, x2)) / (M * x1);
	J[0][1] = (f1(x1, x2 + M * x2) - f1(x1, x2)) / (M * x2);
	J[1][0] = (f2(x1 + M * x1, x2) - f2(x1, x2)) / (M * x1);
	J[1][1] = (f2(x1, x2 + M * x2) - f2(x1, x2)) / (M * x2);
}

void Newton(double* dX,double x1, double x2, int NIT, double E1, double E2)
{	
	double* F = new double[2];
	double** J = new double* [2];
	for (int i = 0; i < 2; i++)
	{
		J[i] = new double[2];
	}
	int k = 1;
	cout << left << setw(4) << "k" << setw(20) << "d1" << setw(20) << "d2" << setw(20) << "x1" << setw(20) << "x2" << endl;
	double x1k, x2k;
	double d1, d2;
	double tmp;
	do {
		F[0] = -f1(x1, x2);
		F[1] = -f2(x1, x2);
		Jcobian(x1, x2, J);

		Gauss(J, F,dX,2);
		x1k = x1 + dX[0];
		x2k = x2 + dX[1];
		d1 = abs(f1(x1k, x2k));
		tmp = abs(f2(x1k, x2k));
		if (tmp > d1)
		{
			d1 = tmp;
		}
		d2 = abs(x1k - x1) / (x1k >= 1 ? x1k : 1);
		tmp = abs(x2k - x2) / (x2k >= 1 ? x2k : 1);
		if (tmp > d2)
		{
			d2 = tmp;
		}
		x1 = x1k;
		x2 = x2k;
		cout << left << setw(4) << k << setw(20) << d1 << setw(20) << d2 << setw(20) << x1 << setw(20) << x2 << endl;
		if (k >= NIT)
		{
			cout << "IER=2\n";
			exit(2);
		}
		k++;
	} while (d1 > E1 || d2 > E2);
	dX[0] = x1;
	dX[1] = x2;
	delete[] F;
	for (int i = 0; i < 2; i++) {
		delete J[i];
	}
}

void Newton_M( double* dX, double x1, double x2, int NIT, double E1, double E2, double M)
{	
	double* F = new double[2];
	double** J = new double* [2];
	for (int i = 0; i < 2; i++)
	{
		J[i] = new double[2];
	}
	int k = 1;
	cout << left << setw(4) << "k" << setw(20) << "d1" << setw(20) << "d2" << setw(20) << "x1" << setw(20) << "x2" << endl;
	double x1k, x2k;
	double d1, d2;
	double tmp;
	do {
		F[0] = -f1(x1, x2);
		F[1] = -f2(x1, x2);
		Jcobian_M(x1, x2, J, M);

		Gauss(J, F, dX,2);
		x1k = x1 + dX[0];
		x2k = x2 + dX[1];
		d1 = abs(f1(x1k, x2k));
		tmp = abs(f2(x1k, x2k));
		if (tmp > d1)
		{
			d1 = tmp;
		}
		d2 = abs(x1k - x1) / (x1k >= 1 ? x1k : 1);
		tmp = abs(x2k - x2) / (x2k >= 1 ? x2k : 1);
		if (tmp > d2)
		{
			d2 = tmp;
		}
		x1 = x1k;
		x2 = x2k;
		cout << left << setw(4) << k << setw(20) << d1 << setw(20) << d2 << setw(20) << x1 << setw(20) << x2 << endl;
		if (k >= NIT)
		{
			cout << "IER=2\n";
			exit(2);
		}
		k++;
	} while (d1 > E1 || d2 > E2);
	dX[0] = x1;
	dX[1] = x2;
	delete[] F;
	for (int i = 0; i < 2; i++) {
		delete J[i];
	}
}

int main()
{
	int NIT = 100;
	int n = 2;
	double* dX = new double[2];
	
	double* x = new double[2];
	double* x1 = new double[2];
	double* x3 = new double[2];
	double* x4 = new double[2];
	double* x5 = new double[2];
	cout << "\nаналитический метод для точки (1, 1): \n";
	Newton( dX, 1, 1, NIT, E1, E2);
	cout << "\nчерез матрицу якоби для точки (1, 1) и M=0.01: \n";
	Newton_M(dX, 1, 1, NIT, E1, E2, 0.01);
	cout << "\nчерез матрицу якоби для точки (2, 1.5) и M=0.05: \n";
	Newton_M(dX, 2, 1.5, NIT, E1, E2, 0.05);
	cout << "\nчерез матрицу якоби для точки (-3, -1,5) и M=0.1: \n";
	Newton_M(dX, -3, -1.5, NIT, E1, E2, 0.1);
	delete[] dX;
	delete[] x;
	delete[] x1;
	delete[] x3;
	delete[] x4;
	delete[] x5;

	return 0;
}



