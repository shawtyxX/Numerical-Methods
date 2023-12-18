#include "ExplicitEulerMethod.h"
#include "ImplicitEulerMethod.h"
#include "cmath"

int main()
{
	double firstLocalError = pow(10, -3), secondLocalError = pow(10, -5);
	vector<double> parameterBoundariesT = { 0,1 }; // 0 <= t <= 1
	double parameterOmega = 30; // 25 <= omega <= 48
	vector<double> valueUAtZeroExp = { 0,-0.412 }; // u1 and u2
	vector<double> rootsOfFirstEquation = SolveExplicitEulerMethod(valueUAtZeroExp, parameterBoundariesT, parameterOmega, firstLocalError);
	cout << "Solving the equation using explicit Euler method: ";
	for (int i = 0; i < rootsOfFirstEquation.size(); i++)
		cout  << rootsOfFirstEquation[i] << " ";
	vector<double> valueUAtZeroImpl = { 10,22,9 };
	vector<double> vectorOfLambda = {-10,-10,-10};
	vector<vector<double>> A = FillMatrix(vectorOfLambda);
	vector<double> b = FillVector(vectorOfLambda);
	vector<double>rootsOfSecondEquation = SolveImplicitEulerMethod(A,b,valueUAtZeroImpl, parameterBoundariesT,firstLocalError);
	cout << "\n\nSolving the equation using implicit Euler method: ";
	for (int i = 0; i < rootsOfSecondEquation.size(); i++)
		cout << rootsOfSecondEquation[i] << " ";
	return 0;
}