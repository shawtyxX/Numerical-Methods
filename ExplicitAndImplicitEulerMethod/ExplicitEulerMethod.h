#pragma once //используется директива #pragma once, которая указывает препроцессору, что файл подключается только один раз.

#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

vector<double> SolveExplicitEulerMethod(vector<double>, vector<double>, double, double);
double u0dt(vector<double>, double);
double u1dt(vector<double>, double, double);
vector<double> CalculateDifferentialEquations(vector<double>, double, double);
double CalculateIntegrationStep(vector<double>, double, double);
vector<double> CalculateSolutionOfEquations(vector<double>, vector<double>, double);
