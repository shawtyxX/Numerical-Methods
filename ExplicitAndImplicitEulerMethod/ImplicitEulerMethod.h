#pragma once //используется директива #pragma once, которая указывает препроцессору, что файл подключается только один раз.

#include <iostream>
#include <vector>
#include "NewtonsMethod.h"

using namespace std;

vector<vector<double>> FillMatrix(vector<double> );
vector<double> FillVector(vector<double>);
vector<double> CalculateLocalError(vector<double> , vector<double> , vector<double> , double , double );
bool CheckValidityForError(vector<double> , double );
double CalculateNextStepOfIntegration(vector<double>, double , double );
vector<double> SolveImplicitEulerMethod(vector<vector<double>> , vector<double> ,vector<double> , vector<double> ,double);
