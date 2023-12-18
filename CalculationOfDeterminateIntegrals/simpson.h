#pragma once

double simpson(double a, double b, int n, double (*F)(double));
double simpson_wrapped(double a, double b, int n, double (*F)(double));