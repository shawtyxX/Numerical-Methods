#include <iostream>
#include <cmath>
#include "trapezoid.h"
#include "simpson.h"
#include "simpson2d.h"


//TUNABLES
#define A 0         // START POINT
#define B M_PI/2    // END POINT
#define N 50        // FRACTION

#define A2d 0       // START POINT FOR X IN SIMPSON 2D
#define B2d M_PI/2  // END POINT FOR X IN SIMPSON 2D
#define C2d 0       // START POINT FOR Y IN SIMPSON 2D
#define D2d M_PI/4  // END POINT FOR Y IN SIMPSON 2D
#define N2d 100      // FRACTION FOR X IN SIMPSON 2D
#define M2d 100     // FRACTION FOR Y IN SIMPSON 2D

using namespace std;


double func(double x){
    return (pow(4-pow(sin(x),2),0.5))/2;
}

double func2(double x, double y){
    return sin(x + y);
}

int main(){
    cout << "Bubna lisaya" << std::endl;
    cout << "FUNCTION: F(x) = ((4-(sinx)^2)^0.5)/2" << endl << endl;
    double trapezoid_result = trapezoid_wrapped(A, B, N, func);
    cout << "TRAPEZOID: " << trapezoid_result << endl;

    double simpson_result = simpson_wrapped(A, B, N, func);
    cout << "SIMPSON: " << simpson_result << endl;

    double simpson2d_result = simpson2d(A2d, B2d, C2d, D2d, N2d, M2d, func2);
    cout << "SIMPSON 2D: " << simpson2d_result << endl;

    return 0;
}