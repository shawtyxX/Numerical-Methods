#include "simpson.h"
#include <cmath>
#define eps 1e-6

// double simpson_middle_part(double ){
//     double result = 
// }

double simpson(double a, double b, int m, double (*F)(double)) {
    int n = 2*m;
    const double width = (b-a)/n;

    double simpson_integral = F(a);
    
    double sum1 = 0;

    // calculate sum 1
    for(int i = 1; i <= m; i++){
        double xx = a + width * (2*i-1);
        sum1 += 4 * F(xx);
    }
    
    //calculate sum 2
    double sum2 = 0;
    for(int i = 1; i < m; i++){
        double xx = a + width * 2*i;
        sum2 += 2 * F(xx);
    }

    simpson_integral += sum1;
    simpson_integral += sum2;

    //last part
    double xlast = a + width*n;
    simpson_integral += F(xlast);

    simpson_integral *= width/3;

    return simpson_integral;
}

double simpson_wrapped(double a, double b, int m, double (*F)(double)){
    double base_m = m;
    double m_2 = base_m*2;

    double simpson1 = simpson(a, b, base_m, F);
    double simpson2 = simpson(a, b, m_2, F);

    while(!(std::abs(simpson1-simpson2) <= eps)){
        base_m = m_2;
        m_2 *= 2;
        simpson1 = simpson(a, b, base_m, F);
        simpson2 = simpson(a, b, m_2, F);
    }

    return simpson2;
}
