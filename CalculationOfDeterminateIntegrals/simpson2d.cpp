#include "simpson2d.h"

double simpson2d(double a, double b, double c, double d, int n, int m, double (*F)(double, double)){
    const double width_x = (b-a)/(2*n);
    const double width_y = (d-c)/(2*m);

    double simpson2d_integral = 0;

    for(int i = 0; i < n; i++){
        // double inner_result = 0;
        for(int j = 0; j < m; j++){
            // calculate f2i2j part
            double x1 = a + width_x*2*i;
            double y1 = c + width_y*2*j;
            double f1 = F(x1, y1);

            //calculate F 2i+1, 2j part
            double x2 = a + width_x*(2*i+1);
            double y2 = c + width_y*(2*j);
            double f2 = 4*F(x2, y2);

            //calculate F 2i+2, j+2 part
            double x3 = a + width_x*(2*i+2);
            double y3 = c + width_y*(2*j);
            double f3 = F(x3, y3);

            //calculate F 2i, 2j+1 part
            double x4 = a + width_x*(2*i);
            double y4 = c + width_y*(j+1);
            double f4 = 4*F(x4, y4);

            //calculate F 2i+1, 2j+1 part
            double x5 = a + width_x*(2*i+1);
            double y5 = c + width_y*(2*j+1);
            double f5 = 16*F(x5, y5);

            //calculate F 2i+2, 2j+1 part
            double x6 = a + width_x*(2*i+2);
            double y6 = c + width_y*(2*j+1);
            double f6 = 4*F(x6, y6);

            //calculate F 2i, 2j+2 part
            double x7 = a + width_x*(2*i);
            double y7 = c + width_y*(2*j+2);
            double f7 = F(x7, y7);

            //calculate F 2i+1, 2j+2 part
            double x8 = a + width_x*(2*i+1);
            double y8 = c + width_y*(2*j+2);
            double f8 = 4*F(x8, y8);

            //calculate F 2i+2, 2j+2 part
            double x9 = a + width_x*(2*i+2);
            double y9 = c + width_y*(2*j+2);
            double f9 = F(x9, y9);

            double inner_result = f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9;
            simpson2d_integral += inner_result;
        }
        // simpson2d_integral += inner_result;
    }

    simpson2d_integral *= width_x*width_y/9;

    return simpson2d_integral;
}