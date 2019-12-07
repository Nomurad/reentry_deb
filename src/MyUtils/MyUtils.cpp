#include "MyUtils.hpp"

double MyUtils::Math::rad2deg(double radian){
    return radian*180.0/M_PI;
}

double MyUtils::Math::deg2rad(double degree){
    return degree*M_PI/180.0;
}

double MyUtils::Integrate::Trapezoidal(int (*f)(int), int n, double st, double fin){
    double h = (fin - st)/double(n);
    double x, s, k;

    x = st;
    s = 0;

    for(k=0;k<n;k++){
        x = x + h;
        s = f(x);
    }
    s = h*((f(st) + f(fin)) / 2 + s);

    return s;
}