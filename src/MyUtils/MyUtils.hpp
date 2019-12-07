#ifndef MYUTILS_HPP
#define MYUTILS_HPP

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <cmath>
#include <limits>
#include <boost/math/special_functions/expint.hpp>

class MyUtils{
    private:
        int para;


    public:
        class Integrate{
            public:
            double Trapezoidal(int (*f)(int), int n, double st, double fin);
            double Simpson(int n, double st, double fin);
        };

        class Math{
            public:
                double rad2deg(double radian);
                double deg2rad(double degree);
                Integrate integrate;
        };
        // Math math;

};

#endif //MYUTILS_HPP