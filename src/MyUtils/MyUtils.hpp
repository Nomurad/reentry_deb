#ifndef MYUTILS_HPP
#define MYUTILS_HPP

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <cmath>
#include <limits>

class MyUtils{
    private:
        int para;


    public:
        class Math{
            public:
                double rad2deg(double radian);
                double deg2rad(double degree);
        };
        // Math math;

};

#endif //MYUTILS_HPP