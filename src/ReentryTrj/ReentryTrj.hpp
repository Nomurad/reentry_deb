#ifndef REENTRYTRJ_HPP
#define REENTRYTRJ_HPP

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <iostream>
#include <math.h>

#include "../MyUtils/MyUtils.hpp"

class DebrisOrbit{
    private:
        int para1 = 0;

        const double r_e = 6378.0; //地球平均半径
        const double g = 9.8*0.001; //[km/s]
        const double mu = 3.986004418*pow(10.0, 5.0); //[km^3/s^2];
        // double mu = 398600; //[km^3/s^2];

        double sma = 6829.677; //[km]
        double ecc = 0.001393; //[-]
        double inc = 96.980; //[deg]
        double raan = 312.158; //[deg]
        double arg_perigee = 35.169; //[deg]

        double Vc = 10.0; //[km/s]
        double V_hat; //non-Dementionnal verosity
        double gamma_c = 0.0;

        double r[2]; //軌道半径(r_D, r_E)
        double V[2]; //verosity
        double gamma[2]; //飛行経路角
        double dV[2]; //ΔV
        double beta; //水平方向に対するΔVの方向

        MyUtils myUtils;
        MyUtils::Math math;
        double deg2rad(double x){return math.deg2rad(x);}
        double rad2deg(double x){return math.rad2deg(x);}

    public:
        DebrisOrbit(double dv, 
                    double beta_d, 
                    double start_heigh,
                    double h_interface);

        void initial_calc();

        double calc_interface_conditions();
};

#endif //REENTRYTRJ_HPP