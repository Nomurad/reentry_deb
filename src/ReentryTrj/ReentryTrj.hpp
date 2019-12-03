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
        const double r_e = 6378.0; //地球平均半径
        const double g = 9.8*0.001; //重力加速度[km/s]
        const double mu = 3.986004418*pow(10.0, 5.0); //重力定数(地球)[km^3/s^2];
        // double mu = 398600; //[km^3/s^2];

        double sma = 6829.677; //軌道長半径[km]
        double ecc = 0.001393; //離心率[-]
        double inc = 96.980; //軌道傾斜角[deg]
        double raan = 312.158; //昇交点赤経[deg]
        double arg_perigee = 35.169; //近地点引数[deg]

        double Vc = sqrt( (g*r_e) ); //10.0; //[km/s]
        double V_hat; //non-Dementionnal verosity
        double gamma_c = 0.0;

        double r[2]; //軌道半径(r_D, r_E)
        double V[2]; //verosity
        double gamma[2]; //飛行経路角
        double dV[2]; //ΔV
        double beta; //水平方向に対するΔVの方向

        double m = 100.0;    //質量 m[kg]
        double A = 10.0;     //基準面積 A[m^2]
        double C_D = 1.0;    //抵抗係数 C_D[-]
        double C_B = m/(C_D*A);    //弾道係数 C_B[kg/m^2] ((m/C_D*A)
        

        // from "Dynamics of Atmospheric Re-Entry",p38
        const double rho_atmosphere = 1.752;    //[kg/m^3]
        const double ScaleHight = 6.7;    //[km]

        MyUtils myUtils;
        MyUtils::Math math;
        double deg2rad(double x){return math.deg2rad(x);}
        double rad2deg(double x){return math.rad2deg(x);}

    public:
        double rho_atomosphere;
        double ScaleHeight;

        DebrisOrbit();

        void initial_calc();

        double calc_interface_conditions(double dv, 
                                        double beta_d, 
                                        double start_heigh,
                                        double h_interface);

        double calc_1dim_order_solution();

        // SpacePlane's Dynamics (coordinate system:ECI)
        double theta_dot();
        double phi_dot();
        double r_dot();
        double V_dot();
        double gamma_dot();
        double kai_dot();

        double atmos_density(double height);

};

#endif //REENTRYTRJ_HPP