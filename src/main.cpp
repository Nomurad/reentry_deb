#include <iostream>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

#include "MyUtils/MyUtils.hpp"
#include "ReentryTrj/ReentryTrj.hpp"

class ClassTest{
    private:
        std::string name;

    public:
        ClassTest(std::string s){
            name = s;
        }

        void output(){
            std::cout<<"hallo, "<<name<<std::endl;
            std::cout<<""<<std::endl;
        }

};


// class DebrisOrbit{
//     private:
//         int para1 = 0;

//         const double r_e = 6378.0; //地球平均半径
//         const double g = 9.8*0.001; //[km/s]
//         const double mu = 3.986004418*pow(10.0, 5.0); //[km^3/s^2];
//         // double mu = 398600; //[km^3/s^2];

//         double sma = 6829.677; //[km]
//         double ecc = 0.001393; //[-]
//         double inc = 96.980; //[deg]
//         double raan = 312.158; //[deg]
//         double arg_perigee = 35.169; //[deg]

//         double Vc = 10.0; //[km/s]
//         double V_hat; //non-Dementionnal verosity
//         double gamma_c = 0.0;

//         double r[2]; //軌道半径(r_D, r_E)
//         double V[2]; //verosity
//         double gamma[2]; //飛行経路角
//         double dV[2]; //ΔV
//         double beta; //水平方向に対するΔVの方向

//         MyUtils myUtils;
//         double deg2rad(double x){return myUtils.math.deg2rad(x);}
//         double rad2deg(double x){return myUtils.math.rad2deg(x);}

//     public:
//         DebrisOrbit(double dv = 4.0, 
//                     double beta_d = (M_PI/6.0), 
//                     double start_height = 200,
//                     double h_interface = 120.0){
//             dV[0] = dv;
//             beta = beta_d;
            
//             r[0] = r_e + 206.18;
//             // r[0] = 206.18;
//             r[1] = r_e + h_interface; //[km] 地球半径＋interface高度

//             std::cout<<"dv "<<dv<<std::endl;
//             std::cout<<"beta_d "<<beta_d*180/M_PI<<std::endl;
//             std::cout<<"interface height "<<h_interface<<"\n"<<std::endl;
            
//             initial_calc();
//         }

//         void initial_calc(){
//             Vc = sqrt(g*pow(r_e,2.0)/r[0]);
//             // Vc = sqrt(g*r[0]);
//             return;
//         }

//         double calc_interface_conditions(){
            
//             V[0] = sqrt( pow(Vc,2.0) + pow(dV[0],2.0) - 2.0*Vc*dV[0]*cos(beta) );
//             gamma[0] = acos( (( Vc - dV[0]*cos(beta) ) / V[0] ) );
            
//             V[1] = sqrt( pow(V[0], 2.0) - (2.0*mu/r[0])+(2.0*mu/r[1]) );
//             gamma[1] = acos( (r[0]*V[0]*cos(gamma[0]))/(r[1]*V[1]) );
            
//             std::cout<< "V_C = " << Vc <<std::endl;
//             std::cout<< "V_D = " << V[0] <<std::endl;
//             std::cout<< "Gamma_D = " << (gamma[0]*180.0/M_PI) <<std::endl;
//             std::cout<< "V_E = " << V[1] <<std::endl;
//             std::cout<< "Gamma_E = " << rad2deg(gamma[1]) <<std::endl;
//             std::cout<< "V_hat = " << V[1]/Vc <<std::endl;
//             // std::cout<< "Gamma_E = " <<  gamma[1] << "\n" <<std::endl;

//             // std::cout<<"V_E = "<< (r[0]*V[0]*cos(2.10305*M_PI/180.0)/(r[1]*cos(-1.6*M_PI/180.0))) <<std::endl;

//             return 0.0;
//         }

//         void print_paras(){
//             std::cout<< sma <<std::endl;
//             std::cout<< ecc <<std::endl;
//             std::cout<< inc <<std::endl;
//             std::cout<< raan <<std::endl;
//             std::cout<< arg_perigee <<std::endl;

//         }
// };

int main(){
    // std::cout<<"hello world!!"<<std::endl;

    ClassTest classtest("world !!");
    classtest.output();

    double beta = (34.0*M_PI/180.0);
    DebrisOrbit debrisOrbit(0.15, beta, 206.18, 122.0);
    debrisOrbit.calc_interface_conditions();

    return 0;
}