#include <iostream>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

class Classtest{
    private:
        std::string name;

    public:
        Classtest(std::string s){
            name = s;
        }

        void output(){
            std::cout<<"hallo, "<<name<<std::endl;
            std::cout<<""<<std::endl;
        }

};

class DebrisOrbit{
    private:
        int para1 = 0;

        double g = 9.8*0.001; //[km/s]
        double mu = 3.986004418*pow(10.0, 5.0); //[km^3/s^2];

        double sma = 6829.677; //[km]
        double ecc = 0.001393; //[-]
        double inc = 96.980; //[deg]
        double raan = 312.158; //[deg]
        double arg_perigee = 35.169; //[deg]

        double Vc = 10.0; //[km/s]
        double gamma_c = 0.0;

        double r[2]; //軌道半径(r_D, r_E)
        double V[2]; //verosity
        double gamma[2]; //飛行経路角
        double dV[2]; //ΔV
        double beta; //水平方向に対するΔVの方向

    public:
        DebrisOrbit(double dv = 4.0, 
                    double beta_d = (M_PI/6.0), 
                    double h_interface = 120.0){
            dV[0] = dv;
            beta = beta_d;
            
            r[0] = 6371.0 + 206.18;
            r[1] = 6371.0 + h_interface; //[km] 地球半径＋interface高度

            std::cout<<"dv "<<dv<<std::endl;
            std::cout<<"beta_d "<<beta_d*180/M_PI<<std::endl;
            std::cout<<"interface height "<<h_interface<<"\n"<<std::endl;
        }

        double calc_interface_conditions(){
            
            Vc = sqrt(g*r[0]);
            V[0] = sqrt( pow(Vc,2.0) + pow(dV[0],2.0) - 2.0*Vc*dV[0]*cos(beta) );
            gamma[0] = acos( (( Vc - dV[0]*cos(beta) ) / V[0] ) );
            
            V[1] = sqrt( pow(V[0], 2.0) - (2.0*mu/r[0])+(2.0*mu/r[1]) );
            gamma[1] = acos( (r[0]*V[0]*cos(gamma[0]))/(r[1]*V[1]) );
            
            std::cout<< "V_C = " << Vc <<std::endl;
            std::cout<< "V_D = " << V[0] <<std::endl;
            std::cout<< "Gamma_D = " << (gamma[0]*180.0/M_PI) <<std::endl;
            std::cout<< "V_E = " << V[1] <<std::endl;
            // std::cout<< "Gamma_E = " <<  (gamma[1]*180.0/M_PI) <<std::endl;
            std::cout<< "Gamma_E = " <<  gamma[1] <<std::endl;

            return 0.0;
        }

        void print_paras(){
            std::cout<< sma <<std::endl;
            std::cout<< ecc <<std::endl;
            std::cout<< inc <<std::endl;
            std::cout<< raan <<std::endl;
            std::cout<< arg_perigee <<std::endl;

        }
};

int main(){
    // std::cout<<"hello world!!"<<std::endl;

    Classtest classtest("world !!");
    classtest.output();

    double beta = (34.0*M_PI/180.0);
    DebrisOrbit debrisOrbit(0.5, beta, 122.0);
    debrisOrbit.calc_interface_conditions();

    return 0;
}