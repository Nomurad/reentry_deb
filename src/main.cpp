#include <iostream>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

#include "MyUtils/MyUtils.hpp"
#include "ReentryTrj/ReentryTrj.hpp"

using std::cout;
using std::endl;

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


int main(){
    int i;
    int iter_max = 100;
    double h_E;
    std::vector<double> h_list, gamma_list, V_list;
    // std::cout<<"hello world!!"<<std::endl;
    ClassTest classtest("world !!");
    classtest.output();

    MyUtils::Math math_s;
    double beta = math_s.deg2rad(34.0);
    DebrisOrbit deb;
    cout<<"ballistic = " <<deb.get_ballistic()<<endl;
    deb.calc_interface_conditions(0.15, beta, 206.18, 122.0);

    deb.r = deb.r_E;
    h_E = deb.r_E- deb.r_e;

    cout<<"r = "<< deb.r << " "<< (deb.r-deb.r_e)*(11.0/iter_max)<< endl;
    for(i=0; i<iter_max+1; i++){
        deb.gamma = deb.calc_flight_path_angle(deb.r);
        deb.V = deb.calc_velocity(deb.r, deb.gamma);
        deb.r = deb.r_E - h_E*double(i)/double(iter_max);
        
        // cout<<i<<":"<<"h = "<<(deb.r-deb.r_e)<<", r = "<< deb.r 
        // <<", gamma = "<< math_s.rad2deg(deb.gamma) << endl;

        h_list.push_back(deb.r-deb.r_e);
        gamma_list.push_back(math_s.rad2deg(deb.gamma));
        V_list.push_back(deb.V);
    }

    // output csv & display
    std::string filename = "reentry_trj.csv";
    std::ofstream writing_file;
    writing_file.open(filename, std::ios::out);
    for(i=0;i<iter_max+1;i++){
        cout<<std::setw(4)<<i <<": h = "<<std::setw(10) <<h_list[i]
            <<", gamma = "<<std::setw(10) <<gamma_list[i]
            <<", V = " <<std::setw(10) <<V_list[i] << endl;
        writing_file << i << "," << h_list[i] << "," 
            << gamma_list[i] << "," << V_list[i] <<endl;
    }
    writing_file.close();

    cout<<"flight span = "<<deb.calc_reentry_span(deb.r_e)<<endl;

    return 0;
}