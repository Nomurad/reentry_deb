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


int main(){
    // std::cout<<"hello world!!"<<std::endl;

    ClassTest classtest("world !!");
    classtest.output();

    MyUtils::Math math_s;
    double beta = math_s.deg2rad(34.0);
    DebrisOrbit debrisOrbit(0.15, beta, 206.18, 122.0);
    debrisOrbit.calc_interface_conditions();

    return 0;
}