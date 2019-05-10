#include "MyUtils.hpp"

double MyUtils::Math::rad2deg(double radian){
    return radian*180.0/M_PI;
}

double MyUtils::Math::deg2rad(double degree){
    return degree*M_PI/180.0;
}