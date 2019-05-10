#include "ReentryTrj.hpp"

// DebrisOrbit::DebrisOrbit(double dv = 4.0, 
//                     double beta_d = (M_PI/6.0), 
//                     double start_height = 200,
//                     double h_interface = 120.0){
//         dV[0] = dv;
//         beta = beta_d;
        
//         r[0] = r_e + 206.18;
//         // r[0] = 206.18;
//         r[1] = r_e + h_interface; //[km] 地球半径＋interface高度

//         std::cout<<"dv "<<dv<<std::endl;
//         std::cout<<"beta_d "<<beta_d*180/M_PI<<std::endl;
//         std::cout<<"interface height "<<h_interface<<"\n"<<std::endl;
        
//         initial_calc();
// }


void DebrisOrbit::initial_calc(){
    Vc = sqrt(g*pow(r_e,2.0)/r[0]);
    // Vc = sqrt(g*r[0]);
    return;
}

double DebrisOrbit::calc_interface_conditions(){
    
    V[0] = sqrt( pow(Vc,2.0) + pow(dV[0],2.0) - 2.0*Vc*dV[0]*cos(beta) );
    gamma[0] = acos( (( Vc - dV[0]*cos(beta) ) / V[0] ) );
    
    V[1] = sqrt( pow(V[0], 2.0) - (2.0*mu/r[0])+(2.0*mu/r[1]) );
    gamma[1] = acos( (r[0]*V[0]*cos(gamma[0]))/(r[1]*V[1]) );
    
    std::cout<< "V_C = " << Vc <<std::endl;
    std::cout<< "V_D = " << V[0] <<std::endl;
    std::cout<< "Gamma_D = " << (gamma[0]*180.0/M_PI) <<std::endl;
    std::cout<< "V_E = " << V[1] <<std::endl;
    std::cout<< "Gamma_E = " << rad2deg(gamma[1]) <<std::endl;
    std::cout<< "V_hat = " << V[1]/Vc <<std::endl;
    // std::cout<< "Gamma_E = " <<  gamma[1] << "\n" <<std::endl;

    // std::cout<<"V_E = "<< (r[0]*V[0]*cos(2.10305*M_PI/180.0)/(r[1]*cos(-1.6*M_PI/180.0))) <<std::endl;

    return 0.0;
}