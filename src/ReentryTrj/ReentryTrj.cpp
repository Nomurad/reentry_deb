#include "ReentryTrj.hpp"
#include "../MyUtils/MyUtils.hpp"

using std::cout;
using std::cin;
using std::endl;

DebrisOrbit::DebrisOrbit(){
    initial_calc();
}


void DebrisOrbit::initial_calc(){
    // Vc = sqrt(g*r[0]);
    Vc = sqrt(g*pow(r_e,2.0)/r[0]);
    return;
}

double DebrisOrbit::atmos_density(double height){
    double rho;
    const double rho_0 = rho_atmosphere;
    const double H = ScaleHeight;

    rho = rho_0 * exp(-height/H);

    return rho;
}

double DebrisOrbit::calc_interface_conditions(double dv = 4.0, 
                                            double beta_d = (M_PI/6.0), 
                                            double start_height = 200,
                                            double h_interface = 120.0){
    
    MyUtils::Math math_s;
    dV[0] = dv;
    beta = beta_d;
    
    r[0] = r_e + 206.18;
    // r[0] = 206.18;
    r[1] = r_e + h_interface; //[km] 地球半径＋interface高度

    cout<<"μ "<< mu << " | "<<g*r_e*r_e << endl;
    cout<<"dv "<< dv <<endl;
    cout<<"beta_d "<< beta*180.0/M_PI <<endl;
    cout<<"beta(rad) "<< beta <<endl;
    cout<<"interface height "<< h_interface <<"\n"<<endl;
    cout<<"Vc(init) = "<< Vc << endl;
    
    initial_calc();
    
    V[0] = sqrt( pow(Vc,2.0) + pow(dV[0],2.0) - 2.0*Vc*dV[0]*cos(beta) );
    gamma[0] = acos( (( Vc - dV[0]*cos(beta) ) / V[0] ) );
    
    V[1] = sqrt( pow(V[0], 2.0) - (2.0*mu/r[0])+(2.0*mu/r[1]) );
    gamma[1] = acos( (r[0]*V[0]*cos(gamma[0]))/(r[1]*V[1]) );
    gamma[1] = deg2rad(-1.6);
    V[1] = (r[0]*V[0]*cos(gamma[0]))/(r[1]*cos(gamma[1]));
    

    /***** output *****/
    cout<< "V_C = " << Vc <<endl;
    cout<< "V_D = " << V[0] <<endl;
    cout<< "Delta_V = " << dV[0] <<endl;
    cout<< "Gamma_D = " << (gamma[0]*180.0/M_PI) << " | " << gamma[0] <<endl;
    cout<< "V_E = " << V[1] <<endl;
    cout<< "Gamma_E = " << rad2deg(gamma[1]) <<endl;
    // cout<< "V_hat = " << V[1]/sqrt(g*pow(r_e,2.0)/r[1]) <<endl;
    cout<< "V_hat = " << V[1]/sqrt(g*pow(r_e,2.0)/r_e) <<endl;
    // cout<< "Gamma_E = " <<  gamma[1] << "\n" <<endl;
    cout<< "g = " << mu/(r_e*r_e) << endl;

    // cout<<"V_E = "<< (r[0]*V[0]*cos(2.10305*M_PI/180.0)/(r[1]*cos(-1.6*M_PI/180.0))) <<endl;

    //指数積分の計算
    auto inf = std::numeric_limits<double>::infinity();
    cout<<std::expint(-inf)<<endl;

    return 0.0;
}

double DebrisOrbit::calc_1dim_order_solution(){
    double trj;
    


    return trj;
}