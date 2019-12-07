#include "ReentryTrj.hpp"
#include "../MyUtils/MyUtils.hpp"

using std::cout;
using std::cin;
using std::endl;

DebrisOrbit::DebrisOrbit(){
    initial_calc();
}


void DebrisOrbit::initial_calc(){
    // Vc = sqrt(g*r_D);
    Vc = sqrt(g*pow(r_e,2.0)/r_D);
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
    
    r_D = r_e + 206.18;
    // r_D = 206.18;
    r_E = r_e + h_interface; //[km] 地球半径＋interface高度

    cout<<"μ "<< mu << " | "<<g*r_e*r_e << endl;
    cout<<"dv "<< dv <<endl;
    cout<<"beta_d "<< beta*180.0/M_PI <<endl;
    cout<<"beta(rad) "<< beta <<endl;
    cout<<"interface height "<< h_interface <<"\n"<<endl;
    cout<<"Vc(init) = "<< Vc << endl;
    
    initial_calc();
    
    V_D = sqrt( pow(Vc,2.0) + pow(dV[0],2.0) - 2.0*Vc*dV[0]*cos(beta) );
    gamma_D = acos( (( Vc - dV[0]*cos(beta) ) / V_D ) );
    
    V_E = sqrt( pow(V_D, 2.0) - (2.0*mu/r_D)+(2.0*mu/r_E) );
    gamma_E = acos( (r_D*V_D*cos(gamma_D))/(r_E*V_E) );
    gamma_E = deg2rad(-1.6);
    V_E = (r_D*V_D*cos(gamma_D))/(r_E*cos(gamma_E));
    
    this->rho_E = atmos_density(h_interface);

    /***** output *****/
    cout<< "V_C = " << Vc <<endl;
    cout<< "V_D = " << V_D <<endl;
    cout<< "Delta_V = " << dV[0] <<endl;
    cout<< "Gamma_D = " << (gamma_D*180.0/M_PI) << " | " << gamma_D <<endl;
    cout<< "V_E = " << V_E <<endl;
    cout<< "Gamma_E = " << rad2deg(gamma_E) <<endl;
    // cout<< "V_hat = " << V_E/sqrt(g*pow(r_e,2.0)/r_E) <<endl;
    cout<< "V_hat = " << V_E/sqrt(g*pow(r_e,2.0)/r_e) <<endl;
    // cout<< "Gamma_E = " <<  gamma_E << "\n" <<endl;
    cout<< "g = " << mu/(r_e*r_e) << endl;

    // cout<<"V_E = "<< (r_D*V_D*cos(2.10305*M_PI/180.0)/(r_E*cos(-1.6*M_PI/180.0))) <<endl;

    //指数積分の計算
    // auto inf = std::numeric_limits<double>::infinity();
    // cout<<std::expint(-inf)<<endl;

    return 0.0;
}

double DebrisOrbit::calc_flight_path_angle(double r){
    double Beta = 1.0/ScaleHeight;
    double h = r - r_e;
    double h_E = r_E - r_e;
    double b = -1.0/(C_B*beta*sin(gamma_E));
    double rho = atmos_density(h);
    double rho_E = atmos_density(h_E);
    
    // cout<<"scaleheight^-1 = "<< 1.0/6.7 <<endl;
    double ei = b*rho;
    double Ei;
    // cout<<"b = "<< b <<", rho = "<< rho <<", Ei="<<boost::math::expint(ei)<< endl;
    try{
        Ei = boost::math::expint(ei) - boost::math::expint(b*rho_E);
        this->rho = rho;
    }
    catch(...){
        ei = this->rho*b;
        Ei = boost::math::expint(ei) - boost::math::expint(b*rho_E);
        // cout <<"skip:"<<ei<<", "<<Ei << endl;
    }

    double bunbo = 1.0 + (2.0*cos(gamma_E)*log(r/r_E)) 
                    + (2*mu)*Ei/(Beta*V_E*V_E*r_E*r_E);
    
    double cos_gamma = cos(gamma_E)/sqrt(bunbo);

    return acos(cos_gamma);
}

double DebrisOrbit::calc_velocity(double r, double gamma){
    double V;

    double Beta = 1.0/ScaleHeight;
    double b = -1.0/(C_B*beta*sin(gamma_E));
    double h = r - r_e;
    double rho = atmos_density(h);
    double rho_E = atmos_density(r_E - r_e);

    double exp_in = -b*(rho-rho_E)/2.0;
    V = V_E*(r_E*cos(gamma_E) / (r*cos(gamma)))*exp(exp_in);

    return V;
}

double DebrisOrbit::Integrand(double r){
    double gamma = calc_flight_path_angle(r);
    double V = calc_velocity(r, gamma);
    return 1.0/(V*sin(gamma));
}

double DebrisOrbit::calc_reentry_span(double r){
    double t_span;
    double n = 1000;
    double st = r_E;
    double fin = r;
    t_span = math.integrate.Trapezoidal(Integrand, n, st, fin);
    
    return t_span;
}

double DebrisOrbit::calc_1dim_order_solution(){
    double trj;
    


    return trj;
}