/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "deom.hpp"

static const int maxIter = 100000000;

// return (x^n/n!)e^{-x}
static double poisson_gauss (const double x, const double lambda) {
    double result = exp(-(x-lambda)*(x-lambda)/(2.0*lambda))/sqrt(2*deom_pi*lambda);
    return result;
}

// return (x^n/n!)e^{-y}
static double poisson_gauss2 (const double x, const double lambda, const double y) {
    double result = exp(-(x-lambda)*(x-lambda)/(2.0*lambda)-y)/sqrt(2*deom_pi*lambda);
    return result;
}

static void collapse (deom& d, cx_cube& ddos, const double I, const double kappa, const double I0, const double tau) {
    double P11 = poisson_gauss(I*tau,(1.0+kappa)*(1.0+kappa)*I0*tau);
    double P22 = poisson_gauss(I*tau,(1.0-kappa)*(1.0-kappa)*I0*tau);
    double P12 = poisson_gauss2(I*tau,(1.0-kappa*kappa)*I0*tau,2*kappa*kappa*I0*tau);
    double pIv = 1.0/(P11*real(ddos(0,0,0))+P22*real(ddos(1,1,0)));
    // update density operators after weak measure
    for (int iddo=0; iddo<d.nddo; ++iddo) {
        ddos(0,0,iddo) *= pIv*P11;
        ddos(1,0,iddo) *= pIv*P12;
        ddos(0,1,iddo) *= pIv*P12;
        ddos(1,1,iddo) *= pIv*P22;
    }     
}

static double measure (const cx_mat& rho, const double I, const double kappa, const double I0, const double tau) {
    double P11 = poisson_gauss(I*tau,(1.0+kappa)*(1.0+kappa)*I0*tau);
    double P22 = poisson_gauss(I*tau,(1.0-kappa)*(1.0-kappa)*I0*tau);
    double pnt = P11*real(rho(0,0))+P22*real(rho(1,1));
    return pnt;
}

int main () {

    ifstream jsonFile("input.json");
    stringstream strStream;
    strStream << jsonFile.rdbuf();
    string jsonStr = strStream.str();
    string err;

    const Json json = Json::parse(jsonStr,err);
    if (!err.empty()) {
        printf ("Error in parsing input file: %s\n", err.c_str());
        return 0;
    }

    deom d1(json["deom"]);

    const int inistate = json["rhot"]["inistate"].int_value();
    const int    nt = json["rhot"]["nt"].int_value();
    const int    nk = json["rhot"]["nk"].int_value();
    const double dt = json["rhot"]["dt"].number_value();
    const double staticErr = json["rhot"]["staticErr"].number_value();
    const double kappa = json["rhot"]["kappa"].number_value();
    const double tau = json["rhot"]["tau"].number_value();
    const double I0 = json["rhot"]["I0"].number_value();
    const double dI = json["rhot"]["dI"].number_value();

    cx_cube rho1 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);
    cx_cube rho2 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);
    rho1(inistate,inistate,0) = 1.0;

    d1.equilibrium(rho1, dt, staticErr, 1, "SCI2");

    for (int n1=0; n1<=10; ++n1) {
        double I1 = dI*n1;
        // collapse first time
        deom d2(d1);
        rho2.slices(0,d2.nddo-1) = rho1.slices(0,d2.nddo-1);
        collapse (d2, rho2, I1, kappa, I0, tau);
        // save to file
        double t = 0.0;
        stringstream ss;
        ss << "pn_" << n1 << ".dat";
        FILE *fpn = fopen(ss.str().c_str(),"w");
        for (int it=0; it<nt; ++it) {
            for (int jt=0; jt<nk; ++jt) {
                d2.rk4 (rho2,t,dt);
                t += dt;
            }
            // collapse first time
            fprintf (fpn, "%16.6e", t);
            for (int n2=0; n2<=10; ++n2) {
                double I2 = dI*n2;
                double pn2 = measure (rho2.slice(0), I2, kappa, I0, tau);
                fprintf (fpn, "%16.6e", pn2);
            }
            fprintf (fpn, "\n");
        }
        fclose (fpn);
    }

    return 0;
}
