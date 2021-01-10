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

    const int    nt = json["rhot"]["nt"].int_value();
    const double dt = json["rhot"]["dt"].number_value();
    const double staticErr = json["rhot"]["staticErr"].number_value();
    const int    nk = json["rhot"]["nk"].int_value();
    const string sdip1File = json["rhot"]["sdip1File"].string_value();
    const string sdip2File = json["rhot"]["sdip2File"].string_value();
    mat  sdip1;
    if (sdip1.load (sdip1File, arma_ascii)) {
        sdip1.print(sdip1File);
    } else {
        printf("Fail to load sdip1!\n");
    }
    mat  sdip2;
    if (sdip2.load (sdip2File, arma_ascii)) {
        sdip2.print(sdip2File);
    } else {
        printf("Fail to load sdip2!\n");
    }

    cx_cube rho0 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);
    cx_cube rho1 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);

    const mat& exph= expmat(-real(d1.ham1)/d1.temperature);
    rho0.slice(0).set_real(exph/trace(exph));
    d1.equilibrium (rho0,dt,staticErr,nk);

    printf("Equilibrium pop:");
    for (int i=0; i<d1.nsys; ++i) {
        printf("%16.6e", real(rho0(i,i,0)));
    }
    printf("\n");

    cx_vec ft = zeros<cx_vec>(nt); 

    d1.oprAct(rho1,sdip1,rho0,'c');
    FILE *log_t = fopen("1d-corr-t.dat","w");
    for (int it=0; it<nt; ++it) {
        double t = it*dt;
        ft(it) = deom_ci*d1.Trace(sdip2,rho1);
        printf ("In 1d-correlation: it=%d, nddo=%d, lddo=%d\n", it, d1.nddo, d1.lddo);
        fprintf(log_t, "%16.6e%16.6e%16.6e\n", t/deom_fs2unit, real(ft(it)), imag(ft(it)));
        d1.rk4 (rho1,t,dt);
    }
    fclose(log_t);

    // 1D FFT
    ft(0) *= 0.5;
    const double dw = 2.0*deom_pi/(nt*dt);
    const cx_vec& fw = ifft(ft)*nt*dt;
    
    FILE *log_w = fopen("1d-corr-w.dat","w");
    for (int iw=nt/2; iw<nt; ++iw) {
        double w = (iw-nt)*dw/deom_cm2unit;
        fprintf(log_w, "%16.6e%16.6e%16.6e\n", w, real(fw(iw)), imag(fw(iw)));
    }
    for (int iw=0; iw<nt/2; ++iw) {
        double w = iw*dw/deom_cm2unit;
        fprintf(log_w, "%16.6e%16.6e%16.6e\n", w, real(fw(iw)), imag(fw(iw)));
    }
    fclose(log_w);

    return 0;
}
