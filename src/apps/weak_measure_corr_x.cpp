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

static const int ntMax = 100000000;

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
    const double I0 = json["rhot"]["I0"].number_value();

    cx_cube rho1 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);
    rho1(inistate,inistate,0) = 1.0;

    d1.equilibrium(rho1, dt, staticErr, 1, "SCI2");

    // D_- Action
    // x -- detecting
    const double nDOTrho = 2*real(rho1(1,0,0));
    for (int iddo=0; iddo<d1.nddo; ++iddo) {
        // x -- detecting
        const double nDOTpsi = 2*real(rho1(1,0,iddo));
        const double bar_psi = real(rho1(0,0,iddo)+rho1(1,1,iddo));
        // Dbar -- Action
        const double tmp0 = nDOTpsi-nDOTrho*bar_psi;
        // Dvec -- Action
        // x -- detecting
        const double tmp1 = bar_psi-nDOTrho*nDOTpsi+nDOTrho*(2*real(rho1(1,0,iddo))-nDOTpsi);
        const double tmp2 = nDOTrho*2*imag(rho1(1,0,iddo));
        const double tmp3 = nDOTrho*(real(rho1(0,0,iddo)-rho1(1,1,iddo)));
        rho1(0,0,iddo) = 0.5*(tmp0+tmp3);
        rho1(0,1,iddo) = 0.5*(tmp1-deom_ci*tmp2);
        rho1(1,0,iddo) = 0.5*(tmp1+deom_ci*tmp2);
        rho1(1,1,iddo) = 0.5*(tmp0-tmp3);
    }

    const double dt_fft = nk*dt;
    cx_vec ft_xx = zeros<cx_vec>(nt); 
    cx_vec ft_yx = zeros<cx_vec>(nt); 
    cx_vec ft_zx = zeros<cx_vec>(nt); 

    // Tr[D+ G(t)\rho(0;D-)]
    FILE *log_t = fopen("weak_measure_corr.t","w");
    for (int it=0; it<nt; ++it) {
        double t = it*dt_fft;
        // xyz -- detecting
        const double cxx = 2*real(rho1(1,0,0));
        const double cyx = 2*imag(rho1(1,0,0));
        const double czx = real(rho1(0,0,0)-rho1(1,1,0));
        ft_xx(it) = (2*kappa*I0)*(2*kappa*I0)*cxx;
        ft_yx(it) = (2*kappa*I0)*(2*kappa*I0)*cyx;
        ft_zx(it) = (2*kappa*I0)*(2*kappa*I0)*czx;
        printf ("Weak-measure-correlation: it=%d, nddo=%d, lddo=%d\n", it, d1.nddo, d1.lddo);
        fprintf(log_t, "%16.6e%16.6e%16.6e%16.6e\n", t/deom_fs2unit, real(ft_xx(it)), real(ft_yx(it)), real(ft_zx(it)));
        for (int kt=0; kt<nk; ++kt) {
            d1.rk4 (rho1,t,dt);
            t += dt;
        }
    }
    fclose(log_t);

    // 1D FFT
    ft_xx(0) *= 0.5;
    ft_yx(0) *= 0.5;
    ft_zx(0) *= 0.5;
    const double dw = 2.0*deom_pi/(nt*dt_fft);
    const cx_vec& fw_xx = ifft(ft_xx)*nt*dt_fft;
    const cx_vec& fw_yx = ifft(ft_yx)*nt*dt_fft;
    const cx_vec& fw_zx = ifft(ft_zx)*nt*dt_fft;
    
    FILE *log_w = fopen("weak_measure_corr.w","w");
    for (int iw=nt/2; iw<nt; ++iw) {
        double w = (iw-nt)*dw/deom_cm2unit;
        fprintf(log_w, "%16.6e%16.6e%16.6e%16.6e\n", w, real(fw_xx(iw)), real(fw_yx(iw)), real(fw_zx(iw)));
    }
    for (int iw=0; iw<nt/2; ++iw) {
        double w = iw*dw/deom_cm2unit;
        fprintf(log_w, "%16.6e%16.6e%16.6e%16.6e\n", w, real(fw_xx(iw)), real(fw_yx(iw)), real(fw_zx(iw)));
    }
    fclose(log_w);

    return 0;
}
