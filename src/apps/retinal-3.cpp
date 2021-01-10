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

    deom d(json["deom"]);

    const int    ti = json["retinal"]["ti"].number_value();
    const int    tf = json["retinal"]["tf"].number_value();
    const double dt = json["retinal"]["dt"].number_value();
    const int    nk = json["retinal"]["nk"].int_value();
    const int    ntor = json["retinal"]["ntor"].int_value();
    const int    nvib = json["retinal"]["nvib"].int_value();
    const double staticErr = json["retinal"]["staticErr"].number_value();
    const string rho0File = json["retinal"]["rho0File"].string_value();
    const string wfunFile = json["retinal"]["wfunFile"].string_value();
    const string sdipFile = json["retinal"]["sdipFile"].string_value();

    cx_mat rho0;
    if (rho0.load(rho0File,arma_ascii)) {
        rho0.print("Initial state!");
    } else {
        printf("Error in loading rho0!\n");
        return 0;
    }

    mat wfun;
    if (wfun.load(wfunFile,arma_ascii)) {
        // wfun.print("wfun!");
    } else {
        printf("Error in loading wfun!\n");
        return 0;
    }

    mat sdip;
    if (sdip.load(sdipFile,arma_ascii)) {
        // sdip.print("Dipole moment operator!");
    } else {
        printf("Error in loading sdip!\n");
        return 0;
    }

    pulse p(json["retinal"]["pulse"]);

    cx_cube ddos = zeros<cx_cube>(d.nsys,d.nsys,d.nmax);
    ddos.slice(0) = rho0;
    // propagation to equilibrium
    const mat& exph= expmat(-real(d.ham1)/d.temperature);
    ddos.slice(0).set_real(exph/trace(exph));
    d.equilibrium(ddos,dt,staticErr,nk);

    const int nt_i = ceil(abs(ti)/dt);
    const int nt_f = ceil(abs(tf)/dt);
    const int nt = nt_i+nt_f;

    double t = -(nt_i-1)*dt;
    p.turnon();

    // rho propagation from ti to tf
    FILE *ftor0 = fopen("prop-tor0.dat","w");
    FILE *ftor1 = fopen("prop-tor1.dat","w");
    FILE *fisom = fopen("prop-isom.dat","w");
    const int nblk = ntor*nvib;
    mat rhot;
    for (int it=0; it<nt; ++it) {
        printf ("Propagation %5.1f%%: nddo=%6d, lddo=%3d\n",
                100*it/static_cast<double>(nt), d.nddo, d.lddo);
        if (it%nk == 0) {
            // Cis/Trans
            double pcis_0 = 0.0;
            double pcis_1 = 0.0;
            double ptrans_0 = 0.0;
            double ptrans_1 = 0.0;
            // Write down torsional density
            fprintf (ftor0, "%16.6e", t/deom_fs2unit);
            fprintf (ftor1, "%16.6e", t/deom_fs2unit);
            for (int i=0; i<ntor; ++i) {
                double prob_theta0 = 0.0;
                double prob_theta1 = 0.0;
                for (int j=0; j<nvib; ++j) {
                    const int ii = i*nvib+j;
                    const int jj = i*nvib+j+nblk;
                    prob_theta0 += as_scalar(wfun.row(ii)*rhot*wfun.row(ii).t());
                    prob_theta1 += as_scalar(wfun.row(jj)*rhot*wfun.row(jj).t());
                }
                if (i<ntor/4 || i>(ntor/4)*3) {
                    pcis_0 += prob_theta0;
                    pcis_1 += prob_theta1;
                } else {
                    ptrans_0 += prob_theta0;
                    ptrans_1 += prob_theta1;
                }
                fprintf (ftor0, "%16.6e", prob_theta0);
                fprintf (ftor1, "%16.6e", prob_theta1);
            }
            fprintf (ftor0, "\n");
            fprintf (ftor1, "\n");
            fprintf (fisom, "%16.6e%16.6e%16.6e%16.6e%16.6e\n", t/deom_fs2unit, pcis_0, pcis_1, ptrans_0, ptrans_1);
            fflush (ftor0);
            fflush (ftor1);
            fflush (fisom);
        }
        d.rk4 (ddos,t,dt,sdip,p);
        t += dt;
        if (t >= abs(ti)) {
            p.turnoff();
        }
    }
    fclose (fisom);
    fclose (ftor1);
    fclose (ftor0);

    return 0;
}
