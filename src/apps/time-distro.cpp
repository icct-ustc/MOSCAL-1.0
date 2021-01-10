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

    const int inistate = json["rhot"]["inistate"].int_value();
    const int obsstate = json["rhot"]["obsstate"].int_value();
    const int    nt = json["rhot"]["nt"].int_value();
    const double tau = json["rhot"]["tau"].number_value();
    const double dt = json["rhot"]["dt"].number_value();
    const double staticErr = json["rhot"]["staticErr"].number_value();
    const int    nk = json["rhot"]["nk"].int_value();

    int ntau = floor(tau/dt);
    double dt_res = tau-dt*ntau;

    cx_cube rho1 = zeros<cx_cube>(d.nsys,d.nsys,d.nmax);
    rho1(inistate,inistate,0) = 1.0;

    // d.equilibrium(rho1, dt, staticErr, nk, "SCI3");

    int iter = 0;
    int nddo_backup(d.nddo);
    hkey keys_backup(d.keys);
    cx_cube ddos_backup = zeros<cx_cube>(size(rho1));
    vec pop = zeros<vec>(d.nsys);
    while (iter<maxIter) {
        double t = 0;
        pop = real(rho1.slice(0).diag());
        pop.print("pop");
        // von Neumann wavefunction collapse
        for (int iddo=0; iddo<d.nddo; ++iddo) {
            for (int i=0; i<d.nsys-1; ++i) {
                for (int j=i+1; j<d.nsys; ++j) {
                    rho1(i,j,iddo) = rho1(j,i,iddo) = 0;
                }
            }
        }
        // Check convergence
        if (d.notconverged(nddo_backup,keys_backup,ddos_backup,rho1,staticErr)) {
            nddo_backup = d.nddo;
            keys_backup.rows(0,d.nddo-1) = d.keys.rows(0,d.nddo-1);
            ddos_backup.head_slices(d.nddo) = rho1.head_slices(d.nddo);
        } else {
            break;
        }
        // propagation
        for (int itau=0; itau<ntau; ++itau) {
            d.rk4 (rho1,t,dt);
            t += dt;
        }
        d.rk4 (rho1,t,dt_res);
        iter += 1;
    } 

    if (iter == maxIter) {
        printf ("Error! Failed to converge!\n");
    }

    for (int iddo=0; iddo<d.nddo; ++iddo) {
        cx_double tmp = rho1(obsstate,obsstate,iddo);
        rho1.slice(iddo).zeros();
        rho1(obsstate,obsstate,iddo) = tmp;
    }
    double pa0 = real(rho1(obsstate,obsstate,0));

    FILE *frho = fopen("prop-rho.dat","w");
    for (int it=0; it<nt; ++it) {
        double t = it*tau;
        if (it%nk == 0) {
            printf ("Propagation %5.1f%%: nddo=%6d, lddo=%3d\n",
                100*it/static_cast<double>(nt), d.nddo, d.lddo);
        }
        // Ua-propagate
        for (int itau=0; itau<ntau; ++itau) {
            d.rk4 (rho1,t,dt);
            t += dt;
        }
        d.rk4 (rho1,t,dt_res);
        // sum-over-b
        double pa1 = real(trace(rho1.slice(0)))-real(rho1(obsstate,obsstate,0));
        fprintf (frho, "%16.6e%16.6e%16.6e\n", t/deom_fs2unit, pa1, pa1/pa0);
        // Ua-measure
        for (int iddo=0; iddo<d.nddo; ++iddo) {
            cx_double tmp = rho1(obsstate,obsstate,iddo);
            rho1.slice(iddo).zeros();
            rho1(obsstate,obsstate,iddo) = tmp;
        }
    }
    fclose (frho);

    return 0;
}
