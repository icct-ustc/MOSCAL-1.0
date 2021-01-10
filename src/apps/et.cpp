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

    const int inistate = json["rhot"]["inistate"].int_value();
    const int    nt = json["rhot"]["nt"].int_value();
    const int    nk = json["rhot"]["nk"].int_value();
    const double dt = json["rhot"]["dt"].number_value();
    const double staticErr = json["rhot"]["staticErr"].number_value();

    cx_cube ddos = zeros<cx_cube>(d.nsys,d.nsys,d.nmax);
    ddos(inistate,inistate,0) = 1.0;

    cx_mat hams(d.ham1);
    d.ham1(0,1) = d.ham1(1,0) = 0.0;
    d.equilibrium(ddos, dt, staticErr, nk, "SCI2");
    d.ham1 = hams;

    FILE *frho = fopen("prop-rho.dat","w");
    for (int it=0; it<nt; ++it) {
        const double t = it*dt;
        printf ("Propagation %5.1f%%: nddo=%6d, lddo=%3d\n",
                100*it/static_cast<double>(nt), d.nddo, d.lddo);
        if (it%nk == 0) {
            fprintf (frho, "%16.6e", t/deom_fs2unit);
            for (int i=0; i<d.nsys; ++i) {
                fprintf (frho, "%20.10e", real(ddos(i,i,0)));
            }
            fprintf (frho, "\n");
        }
        d.rk4 (ddos,t,dt);
    }
    fclose (frho);

    return 0;
}
