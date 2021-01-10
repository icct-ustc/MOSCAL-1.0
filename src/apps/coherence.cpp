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

    if (d.flag_rwa != 0) {
        printf("Rotating Wave Approximation\n");
    }
    const int inistate = json["rhot"]["inistate"].int_value();
    const int    nt = json["rhot"]["nt"].int_value();
    const double dt = json["rhot"]["dt"].number_value();
    const int    nk = json["rhot"]["nk"].int_value();

    cx_cube ddos = zeros<cx_cube>(d.nsys,d.nsys,d.nmax);
    ddos(inistate,inistate,0) = 1.0;

    vec eval = zeros<vec>(d.nsys);
    cx_mat evec = zeros<cx_mat>(d.nsys,d.nsys);
    FILE *fpop = fopen("eval.dat","w");
    FILE *fcoh = fopen("evec.dat","w");
    for (int it=0; it<nt; ++it) {
        const double t = it*dt;
        printf ("Propagation %5.1f%%: nddo=%6d, lddo=%3d\n",
                100*it/static_cast<double>(nt), d.nddo, d.lddo);
        if (it%nk == 0) {
            fprintf (fpop, "%16.6e", t/deom_fs2unit);
            fprintf (fcoh, "%16.6e", t/deom_fs2unit);
            eig_sym(eval, evec, ddos.slice(0));
            for (int i=0; i<d.nsys; ++i) {
                fprintf (fpop, "%16.6e", eval(i));
                for (int j=0; j<d.nsys; ++j) {
                    fprintf (fcoh, "%16.6e%16.6e", real(evec(i,j)), imag(evec(i,j)));
                }
            }
            fprintf (fpop, "\n");
            fprintf (fcoh, "\n");
        }
        d.rk4 (ddos,t,dt);
    }
    fclose (fcoh);
    fclose (fpop);

    return 0;
}
