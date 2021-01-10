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
    const double dt = json["rhot"]["dt"].number_value();
    const int    nk = json["rhot"]["nk"].int_value();

    cx_cube ddos = zeros<cx_cube>(d.nsys,d.nsys,d.nmax);
    ddos(inistate,inistate,0) = 1.0;
    // const mat& exph= expmat(-real(d.ham1)/d.temperature);
    // ddos.slice(0).set_real(exph/trace(exph));
    // ddos.slice(0) = 0.5*(ddos.slice(0)+ddos.slice(0).t());

    FILE *frho = fopen("bloch.dat","w");
    for (int it=0; it<nt; ++it) {
        const double t = it*dt;
        printf ("Propagation %5.1f%%: nddo=%6d, lddo=%3d\n",
                100*it/static_cast<double>(nt), d.nddo, d.lddo);
        if (it%nk == 0) {
            const double x = 2*real(ddos(1,0,0));
            const double y = 2*imag(ddos(1,0,0));
            const double z = real(ddos(0,0,0)-ddos(1,1,0));
            fprintf (frho, "%16.6e%24.16e%24.16e%24.16e\n", t/deom_fs2unit, x, y, z);
        }
        d.rk4 (ddos,t,dt);
    }
    fclose (frho);

    return 0;
}
