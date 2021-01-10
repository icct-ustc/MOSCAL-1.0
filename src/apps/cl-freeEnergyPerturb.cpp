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

    cx_cube rho0 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);

    const mat& exph= expmat(-real(d1.ham1)/d1.temperature);
    rho0.slice(0).set_real(exph/trace(exph));
    d1.equilibrium (rho0,dt,staticErr,nk);

    printf("Equilibrium pop:");
    for (int i=0; i<d1.nsys; ++i) {
        printf("%16.6e", real(rho0(i,i,0)));
    }
    printf("\n");

    double db = 1.0/(d1.temperature*nt);
    for (int it=0; it<nt; ++it) {
        double t = it*db;
        printf ("In 1d-correlation: it=%d, nddo=%d, lddo=%d\n", it, d1.nddo, d1.lddo);
        rho0.slice(0).print("rho0");
        d1.rk4 (rho0,t,db,'c');
    }
    double Acl = d1.temperature*log(real(trace(rho0.slice(0))));
    FILE *fs = fopen("Acl.dat","w");
    fprintf(fs, "%16.6e%16.6e\n", d1.temperature, Acl);
    fclose(fs);

    return 0;
}
