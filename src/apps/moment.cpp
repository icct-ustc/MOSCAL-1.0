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

    const int nlam = json["rhot"]["nlam"].int_value();
    const int nk = json["rhot"]["nk"].int_value();
    const double dt = json["rhot"]["dt"].number_value();
    const double staticErr = json["rhot"]["staticErr"].number_value();

    mat  sdip = zeros<mat>(d.nsys,d.nsys);
    vec  bdip = ones<vec>(d.nind);
    cube pdip = zeros<cube>(d.nsys,d.nsys,d.nmod);

    cx_cube ddos = zeros<cx_cube>(d.nsys,d.nsys,d.nmax);
    cx_cube rho0 = zeros<cx_cube>(d.nsys,d.nsys,d.nmax);
    cx_cube rho1 = zeros<cx_cube>(d.nsys,d.nsys,d.nmax);
    mat eham = expmat_sym(-real(d.ham1)/d.temperature);
    ddos.slice(0) = eham/trace(eham)*deom_c1;
    // solve eq
    d.equilibrium(ddos, dt, staticErr, nk, "SCI2");

    FILE *fs = fopen("moment.dat","w");
    for (int m=0; m<d.nmod; ++m) {
        deom d1(d);
        pdip.zeros();
        pdip.slice(m) = eye<mat>(d1.nsys,d1.nsys);
        rho0.head_slices(d1.nddo) = ddos.head_slices(d1.nddo);
        fprintf (fs, "%d",m);
        for (int o=0; o<4; ++o) {
            d1.oprAct (rho1, sdip, pdip, bdip, rho0, 'l');
            rho0.head_slices(d1.nddo) = rho1.head_slices(d1.nddo);
            fprintf (fs, "%16.6e", real(trace(rho0.slice(0))));
        }
        fprintf (fs, "\n");
    }
    fclose(fs);

    return 0;
}
