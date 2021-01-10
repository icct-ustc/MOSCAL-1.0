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

    const int nalp = json["rhot"]["nalp"].int_value();
    const int nk = json["rhot"]["nk"].int_value();
    const double dt = json["rhot"]["dt"].number_value();
    const double staticErr = json["rhot"]["staticErr"].number_value();

    cx_vec etal = d.coef_lft;
    cx_vec etar = d.coef_rht;
    vec etaa = d.coef_abs;
    vec delt = d.delt_res;

    double Iacc = 0.0;

    cx_cube ddos = zeros<cx_cube>(d.nsys,d.nsys,d.nmax);
    mat eham = expmat_sym(-real(d.ham1)/d.temperature);
    ddos.slice(0) = eham/trace(eham)*deom_c1;
    FILE *f = fopen("thermo.real","w");
    {
        double Svn = d.entropy(ddos);
        double Umi = d.Umicro(ddos);
        double Ami = d.Amicro(ddos);
        double Uma = d.Umacro(ddos);
        double Ama = d.Amacro(ddos);
        fprintf(f,"%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e\n", 0.0, Svn, Umi, Ami, Uma, Ama, 0.0, 1.0);
    }
    for (int ialp=1; ialp<=nalp; ++ialp) {
        double scaling = 1.0*ialp/nalp;
        // update coefs
        d.coef_lft = etal*scaling*scaling;
        d.coef_rht = etar*scaling*scaling;
        d.coef_abs = etaa*scaling*scaling;
        d.delt_res = delt*scaling*scaling;
        // solve eq
        d.equilibrium(ddos, dt, staticErr, nk, "SCI2");
        // evaluate thermodynamic quantities
        double Svn = d.entropy(ddos);
        double Umi = d.Umicro(ddos);
        double Ami = d.Amicro(ddos);
        double Uma = d.Umacro(ddos);
        double Ama = d.Amacro(ddos);
        double hsb = Uma-Umi;
        Iacc += hsb/scaling*0.5;
        double Aint = Iacc/nalp;
        double Zint = exp(-Aint/d.temperature);
        fprintf(f,"%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e\n", 1.0*ialp/nalp, Svn, Umi, Ami, Uma, Ama, Aint, Zint);
        Iacc += hsb/scaling*0.5;
    }
    fclose(f);

    return 0;
}
