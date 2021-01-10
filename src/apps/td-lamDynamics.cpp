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

    deom d0(json["deom"]);

    const int nlam = json["rhot"]["nlam"].int_value();
    const int nt = json["rhot"]["nt"].int_value();
    const int nk = json["rhot"]["nk"].int_value();
    const double dt = json["rhot"]["dt"].number_value();

    cx_vec etal = d0.coef_lft;
    cx_vec etar = d0.coef_rht;
    vec etaa = d0.coef_abs;
    vec delt = d0.delt_res;

    vec Ft = zeros<vec>(nt);
    vec Gt = zeros<vec>(nt);

    cx_cube ddo0 = zeros<cx_cube>(d0.nsys,d0.nsys,d0.nmax);
    mat eham = expmat_sym(-real(d0.ham1)/d0.temperature);
    ddo0.slice(0) = eham/trace(eham)*deom_c1;

    FILE *f = fopen("thermo.0","w");
    {
        double Svn = d0.entropy(ddo0);
        double Umi = d0.Umicro(ddo0);
        double Ami = d0.Amicro(ddo0);
        double Uma = d0.Umacro(ddo0);
        double Ama = d0.Amacro(ddo0);
        fprintf(f,"%d%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e\n", 0, 0.0, Svn, Umi, Ami, Uma, Ama);
    }
    fclose(f);

    for (int ilam=1; ilam<=nlam; ++ilam) {
        deom d(json["deom"]);
        // update coefs
        d.coef_lft = etal*ilam/nlam;
        d.coef_rht = etar*ilam/nlam;
        d.coef_abs = etaa*ilam/nlam;
        d.delt_res = delt*ilam/nlam;
        // propagation for each lam
        //
        cx_cube ddos(ddo0);
        
        stringstream ss;
        ss << "thermo." << ilam;
        FILE *fAS = fopen(ss.str().c_str(),"w");
        for (int it=0; it<nt; ++it) {
            const double t = it*dt;
            if (it%nk == 0) {
                double Svn = d.entropy(ddos);
                double Umi = d.Umicro(ddos);
                double Ami = d.Amicro(ddos);
                double Uma = d.Umacro(ddos);
                double Ama = d.Amacro(ddos);
                double hsb = Uma-Umi;
                Ft(it) += hsb/(2*ilam);
                //Gt(it) += Uma/(2*ilam);//newadd
                fprintf (fAS, "%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e\n", t/deom_fs2unit, Ft(it), Svn, Umi, Ami, Uma, Ama);
                //fprintf (fAS, "%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e\n", t/deom_fs2unit, Gt(it),Ft(it), Svn, Umi, Ami, Uma, Ama);
            }
            d.rk4 (ddos,t,dt);
        }
        fclose(fAS);
    }

    return 0;
}
