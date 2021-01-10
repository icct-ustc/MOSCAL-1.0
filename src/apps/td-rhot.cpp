/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "deom.hpp"

static double IPR (const cx_mat& rho) {
    const double tmp1 = accu(abs(rho));
    const double tmp2 = accu(pow(abs(rho),2));
    return tmp1*tmp1/(rho.n_rows*tmp2);
}

static void td_rhot (deom& d, const double dt, const double ti, const double tf, const int nk, const mat& sdip, const cube& pdip, const vec& bdip, pulse& p) {

    const int nsys = d.nsys;
    const int nt_i = ceil(abs(ti)/dt);
    const int nt_f = ceil(abs(tf)/dt);
    const int nt = nt_i+nt_f;

    double t = -(nt_i-1)*dt;
    cx_cube rho1 = zeros<cx_cube>(nsys,nsys,d.nmax);
    rho1(0,0,0) = 1.0;   
    p.turnon();

    // rho propagation from ti to tf
    FILE *frho = fopen("td-rho.dat","w");
    FILE *fipr = fopen("td-ipr.dat","w");
    for (int it=0; it<nt; ++it) {

        fprintf(frho,"%16.6e",t/deom_fs2unit);
        for (int i=1; i<nsys; ++i) {
            fprintf(frho,"%16.6e",real(rho1(i,i,0)));
        }
        fprintf(frho,"\n");

        const cx_mat& rho0 = rho1.slice(0).submat(1,1,nsys-1,nsys-1);
        fprintf(fipr,"%16.6e%16.6e\n",t/deom_fs2unit,IPR(rho0));

        if (it%nk == 0) {
            printf ("%s: %5.1f%%, nddo = %6d, lddo = %2d\n", "td-proprho",
                    100*it/static_cast<double>(nt), d.nddo, d.lddo);
        }
        d.rk4 (rho1,t,dt,sdip,pdip,bdip,p);
        t += dt;
    }
    fclose(fipr);
    fclose(frho);
}

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

    const double dt = json["td-rhot"]["dt"].number_value();
    const double ti = json["td-rhot"]["ti"].number_value();
    const double tf = json["td-rhot"]["tf"].number_value();
    const int nk = json["td-rhot"]["nk"].int_value();
    const string sdipFile = json["td-rhot"]["sdipFile"].string_value();
    const string pdipFile = json["td-rhot"]["pdipFile"].string_value();
    const string bdipFile = json["td-rhot"]["bdipFile"].string_value();

    pulse p(json["td-rhot"]["pulse"]);

    mat  sdip;
    cube pdip;
    vec  bdip;

    if (sdip.load (sdipFile, arma_ascii)) {
        sdip.print(sdipFile);
    } else {
        printf("Fail to load sdip");
    }

    if (pdip.load (pdipFile, arma_ascii)) {
        pdip.print(pdipFile);
    } else {
        printf("Fail to load pdip");
    }

    if (bdip.load (bdipFile, arma_ascii)) {
        bdip.print(bdipFile);
    } else {
        printf("Fail to load bdip");
    }

    td_rhot (d, dt, ti, tf, nk, sdip, pdip, bdip, p);

    return 0;
}
