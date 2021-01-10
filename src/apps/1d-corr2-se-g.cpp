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
#include "deom2.hpp"

static void corr1st2(const double w_max, const int nt, const double dt,
              const double staticErr, const int nk,
              const mat& sdip, const cube& pdip, const vec& bdip,
              const char& sch_hei, const syst& s, const bath& b, const hidx& h) {

    const double dw1 = w_max/nt;
    const double dt1 = 2.0*deom_pi/w_max;
    const int    mt  = floor(dt1/dt);
    const double dt1_res = dt1-dt*mt; 

    deom2 d1(s,b,h);

    cx_cube rho_t0 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);
    rho_t0(1,1,0) = deom_c1;
    d1.equilibrium (rho_t0,dt,staticErr,nk,"Prop");

    cx_cube rho_t1 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);
    d1.oprAct(rho_t1,sdip,pdip,bdip,rho_t0,'l');

    cx_vec ft = zeros<cx_vec>(nt); 
    
    if (sch_hei == 's') { // sch-picture
        for (int it=0; it<nt; ++it) {
            ft(it) = deom_ci*conj(d1.Trace(sdip,pdip,bdip,rho_t1));
            printf ("In sch-pic: it=%d, nddo=%d, lddo=%d\n", it, d1.nddo, d1.lddo);
            double t1 = it*dt1;
            for (int jt=0; jt<mt; ++jt) {
                d1.rk4 (rho_t1,t1,dt);
                t1 += dt;
            }
            d1.rk4 (rho_t1,t1,dt1_res);
        }
    } else if (sch_hei == 'h') { // hei-picture
        deom2 d2(s,b,h);
        cx_cube oprs = zeros<cx_cube>(d2.nsys,d2.nsys,d2.nmax);
        d2.iniHei(oprs,sdip,pdip,bdip);
        for (int it=0; it<nt; ++it) {
            cx_double ctmp = 0.0;   
            for (int iado=0; iado<d1.nddo; ++iado) {
                TrieNode* p = d2.tree.find(d1.keys(iado).key);
                if (p && p->rank>=0) {
                    const int jado = p->rank;
                    ctmp += trace(oprs.slice(jado)*rho_t1.slice(iado));
                }
            }
            ft(it) = deom_ci*conj(ctmp); 
            printf ("In hei-pic: it=%d, nddo=%d, lddo=%d\n", it, d2.nddo, d2.lddo);
            double t1 = it*dt1;
            for (int jt=0; jt<mt; ++jt) {
                d2.rk4 (oprs,t1,dt,'h');
                t1 += dt;
            }
            d2.rk4 (oprs,t1,dt1_res,'h');
        }
    }

    // 1D FFT
    ft(0) *= 0.5;
    const cx_vec& fw = ifft(ft)*nt*dt1;

    // write time-domain signal
    FILE *f1 = fopen("corr.t","w");
    for (int i=0; i<nt; ++i) {
        fprintf (f1, "%16.6e%16.6e%16.6e\n", dt1*i/deom_fs2unit, real(ft(i)), imag(ft(i)));
    }
    fclose(f1);
    
    // write frequency-domain signal
    FILE *f2 = fopen("corr.w","w");
    for (int i=nt/2; i<nt; ++i) {
        double w = (i-nt)*dw1/deom_cm2unit;
        fprintf(f2, "%16.6e%16.6e%16.6e\n", w, real(fw(i)), imag(fw(i)));
    }
    for (int i=0; i<nt/2; ++i) {
        double w = i*dw1/deom_cm2unit;
        fprintf(f2, "%16.6e%16.6e%16.6e\n", w, real(fw(i)), imag(fw(i)));
    }
    fclose(f2);
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

    syst s(json["deom"]["syst"]);
    bath b(json["deom"]["bath"]);
    hidx h(json["deom"]["hidx"]);

    const double w1_max = json["spec"]["w1max"].number_value();              
    const int    nt1 = json["spec"]["nt1"].int_value();
    const double dt = json["spec"]["dt"].number_value();
    const double staticErr = json["spec"]["staticErr"].number_value();
    const int    nk = json["spec"]["nk"].int_value();
    const string sch_hei = json["spec"]["sch_hei"].string_value();
    const string sdipFile = json["spec"]["sdipFile"].string_value();
    const string pdipFile = json["spec"]["pdipFile"].string_value();
    const string bdipFile = json["spec"]["bdipFile"].string_value();
    mat  sdip;
    cube pdip;
    vec  bdip;
    if (sdip.load (sdipFile, arma_ascii)) {
        sdip.print(sdipFile);
    } else {
        printf("Fail to load sdip!\n");
    }
    if (pdip.load (pdipFile, arma_ascii)) {
        pdip.print(pdipFile);
    } else {
        printf("Fail to load pdip!\n");
    }
    if (bdip.load (bdipFile, arma_ascii)) {
        bdip.print(bdipFile);
    } else {
        printf("Fail to load bdip!\n");
    }

    corr1st2 (w1_max, nt1, dt, staticErr, nk, sdip, pdip, bdip, sch_hei[0], s, b, h);

    return 0;
}
