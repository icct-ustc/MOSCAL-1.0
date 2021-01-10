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
//#include "deom.hpp"
#include "blockdeom.hpp"

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

    const int    nt = json["spec"]["nt"].int_value();
    const double dt = json["spec"]["dt"].number_value();
    const double tf = json["spec"]["tf"].number_value();
    const string sch_hei  = json["spec"]["sch_hei"].string_value();
    const string digeFile = json["spec"]["digeFile"].string_value();
    const string diefFile = json["spec"]["diefFile"].string_value();
    const string difeFile = json["spec"]["difeFile"].string_value();
    const string diegFile = json["spec"]["diegFile"].string_value();
    const string ham0File = json["spec"]["ham0File"].string_value();
    const string ham1File = json["spec"]["ham1File"].string_value();
    const string ham2File = json["spec"]["ham2File"].string_value();
    const string qmd0File = json["spec"]["qmd0File"].string_value();
    const string qmd1File = json["spec"]["qmd1File"].string_value();
    const string qmd2File = json["spec"]["qmd2File"].string_value();
    const string rho0File = json["spec"]["rho0File"].string_value();
    mat  dip_ge, dip_ef, dip_fe, dip_eg;
    cx_mat ham0, ham1, ham2;
    cx_cube qmd0, qmd1, qmd2;

    if (dip_ge.load (digeFile, arma_ascii)) {
        dip_ge.print(digeFile);
    } else {
        printf("Fail to load dip_ge!\n");
    }
    if (dip_ef.load (diefFile, arma_ascii)) {
        dip_ef.print(diefFile);
    } else {
        printf("Fail to load dip_ef!\n");
    }
    if (dip_fe.load (difeFile, arma_ascii)) {
        dip_fe.print(difeFile);
    } else {
        printf("Fail to load dip_fe!\n");
    }
    if (dip_eg.load (diegFile, arma_ascii)) {
        dip_eg.print(diegFile);
    } else {
        printf("Fail to load dip_eg!\n");
    }
    if (ham0.load (ham0File, arma_ascii)) {
        ham0.print(ham0File);
    } else {
        printf("Fail to load ham0!\n");
    }
    if (ham1.load (ham1File, arma_ascii)) {
        ham1.print(ham1File);
    } else {
        printf("Fail to load ham1!\n");
    }
    if (ham2.load (ham2File, arma_ascii)) {
        ham2.print(ham2File);
    } else {
        printf("Fail to load ham2!\n");
    }
    if (qmd0.load (qmd0File, arma_ascii)) {
        qmd0.print(qmd0File);
    } else {
        printf("Fail to load qmd0!\n");
    }
    if (qmd1.load (qmd1File, arma_ascii)) {
        qmd1.print(qmd1File);
    } else {
        printf("Fail to load qmd1!\n");
    }
    if (qmd2.load (qmd2File, arma_ascii)) {
        qmd2.print(qmd2File);
    } else {
        printf("Fail to load qmd2!\n");
    }
    cx_mat rho0;
    if (rho0.load (rho0File, arma_ascii)) {
        rho0.print(rho0File);
    } else {
        printf("Fail to load rho0!\n");
    }

    // respPcc (tf, nt, dt, rho0, dipo_mat, dipoDirs, sch_hei[0], s, b, h);
    
//    corrPcc (tf, nt, dt, rho0, sch_hei[0], 
//             dip_ge, dip_ef, dip_fe, dip_eg, 
//             ham0, ham1, ham2, qmd0, qmd1, qmd2, b, h);

    return 0;
}
