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
#include "blockdeom2.hpp"

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
    const double t2_max = json["spec"]["t2max"].number_value();
    const double w3_max = json["spec"]["w3max"].number_value();
    const int    nt1 = json["spec"]["nt1"].int_value();
    const int    nt2 = json["spec"]["nt2"].int_value();
    const int    nt3 = json["spec"]["nt3"].int_value();
    const double dt = json["spec"]["dt"].number_value();
    const double staticErr = json["spec"]["staticErr"].number_value();
    const int    nk = json["spec"]["nk"].int_value();
    const string sch_hei = json["spec"]["sch_hei"].string_value();
    const string sdipFile = json["spec"]["sdipFile"].string_value();
    mat  sdip;
    if (sdip.load (sdipFile, arma_ascii)) {
        sdip.print(sdipFile);
    } else {
        printf("Fail to load sdip!\n");
    }

    corr3rdPP2_monomer (w1_max, t2_max, w3_max, nt1, nt2, nt3, dt, staticErr, nk, sdip, sch_hei[0], s, b, h);
    
    return 0;
}
