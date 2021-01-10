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

    const int nk = json["rhot"]["nk"].int_value();
    const double dt = json["rhot"]["dt"].number_value();
    const double staticErr = json["rhot"]["staticErr"].number_value();
    const string weightFile = json["rhot"]["weightFile"].string_value();
    const string lambdaFile = json["rhot"]["lambdaFile"].string_value();
    vec weight;
    vec lambda;

    if (weight.load (weightFile, arma_ascii)) {
        weight.print(weightFile);
    } else {
        printf("Fail to load weight!\n");
    }

    if (lambda.load (lambdaFile, arma_ascii)) {
        lambda.print(lambdaFile);
    } else {
        printf("Fail to load lambda!\n");
    }

    cx_vec etal = d.coef_lft;
    cx_vec etar = d.coef_rht;
    vec etaa = d.coef_abs;
    vec delt = d.delt_res;

    double Iacc = 0.0;
    double Ine = 0.0;

    cx_cube ddos = zeros<cx_cube>(d.nsys,d.nsys,d.nmax);
    mat eham = expmat_sym(-real(d.ham1)/d.temperature);
    ddos.slice(0) = eham/trace(eham)*deom_c1;
    FILE *f = fopen("thermo.real","w");
    {
        double Svn = d.entropy(ddos);
        fprintf(f,"%d%20.10e%20.10e%20.10e%20.10e\n",    0, Svn, 0.0, 0.0, 0.0);
    }
    for (unsigned int ilam=0; ilam<weight.n_rows; ++ilam) {
        // update coefs
        d.coef_lft = etal*lambda(ilam);
        d.coef_rht = etar*lambda(ilam);
        d.coef_abs = etaa*lambda(ilam);
        d.delt_res = delt*lambda(ilam);
        // solve eq
        d.equilibrium(ddos, dt, staticErr, nk, "Prop");
        // evaluate thermodynamic quantities
        double Svn = d.entropy(ddos);
        double Umi = d.Umicro(ddos);
        // double Ami = d.Amicro(ddos);
        double Uma = d.Umacro(ddos);
        // double Ama = d.Amacro(ddos);
        double hsb = Uma-Umi;
        Iacc += 0.25*weight(ilam)/lambda(ilam)*hsb;
        //
        double tne = 0.0;
        for (int mp=0; mp<d.nind; ++mp) {
            const int m = d.modLabel(mp);
            ivec key0 = zeros<ivec>(mp+1);
            key0(mp) = 1;
            const double sn = sqrt(d.coef_abs(mp))*real(trace(d.qmd1.slice(m)*ddos.slice(0)));
            TrieNode* p = d.tree.find(key0);
            if (p && p->rank>=0) {
                int loc = p->rank;
                tne += sn*real(trace(ddos.slice(loc)));
            }
        }
        Ine += 0.25*weight(ilam)/lambda(ilam)*tne;
        ////
        fprintf(f,"%d%20.10e%20.10e%20.10e%20.10e\n", ilam, Svn, Iacc, Ine, Iacc-Ine);
    }
    fclose(f);

    return 0;
}
