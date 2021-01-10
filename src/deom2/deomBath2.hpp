/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#ifndef DEOMBATH_H_
#define DEOMBATH_H_

#include <map>
#include <string>
#include "armadillo"
#include "json11.hpp"
#include "deomConst.hpp"

using namespace std;
using namespace arma;
using namespace json11;

class bath {

    public:

    double temperature;
    ivec   modLabel;
    cx_vec coef_lft;
    cx_vec coef_rht;
    vec    coef_abs;
    cx_vec expn_gam;
    vec    delt_res;
    vec    alpha1;
    vec    alpha2;

    bath (const Json& json) {

        const string modeFile = json["modeFile"].string_value();
        const string etalFile = json["etalFile"].string_value();
        const string etarFile = json["etarFile"].string_value();
        const string etaaFile = json["etaaFile"].string_value();
        const string expnFile = json["expnFile"].string_value();
        const string delrFile = json["delrFile"].string_value();
        const string alp1File = json["alp1File"].string_value();
        const string alp2File = json["alp2File"].string_value();

        temperature = json["temp"].number_value();

        printf ("$InitBath\n");
        if (modLabel.load (modeFile, arma_ascii)) {
            modLabel.print("modLabel");
        } else {
            printf ("modLabel is not loaded!\n");
        }
        if (coef_lft.load (etalFile, arma_ascii)) {
            coef_lft.print("coef_lft");
        } else {
            printf ("coef_lft is not loaded!\n");
        }
        if (coef_rht.load (etarFile, arma_ascii)) {
            coef_rht.print("coef_rht");
        } else {
            printf ("coef_rht is not loaded!\n");
        }
        if (coef_abs.load (etaaFile, arma_ascii)) {
            coef_abs.print("coef_abs");
        } else {
            printf ("coef_abs is not loaded!\n");
        }
        if (expn_gam.load (expnFile, arma_ascii)) {
            expn_gam.print("expn_gam");
        } else {
            printf ("expn_gam is not loaded!\n");
        }
        if (delt_res.load (delrFile, arma_ascii)) {
            delt_res.print("delt_res");
        } else {
            printf ("delt_res is not loaded!\n");
        }
        if (alpha1.load (alp1File, arma_ascii)) {
            alpha1.print("alpha1");
        } else {
            printf ("alpha1 is not loaded!\n");
        }
        if (alpha2.load (alp2File, arma_ascii)) {
            alpha2.print("alpha2");
        } else {
            printf ("alpha2 is not loaded!\n");
        }
        printf ("$InitBath\n\n");
    }

    bath (const bath& rhs): 
          temperature(rhs.temperature),
          modLabel(rhs.modLabel),
          coef_lft(rhs.coef_lft), 
          coef_rht(rhs.coef_rht), 
          coef_abs(rhs.coef_abs), 
          expn_gam(rhs.expn_gam), 
          delt_res(rhs.delt_res),
          alpha1(rhs.alpha1), alpha2(rhs.alpha2) {}

    bath (const double temp, const ivec& mlbl, const cx_vec& etal, const cx_vec& etar, const vec& etaa, 
          const cx_vec& expn, const vec& delr, const vec& alp1, const vec& alp2): 
          temperature(temp),
          modLabel(mlbl),
          coef_lft(etal), 
          coef_rht(etar), 
          coef_abs(etaa),
          expn_gam(expn), 
          delt_res(delr),
          alpha1(alp1), alpha2(alp2){}

    bath& operator= (const bath& rhs) {
        if (this != &rhs) {
            temperature = rhs.temperature;
            modLabel = rhs.modLabel;
            coef_lft = rhs.coef_lft;
            coef_rht = rhs.coef_rht;
            coef_abs = rhs.coef_abs;
            expn_gam = rhs.expn_gam;
            delt_res = rhs.delt_res;
            alpha1   = rhs.alpha1;
            alpha2   = rhs.alpha2;
        }
        return *this;
    }

   ~bath () {};

};

#endif
