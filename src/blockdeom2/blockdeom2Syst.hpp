/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#ifndef BLOCKDEOM2SYST_H_
#define BLOCKDEOM2SYST_H_

#include <string>
#include "armadillo"
#include "json11.hpp"
#include "trie.hpp"

using namespace std;
using namespace arma;
using namespace json11;

class syst {
    
    public:
    
    int nsys;
    int nmod;
    cx_mat  ham1;
    cx_cube qmd1;
    int nl;
    int nr;
    cx_mat  haml, hamr;
    cx_cube qmdl, qmdr;
    
    syst (const Json& json) {
        const string hamsFile = json["hamsFile"].string_value();
        const string qmdsFile = json["qmdsFile"].string_value();
        printf ("$InitSyst\n");
        if (ham1.load (hamsFile, arma_ascii)) {
            ham1.print(hamsFile);
        } else {
            printf ("Fail to load ham1!\n");
        }
        if (qmd1.load (qmdsFile, arma_ascii)) {
            qmd1.print(qmdsFile);
        } else {
            printf ("Fail to load qmd1!\n");
        }
        nsys = qmd1.n_rows;
        nmod = qmd1.n_slices;
        printf ("$InitSyst\n\n");
    }
    
    syst (const cx_mat& _h, const cx_cube& _q1): 
         nsys(_q1.n_rows), nmod(_q1.n_slices), ham1(_h), qmd1(_q1) {}

    syst (const syst& _s): nsys(_s.nsys), nmod(_s.nmod), ham1(_s.ham1), qmd1(_s.qmd1) {}

    syst (const cx_mat& hl, const cx_mat hr, const cx_cube& ql, const cx_cube& qr): nmod(ql.n_slices), nl(hl.n_rows), nr(hr.n_rows), haml(hl), hamr(hr), qmdl(ql), qmdr(qr) {}
        
   ~syst () {}
};

#endif
