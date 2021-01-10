/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#ifndef BLOCKDEOM_H_
#define BLOCKDEOM_H_

#include "armadillo"
#include "trie.hpp"

#include "deomConst.hpp"
#include "blockdeom2Syst.hpp"
#include "deomBath2.hpp"
#include "deomHidx.hpp"
#include "deomPulse.hpp"

using namespace std;
using namespace arma;

class blockdeom2: public syst, public bath, public hidx {

    public:

        // Syst and Bath
        int ng;
        int ne;
        int nf;

        cx_cube ddos1;
        cx_cube ddos2;
        cx_cube ddos3;

        blockdeom2 (const syst& s, const bath& b, const hidx& h, const char lc, const char rc): 
            syst (s), bath (b), hidx (h) {
            ng = 1;
            ne = ham1.n_rows;
            nf = ne*(ne-1)/2;

            if (lc == 'g') {
                nl = ng;
                haml = get_h00();
                qmdl = get_q00();
            } else if (lc == 'e') {
                nl = ne;
                haml = get_h11();
                qmdl = get_q11();            
            } else if (lc == 'f') {
                nl = nf;
                haml = get_h22();
                qmdl = get_q22();                        
            }

            if (rc == 'g') {
                nr = ng;
                hamr = get_h00();
                qmdr = get_q00();
            } else if (rc == 'e') {
                nr = ne;
                hamr = get_h11();
                qmdr = get_q11();            
            } else if (rc == 'f') {
                nr = nf;
                hamr = get_h22();
                qmdr = get_q22();                        
            }        

            ddos1.set_size(nl,nr,nmax);
            ddos2.set_size(nl,nr,nmax);
            ddos3.set_size(nl,nr,nmax);
        }

        blockdeom2 (const blockdeom2& rhs, const char lc, const char rc): 
            syst(rhs.ham1,rhs.qmd1), 
            bath(rhs.temperature, rhs.modLabel, rhs.coef_lft, rhs.coef_rht, rhs.coef_abs, rhs.expn_gam, rhs.delt_res, rhs.alpha1, rhs.alpha2), 
            hidx(rhs.nind, rhs.lmax, rhs.nmax, rhs.lddo, rhs.nddo, 
                    rhs.ferr, rhs.keys, rhs.tree, rhs.expn) {
            ng = 1;
            ne = ham1.n_rows;
            nf = ne*(ne-1)/2;

            if (lc == 'g') {
                nl = ng;
                haml = get_h00();
                qmdl = get_q00();
            } else if (lc == 'e') {
                nl = ne;
                haml = get_h11();
                qmdl = get_q11();            
            } else if (lc == 'f') {
                nl = nf;
                haml = get_h22();
                qmdl = get_q22();                        
            }

            if (rc == 'g') {
                nr = ng;
                hamr = get_h00();
                qmdr = get_q00();
            } else if (rc == 'e') {
                nr = ne;
                hamr = get_h11();
                qmdr = get_q11();            
            } else if (rc == 'f') {
                nr = nf;
                hamr = get_h22();
                qmdr = get_q22();                        
            }  

            ddos1.set_size(nl,nr,nmax);
            ddos2.set_size(nl,nr,nmax);
            ddos3.set_size(nl,nr,nmax);
        }

        blockdeom2 (const cx_mat& hl, const cx_mat& hr, const cx_cube& ql, const cx_cube& qr,
                   const bath& b, const hidx& h, const char lc, const char rc): 
            syst (hl,hr,ql,qr), bath (b), hidx (h) {
            ddos1.set_size(nl,nr,nmax);
            ddos2.set_size(nl,nr,nmax);
            ddos3.set_size(nl,nr,nmax);
        }

        blockdeom2 (const blockdeom2& rhs, const cx_mat& hl, const cx_mat& hr, const cx_cube& ql, const cx_cube& qr, const char lc, const char rc): 
            syst(hl,hr,ql,qr), 
            bath(rhs.temperature, rhs.modLabel, rhs.coef_lft, rhs.coef_rht, rhs.coef_abs, rhs.expn_gam, rhs.delt_res, rhs.alpha1, rhs.alpha2), 
            hidx(rhs.nind, rhs.lmax, rhs.nmax, rhs.lddo, rhs.nddo, 
                    rhs.ferr, rhs.keys, rhs.tree, rhs.expn) {
            ddos1.set_size(nl,nr,nmax);
            ddos2.set_size(nl,nr,nmax);
            ddos3.set_size(nl,nr,nmax);
        }

        ~blockdeom2 () {}

        cx_double Trace (const mat& sdip, const cx_cube& ddos) const {
            cx_double result = trace(sdip*ddos.slice(0));
            return result;
        }

        void oprAct (cx_cube& ddos, const mat& sdip, const cx_cube& rhos, const char lrc='l');

        void iniSch (cx_cube& ddos, const mat& sdip);

        void iniHei (cx_cube& oprs, const mat& sdip);

        void remSch (cx_cube& d_ddos, const cx_cube& ddos, const double t);

        void remHei (cx_cube& d_ddos, const cx_cube& ddos, const double t);

        void rem (cx_cube& d_ddos, const cx_cube& ddos, const double t, const char sch_hei='s') {
            if (sch_hei == 's') {
                remSch (d_ddos, ddos, t);
            } else if (sch_hei == 'h') {
                remHei (d_ddos, ddos, t);
            } else {
                printf("sch_hei is invalid!\n");
            }
        }

        void propagation (cx_cube& ddos, const double dt=0.005, const int nt=1000, const int nk=10);

        inline bool is_valid (const cx_mat& ddo) const {
            return any(abs(vectorise(ddo))>ferr);
        }

        void filter (cx_cube& ddos) {
            int n = 1;
            int l = 0;
            for (int iddo=1; iddo<nddo; ++iddo) {
                TrieNode* p = tree.find(keys(iddo).key);
                if (is_valid(ddos.slice(iddo))) {
                    if (n != iddo) {
                        p->rank = n;
                        keys(n) = keys(iddo);
                        ddos.slice(n) = ddos.slice(iddo);
                    }
                    l = l>(p->tier)?l:(p->tier);
                    ++n;
                } else {
                    p->rank = -9527;
                }
            } 
            lddo = l;
            nddo = n;
        }

        bool notconverged (const int& nddo_backup, const hkey& keys_backup, 
                        const cx_cube& ddos_backup, const cx_cube& ddos, const double& tol) {
            for (int iddo=0; iddo<nddo_backup; ++iddo) {
                const auto *nod = tree.find(keys_backup(iddo).key);
                double maxDiff = 0.0;
                if (nod != NULL && nod->rank>=0) {
                    maxDiff = max(abs(vectorise(ddos_backup.slice(iddo)-ddos.slice(nod->rank))));
                } else {
                    maxDiff = max(abs(vectorise(ddos_backup.slice(iddo))));
                }
                if (maxDiff > tol) {
                    printf("Error bigger than %16.6e, while tol=%16.6e\n", maxDiff, tol);
                    return true;
                }
            }
            return false;
        }

        template<typename... Tc>
            void rk4 (cx_cube& ddos, const double t, const double dt, const Tc&... args) {

                const double dt2 = dt*0.5;
                const double dt6 = dt/6.0;

                // K1
                const int nddo0 = nddo;
                rem (ddos1,ddos,t,args...);
                ddos3.slices(0,nddo0-1) = ddos.slices(0,nddo0-1)+ddos1.slices(0,nddo0-1)*dt2;
                if (nddo > nddo0) {
                    ddos3.slices(nddo0,nddo-1) = ddos1.slices(nddo0,nddo-1)*dt2;
                }
                // K2
                const int nddo1 = nddo;
                rem (ddos2,ddos3,t+0.5*dt,args...);
                ddos1.slices(0,nddo1-1) += ddos2.slices(0,nddo1-1)*2.0;
                if (nddo > nddo1) {
                    ddos1.slices(nddo1,nddo-1) = ddos2.slices(nddo1,nddo-1)*2.0;
                }
                ddos3.slices(0,nddo0-1) = ddos.slices(0,nddo0-1)+ddos2.slices(0,nddo0-1)*dt2;
                if (nddo > nddo0) {
                    ddos3.slices(nddo0,nddo-1) = ddos2.slices(nddo0,nddo-1)*dt2;
                }
                // K3
                const int nddo2 = nddo;
                rem (ddos2,ddos3,t+0.5*dt,args...);
                ddos1.slices(0,nddo2-1) += ddos2.slices(0,nddo2-1)*2.0;
                if (nddo > nddo2) {
                    ddos1.slices(nddo2,nddo-1) = ddos2.slices(nddo2,nddo-1)*2.0;
                }
                ddos3.slices(0,nddo0-1) = ddos.slices(0,nddo0-1)+ddos2.slices(0,nddo0-1)*dt;
                if (nddo > nddo0) {
                    ddos3.slices(nddo0,nddo-1) = ddos2.slices(nddo0,nddo-1)*dt;
                }
                // K4
                const int nddo3 = nddo;
                rem (ddos2,ddos3,t+dt,args...);
                ddos1.slices(0,nddo3-1) += ddos2.slices(0,nddo3-1);
                if (nddo > nddo3) {
                    ddos1.slices(nddo3,nddo-1) = ddos2.slices(nddo3,nddo-1);
                }
                ddos.slices(0,nddo0-1) += ddos1.slices(0,nddo0-1)*dt6;
                if (nddo > nddo0) {
                    ddos.slices(nddo0,nddo-1) = ddos1.slices(nddo0,nddo-1)*dt6;
                }
                filter (ddos);
            }

        /**
         * -----------------------------------------------
         * h1+h2
         *       h1+h3
         *             h1+h4
         *                   h2+h3
         *                         h2+h4
         *                               h3+h4
         * -----------------------------------------------
         * f(i,j) = (n-1+n-i-1)/2*(i+1)-(n-i-1)+j-i-1; i<j
         * -----------------------------------------------
         * f(0,3) = 2;
         * f(1,2) = 3;
         * f(2,3) = 5; 
         */   
        inline int get_fno (const int i, const int j) const {
            if (i < j) {
                return (2*ne-2-i)*(i+1)/2-ne+j;
            } else if (j < i) {
                return (2*ne-2-j)*(j+1)/2-ne+i;
            } else {
                return -1;
            }
        }

        cx_mat get_h00 () const {
            return zeros<cx_mat>(1,1); 
        }

        cx_mat get_h11 () const {
            return ham1;
        }

        cx_mat get_h22 () const {
            cx_mat h2 = zeros<cx_mat>(nf,nf);
            for (int i=0; i<ne-1; ++i) {
                for (int j=i+1; j<ne; ++j) {
                    const int ij = get_fno(i,j);
                    for (int k=0; k<ne-1; ++k) {
                        for (int l=k+1; l<ne; ++l) {
                            const int kl = get_fno(k,l);
                            if (i==k && j==l) {
                                h2(ij,kl) = ham1(i,i)+ham1(j,j);
                            } else if (i==k && j!=l) {
                                h2(ij,kl) = ham1(j,l);
                            } else if (i!=k && j==l) {
                                h2(ij,kl) = ham1(i,k);
                            }
                        }
                    }
                }
            }
            return h2;
        }

        cx_cube get_q00 () const {
            return zeros<cx_cube>(1,1,nmod);
        }

        cx_cube get_q11 () const {
            return qmd1;
        }

        cx_cube get_q22 () const {
            cx_cube q2 = zeros<cx_cube>(nf,nf,nmod);
            for (int m=0; m<nmod; ++m) {
                for (int i=0; i<ne-1; ++i) {
                    for (int j=i+1; j<ne; ++j) {
                        const int ij = get_fno(i,j);
                        for (int k=0; k<ne-1; ++k) {
                            for (int l=k+1; l<ne; ++l) {
                                const int kl = get_fno(k,l);
                                if (i==k && j==l) {
                                    q2(ij,kl,m) = qmd1(i,i,m)+qmd1(j,j,m);
                                } else if (i==k && j!=l) {
                                    q2(ij,kl,m) = qmd1(j,l,m);
                                } else if (i!=k && j==l) {
                                    q2(ij,kl,m) = qmd1(i,k,m);
                                }
                            }
                        }
                    }
                }
            }
            return q2;        
        }

        mat  get_d01 (const  mat& sdip_eg) const {
            mat dipo = sdip_eg.t();
            return dipo;
        }

        cube get_d01 (const cube& pdip_eg) const {
            const int nm = pdip_eg.n_slices;
            cube dipo = zeros<cube>(1,ne,nm);
            for (int m=0; m<nm; ++m) {
                dipo.slice(m) = pdip_eg.slice(m).t();
            }
            return dipo;
        }

        mat  get_d10 (const  mat& sdip_eg) const {
            return sdip_eg;
        }

        cube get_d10 (const cube& pdip_eg) const {
            return pdip_eg;
        }

        /*
         * 1      |2  3  4
         *   2    |1        3  4  
         *     3  |   1     2     4
         *       4|      1     2  3
         * -------------------------
         *        |12
         *        |   13
         *        |      14
         *        |         23
         *        |            24
         *        |               34
         */ 

        mat  get_d12 (const mat& sdip_eg) const {
            mat dipo = zeros<mat>(ne,nf);
            for (int m1=0; m1<ne; ++m1) {
                for (int m2=0; m2<ne; ++m2) {
                    if (m1 != m2) {
                        const int fn = get_fno(m1,m2);
                        dipo(m1,fn) += sdip_eg(m2,0);
                    }
                }
            }
            return dipo;
        }    

        cube get_d12 (const cube& pdip_eg) const {
            const int nm = pdip_eg.n_slices;
            cube dipo = zeros<cube>(ne,nf,nm);
            for (int im=0; im<nm; ++im) {
                for (int m1=0; m1<ne; ++m1) {
                    for (int m2=0; m2<ne; ++m2) {
                        if (m1 != m2) {
                            const int fn = get_fno(m1,m2);
                            dipo(m1,fn,im) += pdip_eg(m2,0,im);
                        }
                    }
                }
            }
            return dipo;
        }    

        mat  get_d21 (const mat& sdip_eg) const {
            mat dipo = zeros<mat>(nf,ne);
            for (int m1=0; m1<ne; ++m1) {
                for (int m2=0; m2<ne; ++m2) {
                    if (m1 != m2) {
                        const int fn = get_fno(m1,m2);
                        dipo(fn,m1) = sdip_eg(m2,0);
                    }
                }
            }
            return dipo;
        }    

        cube get_d21 (const cube& pdip_eg) const {
            const int nm = pdip_eg.n_slices;
            cube dipo = zeros<cube>(nf,ne,nm);
            for (int im=0; im<nm; ++im) {
                for (int m1=0; m1<ne; ++m1) {
                    for (int m2=0; m2<ne; ++m2) {
                        if (m1 != m2) {
                            const int fn = get_fno(m1,m2);
                            dipo(fn,m1,im) = pdip_eg(m2,0,im);
                        }
                    }
                }
            }
            return dipo;
        }    

};


void corr3rdPP2_monomer (const double w1_max, const double t2_max, const double w3_max,
        const int nt1, const int nt2, const int nt3, const double dt,
        const double staticErr, const int nk, const mat& sdip, 
        const char& sch_hei, const syst& s, const bath &b, const hidx& h); 

#endif
