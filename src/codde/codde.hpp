/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#ifndef CODDE_H_
#define CODDE_H_

#include "armadillo"
#include "trie.hpp"

#include "deomConst.hpp"
#include "deomSyst.hpp"
#include "deomBath.hpp"
#include "deomPulse.hpp"

using namespace std;
using namespace arma;

class codde: public syst, public bath {

    public:

        int nind;
        int nddo;

        cx_cube ddos1;
        cx_cube ddos2;
        cx_cube ddos3;

        codde (const Json& json): syst (json["syst"]), bath (json["bath"]) {
            nind = expn_gam.n_rows;
            nddo = 2*nind+1;
            ddos1.set_size(nsys,nsys,nddo);
            ddos2.set_size(nsys,nsys,nddo);
            ddos3.set_size(nsys,nsys,nddo);
        }

        codde (const syst& s, const bath& b): syst (s), bath (b) {
            nind = expn_gam.n_rows;
            nddo = 2*nind+1;
            ddos1.set_size(nsys,nsys,nmax);
            ddos2.set_size(nsys,nsys,nmax);
            ddos3.set_size(nsys,nsys,nmax);
        }

        codde (const codde& rhs): syst(rhs.ham1,rhs.qmd1), 
            nind = expn_gam.n_rows;
            nddo = 2*nind+1;
            bath(rhs.temperature, rhs.modLabel, rhs.coef_lft, rhs.coef_rht, rhs.coef_abs, rhs.expn_gam, rhs.delt_res) {
            ddos1.set_size(nsys,nsys,nmax);
            ddos2.set_size(nsys,nsys,nmax);
            ddos3.set_size(nsys,nsys,nmax);
        }

        ~codde () {}

        void oprAct (cx_cube& d_ddos, const mat& sdip, const cx_cube& ddos, const char lrc='l');

        void iniHei (cx_cube& d_ddos, const mat& sdip);

        void remSch (cx_cube& d_ddos, const cx_cube& ddos, const double t);

        void remHei (cx_cube& d_ddos, const cx_cube& ddos, const double t);

        void rem (cx_cube& d_ddos, const cx_cube& ddos, const double t, const int& projection);

        void rem (cx_cube& d_ddos, const cx_cube& ddos, const double t, const ivec& projection);

        void rem (cx_cube& d_ddos, const cx_cube& ddos, const double t, const mat& sdip, const pulse& p);

        void rem (cx_cube& d_ddos, const cx_cube& ddos, const double t, const char sch_hei = 's') {
            if (sch_hei == 's') {
                remSch (d_ddos, ddos, t);
            } else if (sch_hei == 'h') {
                remHei (d_ddos, ddos, t);
            } else {
                printf("sch_hei is invalid!\n");
            }
        }

        void propagation (cx_cube& ddos, const double dt=0.005, const int nt=1000, const int nk=10);

        void equilibrium (cx_cube& ddos, const double dt=0.005, const double err=2.e-8, const int nk=10);

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

        cx_double Trace (const mat& sdip, const cx_cube& ddos) const {
            return trace(sdip*ddos.slice(0));
        }

};


void resp1st (const double w_max, const int nt, const double dt,
              const double staticErr, const int nk,
              const mat& sdip, 
              const char& sch_hei, const syst& s, const bath& b);

void resp2nd (const double w1_max, const double w2_max, const int nt1, const int nt2, const double dt,
              const double staticErr, const int nk,
              const mat& sdip, const cube& pdip, const vec& bdip,
              const char& sch_hei, const syst& s, const bath& b, const hidx& h);

void resp3rd (const double w1_max, const double t2_max, const double w3_max, const int nt1, const int nt2, const int nt3, const double dt,
              const double staticErr, const int nk,
              const mat& sdip, const cube& pdip, const vec& bdip,
              const char& sch_hei, const syst& s, const bath& b, const hidx& h);
              
#endif
