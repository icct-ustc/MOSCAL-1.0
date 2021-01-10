/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cmath>
#include <sstream>
#include "blockdeom.hpp"

void calcPcc (const double t1, const double w1, const double st1, const double sw1,
              const double t2, const double w2, const double st2, const double sw2,
              const int nt, const int nw, const double dt, const double dw,
              const cx_mat& rho0, 
              const mat& dip_ge, const mat& dip_ef, const mat& dip_fe, const mat& dip_eg,
              const cx_mat& ham0, const cx_mat& ham1, const cx_mat& ham2,
              const cx_cube& qmd0, const cx_cube& qmd1, const cx_cube& qmd2,
              const bath& b, const hidx& h) {

    double dt_int = (t1-t2)/nt;
    int    mt = floor(dt_int/dt);
    double dt_res = dt_int-mt*dt;

    mat pccAD = zeros<mat>(nw,nw);
    mat pccBE = zeros<mat>(nw,nw);
    mat pccCF = zeros<mat>(nw,nw);

    /*
     * Heisenberg picture 
     */
    for (int iw1=0; iw1<nw; ++iw1) {
        // 1st
        blockdeom d3(ham1,ham0,qmd1,qmd0,b,h,'e','g');
        cx_cube opr_t3 = zeros<cx_cube>(d3.nl,d3.nr,d3.nmax);
        d3.iniHei(opr_t3,dip_eg);
        double w1 = (iw1-nw/2)*dw;
        d3.freqSolver(opr_s3,opr_t3,deom_ci*w1+st1+sw1,'h');
        // 2nd 
        blockdeom d2(d3,ham1,ham1,qmd1,qmd1,'e','e');
        cx_cube opr_t2 = zeros<cx_cube>(d2.nl,d2.nl,d2.nmax);
        cx_cube opr_s2 = zeros<cx_cube>(d2.nl,d2.nl,d2.nmax);
        d2.oprAct(opr_s3,dip_ge,opr_t2,'r'); 
        d2.freqSolver(opr_s2,opr_t2,2*st1,'h');
        // A & B
        blockdeom dAB(d2,ham2,ham1,qmd2,qmd1,'f','e');
        cx_cube opr_AB = zeros<cx_cube>(dAB.nl,dAB.nl,dAB.nmax);
        dAB.oprAct(opr_s2,dip_fe,opr_AB,'l'); 
        for (int iw2=0; iw2<nw; ++iw2) {
            blockdeom dAB_1(dAB,ham2,ham1,qmd2,qmd1,'f','e');
            cx_cube opr_AB_1 = zeros<cx_cube>(dAB_1.nl,dAB_1.nr,dAB_1.nmax);
            double w2 = (iw2-nw/2)*dw;
            dAB_1.freqSolver(opr_AB_1,opr_AB,deom_ci*w1+2*st1+st2+sw2,'h');
            field<ivec> keys(dAB_1.nddo);
            for (int i=0; i<dAB_1.nddo; ++i) {
                keys(i) = dAB_1.keys(i).key;
            }
            const cx_cube& oprs = opr_AB_1.slices(0,dAB_1.nddo-1);
            stringstream ss1, ss2;
            ss1 << "key_iw1_" << iw1 << "_iw2_" << iw2 << ".AB";
            ss2 << "opr_iw1_" << iw1 << "_iw2_" << iw2 << ".AB";
            keys.save(ss1.str(),arma_binary);
            oprs.save(ss2.str(),arma_binary);
        }
        // D & E
        blockdeom dDE(d2,ham1,ham2,qmd1,qmd2,'e','f');
        cx_cube opr_DE = zeros<cx_cube>(dDE.nl,dDE.nl,dDE.nmax);
        dDE.oprAct(opr_s2,dip_ef,opr_DE,'r'); 
        // D
        for (int iw2=0; iw2<nw; ++iw2) {
            blockdeom dD_1(dD,ham1,ham2,qmd1,qmd2,'e','f');
            cx_cube opr_D_1 = zeros<cx_cube>(dD_1.nl,dD_1.nr,dD_1.nmax);
            double w2 = (iw2-nw/2)*dw;
            dD_1.freqSolver(opr_D_1,opr_DE,-deom_ci*w2+2*st1+st2+sw2,'h');
            field<ivec> keys(D_1.nddo);
            for (int i=0; i<dD_1.nddo; ++i) {
                keys(i) = dD_1.keys(i).key;
            }
            const cx_cube& oprs = opr_D_1.slices(0,dD_1.nddo-1);
            stringstream ss1, ss2;
            ss1 << "key_iw1_" << iw1 << "_iw2_" << iw2 << ".D";
            ss2 << "opr_iw1_" << iw1 << "_iw2_" << iw2 << ".D";
            keys.save(ss1.str(),arma_binary);
            oprs.save(ss2.str(),arma_binary);
        }
        // E
        for (int it=0; it<nt; ++it) {
            field<ivec> keys(dDE.nddo);
            for (int i=0; i<dDE.nddo; ++i) {
                keys(i) = dDE.keys(i).key;
            }
            const cx_cube& oprs = opr_DE.slices(0,dDE.nddo-1);
            stringstream ss1, ss2;
            ss1 << "key_iw1_" << iw1 << "_it1_" << it << ".E";
            ss2 << "opr_iw1_" << iw1 << "_it1_" << it << ".E";
            keys.save(ss1.str(),arma_binary);
            oprs.save(ss2.str(),arma_binary);
            for (int jt=0; jt<mt; ++jt) {
                dDE.rk4(opr_DE,t,dt,'h');
            }
            dDE.rk4(opr_DE,t,dt_res,'h');
        }
        // C & F
        for (int it=0; it<nt; ++it) {
            field<ivec> keys(d2.nddo);
            for (int i=0; i<d2.nddo; ++i) {
                keys(i) = d2.keys(i).key;
            }
            const cx_cube& oprs = opr_s2.slices(0,d2.nddo-1);
            stringstream ss1, ss2;
            ss1 << "key_iw1_" << iw1 << "_it1_" << it << ".CF";
            ss2 << "opr_iw1_" << iw1 << "_it1_" << it << ".CF";
            keys.save(ss1.str(),arma_binary);
            oprs.save(ss2.str(),arma_binary);
            for (int jt=0; jt<mt; ++jt) {
                d2.rk4(opr_s2,t,dt,'h');
            }
            d2.rk4(opr_s2,t,dt_res,'h');
        }
    }

    /*
     * Schrodinger picture
     */
    blockdeom d0_A(ham2,ham2,qmd2,qmd2,b,h,'f','f');
    cx_cube rho_t0_A = zeros<cx_cube>(d0_A.nl,d0_A.nr,d0_A.nmax);
    rho_t0_A.slice(0) = rho0;

    // propagate t2
    for (double ta=0.0; ta<t2; ta+=dt) {
        d0_A.rk4 (rho_t0_A,ta,dt);
    } 
    
    // prepate D
    blockdeom d0_D(d0_A,ham2,ham2,qmd2,qmd2,'f','f');
    cx_cube  rho_s0_D = zeros<cx_cube>(d0_D.nl,d0_D.nl,d0_D.nmax);
    d0_D.freqSolver(rho_s0_D,rho_t0,2*(st1+st2),'s');
    blockdeom d1_D(d0_D,ham2,ham1,qmd2,qmd1,'f','e');
    cx_cube  rho_t1_D = zeros<cx_cube>(d1_D.nl,d1_D.nl,d1_D.nmax);
    d1_D.oprAct(rho_s0_D,dip_fe,rho_t1_D,'r'); 

    // Loop ta
    double ta = 0.0;
    for (it1=0; it1<nt; ++it1) {
        // B & C
        blockdeom d1_BC(d0_A,ham1,ham2,qmd1,qmd2,'e','f');
        cx_cube  rho_t1_BC = zeros<cx_cube>(d1_BC.nl,d1_BC.nl,d1_BC.nmax);
        d1_BC.oprAct(rho_t0_A,dip_ef,rho_t1_BC,'l'); 
        // E
        blockdeom d0_E(d0_A,ham2,ham2,qmd2,qmd2,'f','f');
        cx_cube  rho_s0_E = zeros<cx_cube>(d0_E.nl,d0_E.nl,d0_E.nmax);
        d0_E.freqSolver(rho_s0_E,rho_t0,2*(st1+st2),'s');
        // F
        blockdeom d1_F(d0_A,ham2,ham1,qmd2,qmd1,'f','e');
        cx_cube  rho_t1_F = zeros<cx_cube>(d1_F.nl,d1_F.nl,d1_F.nmax);
        d1_F.oprAct(rho_t0_A,dip_fe,rho_t1_F,'r'); 
        // Loop tb
        double tb = 0.0;
        for (it2=0; it2<nt-it1; ++it2) {
            // C & F
            int it3 = nt-it1-it2;
            for (iw1=0; iw1<nw; ++iw1) {
                 stringstream ss1, ss2;
                 field<ivec>  keys;
                 cx_cube oprs;
                 //
                 ss1 << "key_w1_" << iw1 << "_it_" << it3 << ".CF";
                 ss2 << "opr_w1_" << iw1 << "_it_" << it3 << ".CF";
                 keys.load(ss1.str(),arma_binary);
                 oprs.load(ss2.str(),arma_binary);
                 cx_double r1=0, r2=0;
                 for (int iado=0; iado<keys.n_rows; ++iado) {
                     auto p1 = d1_BC.tree.find(keys(iado));
                     if (p1 && p1->rank>=0) {
                         const int jado = p1->rank;
                         r1 += trace(oprs.slice(iado)*rho_t1_BC.slice(jado)*dip_fe);
                     }
                     A
                     auto p2 = d1_F.tree.find(keys(iado));
                     if (p2 && p2->rank>=0) {
                         const int jado = p2->rank;
                         r2 += trace(oprs.slice(iado)*dip_ef*rho_t1_F.slice(jado));
                     }
                 }
                 double tmp = -(st2+sw2)*tb-2*st2*ta;
                 for (int iw2=0; iw2<nw; ++iw2) {
                     double w2 = (iw2-nw/2)*dw;
                     pccCF(iw1,iw2) += real(r1*exp(tmp-deom_ci*w2*tb)+r2*exp(tmp+deom_ci*w2*tb))
                 }
            }
            double tb = it2*dt_int;
            for (int jt=0; jt<mt; ++jt) {
                d1_BC.rk4(rho_t1_BC,tb,dt,'s');
                d1_F.rk4(rho_t1_F,tb,dt,'s');
                tb += dt;
            }
            d1_BC.rk4(rho_t1_BC,tb,dt_res,'s');
            d1_F.rk4(rho_t1_F,tb,dt_res,'s');
            tb += dt_res;
        }
        // B & E
        for (iw1=0; iw1<nw; ++iw1) {
            for (iw2=0; iw2<nw; ++iw2) {
                 double w2 = (iw2-nw/2)*dw;
                 stringstream ss1, ss2;
                 field<ivec>  keys;
                 cx_cube oprs;
                 //
                 ss1 << "key_w1_" << iw1 << "_w2_" << iw2 << ".AB";
                 ss2 << "opr_w1_" << iw1 << "_w2_" << iw2 << ".AB";
                 keys.load(ss1.str(),arma_binary);
                 oprs.load(ss2.str(),arma_binary);
                 cx_double r1=0;
                 for (int iado=0; iado<keys.n_rows; ++iado) {
                     auto p1 = d1_BC.tree.find(keys(iado));
                     if (p1 && p1->rank>=0) {
                         const int jado = p1->rank;
                         r1 += trace(oprs.slice(iado)*rho_t1_BC.slice(jado));
                     }
                 }
                 pccBE(iw1,iw2) += real(r1*exp(-(deom_ci*w2+st2+sw2)*(t1-t2-ta)-2*st2*ta))
            }
            stringstream ss1, ss2;
            field<ivec>  keys;
            cx_cube oprs;
            //
            ss1 << "key_w1_" << iw1 << "_it_" << it1 << ".E";
            ss2 << "opr_w1_" << iw1 << "_it_" << it1 << ".E";
            keys.load(ss1.str(),arma_binary);
            oprs.load(ss2.str(),arma_binary);
            cx_double r1=0;
            for (int iado=0; iado<keys.n_rows; ++iado) {
                auto p1 = d1_E.tree.find(keys(iado));
                if (p1 && p1->rank>=0) {
                    const int jado = p1->rank;
                    r1 += trace(oprs.slice(iado)*rho_t1_E.slice(jado)*dip_fe);
                }
            }
            for (iw2=0; iw2<nw; ++iw2) {
                 double w2 = (iw2-nw/2)*dw;
                 pccBE(iw1,iw2) += real(r1*exp(-(deom_ci*w2+st2-sw2)*ta));
            }
        }
        // prop A & D
        ta = it1*dt_int;
        for (int jt=0; jt<mt; ++jt) {
            d0_A.rk4(rho_t0_A,ta,dt,'s');
            d1_D.rk4(rho_t1_D,ta,dt,'s');
            ta += dt;
        }
        d0_A.rk4(rho_t0_A,ta,dt_res,'s');
        d1_D.rk4(rho_t1_D,ta,dt_res,'s');
        ta += dt_res;
    }

    // A
    cx_cube rho_s1_A = zeros<cx_cube>(d0_A.nl,d0_A.nr,d0_A.nmax);
    d0_A.freqSolver(rho_s0_A,rho_t0_A,2*(st1+st2),'s');
    // A & D
    for (iw1=0; iw1<nw; ++iw1) {
        for (iw2=0; iw2<nw; ++iw2) {
             stringstream ss1, ss2;
             field<ivec>  keys;
             cx_cube oprs;
             //
             ss1 << "key_w1_" << iw1 << "_w2_" << iw2 << ".AB";
             ss2 << "opr_w1_" << iw1 << "_w2_" << iw2 << ".AB";
             keys.load(ss1.str(),arma_binary);
             oprs.load(ss2.str(),arma_binary);
             cx_double r1=0;
             for (int iado=0; iado<keys.n_rows; ++iado) {
                 auto p1 = d0_A.tree.find(keys(iado));
                 if (p1 && p1->rank>=0) {
                     const int jado = p1->rank;
                     r1 += trace(oprs.slice(iado)*dip_ef*rho_s1_A.slice(jado));
                 }
             }
             //
             ss1 << "key_w1_" << iw1 << "_w2_" << iw2 << ".D";
             ss2 << "opr_w1_" << iw1 << "_w2_" << iw2 << ".D";
             keys.load(ss1.str(),arma_binary);
             oprs.load(ss2.str(),arma_binary);
             cx_double r2=0;
             for (int iado=0; iado<keys.n_rows; ++iado) {
                 auto p1 = d0_A.tree.find(keys(iado));
                 if (p1 && p1->rank>=0) {
                     const int jado = p1->rank;
                     r2 += trace(oprs.slice(iado)*rho_t1_D.slice(jado));
                 }
             }
             pccAD(iw1,iw2) = real(r1*exp(-2*st2*(t1-t2))+r2*exp((deom_ci*w2-st2-sw2)*(t1-t2)));
        }
    }

    // clean tmp files
    for (int iw1=0; iw1<nw; ++iw1) {
        for (int it1=0; it1<nt; ++it1) {
            stringstream ss1, ss2;
            ss1 << "key_iw1_" << iw1 << "_it1_" << it1 << ".E";
            ss2 << "opr_iw1_" << iw1 << "_it1_" << it1 << ".E";
            remove (ss1.str().c_str());
            remove (ss2.str().c_str());
            ss1 << "key_iw1_" << iw1 << "_it1_" << it1 << ".CF";
            ss2 << "opr_iw1_" << iw1 << "_it1_" << it1 << ".CF";
            remove (ss1.str().c_str());
            remove (ss2.str().c_str());
        }
        for (int iw2=0; iw2<nw; ++iw2) {
            stringstream ss1, ss2;
            ss1 << "key_iw1_" << iw1 << "_iw2_" << iw2 << ".AB";
            ss2 << "opr_iw1_" << iw1 << "_iw2_" << iw2 << ".AB";
            remove (ss1.str().c_str());
            remove (ss2.str().c_str());
            ss1 << "key_iw1_" << iw1 << "_iw2_" << iw2 << ".D";
            ss2 << "opr_iw1_" << iw1 << "_iw2_" << iw2 << ".D";
            remove (ss1.str().c_str());
            remove (ss2.str().c_str());
        }
    }
    
    // write to file
    mat pcc = pccAD/(dt_int*dt_int)+pccBE/dt_int+pccCF;
    pcc.save("pcc.mat",raw_ascii);
}
