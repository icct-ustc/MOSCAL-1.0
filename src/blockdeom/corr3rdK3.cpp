/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */

#include <cmath>
#include <sstream>
#include "blockdeom.hpp"

void corr3rdK3 (const double t1_max, const double w2_max, const double w3_max,
              const int nt1, const int nt2, const int nt3, const double dt,
              const double staticErr, const int nk,
              const mat& sdip, const cube& pdip, const vec& bdip,
              const char& sch_hei, const syst& s, const bath &b, const hidx& h) {

    const double dw2 = w2_max/nt2;
    const double dw3 = w3_max/nt3;
    const double dt1 = t1_max/nt1;
    const double dt2 = 2.0*deom_pi/w2_max;
    const double dt3 = 2.0*deom_pi/w3_max;
    const int    mt1 = floor(dt1/dt);
    const int    mt2 = floor(dt2/dt);
    const int    mt3 = floor(dt3/dt);
    const double dt1_res = dt1-dt*mt1;    
    const double dt2_res = dt2-dt*mt2;    
    const double dt3_res = dt3-dt*mt3; 
      
    cx_cube R7 = zeros<cx_cube>(nt3,nt2,nt1);
    cx_cube R8 = zeros<cx_cube>(nt3,nt2,nt1);

    blockdeom d1_eg(s,b,h,'e','g');
        
    const  mat& sdip_eg = d1_eg.get_d10(sdip);
    const  mat& sdip_ge = d1_eg.get_d01(sdip);
    const  mat& sdip_ef = d1_eg.get_d12(sdip);
    const  mat& sdip_fe = d1_eg.get_d21(sdip);
    const cube& pdip_eg = d1_eg.get_d10(pdip);
    const cube& pdip_ge = d1_eg.get_d01(pdip);
    const cube& pdip_ef = d1_eg.get_d12(pdip);
    const cube& pdip_fe = d1_eg.get_d21(pdip);

    cx_cube rho_t1_eg = zeros<cx_cube>(d1_eg.nl,d1_eg.nr,d1_eg.nmax);
    d1_eg.iniSch(rho_t1_eg,sdip_eg,pdip_eg,bdip);

    if (sch_hei == 's') { // sch-pic
        
        for (int it1=0; it1<nt1; ++it1) {
            
            blockdeom d2_fg(d1_eg,'f','g');

            cx_cube rho_t2_fg = zeros<cx_cube>(d2_fg.nl,d2_fg.nr,d2_fg.nmax);

            d2_fg.oprAct(rho_t2_fg,sdip_fe,pdip_fe,bdip,rho_t1_eg);

            for (int it2=0; it2<nt2; ++it2) {

                blockdeom d3_r7(d2_fg,'f','e');
                blockdeom d3_r8(d2_fg,'e','g');

                cx_cube rho_t3_r7 = zeros<cx_cube>(d3_r7.nl,d3_r7.nr,d3_r7.nmax);
                cx_cube rho_t3_r8 = zeros<cx_cube>(d3_r8.nl,d3_r8.nr,d3_r8.nmax);

                d3_r7.oprAct(rho_t3_r7,sdip_ge,pdip_ge,bdip,rho_t2_fg,'r');
                d3_r8.oprAct(rho_t3_r8,sdip_ef,pdip_ef,bdip,rho_t2_fg,'l');

                for (int it3=0; it3<nt3; ++it3) {

                    R7(it3,it2,it1) =-d3_r7.Trace(sdip_ef,pdip_ef,bdip,rho_t3_r7);
                    R8(it3,it2,it1) = d3_r8.Trace(sdip_ge,pdip_ge,bdip,rho_t3_r8);

                    printf("In sch-pic:R7 it3=%d, nddo=%d, lddo=%d\n", it3, d3_r7.nddo, d3_r7.lddo);
                    printf("In sch-pic:R8 it3=%d, nddo=%d, lddo=%d\n", it3, d3_r8.nddo, d3_r8.lddo);

                    double t3 = it3*dt3;
                    for (int jt3=0; jt3<mt3; ++jt3) {
                        d3_r7.rk4 (rho_t3_r7,t3,dt);
                        d3_r8.rk4 (rho_t3_r8,t3,dt);
                        t3 += dt;
                    }
                    d3_r7.rk4 (rho_t3_r7,t3,dt3_res);
                    d3_r8.rk4 (rho_t3_r8,t3,dt3_res);
                    t3 += dt3_res;
                }

                printf("In sch-pic:ee it2=%d, nddo=%d, lddo=%d\n", it2, d2_fg.nddo, d2_fg.lddo);
                double t2 = it2*dt2;                    
                for (int jt2=0; jt2<mt2; ++jt2) {
                    d2_fg.rk4 (rho_t2_fg,t2,dt);
                    t2 += dt;
                }
                d2_fg.rk4 (rho_t2_fg,t2,dt2_res);
                t2 += dt2_res;
            }
            
            printf("In sch-pic:eg it1=%d, nddo=%d, lddo=%d\n", it1, d1_eg.nddo, d1_eg.lddo);
            double t1 = it1*dt1;
            for (int jt1=0; jt1<mt1; ++jt1) {
                d1_eg.rk4 (rho_t1_eg,t1,dt);
                t1 += dt;
            }
            d1_eg.rk4 (rho_t1_eg,t1,dt1_res);
            t1 += dt1_res;
        }
        
    } else if (sch_hei == 'h') { // hei-pic
    
        blockdeom d3_ge(s,b,h,'g','e');
        blockdeom d3_ef(s,b,h,'e','f');

        cx_cube opr_t3_ge = zeros<cx_cube>(d3_ge.nl,d3_ge.nr,d3_ge.nmax);
        cx_cube opr_t3_ef = zeros<cx_cube>(d3_ef.nl,d3_ef.nr,d3_ef.nmax);

        d3_ge.iniHei(opr_t3_ge,sdip_ge,pdip_ge,bdip);
        d3_ef.iniHei(opr_t3_ef,sdip_ef,pdip_ef,bdip);

        for (int it3=0; it3<nt3; ++it3) {
            
            field<ivec> keys_ge(d3_ge.nddo);
            for (int i=0; i<d3_ge.nddo; ++i) {
                keys_ge(i) = d3_ge.keys(i).key;
            }
            const cx_cube& oprs_ge = opr_t3_ge.slices(0,d3_ge.nddo-1);
            field<ivec> keys_ef(d3_ef.nddo);
            for (int i=0; i<d3_ef.nddo; ++i) {
                keys_ef(i) = d3_ef.keys(i).key;
            }
            const cx_cube& oprs_ef = opr_t3_ef.slices(0,d3_ef.nddo-1);
            stringstream ss1_ge, ss2_ge, ss1_ef, ss2_ef;
            ss1_ge << "key_t3_ge_" << it3 << ".bin";
            ss2_ge << "opr_t3_ge_" << it3 << ".bin";
            ss1_ef << "key_t3_ef_" << it3 << ".bin";
            ss2_ef << "opr_t3_ef_" << it3 << ".bin";
            keys_ge.save(ss1_ge.str(),arma_binary);
            oprs_ge.save(ss2_ge.str(),arma_binary);
            keys_ef.save(ss1_ef.str(),arma_binary);
            oprs_ef.save(ss2_ef.str(),arma_binary);

            printf("In hei-pic:ge it3=%d, nddo=%d, lddo=%d\n", it3, d3_ge.nddo, d3_ge.lddo);
            printf("In hei-pic:ef it3=%d, nddo=%d, lddo=%d\n", it3, d3_ef.nddo, d3_ef.lddo);

            double t3 = it3*dt3;
            for (int jt3=0; jt3<mt3; ++jt3) {
                d3_ge.rk4 (opr_t3_ge,t3,dt,'h');
                d3_ef.rk4 (opr_t3_ef,t3,dt,'h');
                t3 += dt;
            }
            d3_ge.rk4 (opr_t3_ge,t3,dt3_res,'h');
            d3_ef.rk4 (opr_t3_ef,t3,dt3_res,'h');
        }
        
        for (int it1=0; it1<nt1; ++it1) {

            blockdeom d2_fg(d1_eg,'f','g');

            cx_cube rho_t2_fg = zeros<cx_cube>(d2_fg.nl,d2_fg.nr,d2_fg.nmax);

            d2_fg.oprAct(rho_t2_fg,sdip_fe,pdip_fe,bdip,rho_t1_eg);

            for (int it2=0; it2<nt2; ++it2) {

                blockdeom d3_r7(d2_fg,'f','e');
                blockdeom d3_r8(d2_fg,'e','g');

                cx_cube rho_t3_r7 = zeros<cx_cube>(d3_r7.nl,d3_r7.nr,d3_r7.nmax);
                cx_cube rho_t3_r8 = zeros<cx_cube>(d3_r8.nl,d3_r8.nr,d3_r8.nmax);

                d3_r7.oprAct(rho_t3_r7,sdip_ge,pdip_ge,bdip,rho_t2_fg,'r');
                d3_r8.oprAct(rho_t3_r8,sdip_ef,pdip_ef,bdip,rho_t2_fg,'l');

                for (int it3=0; it3<nt3; ++it3) {
                
                    stringstream ss1_ge, ss2_ge;
                    field<ivec>  keys_ge;
                    cx_cube oprs_ge;
                    ss1_ge << "key_t3_ge_" << it3 << ".bin";
                    ss2_ge << "opr_t3_ge_" << it3 << ".bin";
                    keys_ge.load(ss1_ge.str(),arma_binary);
                    oprs_ge.load(ss2_ge.str(),arma_binary);

                    stringstream ss1_ef, ss2_ef;
                    field<ivec>  keys_ef;
                    cx_cube oprs_ef;
                    ss1_ef << "key_t3_ef_" << it3 << ".bin";
                    ss2_ef << "opr_t3_ef_" << it3 << ".bin";
                    keys_ef.load(ss1_ef.str(),arma_binary);
                    oprs_ef.load(ss2_ef.str(),arma_binary);

                    const int nddo_ge = keys_ge.n_rows, nddo_ef = keys_ef.n_rows;
                    cx_double r7=0, r8=0;

                    for (int iado=0; iado<nddo_ge; ++iado) {
                        auto p = d3_r8.tree.find(keys_ge(iado));
                        if (p && p->rank>=0) {
                            const int jado = p->rank;
                            r8 += trace(oprs_ge.slice(iado)*rho_t3_r8.slice(jado));
                        }
                    }

                    for (int iado=0; iado<nddo_ef; ++iado) {
                        auto p = d3_r7.tree.find(keys_ef(iado));
                        if (p && p->rank>=0) {
                            const int jado = p->rank;
                            r7 -= trace(oprs_ef.slice(iado)*rho_t3_r7.slice(jado));
                        }
                    }

                    R7(it3,it2,it1) = r7;
                    R8(it3,it2,it1) = r8;
                }

                printf("In sch-pic:fg it2=%d, nddo=%d, lddo=%d\n", it2, d2_fg.nddo, d2_fg.lddo);
                double t2 = it2*dt2;                    
                for (int jt2=0; jt2<mt2; ++jt2) {
                    d2_fg.rk4 (rho_t2_fg,t2,dt);
                    t2 += dt;
                }
                d2_fg.rk4 (rho_t2_fg,t2,dt2_res);
                t2 += dt2_res;
            }
            
            printf("In sch-pic: it1=%d, nddo=%d, lddo=%d\n", it1, d1_eg.nddo, d1_eg.lddo);
            double t1 = it1*dt1;
            for (int jt1=0; jt1<mt1; ++jt1) {
                d1_eg.rk4 (rho_t1_eg,t1,dt);
                t1 += dt;
            }
            d1_eg.rk4 (rho_t1_eg,t1,dt1_res);
            t1 += dt1_res;
        }
       
        // clean .bin files
        for (int it3=0; it3<nt3; ++it3) {
            stringstream ss1_ge, ss2_ge, ss1_ef, ss2_ef;
            ss1_ge << "key_t3_ge_" << it3 << ".bin";
            ss2_ge << "opr_t3_ge_" << it3 << ".bin";
            ss1_ef << "key_t3_ef_" << it3 << ".bin";
            ss2_ef << "opr_t3_ef_" << it3 << ".bin";
            remove (ss1_ge.str().c_str());
            remove (ss2_ge.str().c_str());
            remove (ss1_ef.str().c_str());
            remove (ss2_ef.str().c_str());
        }
    }

    // write time-domain signal
    const vec& ft_t1 = linspace(0.0,dt1*(nt1-1)/deom_fs2unit,nt1);
    const vec& ft_t2 = linspace(0.0,dt2*(nt2-1)/deom_fs2unit,nt2);
    const vec& ft_t3 = linspace(0.0,dt3*(nt3-1)/deom_fs2unit,nt3);
    ft_t1.save("pp-t1",raw_ascii);
    ft_t2.save("pp-t2",raw_ascii);
    ft_t3.save("pp-t3",raw_ascii);

    // write freq-domain signal
    const vec& fw_w2 = linspace(dw2*(-nt2/2)/deom_cm2unit,dw2*(nt2/2-1)/deom_cm2unit,nt2);
    const vec& fw_w3 = linspace(dw3*(-nt3/2)/deom_cm2unit,dw3*(nt3/2-1)/deom_cm2unit,nt3);
    fw_w2.save("pp-w2",raw_ascii);
    fw_w3.save("pp-w3",raw_ascii);

    for (int it1=0; it1<nt1; ++it1) {

        stringstream ss_k3_time, ss_k3_freq;
        ss_k3_time << "K3_t1_" << it1 << ".t";
        ss_k3_freq << "K3_t1_" << it1 << ".w";

        cx_mat ft_k3 = R7.slice(it1)+R8.slice(it1);
        ft_k3.save(ss_k3_time.str(),raw_ascii);
    
        ft_k3.row(0) *= 0.5;
        ft_k3.col(0) *= 0.5;
        const mat& fw_k3 = real(ifft2(ft_k3))*nt2*nt3*dt2*dt3;
        const mat& s3mat = shift(fw_k3,nt3/2,0);
        const mat& out3 = shift(s3mat,nt2/2,1);
        out3.save(ss_k3_freq.str(),raw_ascii);
    }
}
