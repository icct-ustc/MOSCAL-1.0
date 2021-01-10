/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */

#include <cmath>
#include <sstream>
#include "blockdeom.hpp"

void corr1stLA (const double w1_max, const int nt1, const double dt,
              const double staticErr, const int nk,
              const mat& sdip, const cube& pdip, const vec& bdip,
              const char& sch_hei, const syst& s, const bath &b, const hidx& h) {

    const double dw1 = w1_max/nt1;
    const double dt1 = 2.0*deom_pi/w1_max;
    const int    mt1 = floor(dt1/dt);
    const double dt1_res = dt1-dt*mt1;    
    
    cx_vec ft = zeros<cx_vec>(nt1);  

    blockdeom d1_eg(s,b,h,'e','g');
        
    const mat&  sdip_eg = d1_eg.get_d10(sdip);
    const mat&  sdip_ge = d1_eg.get_d01(sdip);
    const cube& pdip_eg = d1_eg.get_d10(pdip);
    const cube& pdip_ge = d1_eg.get_d01(pdip);

    cx_cube rho_t1_eg = zeros<cx_cube>(d1_eg.nl,d1_eg.nr,d1_eg.nmax);
    d1_eg.iniSch(rho_t1_eg,sdip_eg,pdip_eg,bdip,'l');

    if (sch_hei == 's') { // sch-pic
        
        for (int it1=0; it1<nt1; ++it1) {
  
            ft(it1) = d1_eg.Trace(sdip_ge,pdip_ge,bdip,rho_t1_eg);
            
            printf("In sch-pic:eg it1=%d, nddo=%d, lddo=%d\n", it1, d1_eg.nddo, d1_eg.lddo);
            
            double t1 = it1*dt1;
            for (int jt1=0; jt1<mt1; ++jt1) {
                d1_eg.rk4 (rho_t1_eg,t1,dt,'s');
                t1 += dt;
            }
            d1_eg.rk4 (rho_t1_eg,t1,dt1_res,'s');
            t1 += dt1_res;
        }
        
    } else if (sch_hei == 'h') { // hei-pic
    
        blockdeom d1_ge(s,b,h,'g','e');

        cx_cube opr_t1_ge = zeros<cx_cube>(d1_ge.nl,d1_ge.nr,d1_ge.nmax);

        d1_ge.iniHei(opr_t1_ge,sdip_ge,pdip_ge,bdip);

        for (int it1=0; it1<nt1; ++it1) {
            
            cx_double ctmp = 0;
            for (int iado=0; iado<d1_eg.nddo; ++iado) {
                auto p = d1_ge.tree.find(d1_eg.keys(iado).key);
                if (p && p->rank>=0) {
                    const int jado = p->rank;
                    ctmp += trace(opr_t1_ge.slice(jado)*rho_t1_eg.slice(iado));
                }
            }
            
            ft(it1) = ctmp;
            
            printf("In hei-pic:ge it1=%d, nddo=%d, lddo=%d\n", it1, d1_ge.nddo, d1_ge.lddo);

            double t1 = it1*dt1;
            for (int jt1=0; jt1<mt1; ++jt1) {
                d1_ge.rk4 (opr_t1_ge,t1,dt,'h');
                t1 += dt;
            }
            d1_ge.rk4 (opr_t1_ge,t1,dt1_res,'h');
        }       
    }

    // write time-domain signal
    const vec& ft_t1 = linspace(0.0,dt1*(nt1-1)/deom_fs2unit,nt1);
    ft_t1.save("la-t1",raw_ascii);
    const vec& ft_re = real(ft);
    const vec& ft_im = imag(ft);
    ft_re.save("la-re.t1",raw_ascii);
    ft_im.save("la-im.t1",raw_ascii);
    
    // write freq-domain signal
    const vec& fw_w1 = linspace(dw1*(-nt1/2)/deom_cm2unit,dw1*(nt1/2-1)/deom_cm2unit,nt1);
    fw_w1.save("la-w1",raw_ascii);

    ft(0) *= 0.5;
    const cx_vec& fw = shift(ifft(ft),nt1/2)*nt1*dt1;
    const vec& fw_re = real(fw);
    const vec& fw_im = imag(fw);
    fw_re.save("la-re.w1",raw_ascii);
    fw_im.save("la-im.w1",raw_ascii);
}
