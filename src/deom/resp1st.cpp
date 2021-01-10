/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cmath>
#include <sstream>
#include "deom.hpp"

void resp1st (const double w_max, const int nt, const double dt,
              const double staticErr, const int nk,
              const mat& sdip, const cube& pdip, const vec& bdip,
              const char& sch_hei, const syst& s, const bath& b, const hidx& h) {

    const double dw1 = w_max/nt;
    const double dt1 = 2.0*deom_pi/w_max;
    const int    mt  = floor(dt1/dt);
    const double dt1_res = dt1-dt*mt; 

    deom d1(s,b,h);

    cx_cube rho_t0 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);

    const mat& exph= expmat(-real(d1.ham1)/d1.temperature);
    rho_t0.slice(0).set_real(exph/trace(exph));
    d1.equilibrium (rho_t0,dt,staticErr,nk);

    cx_cube rho_t1 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);
    d1.oprAct(rho_t1,sdip,pdip,bdip,rho_t0,'c');

    cx_vec ft = zeros<cx_vec>(nt); 
    
    if (sch_hei == 's') { // sch-picture

        for (int it=0; it<nt; ++it) {

            ft(it) = deom_ci*d1.Trace(sdip,pdip,bdip,rho_t1);

            printf ("In sch-pic: it=%d, nddo=%d, lddo=%d\n", it, d1.nddo, d1.lddo);
            double t1 = it*dt1;
            for (int jt=0; jt<mt; ++jt) {
                d1.rk4 (rho_t1,t1,dt);
                t1 += dt;
            }
            d1.rk4 (rho_t1,t1,dt1_res);
        }

    } else if (sch_hei == 'h') { // hei-picture

        deom d2(s,b,h);
        cx_cube oprs = zeros<cx_cube>(d2.nsys,d2.nsys,d2.nmax);
        d2.iniHei(oprs,sdip,pdip,bdip);

        for (int it=0; it<nt; ++it) {

            cx_double ctmp = 0.0;   
            for (int iado=0; iado<d1.nddo; ++iado) {
                TrieNode* p = d2.tree.find(d1.keys(iado).key);
                if (p && p->rank>=0) {
                    const int jado = p->rank;
                    ctmp += trace(oprs.slice(jado)*rho_t1.slice(iado));
                }
            }
            ft(it) = deom_ci*ctmp; 

            printf ("In hei-pic: it=%d, nddo=%d, lddo=%d\n", it, d2.nddo, d2.lddo);
            double t1 = it*dt1;
            for (int jt=0; jt<mt; ++jt) {
                d2.rk4 (oprs,t1,dt,'h');
                t1 += dt;
            }
            d2.rk4 (oprs,t1,dt1_res,'h');
        }
    }

    // write time-domain signal
    const vec& ft_t1 = linspace(0.0,dt1*(nt-1)/deom_fs2unit,nt);
    const vec& ft_re = real(ft);
    const vec& ft_im = imag(ft);
    ft_t1.save("resp1st.t1",raw_ascii);
    ft_re.save("resp1st_re.t",raw_ascii);
    ft_im.save("resp1st_im.t",raw_ascii);
    
    // 1D FFT
    ft(0) *= 0.5;
    const cx_vec& fw = ifft(ft)*nt*dt1;
    
    // write freq-domain signal
    const vec& fw_w1 = linspace(0.0,dw1*(nt/2-1)/deom_cm2unit,nt/2);
    const vec& fw_re = real(fw.rows(0,nt/2-1));
    const vec& fw_im = imag(fw.rows(0,nt/2-1));
    fw_w1.save("resp1st.w1",raw_ascii);
    fw_re.save("resp1st_re.w",raw_ascii);
    fw_im.save("resp1st_im.w",raw_ascii);
}
