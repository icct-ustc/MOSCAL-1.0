/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "deom.hpp"
#include <cmath>
#include <sstream>


void resp3rd(const double w1_max, const double t2_max, const double w3_max,
             const int nt1, const int nt2, const int nt3, const double dt,
             const double staticErr, const int nk,
             const mat &sdip, const cube &pdip, const vec &bdip,
             const char &sch_hei, const syst &s, const bath &b, const hidx &h) {

  const double dw1 = w1_max / nt1;
  const double dw3 = w3_max / nt3;
  const double dt1 = 2.0 * deom_pi / w1_max;
  const double dt2 = t2_max / nt2;
  const double dt3 = 2.0 * deom_pi / w3_max;
  const int mt1 = floor(dt1 / dt);
  const int mt2 = floor(dt2 / dt);
  const int mt3 = floor(dt3 / dt);
  const double dt1_res = dt1 - dt * mt1;
  const double dt2_res = dt2 - dt * mt2;
  const double dt3_res = dt3 - dt * mt3;

  deom d1(s, b, h);

  cx_cube rho_t0 = zeros<cx_cube>(d1.nsys, d1.nsys, d1.nmax);

  const mat &exph = expmat(-real(d1.ham1) / d1.temperature);
  rho_t0.slice(0).set_real(exph / trace(exph));
  d1.equilibrium(rho_t0, dt, staticErr, nk);

  cx_cube rho_t1 = zeros<cx_cube>(d1.nsys, d1.nsys, d1.nmax);
  d1.oprAct(rho_t1, sdip, pdip, bdip, rho_t0, 'c');

  cx_cube ft = zeros<cx_cube>(nt3, nt1, nt2);

  if (sch_hei == 's') { // sch-pic

    for (int it1 = 0; it1 < nt1; ++it1) {

      deom d2(d1);
      cx_cube rho_t2 = zeros<cx_cube>(d2.nsys, d2.nsys, d2.nmax);
      d2.oprAct(rho_t2, sdip, pdip, bdip, rho_t1, 'c');

      for (int it2 = 0; it2 < nt2; ++it2) {

        deom d3(d2);
        cx_cube rho_t3 = zeros<cx_cube>(d3.nsys, d3.nsys, d3.nmax);
        d3.oprAct(rho_t3, sdip, pdip, bdip, rho_t2, 'c');

        for (int it3 = 0; it3 < nt3; ++it3) {

          ft(it3, it1, it2) = -deom_ci * d3.Trace(sdip, pdip, bdip, rho_t3);

          printf("In sch-pic: it3=%d, nddo=%d, lddo=%d\n", it3, d3.nddo, d3.lddo);
          double t3 = it3 * dt3;
          for (int jt3 = 0; jt3 < mt3; ++jt3) {
            d3.rk4(rho_t3, t3, dt);
            t3 += dt;
          }
          d3.rk4(rho_t3, t3, dt3_res);
          t3 += dt3_res;
        }

        printf("In sch-pic: it2=%d, nddo=%d, lddo=%d\n", it2, d2.nddo, d2.lddo);
        double t2 = it2 * dt2;
        for (int jt2 = 0; jt2 < mt2; ++jt2) {
          d2.rk4(rho_t2, t2, dt);
          t2 += dt;
        }
        d2.rk4(rho_t2, t2, dt2_res);
        t2 += dt2_res;
      }

      printf("In sch-pic: it1=%d, nddo=%d, lddo=%d\n", it1, d1.nddo, d1.lddo);
      double t1 = it1 * dt1;
      for (int jt1 = 0; jt1 < mt1; ++jt1) {
        d1.rk4(rho_t1, t1, dt);
        t1 += dt;
      }
      d1.rk4(rho_t1, t1, dt1_res);
      t1 += dt1_res;
    }

  } else if (sch_hei == 'h') { // hei-pic

    deom d4(s, b, h);
    cx_cube opr_t3 = zeros<cx_cube>(d4.nsys, d4.nsys, d4.nmax);
    d4.iniHei(opr_t3, sdip, pdip, bdip);

    for (int it3 = 0; it3 < nt3; ++it3) {

      field<ivec> keys(d4.nddo);
      for (int i = 0; i < d4.nddo; ++i) {
        keys(i) = d4.keys(i).key;
      }
      const cx_cube &oprs = opr_t3.slices(0, d4.nddo - 1);
      stringstream ss1, ss2;
      ss1 << "key_t3_" << it3 << ".tmp";
      ss2 << "opr_t3_" << it3 << ".tmp";
      keys.save(ss1.str(), arma_binary);
      oprs.save(ss2.str(), arma_binary);

      printf("In hei-pic: it3=%d, nddo=%d, lddo=%d\n", it3, d4.nddo, d4.lddo);

      double t3 = it3 * dt3;
      for (int jt3 = 0; jt3 < mt3; ++jt3) {
        d4.rk4(opr_t3, t3, dt, 'h');
        t3 += dt;
      }
      d4.rk4(opr_t3, t3, dt3_res, 'h');
    }

    for (int it1 = 0; it1 < nt1; ++it1) {

      deom d2(d1);

      cx_cube rho_t2 = zeros<cx_cube>(d2.nsys, d2.nsys, d2.nmax);
      d2.oprAct(rho_t2, sdip, pdip, bdip, rho_t1);

      for (int it2 = 0; it2 < nt2; ++it2) {

        deom d3(d2);
        cx_cube rho_t3 = zeros<cx_cube>(d3.nsys, d3.nsys, d3.nmax);
        d3.oprAct(rho_t3, sdip, pdip, bdip, rho_t2);

        for (int it3 = 0; it3 < nt3; ++it3) {

          stringstream ss1, ss2;
          field<ivec> keys;
          cx_cube oprs;
          ss1 << "key_t3_" << it3 << ".tmp";
          ss2 << "opr_t3_" << it3 << ".tmp";
          keys.load(ss1.str(), arma_binary);
          oprs.load(ss2.str(), arma_binary);

          const int nddo = keys.n_rows;
          cx_double ctmp = 0.0;
          for (int iado = 0; iado < nddo; ++iado) {
            TrieNode *p = d3.tree.find(keys(iado));
            if (p && p->rank >= 0) {
              const int jado = p->rank;
              ctmp += trace(oprs.slice(iado) * rho_t3.slice(jado));
            }
          }
          ft(it3, it1, it2) = -deom_ci * ctmp;
        }

        printf("In sch-pic: it2=%d, nddo=%d, lddo=%d\n", it2, d2.nddo, d2.lddo);
        double t2 = it2 * dt2;
        for (int jt2 = 0; jt2 < mt2; ++jt2) {
          d2.rk4(rho_t2, t2, dt);
          t2 += dt;
        }
        d2.rk4(rho_t2, t2, dt2_res);
        t2 += dt2_res;
      }

      printf("In sch-pic: it1=%d, nddo=%d, lddo=%d\n", it1, d1.nddo, d1.lddo);
      double t1 = it1 * dt1;
      for (int jt1 = 0; jt1 < mt1; ++jt1) {
        d1.rk4(rho_t1, t1, dt);
        t1 += dt;
      }
      d1.rk4(rho_t1, t1, dt1_res);
      t1 += dt1_res;
    }

    // clean tmp files
    for (int it2 = 0; it2 < nt2; ++it2) {
      stringstream ss1, ss2;
      ss1 << "key_t3_" << it2 << ".tmp";
      ss2 << "opr_t3_" << it2 << ".tmp";
      const string &s1 = ss1.str();
      const string &s2 = ss2.str();
      remove(s1.c_str());
      remove(s2.c_str());
    }
  }

  // write time-domain signal
  const vec &ft_t1 = linspace(0.0, dt1 * (nt1 - 1) / deom_fs2unit, nt1);
  const vec &ft_t2 = linspace(0.0, dt2 * (nt2 - 1) / deom_fs2unit, nt2);
  const vec &ft_t3 = linspace(0.0, dt2 * (nt3 - 1) / deom_fs2unit, nt3);
  ft_t1.save("resp3rd.t1", raw_ascii);
  ft_t2.save("resp3rd.t2", raw_ascii);
  ft_t3.save("resp3rd.t3", raw_ascii);

  const vec &fw_w1 = linspace(0.0, dw1 * (nt1 / 2 - 1) / deom_cm2unit, nt1 / 2);
  const vec &fw_w3 = linspace(0.0, dw3 * (nt3 / 2 - 1) / deom_cm2unit, nt3 / 2);
  fw_w1.save("resp3rd.w1", raw_ascii);
  fw_w3.save("resp3rd.w2", raw_ascii);

  // 2D FFT
  for (int it2 = 0; it2 < nt2; ++it2) {

    stringstream ss_re_time, ss_im_time, ss_re_freq, ss_im_freq;
    ss_re_time << "resp3rd_re_t2_" << it2 << ".t";
    ss_im_time << "resp3rd_im_t2_" << it2 << ".t";
    ss_re_freq << "resp3rd_re_t2_" << it2 << ".w";
    ss_im_freq << "resp3rd_im_t2_" << it2 << ".w";

    const mat &ft_re = real(ft.slice(it2));
    const mat &ft_im = imag(ft.slice(it2));
    ft_re.save(ss_re_time.str(), raw_ascii);
    ft_im.save(ss_im_time.str(), raw_ascii);

    ft.slice(it2).row(0) *= 0.5;
    ft.slice(it2).col(0) *= 0.5;
    const cx_mat &fw = ifft2(ft.slice(it2)) * nt1 * nt3 * dt1 * dt3;

    const mat &fw_re = real(fw.submat(0, 0, nt3 / 2 - 1, nt1 / 2 - 1));
    const mat &fw_im = imag(fw.submat(0, 0, nt3 / 2 - 1, nt1 / 2 - 1));
    fw_re.save(ss_re_freq.str(), raw_ascii);
    fw_im.save(ss_im_freq.str(), raw_ascii);
  }
}
