/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */

#include "blockdeom.hpp"
#include <cmath>
#include <sstream>

void corr3rdPP_monomer(const double w1_max, const double t2_max, const double w3_max,
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

  cx_cube R1 = zeros<cx_cube>(nt3, nt1, nt2);
  cx_cube R2 = zeros<cx_cube>(nt3, nt1, nt2);
  cx_cube R3 = zeros<cx_cube>(nt3, nt1, nt2);
  cx_cube R4 = zeros<cx_cube>(nt3, nt1, nt2);

  blockdeom d1_eg(s, b, h, 'e', 'g');
  blockdeom d1_ge(s, b, h, 'g', 'e');

  const mat &sdip_eg = d1_eg.get_d10(sdip);
  const mat &sdip_ge = d1_eg.get_d01(sdip);
  const cube &pdip_eg = d1_eg.get_d10(pdip);
  const cube &pdip_ge = d1_eg.get_d01(pdip);

  cx_cube rho_t1_eg = zeros<cx_cube>(d1_eg.nl, d1_eg.nr, d1_eg.nmax);
  cx_cube rho_t1_ge = zeros<cx_cube>(d1_ge.nl, d1_ge.nr, d1_ge.nmax);
  d1_eg.iniSch(rho_t1_eg, sdip_eg, pdip_eg, bdip, 'l');
  d1_ge.iniSch(rho_t1_ge, sdip_ge, pdip_ge, bdip, 'r');

  if (sch_hei == 's') { // sch-pic

    for (int it1 = 0; it1 < nt1; ++it1) {

      blockdeom d2_ee1(d1_eg, 'e', 'e');
      blockdeom d2_ee2(d1_ge, 'e', 'e');
      blockdeom d2_gg2(d1_ge, 'g', 'g');
      blockdeom d2_gg1(d1_eg, 'g', 'g');

      cx_cube rho_t2_ee1 = zeros<cx_cube>(d2_ee1.nl, d2_ee1.nr, d2_ee1.nmax);
      cx_cube rho_t2_ee2 = zeros<cx_cube>(d2_ee2.nl, d2_ee2.nr, d2_ee2.nmax);
      cx_cube rho_t2_gg2 = zeros<cx_cube>(d2_gg2.nl, d2_gg2.nr, d2_gg2.nmax);
      cx_cube rho_t2_gg1 = zeros<cx_cube>(d2_gg1.nl, d2_gg1.nr, d2_gg1.nmax);

      d2_ee1.oprAct(rho_t2_ee1, sdip_ge, pdip_ge, bdip, rho_t1_eg, 'r');
      d2_ee2.oprAct(rho_t2_ee2, sdip_eg, pdip_eg, bdip, rho_t1_ge, 'l');
      d2_gg2.oprAct(rho_t2_gg2, sdip_eg, pdip_eg, bdip, rho_t1_ge, 'r');
      d2_gg1.oprAct(rho_t2_gg1, sdip_ge, pdip_ge, bdip, rho_t1_eg, 'l');

      for (int it2 = 0; it2 < nt2; ++it2) {

        blockdeom d3_r1(d2_ee1, 'e', 'g');
        blockdeom d3_r2(d2_ee2, 'e', 'g');
        blockdeom d3_r3(d2_gg2, 'e', 'g');
        blockdeom d3_r4(d2_gg1, 'e', 'g');

        cx_cube rho_t3_r1 = zeros<cx_cube>(d3_r1.nl, d3_r1.nr, d3_r1.nmax);
        cx_cube rho_t3_r2 = zeros<cx_cube>(d3_r2.nl, d3_r2.nr, d3_r2.nmax);
        cx_cube rho_t3_r3 = zeros<cx_cube>(d3_r3.nl, d3_r3.nr, d3_r3.nmax);
        cx_cube rho_t3_r4 = zeros<cx_cube>(d3_r4.nl, d3_r4.nr, d3_r4.nmax);

        d3_r1.oprAct(rho_t3_r1, sdip_eg, pdip_eg, bdip, rho_t2_ee1, 'r');
        d3_r2.oprAct(rho_t3_r2, sdip_eg, pdip_eg, bdip, rho_t2_ee2, 'r');
        d3_r3.oprAct(rho_t3_r3, sdip_eg, pdip_eg, bdip, rho_t2_gg2, 'l');
        d3_r4.oprAct(rho_t3_r4, sdip_eg, pdip_eg, bdip, rho_t2_gg1, 'l');

        for (int it3 = 0; it3 < nt3; ++it3) {

          R1(it3, it1, it2) = d3_r1.Trace(sdip_ge, pdip_ge, bdip, rho_t3_r1);
          R2(it3, it1, it2) = d3_r2.Trace(sdip_ge, pdip_ge, bdip, rho_t3_r2);
          R3(it3, it1, it2) = d3_r3.Trace(sdip_ge, pdip_ge, bdip, rho_t3_r3);
          R4(it3, it1, it2) = d3_r4.Trace(sdip_ge, pdip_ge, bdip, rho_t3_r4);

          printf("In sch-pic:R1 it3=%d, nddo=%d, lddo=%d\n", it3, d3_r1.nddo, d3_r1.lddo);
          printf("In sch-pic:R2 it3=%d, nddo=%d, lddo=%d\n", it3, d3_r2.nddo, d3_r2.lddo);
          printf("In sch-pic:R3 it3=%d, nddo=%d, lddo=%d\n", it3, d3_r3.nddo, d3_r3.lddo);
          printf("In sch-pic:R4 it3=%d, nddo=%d, lddo=%d\n", it3, d3_r4.nddo, d3_r4.lddo);

          double t3 = it3 * dt3;
          for (int jt3 = 0; jt3 < mt3; ++jt3) {
            d3_r1.rk4(rho_t3_r1, t3, dt);
            d3_r2.rk4(rho_t3_r2, t3, dt);
            d3_r3.rk4(rho_t3_r3, t3, dt);
            d3_r4.rk4(rho_t3_r4, t3, dt);
            t3 += dt;
          }
          d3_r1.rk4(rho_t3_r1, t3, dt3_res);
          d3_r2.rk4(rho_t3_r2, t3, dt3_res);
          d3_r3.rk4(rho_t3_r3, t3, dt3_res);
          d3_r4.rk4(rho_t3_r4, t3, dt3_res);
          t3 += dt3_res;
        }

        printf("In sch-pic:ee1 it2=%d, nddo=%d, lddo=%d\n", it2, d2_ee1.nddo, d2_ee1.lddo);
        printf("In sch-pic:gg1 it2=%d, nddo=%d, lddo=%d\n", it2, d2_gg1.nddo, d2_gg1.lddo);
        printf("In sch-pic:ee2 it2=%d, nddo=%d, lddo=%d\n", it2, d2_ee2.nddo, d2_ee2.lddo);
        printf("In sch-pic:gg2 it2=%d, nddo=%d, lddo=%d\n", it2, d2_gg2.nddo, d2_gg2.lddo);
        double t2 = it2 * dt2;
        for (int jt2 = 0; jt2 < mt2; ++jt2) {
          d2_ee1.rk4(rho_t2_ee1, t2, dt);
          d2_gg1.rk4(rho_t2_gg1, t2, dt);
          d2_ee2.rk4(rho_t2_ee2, t2, dt);
          d2_gg2.rk4(rho_t2_gg2, t2, dt);
          t2 += dt;
        }
        d2_ee1.rk4(rho_t2_ee1, t2, dt2_res);
        d2_gg1.rk4(rho_t2_gg1, t2, dt2_res);
        d2_ee2.rk4(rho_t2_ee2, t2, dt2_res);
        d2_gg2.rk4(rho_t2_gg2, t2, dt2_res);
        t2 += dt2_res;
      }

      printf("In sch-pic:eg it1=%d, nddo=%d, lddo=%d\n", it1, d1_eg.nddo, d1_eg.lddo);
      printf("In sch-pic:ge it1=%d, nddo=%d, lddo=%d\n", it1, d1_ge.nddo, d1_ge.lddo);
      double t1 = it1 * dt1;
      for (int jt1 = 0; jt1 < mt1; ++jt1) {
        d1_eg.rk4(rho_t1_eg, t1, dt);
        d1_ge.rk4(rho_t1_ge, t1, dt);
        t1 += dt;
      }
      d1_eg.rk4(rho_t1_eg, t1, dt1_res);
      d1_ge.rk4(rho_t1_ge, t1, dt1_res);
      t1 += dt1_res;
    }

  } else if (sch_hei == 'h') { // hei-pic

    blockdeom d3_ge(s, b, h, 'g', 'e');

    cx_cube opr_t3_ge = zeros<cx_cube>(d3_ge.nl, d3_ge.nr, d3_ge.nmax);

    d3_ge.iniHei(opr_t3_ge, sdip_ge, pdip_ge, bdip);

    for (int it3 = 0; it3 < nt3; ++it3) {

      field<ivec> keys_ge(d3_ge.nddo);
      for (int i = 0; i < d3_ge.nddo; ++i) {
        keys_ge(i) = d3_ge.keys(i).key;
      }
      const cx_cube &oprs_ge = opr_t3_ge.slices(0, d3_ge.nddo - 1);
      stringstream ss1_ge, ss2_ge;
      ss1_ge << "key_t3_ge_" << it3 << ".bin";
      ss2_ge << "opr_t3_ge_" << it3 << ".bin";
      keys_ge.save(ss1_ge.str(), arma_binary);
      oprs_ge.save(ss2_ge.str(), arma_binary);

      printf("In hei-pic:ge it3=%d, nddo=%d, lddo=%d\n", it3, d3_ge.nddo, d3_ge.lddo);

      double t3 = it3 * dt3;
      for (int jt3 = 0; jt3 < mt3; ++jt3) {
        d3_ge.rk4(opr_t3_ge, t3, dt, 'h');
        t3 += dt;
      }
      d3_ge.rk4(opr_t3_ge, t3, dt3_res, 'h');
    }

    for (int it1 = 0; it1 < nt1; ++it1) {

      blockdeom d2_ee1(d1_eg, 'e', 'e');
      blockdeom d2_gg1(d1_eg, 'g', 'g');
      blockdeom d2_ee2(d1_ge, 'e', 'e');
      blockdeom d2_gg2(d1_ge, 'g', 'g');

      cx_cube rho_t2_ee1 = zeros<cx_cube>(d2_ee1.nl, d2_ee1.nr, d2_ee1.nmax);
      cx_cube rho_t2_gg1 = zeros<cx_cube>(d2_gg1.nl, d2_gg1.nr, d2_gg1.nmax);
      cx_cube rho_t2_ee2 = zeros<cx_cube>(d2_ee2.nl, d2_ee2.nr, d2_ee2.nmax);
      cx_cube rho_t2_gg2 = zeros<cx_cube>(d2_gg2.nl, d2_gg2.nr, d2_gg2.nmax);

      d2_ee1.oprAct(rho_t2_ee1, sdip_ge, pdip_ge, bdip, rho_t1_eg, 'r');
      d2_gg1.oprAct(rho_t2_gg1, sdip_ge, pdip_ge, bdip, rho_t1_eg, 'l');
      d2_ee2.oprAct(rho_t2_ee2, sdip_eg, pdip_eg, bdip, rho_t1_ge, 'l');
      d2_gg2.oprAct(rho_t2_gg2, sdip_eg, pdip_eg, bdip, rho_t1_ge, 'r');

      for (int it2 = 0; it2 < nt2; ++it2) {

        blockdeom d3_r1(d2_ee1, 'e', 'g');
        blockdeom d3_r2(d2_ee2, 'e', 'g');
        blockdeom d3_r3(d2_gg2, 'e', 'g');
        blockdeom d3_r4(d2_gg1, 'e', 'g');

        cx_cube rho_t3_r1 = zeros<cx_cube>(d3_r1.nl, d3_r1.nr, d3_r1.nmax);
        cx_cube rho_t3_r2 = zeros<cx_cube>(d3_r2.nl, d3_r2.nr, d3_r2.nmax);
        cx_cube rho_t3_r3 = zeros<cx_cube>(d3_r3.nl, d3_r3.nr, d3_r3.nmax);
        cx_cube rho_t3_r4 = zeros<cx_cube>(d3_r4.nl, d3_r4.nr, d3_r4.nmax);

        d3_r1.oprAct(rho_t3_r1, sdip_eg, pdip_eg, bdip, rho_t2_ee1, 'r');
        d3_r2.oprAct(rho_t3_r2, sdip_eg, pdip_eg, bdip, rho_t2_ee2, 'r');
        d3_r3.oprAct(rho_t3_r3, sdip_eg, pdip_eg, bdip, rho_t2_gg2, 'l');
        d3_r4.oprAct(rho_t3_r4, sdip_eg, pdip_eg, bdip, rho_t2_gg1, 'l');

        for (int it3 = 0; it3 < nt3; ++it3) {

          stringstream ss1_ge, ss2_ge;
          field<ivec> keys_ge;
          cx_cube oprs_ge;
          ss1_ge << "key_t3_ge_" << it3 << ".bin";
          ss2_ge << "opr_t3_ge_" << it3 << ".bin";
          keys_ge.load(ss1_ge.str(), arma_binary);
          oprs_ge.load(ss2_ge.str(), arma_binary);

          const int nddo_ge = keys_ge.n_rows;
          cx_double r1 = 0, r2 = 0, r3 = 0, r4 = 0;

          for (int iado = 0; iado < nddo_ge; ++iado) {
            auto p1 = d3_r1.tree.find(keys_ge(iado));
            if (p1 && p1->rank >= 0) {
              const int jado = p1->rank;
              r1 += trace(oprs_ge.slice(iado) * rho_t3_r1.slice(jado));
            }
            auto p2 = d3_r2.tree.find(keys_ge(iado));
            if (p2 && p2->rank >= 0) {
              const int jado = p2->rank;
              r2 += trace(oprs_ge.slice(iado) * rho_t3_r2.slice(jado));
            }
            auto p3 = d3_r3.tree.find(keys_ge(iado));
            if (p3 && p3->rank >= 0) {
              const int jado = p3->rank;
              r3 += trace(oprs_ge.slice(iado) * rho_t3_r3.slice(jado));
            }
            auto p4 = d3_r4.tree.find(keys_ge(iado));
            if (p4 && p4->rank >= 0) {
              const int jado = p4->rank;
              r4 += trace(oprs_ge.slice(iado) * rho_t3_r4.slice(jado));
            }
          }

          R1(it3, it1, it2) = r1;
          R2(it3, it1, it2) = r2;
          R3(it3, it1, it2) = r3;
          R4(it3, it1, it2) = r4;
        }

        printf("In sch-pic:ee1 it2=%d, nddo=%d, lddo=%d\n", it2, d2_ee1.nddo, d2_ee1.lddo);
        printf("In sch-pic:gg1 it2=%d, nddo=%d, lddo=%d\n", it2, d2_gg1.nddo, d2_gg1.lddo);
        printf("In sch-pic:ee2 it2=%d, nddo=%d, lddo=%d\n", it2, d2_ee2.nddo, d2_ee2.lddo);
        printf("In sch-pic:gg2 it2=%d, nddo=%d, lddo=%d\n", it2, d2_gg2.nddo, d2_gg2.lddo);
        double t2 = it2 * dt2;
        for (int jt2 = 0; jt2 < mt2; ++jt2) {
          d2_ee1.rk4(rho_t2_ee1, t2, dt);
          d2_gg1.rk4(rho_t2_gg1, t2, dt);
          d2_ee2.rk4(rho_t2_ee2, t2, dt);
          d2_gg2.rk4(rho_t2_gg2, t2, dt);
          t2 += dt;
        }
        d2_ee1.rk4(rho_t2_ee1, t2, dt2_res);
        d2_gg1.rk4(rho_t2_gg1, t2, dt2_res);
        d2_ee2.rk4(rho_t2_ee2, t2, dt2_res);
        d2_gg2.rk4(rho_t2_gg2, t2, dt2_res);
        t2 += dt2_res;
      }

      printf("In sch-pic:eg it1=%d, nddo=%d, lddo=%d\n", it1, d1_eg.nddo, d1_eg.lddo);
      printf("In sch-pic:ge it1=%d, nddo=%d, lddo=%d\n", it1, d1_ge.nddo, d1_ge.lddo);
      double t1 = it1 * dt1;
      for (int jt1 = 0; jt1 < mt1; ++jt1) {
        d1_eg.rk4(rho_t1_eg, t1, dt);
        d1_ge.rk4(rho_t1_ge, t1, dt);
        t1 += dt;
      }
      d1_eg.rk4(rho_t1_eg, t1, dt1_res);
      d1_ge.rk4(rho_t1_ge, t1, dt1_res);
      t1 += dt1_res;
    }

    // clean .bin files
    for (int it3 = 0; it3 < nt3; ++it3) {
      stringstream ss1_ge, ss2_ge;
      ss1_ge << "key_t3_ge_" << it3 << ".bin";
      ss2_ge << "opr_t3_ge_" << it3 << ".bin";
      remove(ss1_ge.str().c_str());
      remove(ss2_ge.str().c_str());
    }
  }

  // write time-domain signal
  const vec &ft_t1 = linspace(0.0, dt1 * (nt1 - 1) / deom_fs2unit, nt1);
  const vec &ft_t2 = linspace(0.0, dt2 * (nt2 - 1) / deom_fs2unit, nt2);
  const vec &ft_t3 = linspace(0.0, dt3 * (nt3 - 1) / deom_fs2unit, nt3);
  ft_t1.save("pp-t1", raw_ascii);
  ft_t2.save("pp-t2", raw_ascii);
  ft_t3.save("pp-t3", raw_ascii);

  // write freq-domain signal
  const vec &fw_w1 = linspace(dw1 * (-nt1 / 2) / deom_cm2unit, dw1 * (nt1 / 2 - 1) / deom_cm2unit, nt1);
  const vec &fw_w3 = linspace(dw3 * (-nt3 / 2) / deom_cm2unit, dw1 * (nt3 / 2 - 1) / deom_cm2unit, nt3);
  fw_w1.save("pp-w1", raw_ascii);
  fw_w3.save("pp-w3", raw_ascii);

  for (int it2 = 0; it2 < nt2; ++it2) {

    stringstream ss_r1_time, ss_r2_time, ss_r3_time, ss_r4_time;

    ss_r1_time << "R1_t2_" << it2 << ".t";
    ss_r2_time << "R2_t2_" << it2 << ".t";
    ss_r3_time << "R3_t2_" << it2 << ".t";
    ss_r4_time << "R4_t2_" << it2 << ".t";

    R1.slice(it2).save(ss_r1_time.str(), raw_ascii);
    R2.slice(it2).save(ss_r2_time.str(), raw_ascii);
    R3.slice(it2).save(ss_r3_time.str(), raw_ascii);
    R4.slice(it2).save(ss_r4_time.str(), raw_ascii);

    stringstream ss_r1_freq, ss_r2_freq, ss_r3_freq, ss_r4_freq;
    ss_r1_freq << "R1_t2_" << it2 << ".w";
    ss_r2_freq << "R2_t2_" << it2 << ".w";
    ss_r3_freq << "R3_t2_" << it2 << ".w";
    ss_r4_freq << "R4_t2_" << it2 << ".w";

    cx_mat cm1, cm2;
    mat dm1, dm2;

    cm1 = R1.slice(it2);
    cm1.row(0) *= 0.5;
    cm1.col(0) *= 0.5;
    dm1 = real(ifft2(cm1)) * nt1 * nt3 * dt1 * dt3;
    dm2 = shift(shift(dm1, nt3 / 2, 0), nt1 / 2, 1);
    dm2.save(ss_r1_freq.str(), raw_ascii);

    cm1 = R2.slice(it2);
    cm1.row(0) *= 0.5;
    cm1.col(0) *= 0.5;
    cm2 = ifft(cm1) * nt3 * dt3;
    dm1 = real(fft(cm2.st()).st()) * dt1;
    dm2 = shift(shift(dm1, nt3 / 2, 0), nt1 / 2, 1);
    dm2.save(ss_r2_freq.str(), raw_ascii);

    cm1 = R3.slice(it2);
    cm1.row(0) *= 0.5;
    cm1.col(0) *= 0.5;
    cm2 = ifft(cm1) * nt3 * dt3;
    dm1 = real(fft(cm2.st()).st()) * dt1;
    dm2 = shift(shift(dm1, nt3 / 2, 0), nt1 / 2, 1);
    dm2.save(ss_r3_freq.str(), raw_ascii);

    cm1 = R4.slice(it2);
    cm1.row(0) *= 0.5;
    cm1.col(0) *= 0.5;
    dm1 = real(ifft2(cm1)) * nt1 * nt3 * dt1 * dt3;
    dm2 = shift(shift(dm1, nt3 / 2, 0), nt1 / 2, 1);
    dm2.save(ss_r4_freq.str(), raw_ascii);
  }
}
