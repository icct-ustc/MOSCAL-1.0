/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#ifndef DEOM_H_
#define DEOM_H_

#include "armadillo"
#include "trie.hpp"

#include "deomBath.hpp"
#include "deomConst.hpp"
#include "deomHidx.hpp"
#include "deomPulse.hpp"
#include "deomSyst.hpp"


using namespace std;
using namespace arma;

class deom : public syst, public bath, public hidx {

public:
  int flag_rwa;
  cx_cube ddos1;
  cx_cube ddos2;
  cx_cube ddos3;

  deom(const Json &json) : syst(json["syst"]), bath(json["bath"]), hidx(json["hidx"]) {
    flag_rwa = json["rwa"].int_value();
    ddos1.set_size(nsys, nsys, nmax);
    ddos2.set_size(nsys, nsys, nmax);
    ddos3.set_size(nsys, nsys, nmax);
  }

  deom(const syst &s, const bath &b, const hidx &h) : syst(s), bath(b), hidx(h) {
    flag_rwa = 0;
    ddos1.set_size(nsys, nsys, nmax);
    ddos2.set_size(nsys, nsys, nmax);
    ddos3.set_size(nsys, nsys, nmax);
  }

  deom(const deom &rhs) : syst(rhs.ham1, rhs.qmd1),
                          bath(rhs.temperature, rhs.modLabel, rhs.coef_lft, rhs.coef_rht, rhs.coef_abs, rhs.expn_gam, rhs.delt_res),
                          hidx(rhs.nind, rhs.lmax, rhs.nmax, rhs.lddo, rhs.nddo, rhs.ferr, rhs.keys, rhs.tree, rhs.expn) {
    flag_rwa = rhs.flag_rwa;
    ddos1.set_size(nsys, nsys, nmax);
    ddos2.set_size(nsys, nsys, nmax);
    ddos3.set_size(nsys, nsys, nmax);
  }

  ~deom() {}

  cx_mat get_superOpr(const cx_mat &Opr);

  void oprAct(cx_cube &d_ddos, const mat &sdip, const cx_cube &ddos, const char lrc = 'l');

  void oprAct(cx_cube &d_ddos, const mat &sdip, const cube &pdip, const vec &bdip, const cx_cube &ddos, const char lrc = 'l');

  void iniHei(cx_cube &d_ddos, const mat &sdip);

  void iniHei(cx_cube &d_ddos, const mat &sdip, const cube &pdip, const vec &bdip);

  void rem_diag(cx_cube &d_ddos, const cx_cube &ddos, const double t);

  void rem_offdiag(cx_cube &d_ddos, const cx_cube &ddos, const double t);

  void remSch(cx_cube &d_ddos, const cx_cube &ddos, const double t);

  void remHei(cx_cube &d_ddos, const cx_cube &ddos, const double t);

  void remHSB(cx_cube &d_ddos, const cx_cube &ddos, const double t);

  void rem_freeze(cx_cube &d_ddos, const cx_cube &ddos, const double t);

  void rem(cx_cube &d_ddos, const cx_cube &ddos, const double t, const int &projection);

  void rem(cx_cube &d_ddos, const cx_cube &ddos, const double t, const ivec &projection);

  void rem(cx_cube &d_ddos, const cx_cube &ddos, const double t, const mat &sdip, const pulse &p, const mat &damp = zeros<mat>(0));

  void rem(cx_cube &d_ddos, const cx_cube &ddos, const double t, const mat &sdip, const cube &pdip, const vec &bdip, const pulse &p);

  void rem(cx_cube &d_ddos, const cx_cube &ddos, const double t, const mat &sdip, const cube &pdip, const vec &bdip, const pulse &p, const ivec &projection);

  void rem(cx_cube &d_ddos, const cx_cube &ddos, const double t, const char sch_hei = 's') {
    if (sch_hei == 's') {
      remSch(d_ddos, ddos, t);
    } else if (sch_hei == 'h') {
      remHei(d_ddos, ddos, t);
    } else if (sch_hei == 'c') { // exp(beta H_SB)
      remHSB(d_ddos, ddos, t);
    } else {
      printf("sch_hei is invalid!\n");
    }
  }

  void rem(cx_cube &d_ddos, const cx_cube &ddos, const double t, const double kappa, const double I0);

  void propagation(cx_cube &ddos, const double dt = 0.005, const int nt = 1000, const int nk = 10);

  void equilibrium(cx_cube &ddos, const double dt = 0.005, const double err = 2.0e-8, const int nk = 10, const string &method = string("SCI2")) {
    if (method == "SCI1") {
      EqSolverSCI1(ddos, dt, err, nk);
    } else if (method == "SCI2") {
      EqSolverSCI2(ddos, dt, err, nk);
    } else if (method == "SCI3") {
      EqSolverSCI3(ddos, dt, err, nk);
    } else if (method == "Prop") {
      EqSolverProp(ddos, dt, err, nk);
    } else if (method == "Krylov") {
      EqSolverKrylov(ddos, dt, err, nk);
    } else {
      printf("Wrong method!\n");
    }
  }

  void EqSolverSCI1(cx_cube &ddos, const double dt = 0.005, const double err = 2.0e-8, const int nk = 10);

  void EqSolverSCI2(cx_cube &ddos, const double dt = 0.005, const double err = 2.0e-8, const int nk = 10);

  void EqSolverSCI3(cx_cube &ddos, const double dt = 0.005, const double err = 2.0e-8, const int nk = 10);

  void EqSolverProp(cx_cube &ddos, const double dt = 0.005, const double err = 2.0e-8, const int nk = 10);

  void EqSolverKrylov(cx_cube &ddos, const double dt = 0.005, const double err = 2.0e-8, const int nk = 10);

  void tfqmr(cx_mat &x0, const cx_mat &b0, const cx_double decay, const int iddo, const double tol = 2.0e-10);

  void tfqmr(cx_cube &ddos, const double err = 2.0e-8);

  void bicgstab(cx_cube &ddos, const double err = 2.0e-8);

  inline bool is_valid(const cx_mat &ddo) const {
    return any(abs(vectorise(ddo)) > ferr);
  }

  void filter(cx_cube &ddos) {
    int n = 1;
    int l = 0;
    for (int iddo = 1; iddo < nddo; ++iddo) {
      TrieNode *p = tree.find(keys(iddo).key);
      if (is_valid(ddos.slice(iddo))) {
        if (n != iddo) {
          p->rank = n;
          keys(n) = keys(iddo);
          ddos.slice(n) = ddos.slice(iddo);
        }
        l = l > (p->tier) ? l : (p->tier);
        ++n;
      } else {
        keys(iddo).info = 0;
        p->rank = -9527;
      }
    }
    lddo = l;
    nddo = n;
  }

  bool notconverged(const int &nddo_backup, const hkey &keys_backup,
                    const cx_cube &ddos_backup, const cx_cube &ddos, const double &tol) {
    for (int iddo = 0; iddo < nddo_backup; ++iddo) {
      const auto *nod = tree.find(keys_backup(iddo).key);
      double maxDiff = 0.0;
      if (nod != NULL && nod->rank >= 0) {
        maxDiff = max(abs(vectorise(ddos_backup.slice(iddo) - ddos.slice(nod->rank))));
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

  template <typename... Tc>
  void spo(cx_cube &ddos, const double t, const double dt, const Tc &... args) {

    const int niter = 5;
    const double dt1 = dt / (3.0 * niter);
    const double dt2 = dt / (2.0 * niter);

    for (int k = 0; k < niter; ++k) {

      rem_diag(ddos1, ddos, t, args...);
      ddos.slices(0, nddo - 1) += ddos1.slices(0, nddo - 1) * dt1;

      int nddo1 = nddo;
      rem_offdiag(ddos1, ddos, t, args...);
      ddos.slices(0, nddo1 - 1) += ddos1.slices(0, nddo1 - 1) * dt2;
      if (nddo > nddo1) {
        ddos.slices(nddo1, nddo - 1) = ddos1.slices(nddo1, nddo - 1) * dt2;
      }

      rem_diag(ddos1, ddos, t, args...);
      ddos.slices(0, nddo - 1) += ddos1.slices(0, nddo - 1) * dt1;

      nddo1 = nddo;
      rem_offdiag(ddos1, ddos, t, args...);
      ddos.slices(0, nddo1 - 1) += ddos1.slices(0, nddo1 - 1) * dt2;
      if (nddo > nddo1) {
        ddos.slices(nddo1, nddo - 1) = ddos1.slices(nddo1, nddo - 1) * dt2;
      }

      rem_diag(ddos1, ddos, t, args...);
      ddos.slices(0, nddo - 1) += ddos1.slices(0, nddo - 1) * dt1;

      filter(ddos);
    }
  }

  template <typename... Tc>
  void rk4(cx_cube &ddos, const double t, const double dt, const Tc &... args) {

    const double dt2 = dt * 0.5;
    const double dt6 = dt / 6.0;

    // K1
    const int nddo0 = nddo;
    rem(ddos1, ddos, t, args...);
    ddos3.slices(0, nddo0 - 1) = ddos.slices(0, nddo0 - 1) + ddos1.slices(0, nddo0 - 1) * dt2;
    if (nddo > nddo0) {
      ddos3.slices(nddo0, nddo - 1) = ddos1.slices(nddo0, nddo - 1) * dt2;
    }
    // K2
    const int nddo1 = nddo;
    rem(ddos2, ddos3, t + 0.5 * dt, args...);
    ddos1.slices(0, nddo1 - 1) += ddos2.slices(0, nddo1 - 1) * 2.0;
    if (nddo > nddo1) {
      ddos1.slices(nddo1, nddo - 1) = ddos2.slices(nddo1, nddo - 1) * 2.0;
    }
    ddos3.slices(0, nddo0 - 1) = ddos.slices(0, nddo0 - 1) + ddos2.slices(0, nddo0 - 1) * dt2;
    if (nddo > nddo0) {
      ddos3.slices(nddo0, nddo - 1) = ddos2.slices(nddo0, nddo - 1) * dt2;
    }
    // K3
    const int nddo2 = nddo;
    rem(ddos2, ddos3, t + 0.5 * dt, args...);
    ddos1.slices(0, nddo2 - 1) += ddos2.slices(0, nddo2 - 1) * 2.0;
    if (nddo > nddo2) {
      ddos1.slices(nddo2, nddo - 1) = ddos2.slices(nddo2, nddo - 1) * 2.0;
    }
    ddos3.slices(0, nddo0 - 1) = ddos.slices(0, nddo0 - 1) + ddos2.slices(0, nddo0 - 1) * dt;
    if (nddo > nddo0) {
      ddos3.slices(nddo0, nddo - 1) = ddos2.slices(nddo0, nddo - 1) * dt;
    }
    // K4
    const int nddo3 = nddo;
    rem(ddos2, ddos3, t + dt, args...);
    ddos1.slices(0, nddo3 - 1) += ddos2.slices(0, nddo3 - 1);
    if (nddo > nddo3) {
      ddos1.slices(nddo3, nddo - 1) = ddos2.slices(nddo3, nddo - 1);
    }
    ddos.slices(0, nddo0 - 1) += ddos1.slices(0, nddo0 - 1) * dt6;
    if (nddo > nddo0) {
      ddos.slices(nddo0, nddo - 1) = ddos1.slices(nddo0, nddo - 1) * dt6;
    }
    filter(ddos);
  }

  template <typename... Tc>
  void rk20(cx_cube &ddos, const double t, const double dt, const Tc &... args) {
    // vec cn = {1.0/420, 1.0/209, 1.0/138, 1.0/102, 1.0/80, 1.0/65, 1.0/54, 2.0/91,
    //           3.0/116, 1.0/33, 11.0/310, 1.0/24, 13.0/264, 1.0/17, 1.0/14, 4.0/45,
    //           17.0/148, 3.0/19, 19.0/78, 1.0/2};
    // vec cn = {1.0/4, 1.0/3, 1.0/2, 1.0};
    const int nrk = 5;
    cx_cube *rhot = &ddos;
    for (int k = 0; k < nrk; ++k) {
      const double dtp = dt / (nrk - k);
      const int ntmp = nddo;
      rem(ddos1, *rhot, t, args...);
      rhot = (k < nrk - 1) ? (&ddos2) : (&ddos);
      rhot->head_slices(ntmp) = ddos.head_slices(ntmp) + ddos1.head_slices(ntmp) * dtp;
      if (nddo > ntmp) {
        rhot->slices(ntmp, nddo - 1) = ddos1.slices(ntmp, nddo - 1) * dtp;
      }
    }
    filter(ddos);
  }

  cx_double Trace(const mat &sdip, const cx_cube &ddos) const {
    return trace(sdip * ddos.slice(0));
  }

  cx_double Trace(const mat &sdip, const cube &pdip, const vec &bdip, const cx_cube &ddos) const {
    cx_double result = trace(sdip * ddos.slice(0));
    for (int mp = 0; mp < nind; ++mp) {
      const int m = modLabel(mp);
      ivec key0 = zeros<ivec>(mp + 1);
      key0(mp) = 1;
      const cx_double sn = bdip(mp) * sqrt(coef_abs(mp));
      TrieNode *p = tree.find(key0);
      if (p && p->rank >= 0) {
        int loc = p->rank;
        result += sn * trace(pdip.slice(m) * ddos.slice(loc));
      }
    }
    return result;
  }

  cx_double Trace(const mat &sdip, const cube &pdip, const cx_vec &bdip, const cx_cube &ddos) const {
    cx_double result = trace(sdip * ddos.slice(0));
    for (int mp = 0; mp < nind; ++mp) {
      const int m = modLabel(mp);
      ivec key0 = zeros<ivec>(mp + 1);
      key0(mp) = 1;
      const cx_double sn = bdip(mp) * sqrt(coef_abs(mp));
      TrieNode *p = tree.find(key0);
      if (p && p->rank >= 0) {
        int loc = p->rank;
        result += sn * trace(pdip.slice(m) * ddos.slice(loc));
      }
    }
    return result;
  }

  double entropy(const cx_cube &ddos, const string &type = string("vn")) const {
    vec popt = zeros<vec>(nsys);
    if (type == "vn") {
      eig_sym(popt, ddos.slice(0));
    } else if (type == "sh") {
      vec eval = zeros<vec>(nsys);
      cx_mat evec = zeros<cx_mat>(nsys, nsys);
      eig_sym(eval, evec, ham1);
      for (int a = 0; a < nsys; ++a) {
        popt(a) = real(as_scalar(evec.col(a).st() * ddos.slice(0) * evec.col(a)));
      }
    }
    double s = 0.0;
    for (int i = 0; i < nsys; ++i) {
      if (popt(i) < 0) {
        return -9527;
      } else if (popt(i) > 1.e-20) {
        s -= popt(i) * log(popt(i));
      }
    }
    return s;
  }

  double Umicro(const cx_cube &ddos) const {
    return real(trace(ham1 * ddos.slice(0)));
  }

  double Umacro(const cx_cube &ddos) const {
    cx_double result = trace(ham1 * ddos.slice(0));
    for (int mp = 0; mp < nind; ++mp) {
      const int m = modLabel(mp);
      ivec key0 = zeros<ivec>(mp + 1);
      key0(mp) = 1;
      const double sn = sqrt(coef_abs(mp));
      TrieNode *p = tree.find(key0);
      if (p && p->rank >= 0) {
        int loc = p->rank;
        result += sn * trace(qmd1.slice(m) * ddos.slice(loc));
      }
    }
    return real(result);
  }

  double Amicro(const cx_cube &ddos) const {
    return Umicro(ddos) - temperature * entropy(ddos);
  }

  double Amacro(const cx_cube &ddos) const {
    return Umacro(ddos) - temperature * entropy(ddos);
  }
};

void resp1st(const double w_max, const int nt, const double dt,
             const double staticErr, const int nk,
             const mat &sdip, const cube &pdip, const vec &bdip,
             const char &sch_hei, const syst &s, const bath &b, const hidx &h);

void resp2nd(const double w1_max, const double w2_max, const int nt1, const int nt2, const double dt,
             const double staticErr, const int nk,
             const mat &sdip, const cube &pdip, const vec &bdip,
             const char &sch_hei, const syst &s, const bath &b, const hidx &h);

void resp3rd(const double w1_max, const double t2_max, const double w3_max, const int nt1, const int nt2, const int nt3, const double dt,
             const double staticErr, const int nk,
             const mat &sdip, const cube &pdip, const vec &bdip,
             const char &sch_hei, const syst &s, const bath &b, const hidx &h);

void respPcc(const double tf, const int nt, const double dt,
             const cx_mat &rho0, const cube &dipo_mat, const string &dipo_lrc,
             const char &sch_hei, const syst &s, const bath &b, const hidx &h);

#endif
