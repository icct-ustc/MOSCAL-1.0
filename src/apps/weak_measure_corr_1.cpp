/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "deom.hpp"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

static const int ntMax = 100000000;

int main() {

  ifstream jsonFile("input.json");
  stringstream strStream;
  strStream << jsonFile.rdbuf();
  string jsonStr = strStream.str();
  string err;

  const Json json = Json::parse(jsonStr, err);
  if (!err.empty()) {
    printf("Error in parsing input file: %s\n", err.c_str());
    return 0;
  }

  deom d1(json["deom"]);

  const int inistate = json["rhot"]["inistate"].int_value();
  const int nt = json["rhot"]["nt"].int_value();
  const int nk = json["rhot"]["nk"].int_value();
  const double dt = json["rhot"]["dt"].number_value();
  const double staticErr = json["rhot"]["staticErr"].number_value();
  const double kappa = json["rhot"]["kappa"].number_value();
  const double I0 = json["rhot"]["I0"].number_value();

  cx_cube rho1 = zeros<cx_cube>(d1.nsys, d1.nsys, d1.nmax);
  rho1(inistate, inistate, 0) = 1.0;

  d1.equilibrium(rho1, dt, staticErr, 1, "SCI2");

  int nddo_backup(d1.nddo);
  hkey keys_backup(d1.keys);
  cx_cube ddos_backup(rho1);

  int it = 0;
  FILE *frho = fopen("equilibrium.log", "w");
  do {
    double t = it * dt;
    printf("Propagate to equilibrium: step %d, nddo=%d, lddo=%d\n", it, d1.nddo, d1.lddo);
    if ((it + 1) % nk == 0) {
      fprintf(frho, "%16.6e", t / deom_fs2unit);
      for (int i = 0; i < d1.nsys; ++i) {
        fprintf(frho, "%16.6e", real(rho1(i, i, 0)));
      }
      fprintf(frho, "\n");

      if (d1.notconverged(nddo_backup, keys_backup, ddos_backup, rho1, staticErr)) {
        nddo_backup = d1.nddo;
        keys_backup.rows(0, d1.nddo - 1) = d1.keys.rows(0, d1.nddo - 1);
        ddos_backup.head_slices(d1.nddo) = rho1.head_slices(d1.nddo);
      } else {
        break;
      }
    }
    d1.rk4(rho1, t, dt, kappa, I0);
    it += 1;
  } while (it < ntMax);
  fclose(frho);

  if (it >= ntMax) {
    printf("Fail to reach equilibrium at error %16.6e\n", staticErr);
  } else {
    printf("Equilibrium reached with %d steps!\n", it);
  }

  // D_- Action
  const double nDOTrho = 2 * (real(rho1(1, 0, 0)) + imag(rho1(1, 0, 0))) + real(rho1(0, 0, 0) - rho1(1, 1, 0));
  for (int iddo = 0; iddo < d1.nddo; ++iddo) {
    const double nDOTpsi = 2 * (real(rho1(1, 0, iddo)) + imag(rho1(1, 0, iddo))) + real(rho1(0, 0, iddo) - rho1(1, 1, iddo));
    const double bar_psi = real(rho1(0, 0, iddo) + rho1(1, 1, iddo));
    // Dbar -- Action
    const double tmp0 = nDOTpsi - nDOTrho * bar_psi;
    // Dvec -- Action
    const double tmp1 = bar_psi - nDOTrho * nDOTpsi + nDOTrho * (2 * real(rho1(1, 0, iddo)) - nDOTpsi);
    const double tmp2 = bar_psi - nDOTrho * nDOTpsi + nDOTrho * (2 * imag(rho1(1, 0, iddo)) - nDOTpsi);
    const double tmp3 = bar_psi - nDOTrho * nDOTpsi + nDOTrho * (real(rho1(0, 0, iddo) - rho1(1, 1, iddo)) - nDOTpsi);
    rho1(0, 0, iddo) = 0.5 * (tmp0 + tmp3);
    rho1(0, 1, iddo) = 0.5 * (tmp1 - deom_ci * tmp2);
    rho1(1, 0, iddo) = 0.5 * (tmp1 + deom_ci * tmp2);
    rho1(1, 1, iddo) = 0.5 * (tmp0 - tmp3);
  }

  const double dt_fft = nk * dt;
  cx_vec ft = zeros<cx_vec>(nt);

  // Tr[D+ G(t)\rho(0;D-)]
  FILE *log_t = fopen("weak_measure_corr.t", "w");
  for (int it = 0; it < nt; ++it) {
    double t = it * dt_fft;
    const double nDOTrho = 2 * (real(rho1(1, 0, 0)) + imag(rho1(1, 0, 0))) + real(rho1(0, 0, 0) - rho1(1, 1, 0));
    ft(it) = (2 * kappa * I0) * (2 * kappa * I0) * nDOTrho;
    printf("Weak-measure-correlation: it=%d, nddo=%d, lddo=%d\n", it, d1.nddo, d1.lddo);
    fprintf(log_t, "%16.6e%16.6e%16.6e\n", t / deom_fs2unit, real(ft(it)), imag(ft(it)));
    for (int kt = 0; kt < nk; ++kt) {
      d1.rk4(rho1, t, dt, kappa, I0);
      t += dt;
    }
  }
  fclose(log_t);

  // 1D FFT
  ft(0) *= 0.5;
  const double dw = 2.0 * deom_pi / (nt * dt_fft);
  const cx_vec &fw = ifft(ft) * nt * dt_fft;

  FILE *log_w = fopen("weak_measure_corr.w", "w");
  for (int iw = nt / 2; iw < nt; ++iw) {
    double w = (iw - nt) * dw / deom_cm2unit;
    fprintf(log_w, "%16.6e%16.6e%16.6e\n", w, real(fw(iw)), imag(fw(iw)));
  }
  for (int iw = 0; iw < nt / 2; ++iw) {
    double w = iw * dw / deom_cm2unit;
    fprintf(log_w, "%16.6e%16.6e%16.6e\n", w, real(fw(iw)), imag(fw(iw)));
  }
  fclose(log_w);

  return 0;
}
