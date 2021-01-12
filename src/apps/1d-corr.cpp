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

static void corr(const int &nt, const double &dt, const double &staticErr, const int &nk,
                 const mat &sdip, const cube &pdip, const vec &bdip,
                 const syst &s, const bath &b, const hidx &h) {

  deom d1(s, b, h);

  cx_cube rho0 = zeros<cx_cube>(d1.nsys, d1.nsys, d1.nmax);
  cx_cube rho1 = zeros<cx_cube>(d1.nsys, d1.nsys, d1.nmax);

  const mat &exph = expmat(-real(d1.ham1) / d1.temperature);
  rho0.slice(0).set_real(exph / trace(exph));
  d1.equilibrium(rho0, dt, staticErr, nk);

  printf("Equilibrium pop:");
  for (int i = 0; i < d1.nsys; ++i) {
    printf("%16.6e", real(rho0(i, i, 0)));
  }
  printf("\n");

  cx_vec ft = zeros<cx_vec>(nt);

  // \mu\rho
  d1.oprAct(rho1, sdip, pdip, bdip, rho0, 'l');
  // iTr[\mu G(t)\mu\rho]
  FILE *log_t = fopen("1d-corr-t.dat", "w");
  for (int it = 0; it < nt; ++it) {
    double t = it * dt;
    ft(it) = d1.Trace(sdip, pdip, bdip, rho1);
    printf("In 1d-correlation: it=%d, nddo=%d, lddo=%d\n", it, d1.nddo, d1.lddo);
    fprintf(log_t, "%16.6e%16.6e%16.6e\n", t / deom_fs2unit, real(ft(it)), imag(ft(it)));
    d1.rk4(rho1, t, dt);
  }
  fclose(log_t);

  // 1D FFT
  ft(0) *= 0.5;
  const double dw = 2.0 * deom_pi / (nt * dt);
  const cx_vec &fw = ifft(ft) * nt * dt;

  FILE *log_w = fopen("1d-corr-w.dat", "w");
  for (int iw = nt / 2; iw < nt; ++iw) {
    double w = (iw - nt) * dw / deom_cm2unit;
    fprintf(log_w, "%16.6e%16.6e%16.6e\n", w, real(fw(iw)), imag(fw(iw)));
  }
  for (int iw = 0; iw < nt / 2; ++iw) {
    double w = iw * dw / deom_cm2unit;
    fprintf(log_w, "%16.6e%16.6e%16.6e\n", w, real(fw(iw)), imag(fw(iw)));
  }
  fclose(log_w);

  printf("fw(0) = %16.6e  +%16.6eI\n", real(fw(0)), imag(fw(0)));
}

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

  syst s(json["deom"]["syst"]);
  bath b(json["deom"]["bath"]);
  hidx h(json["deom"]["hidx"]);

  const int nt = json["spec"]["nt"].int_value();
  const double dt = json["spec"]["dt"].number_value();
  const double staticErr = json["spec"]["staticErr"].number_value();
  const int nk = json["spec"]["nk"].int_value();
  const string sdipFile = json["spec"]["sdipFile"].string_value();
  const string pdipFile = json["spec"]["pdipFile"].string_value();
  const string bdipFile = json["spec"]["bdipFile"].string_value();
  mat sdip;
  cube pdip;
  vec bdip;
  if (sdip.load(sdipFile, arma_ascii)) {
    sdip.print(sdipFile);
  } else {
    printf("Fail to load sdip!\n");
  }
  if (pdip.load(pdipFile, arma_ascii)) {
    pdip.print(pdipFile);
  } else {
    printf("Fail to load pdip!\n");
  }
  if (bdip.load(bdipFile, arma_ascii)) {
    bdip.print(bdipFile);
  } else {
    printf("Fail to load bdip!\n");
  }

  corr(nt, dt, staticErr, nk, sdip, pdip, bdip, s, b, h);

  return 0;
}
