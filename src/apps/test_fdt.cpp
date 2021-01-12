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

  const int nt = json["rhot"]["nt"].int_value();
  const double dt = json["rhot"]["dt"].number_value();
  const double staticErr = json["rhot"]["staticErr"].number_value();
  const int nk = json["rhot"]["nk"].int_value();
  const string sdipFile = json["rhot"]["sdipFile"].string_value();
  const string pdipFile = json["rhot"]["pdipFile"].string_value();
  const string bdipFile = json["rhot"]["bdipFile"].string_value();
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

  cx_cube rho0 = zeros<cx_cube>(d1.nsys, d1.nsys, d1.nmax);
  cx_cube rho1 = zeros<cx_cube>(d1.nsys, d1.nsys, d1.nmax);
  cx_cube rho2 = zeros<cx_cube>(d1.nsys, d1.nsys, d1.nmax);

  const mat &exph = expmat(-real(d1.ham1) / d1.temperature);
  rho0.slice(0).set_real(exph / trace(exph));
  d1.equilibrium(rho0, dt, staticErr, nk);

  printf("Equilibrium pop:");
  for (int i = 0; i < d1.nsys; ++i) {
    printf("%16.6e", real(rho0(i, i, 0)));
  }
  printf("\n");

  double qavg = real(d1.Trace(sdip, rho0));
  double favg = real(d1.Trace(zeros<mat>(d1.nsys, d1.nsys), pdip, bdip, rho0));

  cx_vec ft_qq = zeros<cx_vec>(nt);
  cx_vec ft_fq = zeros<cx_vec>(nt);
  d1.oprAct(rho1, sdip, rho0, 'c');
  double qq = real(d1.Trace(sdip, rho2));
  double fq = real(d1.Trace(zeros<mat>(d1.nsys, d1.nsys), pdip, bdip, rho2));
  FILE *fs = fopen("Rt-qq-fq.dat", "w");
  for (int it = 0; it < nt; ++it) {
    double t = it * dt;
    ft_qq(it) = deom_ci * d1.Trace(sdip, rho1);
    ft_fq(it) = deom_ci * d1.Trace(zeros<mat>(d1.nsys, d1.nsys), pdip, bdip, rho1);
    if (it % nk == 0) {
      printf("d1: it=%d, nddo=%d, lddo=%d\n", it, d1.nddo, d1.lddo);
      fprintf(fs, "%16.6e%16.6e%16.6e%16.6e%16.6e\n", t,
              real(ft_qq(it)), imag(ft_qq(it)),
              real(ft_fq(it)), imag(ft_fq(it)));
    }
    d1.rk4(rho1, t, dt);
  }
  fclose(fs);

  // 1D FFT
  ft_qq(0) *= 0.5;
  ft_fq(0) *= 0.5;
  const double dw = 2.0 * deom_pi / (nt * dt);
  const cx_vec &fw_qq = ifft(ft_qq) * nt * dt;
  const cx_vec &fw_fq = ifft(ft_fq) * nt * dt;

  double beta = 1.0 / d1.temperature;
  double qq_int = 0;
  double fq_int = 0;
  FILE *log_w = fopen("S-qq-fq.w", "w");
  for (int iw = nt / 2; iw < nt; ++iw) {
    double w = (iw - nt) * dw / deom_cm2unit;
    fprintf(log_w, "%16.6e%16.6e%16.6e%16.6e%16.6e\n", w,
            real(fw_qq(iw)), imag(fw_qq(iw)),
            real(fw_fq(iw)), imag(fw_fq(iw)));
  }
  for (int iw = 0; iw < nt / 2; ++iw) {
    double w = iw * dw / deom_cm2unit;
    if (iw != 0) {
      qq_int += imag(fw_qq(iw)) / tanh(0.5 * beta * w);
      fq_int += imag(fw_fq(iw)) / tanh(0.5 * beta * w);
    }
    fprintf(log_w, "%16.6e%16.6e%16.6e%16.6e%16.6e\n", w,
            real(fw_qq(iw)), imag(fw_qq(iw)),
            real(fw_fq(iw)), imag(fw_fq(iw)));
  }
  fclose(log_w);
  qq_int *= dw / deom_pi;
  fq_int *= dw / deom_pi;
  qq_int += qavg * qavg;
  fq_int += favg * qavg;

  FILE *f_fdt = fopen("fdt.comp", "w");
  fprintf(f_fdt, "qq = %16.6e\n", qq);
  fprintf(f_fdt, "fq = %16.6e\n", fq);
  fprintf(f_fdt, "qq_int = %16.6e\n", qq_int);
  fprintf(f_fdt, "fq_int = %16.6e\n", fq_int);
  fclose(f_fdt);

  return 0;
}
