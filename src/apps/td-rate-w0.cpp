/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "deom.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

static void td_rate_w0(deom &d, const double dt, const double ti, const double tf, const int nk, const mat &sdip, const cube &pdip, const vec &bdip, pulse &p, const ivec &projection) {

  const int nsys = d.nsys;
  const int nt_i = ceil(abs(ti) / dt);
  const int nt_f = ceil(abs(tf) / dt);

  double t = -(nt_i - 1) * dt;
  cx_cube rho1 = zeros<cx_cube>(nsys, nsys, d.nmax);
  cx_cube rho2 = zeros<cx_cube>(nsys, nsys, d.nmax);
  rho1(0, 0, 0) = 1.0;
  p.turnon();

  // rho propagation from ti to t0
  for (int it = 0; it < nt_i; ++it) {
    d.rk4(rho1, t, dt, sdip, pdip, bdip, p);
    t += dt;
  }

  // projection to coherence subspace
  for (int i = 1; i < nsys; ++i) {
    if (projection(i) != 0) {
      rho1(i, i, 0) = 0.0;
    }
  }

  // w_0 propagation from t0 to tf
  FILE *flog = fopen("td-rate-w0.dat", "w");
  for (int it = 0; it < nt_f; ++it) {
    d.rem(rho2, rho1, t, sdip, pdip, bdip, p);
    fprintf(flog, "%16.6e", t / deom_fs2unit);
    for (int i = 1; i < nsys; ++i) {
      fprintf(flog, "%16.6e%16.6e", real(rho2(i, i, 0)), imag(rho2(i, i, 0)));
    }
    fprintf(flog, "\n");
    if (it % nk == 0) {
      printf("%s: %5.1f%%, nddo = %6d, lddo = %2d\n", "td-rate-w0",
             100 * it / static_cast<double>(nt_f), d.nddo, d.lddo);
    }
    d.rk4(rho1, t, dt, sdip, pdip, bdip, p, projection);
    t += dt;
  }
  fclose(flog);
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

  deom d(json["deom"]);

  const double dt = json["td-rate"]["dt"].number_value();
  const double ti = json["td-rate"]["ti"].number_value();
  const double tf = json["td-rate"]["tf"].number_value();
  const int nk = json["td-rate"]["nk"].int_value();
  const string sdipFile = json["td-rate"]["sdipFile"].string_value();
  const string pdipFile = json["td-rate"]["pdipFile"].string_value();
  const string bdipFile = json["td-rate"]["bdipFile"].string_value();
  const string projFile = json["td-rate"]["projFile"].string_value();

  pulse p(json["td-rate"]["pulse"]);

  mat sdip;
  cube pdip;
  vec bdip;
  ivec proj;

  if (sdip.load(sdipFile, arma_ascii)) {
    sdip.print(sdipFile);
  } else {
    printf("Fail to load sdip");
  }

  if (pdip.load(pdipFile, arma_ascii)) {
    pdip.print(pdipFile);
  } else {
    printf("Fail to load pdip");
  }

  if (bdip.load(bdipFile, arma_ascii)) {
    bdip.print(bdipFile);
  } else {
    printf("Fail to load bdip");
  }

  if (proj.load(projFile, arma_ascii)) {
    proj.print(projFile);
  } else {
    printf("Fail to load proj");
  }

  if (proj.n_rows != sdip.n_rows) {
    printf("Error! Invalid projection operator!\n");
    return 0;
  }

  td_rate_w0(d, dt, ti, tf, nk, sdip, pdip, bdip, p, proj);

  return 0;
}
