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

  //const mat& exph= expmat(-real(d1.ham1)/d1.temperature);
  //rho0.slice(0).set_real(exph/trace(exph));
  rho0(0, 0, 0) = 1;
  //rho0.slice(0)(0,0)=1;

  cx_vec ft = zeros<cx_vec>(nt);

  FILE *log_t = fopen("At1d-corr.dat", "w");
  for (int it = 0; it < nt; ++it) {
    double t = it * dt;
    // cout << "hello\n";exit(0);
    ft(it) = d1.Trace(sdip, pdip, bdip, rho0);
    printf("In 1d-correlation: it=%d, nddo=%d, lddo=%d\n", it, d1.nddo, d1.lddo);
    fprintf(log_t, "%16.6e%16.6e%16.6e\n", t / deom_fs2unit, real(ft(it)), imag(ft(it)));
    d1.rk4(rho0, t, dt);
  }
  fclose(log_t);
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
