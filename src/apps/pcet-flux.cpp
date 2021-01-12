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

  deom d(json["deom"]);

  const int inistate = json["rhot"]["inistate"].int_value();
  const int nt = json["rhot"]["nt"].int_value();
  const double dt = json["rhot"]["dt"].number_value();
  const int nk = json["rhot"]["nk"].int_value();

  cx_cube ddos = zeros<cx_cube>(d.nsys, d.nsys, d.nmax);
  ddos(inistate, inistate, 0) = 1.0;

  const double ve = real(d.ham1(0, 1));
  const double vp = real(d.ham1(0, 2));
  const double vc = real(d.ham1(0, 3));

  double phi0 = 0.0;
  double phi1 = 0.0;
  double phi2 = 0.0;

  double tmp0 = 0.0;
  double tmp1 = 0.0;
  double tmp2 = 0.0;
  FILE *frhot = fopen("prop-rhot.dat", "w");
  FILE *fflux = fopen("prop-flux.dat", "w");
  for (int it = 0; it < nt; ++it) {
    const double t = it * dt;
    if (it % nk == 0) {
      printf("Propagation %5.1f%%: nddo=%6d, lddo=%3d\n",
             100 * it / static_cast<double>(nt), d.nddo, d.lddo);
      fprintf(frhot, "%16.6e", t / deom_fs2unit);
      for (int i = 0; i < d.nsys; ++i) {
        fprintf(frhot, "%20.10e", real(ddos(i, i, 0)));
      }
      fprintf(frhot, "\n");
    }

    tmp0 = 2.0 * vc * imag(ddos(0, 3, 0)) * dt;
    tmp1 = 2.0 * vp * imag(ddos(1, 3, 0)) * dt;
    tmp2 = 2.0 * ve * imag(ddos(2, 3, 0)) * dt;
    phi0 += tmp0;
    phi1 += tmp1;
    phi2 += tmp2;
    fprintf(fflux, "%16.6e%20.10e%20.10e%20.10e\n", t / deom_fs2unit, tmp0, tmp1, tmp2);

    d.rk4(ddos, t, dt);
  }
  fclose(fflux);
  fclose(frhot);

  phi0 -= tmp0 * nt;
  phi1 -= tmp1 * nt;
  phi2 -= tmp2 * nt;

  printf("PCET Flux:\n");
  printf("%16.6e%16.6e%16.6e\n", phi0, phi1, phi2);
  printf("kappa=%16.6e\n", phi0 / (phi1 + phi2));

  return 0;
}
