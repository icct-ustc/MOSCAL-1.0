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


mat rhox2w(const cx_mat &rhox, const double &xi, const double &xf, const int &n) {
  mat rhow = zeros<mat>(n, n);
  const double dx = (xf - xi) / n;
  const double dp = 2 * deom_pi / (xf - xi);
  for (int k = 0; k < n; ++k) {
    const double pk = -deom_pi / dx + k * dp;
    for (int j = 0; j < n; ++j) {
      cx_double tmp = 0;
      for (int jp = 0; jp <= min(j, n - j - 1); ++jp) {
        const double sjp = dx * jp;
        tmp += rhox(j - jp, j + jp) * exp(deom_ci * pk * sjp);
      }
      rhow(k, j) = real(tmp) / n;
    }
  }
  return rhow;
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

  if (d.flag_rwa != 0) {
    printf("Rotating Wave Approximation\n");
  }
  const int inistate = json["rhot"]["inistate"].int_value();
  const int nt = json["rhot"]["nt"].int_value();
  const double dt = json["rhot"]["dt"].number_value();
  const int nk = json["rhot"]["nk"].int_value();
  const double xi = json["rhot"]["xi"].number_value();
  const double xf = json["rhot"]["xf"].number_value();
  const double m0 = json["rhot"]["m0"].number_value();
  const int wig = json["rhot"]["wig"].int_value();
  const string evecFile = json["rhot"]["evecFile"].string_value();
  mat evec;
  if (evec.load(evecFile, arma_ascii)) {
    evec.print(evecFile);
  } else {
    printf("Fail to load evec!\n");
  }

  cx_cube ddos = zeros<cx_cube>(d.nsys, d.nsys, d.nmax);
  ddos(inistate, inistate, 0) = 1.0;

  cx_mat xopr = d.qmd1.slice(0);
  cx_mat popr = (deom_ci * m0) * (d.ham1 * xopr - xopr * d.ham1);

  FILE *fs = fopen("x-p-xx-pp-xp.dat", "w");
  for (int it = 0; it < nt; ++it) {
    const double t = it * dt;
    if (it % nk == 0) {
      printf("Propagation %5.1f%%: nddo=%6d, lddo=%3d\n",
             100 * it / static_cast<double>(nt), d.nddo, d.lddo);

      if (wig) {
        cx_mat rhox = evec * ddos.slice(0) * evec.t();
        mat rhow = rhox2w(rhox, xi, xf, rhox.n_rows);
        string fname = "rho-" + to_string(it) + ".dat";
        rhow.save(fname, raw_ascii);
      }

      double x = real(trace(xopr * ddos.slice(0)));
      double p = real(trace(popr * ddos.slice(0)));
      double xx = real(trace(xopr * xopr * ddos.slice(0)));
      double pp = real(trace(popr * popr * ddos.slice(0)));
      cx_double xp = trace(xopr * popr * ddos.slice(0));
      fprintf(fs, "%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n", t, x, p, xx, pp, real(xp), imag(xp));
    }
    d.rk4(ddos, t, dt);
  }
  fclose(fs);

  return 0;
}
