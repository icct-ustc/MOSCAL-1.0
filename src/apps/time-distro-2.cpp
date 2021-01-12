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

  const int inistate = json["rhot"]["inistate"].int_value();
  const int obstate1 = json["rhot"]["obstate1"].int_value();
  const int obstate2 = json["rhot"]["obstate2"].int_value();
  const int nt = json["rhot"]["nt"].int_value();
  const double tau = json["rhot"]["tau"].number_value();
  const double dt = json["rhot"]["dt"].number_value();
  const double staticErr = json["rhot"]["staticErr"].number_value();
  // const int    nk = json["rhot"]["nk"].int_value();

  int ntau = floor(tau / dt);
  double dt_res = tau - dt * ntau;

  cx_cube rho1 = zeros<cx_cube>(d1.nsys, d1.nsys, d1.nmax);
  cx_cube rho2 = zeros<cx_cube>(d1.nsys, d1.nsys, d1.nmax);
  rho1(inistate, inistate, 0) = 1.0;

  // d1.equilibrium(rho1, dt, staticErr, nk, "SCI3");

  int iter = 0;
  int nddo_backup(d1.nddo);
  hkey keys_backup(d1.keys);
  cx_cube ddos_backup = zeros<cx_cube>(size(rho1));
  vec pop = zeros<vec>(d1.nsys);
  while (iter < 1000000) {
    double t = 0;
    pop = real(rho1.slice(0).diag());
    pop.print("pop");
    // von Neumann wavefunction collapse
    for (int iddo = 0; iddo < d1.nddo; ++iddo) {
      for (int i = 0; i < d1.nsys - 1; ++i) {
        for (int j = i + 1; j < d1.nsys; ++j) {
          rho1(i, j, iddo) = rho1(j, i, iddo) = 0;
        }
      }
    }
    // Check convergence
    if (d1.notconverged(nddo_backup, keys_backup, ddos_backup, rho1, staticErr)) {
      nddo_backup = d1.nddo;
      keys_backup.rows(0, d1.nddo - 1) = d1.keys.rows(0, d1.nddo - 1);
      ddos_backup.head_slices(d1.nddo) = rho1.head_slices(d1.nddo);
    } else {
      break;
    }
    // propagation
    for (int itau = 0; itau < ntau; ++itau) {
      d1.rk4(rho1, t, dt);
      t += dt;
    }
    d1.rk4(rho1, t, dt_res);
    iter += 1;
  }

  for (int iddo = 0; iddo < d1.nddo; ++iddo) {
    cx_double tmp = rho1(obstate1, obstate1, iddo);
    rho1.slice(iddo).zeros();
    rho1(obstate1, obstate1, iddo) = tmp;
  }
  double pa0 = real(rho1(obstate1, obstate1, 0));

  FILE *frho = fopen("prop-rho.dat", "w");
  for (int it1 = 0; it1 < nt; ++it1) {
    double t1 = it1 * tau;
    deom d2(d1);
    rho2.head_slices(d2.nddo) = rho1.head_slices(d2.nddo);
    for (int it2 = 0; it2 < nt; ++it2) {
      double t2 = it2 * tau;
      printf("Propagation %5.1f%%: nddo=%6d, lddo=%3d\n",
             100 * it2 / static_cast<double>(nt), d2.nddo, d2.lddo);
      // Ua2-propagate
      for (int itau = 0; itau < ntau; ++itau) {
        d2.rk4(rho2, t2, dt);
        t2 += dt;
      }
      d2.rk4(rho2, t2, dt_res);
      // sum-over-b
      double pa1a2 = real(trace(rho2.slice(0))) - real(rho2(obstate2, obstate2, 0));
      fprintf(frho, "%16.6e%16.6e%16.6e\n", t1 / deom_fs2unit, t2 / deom_fs2unit, pa1a2 / pa0);
      // Ua2-measure
      for (int iddo = 0; iddo < d2.nddo; ++iddo) {
        cx_double tmp = rho2(obstate2, obstate2, iddo);
        rho2.slice(iddo).zeros();
        rho2(obstate2, obstate2, iddo) = tmp;
      }
    }
    // Ua1-propagate
    for (int itau = 0; itau < ntau; ++itau) {
      d1.rk4(rho1, t1, dt);
      t1 += dt;
    }
    d1.rk4(rho1, t1, dt_res);
    // Ua1-measure
    for (int iddo = 0; iddo < d1.nddo; ++iddo) {
      cx_double tmp = rho1(obstate1, obstate1, iddo);
      rho1.slice(iddo).zeros();
      rho1(obstate1, obstate1, iddo) = tmp;
    }
  }
  fclose(frho);

  return 0;
}
