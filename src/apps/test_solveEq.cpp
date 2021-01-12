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


static double entropy(const cx_mat &rho) {
  vec eval = eig_sym(rho);
  eval.print("eval of rho0");
  if (any(eval < 0)) {
    return -9527;
  } else {
    return -dot(eval, log(eval));
  }
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

  const int nk = json["rhot"]["nk"].int_value();
  const double dt = json["rhot"]["dt"].number_value();
  const double staticErr = json["rhot"]["staticErr"].number_value();

  cx_cube ddos;
  d.equilibrium(ddos, dt, staticErr, nk, "SCI2");

  ddos.slice(0).save("rho0_eq.mat", arma_ascii);

  FILE *f = fopen("entropy.dat", "w");
  fprintf(f, "%16.6e\n", entropy(ddos.slice(0)));
  fclose(f);

  return 0;
}
