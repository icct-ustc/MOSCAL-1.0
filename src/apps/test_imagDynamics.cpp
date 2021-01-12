/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "ideom.hpp"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>


static double entropy(const cx_mat &rho) {
  vec eval = eig_sym(rho);
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

  ideom d(json["deom"]);

  // const int inistate = json["rhot"]["inistate"].int_value();
  const double dt = json["rhot"]["dt"].number_value();
  const double staticErr = json["rhot"]["staticErr"].number_value();
  const int nk = json["rhot"]["nk"].int_value();

  cx_cube ddos = zeros<cx_cube>(d.nsys, d.nsys, d.nmax);

  d.equilibrium(ddos, dt, staticErr, nk);

  double Svn = d.entropy(ddos);
  double Umi = d.Umicro(ddos);
  double Ami = d.Amicro(ddos);
  double Uma = d.Umacro(ddos);
  double Ama = d.Amacro(ddos);
  double Ahyb = d.Umacro(ddos) - d.Umicro(ddos);
  //double zpr = real(trace(expmat_sym(-d.ham1/d.temperature)))*exp((Umi-Uma)/d.temperature);
  double zpr = real(exp(-(Umi - Uma) / d.temperature));
  FILE *f = fopen("thermo.imag", "w");
  fprintf(f, "%16.6e%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e\n", 1.0 / d.temperature, zpr, Ahyb, Umi, Ami, Uma, Ama);
  fclose(f);

  return 0;
}
