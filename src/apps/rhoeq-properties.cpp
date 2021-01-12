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

  const int nk = json["rhot"]["nk"].int_value();
  const double dt = json["rhot"]["dt"].number_value();
  const double staticErr = json["rhot"]["staticErr"].number_value();

  cx_cube ddos;
  d.equilibrium(ddos, dt, staticErr, nk, "SCI2");

  double Svn = d.entropy(ddos);
  double Umi = d.Umicro(ddos);
  double Ami = d.Amicro(ddos);
  double Uma = d.Umacro(ddos);
  double Ama = d.Amacro(ddos);
  double hsb = Uma - Umi;
  FILE *f = fopen("thermo.real", "w");
  fprintf(f, "%16.6e%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e\n", 1.0 / d.temperature, Svn, Umi, Ami, Uma, Ama, hsb);
  fclose(f);

  return 0;
}
