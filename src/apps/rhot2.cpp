/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "deom2.hpp"
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

  deom2 d(json["deom"]);

  const int inistate = json["rhot"]["inistate"].int_value();
  const int nt = json["rhot"]["nt"].int_value();
  const double dt = json["rhot"]["dt"].number_value();
  const int nk = json["rhot"]["nk"].int_value();

  cx_cube ddos = zeros<cx_cube>(d.nsys, d.nsys, d.nmax);
  ddos(inistate, inistate, 0) = 1.0;

  FILE *frho = fopen("proprho.dat", "w");
  FILE *fent = fopen("entropy.dat", "w");
  for (int it = 0; it < nt; ++it) {
    const double t = it * dt;
    printf("Propagation %5.1f%%: nddo=%6d, lddo=%3d\n",
           100 * it / static_cast<double>(nt), d.nddo, d.lddo);
    if (it % nk == 0) {
      fprintf(frho, "%16.6e", t / deom_fs2unit);
      for (int i = 0; i < d.nsys; ++i) {
        fprintf(frho, "%20.10e", real(ddos(i, i, 0)));
      }
      fprintf(frho, "\n");
      double tmp = entropy(ddos.slice(0));
      if (tmp < 0) {
        fprintf(fent, "%16.6e nan\n", t / deom_fs2unit);
      } else {
        fprintf(fent, "%16.6e%16.6e\n", t / deom_fs2unit, tmp);
      }
      d.status();
    }
    d.rk4(ddos, t, dt);
  }
  fclose(fent);
  fclose(frho);

  return 0;
}
