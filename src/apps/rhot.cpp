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

  if (d.flag_rwa != 0) {
    printf("Rotating Wave Approximation\n");
  }
  const int inistate = json["rhot"]["inistate"].int_value();
  const int nt = json["rhot"]["nt"].int_value();
  const double dt = json["rhot"]["dt"].number_value();
  const int nk = json["rhot"]["nk"].int_value();

  cx_cube ddos = zeros<cx_cube>(d.nsys, d.nsys, d.nmax);
  ddos(inistate, inistate, 0) = 1.0;

  FILE *frho = fopen("prop-rho.dat", "w");
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
      double s1 = d.entropy(ddos, "vn");
      double s2 = d.entropy(ddos, "sh");
      fprintf(fent, "%16.6e%16.6e%16.6e\n", t / deom_fs2unit, s1, s2);
    }
    d.rk4(ddos, t, dt);
  }
  fclose(fent);
  fclose(frho);
  //cx_mat rho0 = ddos.slice(0);
  //rho0.print("rho0");
  //return 0;
  // save rho0
  ddos.slice(0).save("rho0.real2", raw_ascii);
  // save rhon
  FILE *fs = fopen("rhonreal2.dat", "w");
  for (int i = 0; i < d.nddo; ++i) {
    for (size_t j = 0; j < d.keys(i).key.n_rows; ++j)
      fprintf(fs, "%d\t", (int)d.keys(i).key(j));
    fprintf(fs, "\n");
    for (int j = 0; j < d.nsys; ++j)
      for (int k = 0; k < d.nsys; ++k)
        fprintf(fs, "(%26.16e,%26.16e)\t", ddos(j, k, i).real(), ddos(j, k, i).imag());
    fprintf(fs, "\n\n");
  }
  fclose(fs);
}
