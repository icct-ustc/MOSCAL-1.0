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

  const int nt = json["retinal"]["nt"].int_value();
  const double dt = json["retinal"]["dt"].number_value();
  const int nk = json["retinal"]["nk"].int_value();
  const int ng = json["retinal"]["ng"].int_value();
  const int ne = json["retinal"]["ne"].int_value();
  const int ntor = json["retinal"]["ntor"].int_value();
  const string rho0File = json["retinal"]["rho0File"].string_value();
  const string wfn0File = json["retinal"]["wfn0File"].string_value();
  const string wfn1File = json["retinal"]["wfn1File"].string_value();

  if (ng + ne != d.nsys) {
    printf("No. of system levels is wrong!\n");
    return 0;
  }

  cx_mat rho0;
  if (rho0.load(rho0File, arma_ascii)) {
    rho0.print("Initial state!");
  } else {
    printf("Error in loading rho0!\n");
    return 0;
  }

  mat wfn0;
  if (wfn0.load(wfn0File, arma_ascii)) {
    wfn0.print("Torsional operator on S0 state!");
  } else {
    printf("Error in loading wfn0!\n");
    return 0;
  }

  mat wfn1;
  if (wfn1.load(wfn1File, arma_ascii)) {
    wfn1.print("Torsional operator on S1 state!");
  } else {
    printf("Error in loading wfn1!\n");
    return 0;
  }

  cx_cube ddos = zeros<cx_cube>(d.nsys, d.nsys, d.nmax);
  ddos.slice(0) = rho0;

  FILE *ftor0 = fopen("prop-tor0.dat", "w");
  FILE *ftor1 = fopen("prop-tor1.dat", "w");
  FILE *fisom = fopen("prop-isom.dat", "w");
  for (int it = 0; it < nt; ++it) {
    const double t = it * dt;
    printf("Propagation %5.1f%%: nddo=%6d, lddo=%3d\n",
           100 * it / static_cast<double>(nt), d.nddo, d.lddo);
    if (it % nk == 0) {
      // Transform into coordinate representation
      cx_mat r0 = ddos.slice(0).submat(span(0, ng - 1), span(0, ng - 1));
      cx_mat r1 = ddos.slice(0).submat(span(ng, d.nsys - 1), span(ng, d.nsys - 1));
      r0 = wfn0 * r0 * wfn0.t();
      r1 = wfn1 * r1 * wfn1.t();
      // Cis/Trans
      double pcis_0 = 0.0;
      double pcis_1 = 0.0;
      double ptrans_0 = 0.0;
      double ptrans_1 = 0.0;
      // Write down torsional density
      fprintf(ftor0, "%16.6e", t / deom_fs2unit);
      fprintf(ftor1, "%16.6e", t / deom_fs2unit);
      for (int i = 0; i < ntor; ++i) {
        double prob_theta0 = real(r0(i, i));
        double prob_theta1 = real(r1(i, i));
        if (i < ntor / 4 || i > (ntor / 4) * 3) {
          pcis_0 += prob_theta0;
          pcis_1 += prob_theta1;
        } else {
          ptrans_0 += prob_theta0;
          ptrans_1 += prob_theta1;
        }
        fprintf(ftor0, "%16.6e", prob_theta0);
        fprintf(ftor1, "%16.6e", prob_theta1);
      }
      fprintf(ftor0, "\n");
      fprintf(ftor1, "\n");
      fprintf(fisom, "%16.6e%16.6e%16.6e%16.6e%16.6e\n", t / deom_fs2unit, pcis_0, pcis_1, ptrans_0, ptrans_1);
    }
    d.rk4(ddos, t, dt);
  }
  fclose(fisom);
  fclose(ftor1);
  fclose(ftor0);

  return 0;
}
