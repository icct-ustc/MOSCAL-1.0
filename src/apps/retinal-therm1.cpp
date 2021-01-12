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
  const int nvib = json["retinal"]["nvib"].int_value();
  const string rho0File = json["retinal"]["rho0File"].string_value();
  const string tup0File = json["retinal"]["tup0File"].string_value();
  const string tup1File = json["retinal"]["tup1File"].string_value();
  const string vt0File = json["retinal"]["vt0File"].string_value();
  const string vv0File = json["retinal"]["vv0File"].string_value();
  const string vt1File = json["retinal"]["vt1File"].string_value();
  const string vv1File = json["retinal"]["vv1File"].string_value();

  cx_mat rho0;
  if (rho0.load(rho0File, arma_ascii)) {
    rho0.print("Initial state!");
  } else {
    printf("Error in loading rho0!\n");
    return 0;
  }

  mat vt0, vv0, vt1, vv1;
  if (vt0.load(vt0File, arma_ascii)) {
    vt0.print("vt0!");
  } else {
    printf("Error in loading vt0!\n");
    return 0;
  }
  if (vv0.load(vv0File, arma_ascii)) {
    vv0.print("vv0!");
  } else {
    printf("Error in loading vv0!\n");
    return 0;
  }
  if (vt1.load(vt1File, arma_ascii)) {
    vt1.print("vt1!");
  } else {
    printf("Error in loading vt1!\n");
    return 0;
  }
  if (vv1.load(vv1File, arma_ascii)) {
    vv1.print("vv1!");
  } else {
    printf("Error in loading vv1!\n");
    return 0;
  }

  imat tup0, tup1;
  if (tup0.load(tup0File, arma_ascii)) {
    tup0.print("tup0!");
  } else {
    printf("Error in loading tup0!\n");
    return 0;
  }
  if (tup1.load(tup1File, arma_ascii)) {
    tup1.print("tup1!");
  } else {
    printf("Error in loading tup1!\n");
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
      // Cis/Trans
      double pcis_0 = 0.0;
      double pcis_1 = 0.0;
      double ptrans_0 = 0.0;
      double ptrans_1 = 0.0;
      // Write down torsional density
      fprintf(ftor0, "%16.6e", t / deom_fs2unit);
      fprintf(ftor1, "%16.6e", t / deom_fs2unit);
      for (int i = 0; i < ntor; ++i) {
        double prob_theta0 = 0.0;
        double prob_theta1 = 0.0;
        for (int j = 0; j < nvib; ++j) {
          // const int ii = i*nvib+j;
          // P_ii =  \sum_{jk} W_{ij}\rho_{jk}V_{ki}
          double tmp0 = 0.0;
          for (int jj = 0; jj < ng; ++jj) {
            for (int kk = 0; kk < ng; ++kk) {
              const int it0 = tup0(jj, 0);
              const int iv0 = tup0(jj, 1);
              const int jt0 = tup0(kk, 0);
              const int jv0 = tup0(kk, 1);
              const double wfnl = vt0(i, it0) * vv0(j, iv0);
              const double wfnr = vt0(i, jt0) * vv0(j, jv0);
              tmp0 += wfnl * wfnr * real(r0(jj, kk));
            }
          }
          double tmp1 = 0.0;
          for (int jj = 0; jj < ne; ++jj) {
            for (int kk = 0; kk < ne; ++kk) {
              const int it1 = tup1(jj, 0);
              const int iv1 = tup1(jj, 1);
              const int jt1 = tup1(kk, 0);
              const int jv1 = tup1(kk, 1);
              const double wfnl = vt1(i, it1) * vv1(j, iv1);
              const double wfnr = vt1(i, jt1) * vv1(j, jv1);
              tmp1 += wfnl * wfnr * real(r1(jj, kk));
            }
          }
          prob_theta0 += tmp0;
          prob_theta1 += tmp1;
        }
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
