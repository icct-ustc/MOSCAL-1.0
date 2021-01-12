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

static void pcet_flux(deom &d, const double dt, const int nt, const int nk, const ivec &projection) {

  const int nsys = d.nsys;
  cx_cube rho1(nsys, nsys, d.nmax);
  cx_cube rho2(nsys, nsys, d.nmax);
  cx_mat fts(nt, nsys);
  mat k0(nsys, nsys);

  for (int si = 0; si < nsys; ++si) {
    if (projection(si) != 0) {

      // Time-domain
      rho1.zeros();
      rho1(si, si, 0) = 1.0;
      d.rem(rho2, rho1, 0.0, projection);
      char filename[64];
      sprintf(filename, "rateKernel_from%d", si);
      FILE *flog = fopen(filename, "w");
      for (int it = 0; it < nt; ++it) {
        d.rem(rho1, rho2, 0.0);
        for (int sf = 0; sf < nsys; ++sf) {
          fts(it, sf) = -rho1(sf, sf, 0);
        }
        //
        double t = it * dt;
        if (it % nk == 0) {
          printf("%s: %5.1f%%, nddo = %6d, lddo = %2d\n", filename,
                 100 * it / static_cast<double>(nt), d.nddo, d.lddo);
        }
        fprintf(flog, "%16.6e", t / deom_fs2unit);
        for (int sf = 0; sf < nsys; ++sf) {
          fprintf(flog, "%20.12e%20.12e", real(fts(it, sf)), imag(fts(it, sf)));
        }
        fprintf(flog, "\n");
        d.rk4(rho2, t, dt, projection);
      }
      fclose(flog);

      // Freq-domain
      double dw = 2. * deom_pi / (nt * dt);
      fts.row(0) *= 0.5;
      cx_mat fws = ifft(fts) * nt * dt;

      for (int sf = 0; sf < nsys; ++sf) {
        // output K(0)
        k0(sf, si) = real(fws(0, sf));
        // output K(s)
        char filename[64];
        sprintf(filename, "rate_from%dto%d.out", si, sf);
        FILE *fout = fopen(filename, "w");
        for (int iw = nt / 2; iw < nt; ++iw) {
          double w = (iw - nt) * dw / deom_cm2unit;
          double re = real(fws(iw, sf));
          double im = imag(fws(iw, sf));
          fprintf(fout, "%16.6e%16.6e%16.6e\n", w, re, im);
        }
        for (int iw = 0; iw < nt / 2; ++iw) {
          double w = iw * dw / deom_cm2unit;
          double re = real(fws(iw, sf));
          double im = imag(fws(iw, sf));
          fprintf(fout, "%16.6e%16.6e%16.6e\n", w, re, im);
        }
        fclose(fout);
      }
    }
  }
  k0.print("RateConstant:");
  k0.save("RateConstant", raw_ascii);

  // equilibrium populations
  mat amat(k0);
  vec pt(d.nsys);
  pt(0) = 1.0;
  for (int i = 0; i < d.nsys; ++i) {
    amat(0, i) += 1;
  }
  vec pf = solve(amat, pt);
  pf.print("Equilibrium populations");

  double phi0 = 0.0;
  pt.zeros();
  pt(0) = 1.0;
  pt -= pf;
  do {
    phi0 += -k0(3, 0) * (pt(0)) + k0(0, 3) * (pt(3));
    // RK4
    vec p1 = -k0 * pt;
    vec p2 = -k0 * (pt + p1 * dt * 0.5);
    p1 += 2 * p2;
    p2 = -k0 * (pt + p2 * dt * 0.5);
    p1 += 2 * p2;
    p2 = -k0 * (pt + p2 * dt);
    p1 += p2;
    pt += p2 * dt / 6.0;
  } while (any(abs(pt) > 1.0e-10));
  phi0 *= dt;

  printf("PCET Flux:\n");
  printf("%16.6e\n", phi0);
  printf("kappa=%16.6e\n", phi0 / (pf(3) - phi0));
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

  const double dt = json["rate"]["dt"].number_value();
  const int nt = json["rate"]["nt"].int_value();
  const int nk = json["rate"]["nk"].int_value();
  const string projFile = json["rate"]["projFile"].string_value();
  ivec proj;
  if (proj.load(projFile, arma_ascii)) {
    proj.print(projFile);
  } else {
    printf("Fail to load proj");
  }

  pcet_flux(d, dt, nt, nk, proj);

  return 0;
}
