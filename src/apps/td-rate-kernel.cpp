/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "deom.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

static void td_rate_kernel(deom &d, const double dt, const double ti, const double tf, const int nk, const mat &sdip, const cube &pdip, const vec &bdip, pulse &p, const ivec &projection) {

  const int nsys = d.nsys;
  const int nt_i = ceil(abs(ti) / dt);
  const int nt_f = ceil(abs(tf) / dt);

  double t = -(nt_i - 1) * dt;
  cx_cube rho1 = zeros<cx_cube>(nsys, nsys, d.nmax);
  cx_cube rho2 = zeros<cx_cube>(nsys, nsys, d.nmax);
  cx_cube rho3 = zeros<cx_cube>(nsys, nsys, d.nmax);
  cx_mat fts(nt_f, nsys);
  mat k0(nsys, nsys);

  p.turnon();

  for (int si = 1; si < nsys; ++si) {
    if (projection(si) != 0) {

      t = 0.0;
      // Time-domain
      rho1.zeros();
      rho1(si, si, 0) = 1.0;
      d.rem(rho2, rho1, 0.0, sdip, pdip, bdip, p, projection);
      char filename[64];
      sprintf(filename, "td-rate-kernel-from%d", si);
      FILE *flog = fopen(filename, "w");
      for (int it = 0; it < nt_f; ++it) {
        d.rem(rho3, rho2, t, sdip, pdip, bdip, p);
        for (int sf = 0; sf < nsys; ++sf) {
          fts(it, sf) = -rho3(sf, sf, 0);
        }
        //
        if (it % nk == 0) {
          printf("%s: %5.1f%%, nddo = %6d, lddo = %2d\n", filename,
                 100 * it / static_cast<double>(nt_f), d.nddo, d.lddo);
        }
        fprintf(flog, "%16.6e", t / deom_fs2unit);
        for (int sf = 0; sf < nsys; ++sf) {
          fprintf(flog, "%16.6e%16.6e", real(fts(it, sf)), imag(fts(it, sf)));
        }
        fprintf(flog, "\n");
        d.rk4(rho2, t, dt, sdip, pdip, bdip, p, projection);
        t += dt;
      }
      fclose(flog);

      // Freq-domain
      double dw = 2. * deom_pi / (nt_f * dt);
      fts.row(0) *= 0.5;
      cx_mat fws = ifft(fts) * nt_f * dt;

      for (int sf = 0; sf < nsys; ++sf) {
        // output K(0)
        k0(sf, si) = real(fws(0, sf));
        // output K(s)
        char filename[64];
        sprintf(filename, "td-rate-const-from%dto%d.out", si, sf);
        FILE *fout = fopen(filename, "w");
        for (int iw = nt_f / 2; iw < nt_f; ++iw) {
          double w = (iw - nt_f) * dw / deom_cm2unit;
          double re = real(fws(iw, sf));
          double im = imag(fws(iw, sf));
          fprintf(fout, "%16.6e%16.6e%16.6e\n", w, re, im);
        }
        for (int iw = 0; iw < nt_f / 2; ++iw) {
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

  const double dt = json["td-rate"]["dt"].number_value();
  const double ti = json["td-rate"]["ti"].number_value();
  const double tf = json["td-rate"]["tf"].number_value();
  const int nk = json["td-rate"]["nk"].int_value();
  const string sdipFile = json["td-rate"]["sdipFile"].string_value();
  const string pdipFile = json["td-rate"]["pdipFile"].string_value();
  const string bdipFile = json["td-rate"]["bdipFile"].string_value();
  const string projFile = json["td-rate"]["projFile"].string_value();

  pulse p(json["td-rate"]["pulse"]);

  mat sdip;
  cube pdip;
  vec bdip;
  ivec proj;

  if (sdip.load(sdipFile, arma_ascii)) {
    sdip.print(sdipFile);
  } else {
    printf("Fail to load sdip");
  }

  if (pdip.load(pdipFile, arma_ascii)) {
    pdip.print(pdipFile);
  } else {
    printf("Fail to load pdip");
  }

  if (bdip.load(bdipFile, arma_ascii)) {
    bdip.print(bdipFile);
  } else {
    printf("Fail to load bdip");
  }

  if (proj.load(projFile, arma_ascii)) {
    proj.print(projFile);
  } else {
    printf("Fail to load proj");
  }

  if (proj.n_rows != sdip.n_rows) {
    printf("Error! Invalid projection operator!\n");
    return 0;
  }

  td_rate_kernel(d, dt, ti, tf, nk, sdip, pdip, bdip, p, proj);

  return 0;
}
