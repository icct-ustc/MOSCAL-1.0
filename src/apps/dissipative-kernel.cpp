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

static void dissipativekernel(deom &d, const double dt, const int nt, const int nk, const int &projection) {

  const int nsys = d.nsys;
  const int nsys2 = nsys * nsys;
  cx_cube rho1(nsys, nsys, d.nmax);
  cx_cube rho2(nsys, nsys, d.nmax);
  cx_mat kernelt(nt, nsys2);
  cx_mat kernel0(nsys2, nsys2);

  // L_{r1r2},{c1c2}
  for (int c1 = 0; c1 < nsys; ++c1) {
    for (int c2 = 0; c2 < nsys; ++c2) {
      const int cc = c1 * nsys + c2;
      // Time-domain
      rho1.zeros();
      rho1(c1, c2, 0) = 1.0;
      d.rem(rho2, rho1, 0.0, projection);
      char filename[64];
      sprintf(filename, "KernelT_from_%d", cc);
      FILE *flog = fopen(filename, "w");
      for (int it = 0; it < nt; ++it) {
        d.rem(rho1, rho2, 0.0);
        kernelt.row(it) = -vectorise(rho1.slice(0), 1);
        //
        double t = it * dt;
        if (it % nk == 0) {
          printf("%s: %5.1f%%, nddo = %6d, lddo = %2d\n", filename,
                 100 * it / static_cast<double>(nt), d.nddo, d.lddo);
        }
        fprintf(flog, "%16.6e", t / deom_fs2unit);
        for (int rr = 0; rr < nsys2; ++rr) {
          fprintf(flog, "%20.12e%20.12e", real(kernelt(it, rr)), imag(kernelt(it, rr)));
        }
        fprintf(flog, "\n");
        d.rk4(rho2, t, dt, projection);
      }
      fclose(flog);

      // Freq-domain
      const double dw = 2. * deom_pi / (nt * dt);
      kernelt.row(0) *= 0.5;
      cx_mat kernelw = ifft(kernelt) * nt * dt;

      for (int rr = 0; rr < nsys2; ++rr) {
        // output K(0)
        kernel0(rr, cc) = kernelw(0, rr);
        // output K(s)
        char filename[64];
        sprintf(filename, "KernelW_from%dto%d.out", cc, rr);
        FILE *fout = fopen(filename, "w");
        for (int iw = nt / 2; iw < nt; ++iw) {
          double w = (iw - nt) * dw / deom_cm2unit;
          double re = real(kernelw(iw, rr));
          double im = imag(kernelw(iw, rr));
          fprintf(fout, "%16.6e%16.6e%16.6e\n", w, re, im);
        }
        for (int iw = 0; iw < nt / 2; ++iw) {
          double w = iw * dw / deom_cm2unit;
          double re = real(kernelw(iw, rr));
          double im = imag(kernelw(iw, rr));
          fprintf(fout, "%16.6e%16.6e%16.6e\n", w, re, im);
        }
        fclose(fout);
      }
    }
  }
  kernel0.print("DissipativeTensor");
  kernel0.save("DissipativeTensor", raw_ascii);
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

  const double dt = json["rhot"]["dt"].number_value();
  const int nt = json["rhot"]["nt"].int_value();
  const int nk = json["rhot"]["nk"].int_value();
  const int projection = 1;
  dissipativekernel(d, dt, nt, nk, projection);

  return 0;
}
