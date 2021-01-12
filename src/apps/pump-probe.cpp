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

static void pump_probe(const double &w1_max, const int &nt1, const double &dt, const double &staticErr, const int &nk,
                       const mat &sdip, const cube &pdip, const vec &bdip,
                       pulse &pump, const vec &tvec,
                       const syst &s, const bath &b, const hidx &h) {

  const double dw1 = w1_max / nt1;
  const double dt1 = 2.0 * deom_pi / w1_max;
  const int mt = floor(dt1 / dt);
  const double dt1_res = dt1 - dt * mt;

  deom d1(s, b, h);

  cx_cube rho0 = zeros<cx_cube>(d1.nsys, d1.nsys, d1.nmax);
  cx_cube rho1 = zeros<cx_cube>(d1.nsys, d1.nsys, d1.nmax);

  const mat &exph = expmat(-real(d1.ham1) / d1.temperature);
  rho0.slice(0).set_real(exph / trace(exph));
  d1.equilibrium(rho0, dt, staticErr, nk);

  int iprob = -1;
  const int nprob = tvec.n_rows;
  const double ti = -pump.sigm * 3;
  const double tf = tvec(nprob - 1);
  const int nt = static_cast<int>((tf - ti + dt) / dt);

  cx_vec ft = zeros<cx_vec>(nt1);

  pump.turnon();
  for (int it = 0; it < nt; ++it) {
    double t = ti + dt * nt;
    if (iprob == -1 || abs(t - tvec(iprob)) < 0.5 * dt) {
      deom d2(d1);
      // \mu\rho
      d2.oprAct(rho1, sdip, pdip, bdip, rho0, 'l');
      // output filename
      stringstream sst, ssw;
      sst << "corr_t_" << iprob << ".dat";
      ssw << "corr_w_" << iprob << ".dat";
      // iTr[\mu G(t)\mu\rho]
      FILE *log_t = fopen(sst.str().c_str(), "w");
      for (int it1 = 0; it < nt1; ++it1) {
        double t1 = it1 * dt1;
        ft(it1) = deom_ci * d2.Trace(sdip, pdip, bdip, rho1);
        printf("In probe %d: it1=%d, nddo=%d, lddo=%d\n", iprob, it1, d2.nddo, d2.lddo);
        fprintf(log_t, "%16.6e%16.6e%16.6e\n", t1 / deom_fs2unit, real(ft(it1)), imag(ft(it1)));
        for (int jt1 = 0; jt1 < mt; ++jt1) {
          d2.rk4(rho1, t1, dt);
          t1 += dt;
        }
        d2.rk4(rho1, t1, dt1_res);
      }
      fclose(log_t);

      // 1D FFT
      ft(0) *= 0.5;
      const cx_vec &fw = ifft(ft) * nt1 * dt1;

      FILE *log_w = fopen(ssw.str().c_str(), "w");
      for (int iw = nt1 / 2; iw < nt1; ++iw) {
        double w = (iw - nt1) * dw1 / deom_cm2unit;
        fprintf(log_w, "%16.6e%16.6e%16.6e\n", w, real(fw(iw)), imag(fw(iw)));
      }
      for (int iw = 0; iw < nt1 / 2; ++iw) {
        double w = iw * dw1 / deom_cm2unit;
        fprintf(log_w, "%16.6e%16.6e%16.6e\n", w, real(fw(iw)), imag(fw(iw)));
      }
      fclose(log_w);

      iprob += 1; // wait for next probe
    }
    d1.rk4(rho0, it, dt);
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

  syst s(json["deom"]["syst"]);
  bath b(json["deom"]["bath"]);
  hidx h(json["deom"]["hidx"]);

  const double w1_max = json["spec"]["w1max"].number_value();
  const int nt1 = json["spec"]["nt1"].int_value();
  const double dt = json["spec"]["dt"].number_value();
  const double staticErr = json["spec"]["staticErr"].number_value();
  const int nk = json["spec"]["nk"].int_value();
  const string sch_hei = json["spec"]["sch_hei"].string_value();
  const string sdipFile = json["spec"]["sdipFile"].string_value();
  const string pdipFile = json["spec"]["pdipFile"].string_value();
  const string bdipFile = json["spec"]["bdipFile"].string_value();
  const string tvecFile = json["spec"]["tvecFile"].string_value();
  mat sdip;
  cube pdip;
  vec bdip;
  vec tvec;
  if (sdip.load(sdipFile, arma_ascii)) {
    sdip.print(sdipFile);
  } else {
    printf("Fail to load sdip!\n");
  }
  if (pdip.load(pdipFile, arma_ascii)) {
    pdip.print(pdipFile);
  } else {
    printf("Fail to load pdip!\n");
  }
  if (bdip.load(bdipFile, arma_ascii)) {
    bdip.print(bdipFile);
  } else {
    printf("Fail to load bdip!\n");
  }
  if (tvec.load(tvecFile, arma_ascii)) {
    tvec.print(tvecFile);
  } else {
    printf("Fail to load tvec!\n");
  }

  pulse pump(json["spec"]["pump"]);

  pump_probe(w1_max, nt1, dt, staticErr, nk, sdip, pdip, bdip, pump, tvec, s, b, h);

  return 0;
}
