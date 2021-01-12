/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 *
 * TFQMR for HEOM steady state
 */
#include "deom.hpp"
#include <cmath>

const int ntMax = 10000000;

static cx_double cdot(const cx_cube &x, const cx_cube &y) {
  const int n1 = x.n_slices, n2 = x.n_rows, n3 = x.n_cols;
  cx_double tmp = 0.0;
  for (int k = 0; k < n1; ++k) {
    for (int i = 0; i < n2; ++i) {
      for (int j = 0; j < n3; ++j) {
        tmp += conj(x(i, j, k)) * y(i, j, k);
      }
    }
  }
  return tmp;
}

static double norm(const cx_cube &x) {
  const int n1 = x.n_slices, n2 = x.n_rows, n3 = x.n_cols;
  double nrm = 0.0;
  for (int k = 0; k < n1; ++k) {
    for (int i = 0; i < n2; ++i) {
      for (int j = 0; j < n3; ++j) {
        double tmp = abs(x(i, j, k));
        nrm += tmp * tmp;
      }
    }
  }
  return sqrt(nrm);
}

void deom::tfqmr(cx_mat &x0, const cx_mat &b0, const cx_double decay, const int iddo, const double tol) {

  cx_mat qrho(size(x0));

  cx_mat rk = b0 - deom_ci * (ham1 * x0 - x0 * ham1) - decay * x0;
  for (int m = 0; m < nmod; ++m) {
    if (abs(delt_res(m) > 1.0e-15)) {
      qrho = qmd1.slice(m) * x0 - x0 * qmd1.slice(m);
      rk -= delt_res(m) * (qmd1.slice(m) * qrho - qrho * qmd1.slice(m));
    }
  }
  if (iddo == 0) {
    rk(0, 0) -= trace(x0);
  }
  cx_mat wk(rk);
  cx_mat uk(rk);

  cx_mat vk = deom_ci * (ham1 * uk - uk * ham1) + decay * uk;
  if (iddo == 0) {
    vk(0, 0) += trace(uk);
  }

  double tau = norm(rk);
  double theta = 0.0;
  cx_double eta = 0.0;
  cx_double alpha = 0.0;

  cx_mat rp = ones<mat>(size(rk)) * deom_c1;
  cx_mat dk = zeros<cx_mat>(size(x0));
  cx_mat zk = zeros<cx_mat>(size(x0));

  cx_double rho0 = cdot(vectorise(rp), vectorise(rk));

  int m = 0;
  int counter = 0;
  // FILE *ferr = fopen("tfqmr-error.dat","w");
  do {
    rk = uk;
    if (m % 2 == 0) {
      alpha = rho0 / cdot(vectorise(rp), vectorise(vk));
      uk -= alpha * vk;
    }
    dk = rk + (theta * theta * eta / alpha) * dk;

    zk = deom_ci * (ham1 * rk - rk * ham1) + decay * rk;
    for (int m = 0; m < nmod; ++m) {
      if (abs(delt_res(m) > 1.0e-15)) {
        qrho = qmd1.slice(m) * rk - rk * qmd1.slice(m);
        zk += delt_res(m) * (qmd1.slice(m) * qrho - qrho * qmd1.slice(m));
      }
    }
    if (iddo == 0) {
      zk(0, 0) += trace(rk);
    }
    counter++;
    wk -= alpha * zk;

    theta = norm(wk) / tau;
    double cm = 1.0 / sqrt(1.0 + theta * theta);
    tau = tau * theta * cm;
    eta = cm * cm * alpha;
    x0 += eta * dk;

    rk = deom_ci * (ham1 * x0 - x0 * ham1) + decay * x0;
    for (int m = 0; m < nmod; ++m) {
      if (abs(delt_res(m) > 1.0e-15)) {
        qrho = qmd1.slice(m) * x0 - x0 * qmd1.slice(m);
        rk += delt_res(m) * (qmd1.slice(m) * qrho - qrho * qmd1.slice(m));
      }
    }
    if (iddo == 0) {
      rk(0, 0) += trace(x0);
    }
    double error = norm(b0 - rk);
    printf("iddo = %d, %d %16.6e\n", iddo, counter, error);
    // fprintf(ferr, "%d %16.6e\n", counter, error);
    if (error < tol)
      break;

    if (m % 2 == 1) {
      cx_double rho1 = rho0;
      rho0 = cdot(vectorise(rp), vectorise(wk));
      cx_double rho2 = rho0 / rho1;
      uk = wk + rho2 * uk;
      vk = rho2 * (zk + rho2 * vk);

      zk = deom_ci * (ham1 * uk - uk * ham1) + decay * uk;
      for (int m = 0; m < nmod; ++m) {
        if (abs(delt_res(m) > 1.0e-15)) {
          qrho = qmd1.slice(m) * uk - uk * qmd1.slice(m);
          zk += delt_res(m) * (qmd1.slice(m) * qrho - qrho * qmd1.slice(m));
        }
      }
      if (iddo == 0) {
        zk(0, 0) += trace(uk);
      }
      counter++;
      vk += zk;
    }
    m++;
  } while (m < ntMax);
  // fclose(ferr);
}

void deom::tfqmr(cx_cube &x0, const double tol) {

  cx_cube rk = zeros<cx_cube>((size(x0)));
  rem_freeze(rk, x0, 0.0);
  rk(0, 0, 0) += trace(x0.slice(0)) - 1.0;
  rk *= -1.0;
  cx_cube wk(rk);
  cx_cube uk(rk);

  cx_cube vk = zeros<cx_cube>((size(x0)));
  rem_freeze(vk, uk, 0.0);
  vk(0, 0, 0) += trace(uk.slice(0));

  double tau = norm(rk);
  double theta = 0.0;
  cx_double eta = 0.0;
  cx_double alpha = 0.0;

  cx_cube rp = ones<cube>(size(rk)) * deom_c1;
  cx_double rho0 = cdot(rp, rk);
  cx_cube dk = zeros<cx_cube>((size(x0)));
  cx_cube zk = zeros<cx_cube>((size(x0)));

  int m = 0;
  int counter = 0;
  FILE *ferr = fopen("tfqmr-error.dat", "w");
  do {
    rk = uk;
    if (m % 2 == 0) {
      alpha = rho0 / cdot(rp, vk);
      uk -= alpha * vk;
    }
    dk = rk + (theta * theta * eta / alpha) * dk;

    rem_freeze(zk, rk, 0.0);
    zk(0, 0, 0) += trace(rk.slice(0));
    counter++;
    wk -= alpha * zk;

    theta = norm(wk) / tau;
    double cm = 1.0 / sqrt(1.0 + theta * theta);
    tau = tau * theta * cm;
    eta = cm * cm * alpha;
    x0 += eta * dk;

    rem_freeze(rk, x0, 0.0);
    double error = norm(rk);
    fprintf(ferr, "%d %16.6e\n", counter, error);
    if (error < tol)
      break;

    if (m % 2 == 1) {
      cx_double rho1 = rho0;
      rho0 = cdot(rp, wk);
      cx_double rho2 = rho0 / rho1;
      uk = wk + rho2 * uk;
      vk = rho2 * (zk + rho2 * vk);
      rem_freeze(zk, uk, 0.0);
      zk(0, 0, 0) += trace(uk.slice(0));
      counter++;
      vk += zk;
    }
    m++;
  } while (m < ntMax);
  fclose(ferr);
}

void deom::bicgstab(cx_cube &x0, const double tol) {

  cx_cube rk = zeros<cx_cube>(size(x0));
  rem_freeze(rk, x0, 0.0);
  rk(0, 0, 0) += trace(x0.slice(0)) - 1.0;
  rk *= -1.0;
  cx_cube rp(rk);
  cx_cube pk(rk);
  cx_cube sk = zeros<cx_cube>(size(x0));
  cx_cube qk = zeros<cx_cube>(size(x0));
  cx_cube tk = zeros<cx_cube>(size(x0));
  cx_cube zk = zeros<cx_cube>(size(x0));

  int counter = 0;
  FILE *ferr = fopen("bicgstab-error.dat", "w");
  do {
    cx_double rho0 = cdot(rp, rk);
    rem_freeze(qk, pk, 0.0);
    qk(0, 0, 0) += trace(pk.slice(0));
    counter++;
    cx_double alpha = rho0 / cdot(rp, qk);
    sk = rk - alpha * qk;
    rem_freeze(tk, sk, 0.0);
    tk(0, 0, 0) += trace(sk.slice(0));
    counter++;
    cx_double omega = cdot(sk, tk) / norm(sk);
    x0 += alpha * pk + omega * sk;
    ///
    rem_freeze(zk, x0, 0.0);
    double error = norm(zk);
    fprintf(ferr, "%d %16.6e\n", counter, error);
    if (error < tol)
      break;
    ///
    rk = sk - omega * tk;
    cx_double theta = cdot(rp, rk) * alpha / (omega * rho0);
    pk = rk + theta * (pk - omega * qk);
  } while (counter < ntMax);
  fclose(ferr);
}
