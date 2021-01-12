/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "blockdeom.hpp"
#include <cstdio>
#include <cstdlib>

static const int ntMax = 100000000;

void blockdeom::EqSolverSCI3(cx_cube &ddos, const double dt, const double err, const int nk) {

  if (nddo != 1) {
    printf("Warning: Hidx is not correctly initialized!\n");
  }
  vec eigh = eig_sym(ham1);
  const double OMG = max(abs(eigh)) * 500;
  printf("OMG=%16.6e\n", OMG / 4.5554927e-6);

  if (ddos.is_empty()) {
    ddos = zeros<cx_cube>(nsys, nsys, nmax);
  } else {
    ddos.zeros();
  }
  mat eham = expmat_sym(-real(ham1) / temperature);
  ddos.slice(0) = eham / trace(eham) * deom_c1;

  cx_mat rho_old(nsys, nsys);
  cx_mat rhob(nsys, nsys);
  cx_mat rhox(nsys, nsys);
  cx_mat rhoy(nsys, nsys);

  double max_diff = 0.0;
  FILE *ferr = fopen("sci-error.dat", "w");
  int iter = 0;
  do {
    max_diff = 0.0;
    for (int iddo = 0; iddo < nddo; ++iddo) {

      hnod &nod = keys(iddo);
      ivec key0(nod.key);
      const int tier = tree.find(key0)->tier;
      const bool flag = is_valid(ddos.slice(iddo));

      rhob.zeros();
      for (int m = 0; m < nmod; ++m) {
        rho_old = qmd1.slice(m) * ddos.slice(iddo) - ddos.slice(iddo) * qmd1.slice(m);
        if (abs(delt_res(m)) > 1.0e-15) {
          rhob -= delt_res(m) * (qmd1.slice(m) * rho_old - rho_old * qmd1.slice(m));
        }
      }
      if (tier < lmax) {
        for (int mp = 0; mp < nind; ++mp) {
          ivec key1 = gen_key(key0, mp, 1);
          const int m = modLabel(mp);
          const int n = key1(mp) - 1;
          const cx_double sn = -deom_ci * sqrt((n + 1) * coef_abs(mp));
          auto *ptr = tree.find(key1);
          if (ptr && (ptr->rank) >= 0) {
            int loc = ptr->rank;
            rhob += sn * (qmd1.slice(m) * ddos.slice(loc) - ddos.slice(loc) * qmd1.slice(m));
          } else if (flag) {
            tree.try_insert(key1, nddo);
            keys(nddo) = hnod(nod.gams + expn(mp), key1);
            ddos.slice(nddo).zeros();
            nddo += 1;
          }
        }
      }
      for (int mp = 0; mp < nind; ++mp) {
        ivec key1 = gen_key(key0, mp, -1);
        if (!key1.is_empty()) {
          const int m = modLabel(mp);
          const int n = (key1.n_rows < (unsigned int)(mp + 1)) ? 1 : (key1(mp) + 1);
          const cx_double sn = -deom_ci * sqrt(n / coef_abs(mp));
          const cx_double cl = sn * coef_lft(mp);
          const cx_double cr = sn * coef_rht(mp);
          auto *ptr = tree.find(key1);
          if (ptr && (ptr->rank) >= 0) {
            int loc = ptr->rank;
            rhob += cl * qmd1.slice(m) * ddos.slice(loc) - cr * ddos.slice(loc) * qmd1.slice(m);
          } else if (flag) {
            tree.try_insert(key1, nddo);
            keys(nddo) = hnod(nod.gams - expn(mp), key1);
            ddos.slice(nddo).zeros();
            nddo += 1;
          }
        }
      }

      rhob += OMG * ddos.slice(iddo);
      // Dyson expansion
      for (int i = 0; i < nsys; ++i) {
        for (int j = 0; j < nsys; ++j) {
          rhob(i, j) /= (deom_ci * (ham1(i, i) - ham1(j, j)) + OMG + nod.gams);
        }
      }
      rhoy = rhob;
      int idyson = 0;
      while (idyson < 50) {
        rhox = rhob;
        for (int i = 0; i < nsys; ++i) {
          for (int j = 0; j < nsys; ++j) {
            cx_double tmp = 0.0;
            for (int k = 0; k < nsys; ++k) {
              if (i != k)
                tmp += deom_ci * ham1(i, k) * rhox(k, j);
              if (k != j)
                tmp -= deom_ci * rhox(i, k) * ham1(k, j);
            }
            rhob(i, j) = tmp / (deom_ci * (ham1(i, i) - ham1(j, j)) + OMG + nod.gams);
          }
        }
        rhoy += (idyson % 2 == 0 ? -1.0 : 1.0) * rhob;
        idyson++;
        if (all(abs(vectorise(rhob)) < err * 1.0e-3) && idyson >= 5)
          break;
      }
      // printf ("Dyson terms: %d\n", idyson);
      // update ddos(iddo)
      rho_old = ddos.slice(iddo);
      if (iddo != 0) {
        ddos.slice(iddo) = rhoy;
      } else {
        ddos.slice(0) = 0.5 * (rhoy + rhoy.t());
      }
      double diff = max(vectorise(abs(rho_old - ddos.slice(iddo))));
      max_diff = max_diff > diff ? max_diff : diff;
    }

    iter += 1;

    filter(ddos);

    printf("SCI step %d, nddo=%d, lddo=%d, max_diff=%16.6e\n", iter, nddo, lddo, max_diff);
    fprintf(ferr, "%d %16.6e \n", iter, max_diff);

    if (max_diff < err)
      break;

  } while (iter < ntMax);
  fclose(ferr);

  // real-time propagation
  // EqSolverProp (ddos, dt, err, nk);

  // save rho0
  ddos.slice(0).save("rho0.mat", raw_ascii);
}
