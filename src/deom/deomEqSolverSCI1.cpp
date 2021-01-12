/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "deom.hpp"
#include <cstdio>
#include <cstdlib>


static const int ntMax = 100000000;

static double norm(const int nddo, const cx_cube &x) {
  const int n1 = nddo, n2 = x.n_rows, n3 = x.n_cols;
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

void deom::EqSolverSCI1(cx_cube &ddos, const double dt, const double err, const int nk) {

  if (nddo != 1) {
    printf("Warning: Hidx is not correctly initialized!\n");
  }
  vec eigh = eig_sym(ham1);
  const double OMG = max(abs(eigh)) * 1.01;
  printf("OMG=%16.6e\n", OMG / 4.5554927e-6);

  if (ddos.is_empty()) {
    ddos = zeros<cx_cube>(nsys, nsys, nmax);
  } else {
    ddos.zeros();
  }
  mat eham = expmat_sym(-real(ham1) / temperature);
  ddos.slice(0) = eham / trace(eham) * deom_c1;

  cx_mat iLs = deom_ci * get_superOpr(ham1);
  for (int m = 0; m < nmod; ++m) {
    if (abs(delt_res(m)) > 1.0e-15) {
      cx_mat Qs = get_superOpr(qmd1.slice(m));
      iLs += delt_res(m) * Qs * Qs;
    }
  }

  cx_cube ddos_diff = zeros<cx_cube>(size(ddos));
  cx_mat rho_old(nsys, nsys);
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

      if (nod.info == 0 || iter % 5 == 0) {
        cx_mat rhob = zeros<cx_mat>(nsys, nsys);
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

        rho_old = ddos.slice(iddo);
        cx_mat iL(iLs);
        if (iddo != 0) {
          iL.diag() += nod.gams + OMG;
          ddos.slice(iddo) = reshape(solve(iL, vectorise(rhob.st())), nsys, nsys).st();
        } else {
          for (int is = 0; is < nsys; ++is) {
            iL(0, is * nsys + is) += 1.0;
          }
          iL.diag() += OMG;
          rhob(0, 0) += 1.0;
          cx_vec rx = solve(iL, vectorise(rhob.st()));
          for (int is = 0; is < nsys; ++is) {
            ddos(is, is, 0) = real(rx(is * nsys + is));
            for (int js = is + 1; js < nsys; ++js) {
              cx_double tmp1 = rx(is * nsys + js);
              cx_double tmp2 = rx(js * nsys + is);
              ddos(is, js, 0) = 0.5 * real(tmp1 + tmp2) + 0.5 * deom_ci * imag(tmp1 - tmp2);
              ddos(js, is, 0) = 0.5 * real(tmp1 + tmp2) - 0.5 * deom_ci * imag(tmp1 - tmp2);
            }
          }
        }
        double diff = max(vectorise(abs(rho_old - ddos.slice(iddo))));
        if (diff < 0.1 * err) {
          nod.info = 1;
        }
        max_diff = max_diff > diff ? max_diff : diff;
      }
    }

    iter += 1;

    filter(ddos);

    rem_freeze(ddos_diff, ddos, 0.0);
    double error = norm(nddo, ddos_diff); ///fmax(ddos_diff)*41.34902;

    printf("SCI step %d, nddo=%d, lddo=%d, max_diff=%16.6e\n", iter, nddo, lddo, max_diff);
    fprintf(ferr, "%d %16.6e %16.6e\n", iter, error, max_diff);

    if (max_diff < err)
      break;

  } while (iter < ntMax);
  fclose(ferr);

  //int nddo_backup(nddo);
  //hkey keys_backup(keys);
  //cx_cube ddos_backup(ddos);

  int it = 0;
  //FILE *frho = fopen("equilibrium.log","w");
  //do {
  //    double t = it*dt;
  //    printf("Propagate to equilibrium: step %d, nddo=%d, lddo=%d\n", it, nddo, lddo);
  //    if ((it+1)%nk == 0) {
  //        fprintf (frho, "%16.6e", t/deom_fs2unit);
  //        for (int i=0; i<nsys; ++i) {
  //            fprintf(frho, "%16.6e", real(ddos(i,i,0)));
  //        }
  //        fprintf(frho, "\n");

  //        if (notconverged(nddo_backup,keys_backup,ddos_backup,ddos,err)) {
  //            nddo_backup = nddo;
  //            keys_backup.rows(0,nddo-1) = keys.rows(0,nddo-1);
  //            ddos_backup.head_slices(nddo) = ddos.head_slices(nddo);
  //        } else {
  //            break;
  //        }
  //    }
  //    rk4 (ddos,t,dt);
  //    it += 1;
  //} while (it < ntMax);
  //fclose(frho);

  if (it >= ntMax) {
    printf("Fail to reach equilibrium at error %16.6e\n", err);
  } else {
    printf("Equilibrium reached with %d steps!\n", it);
  }
}
