/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "blockdeom.hpp"
#include <cmath>

void blockdeom::remImag(cx_cube &dtotal, const cx_cube &total, const double t) {

  const int nsav = nddo;
  dtotal.slices(0, nddo - 1).zeros();

  cx_cube qddo(size(qmd1));

  for (int iado = 0; iado < nsav; ++iado) {
    const cx_mat &ado = total.slice(iado);
    if (iado == 0 || is_valid(ado)) {
      const hnod &nod = keys(iado);
      ivec key0(nod.key);
      int tier = tree.find(key0)->tier;

      dtotal.slice(iado) += -ham1 * ado + deom_ci * nod.gams * ado;
      for (int m = 0; m < nmod; ++m) {
        qddo.slice(m) = qmd1.slice(m) * ado;
        //if (abs(delt_res(m)) > 1.0e-15) {
        //    dtotal.slice(iado) += deom_ci*delt_res(m)*qmd1.slice(m)*qddo.slice(m);
        //}
      }

      if (tier < lmax) {
        for (int mp = 0; mp < nind; ++mp) {
          ivec key1 = gen_key(key0, mp, 1);
          const int m = modLabel(mp);
          const int n = key1(mp) - 1;
          const cx_double sn = -sqrt((n + 1) / coef_abs(mp)) * coef_lft(mp);
          if (!tree.try_insert(key1, nddo)) {
            int loc = tree.find(key1)->rank;
            dtotal.slice(loc) += sn * qddo.slice(m);
          } else {
            keys(nddo) = hnod(nod.gams + expn(mp), key1);
            dtotal.slice(nddo) = sn * qddo.slice(m);
            nddo += 1;
          }
        }
      }

      for (int mp = 0; mp < nind; ++mp) {
        ivec key1 = gen_key(key0, mp, -1);
        if (!key1.is_empty()) {
          const int m = modLabel(mp);
          const int n = (key1.n_rows < (unsigned int)(mp + 1)) ? 1 : (key1(mp) + 1);
          const cx_double sn = -sqrt(n * coef_abs(mp));
          if (!tree.try_insert(key1, nddo)) {
            int loc = tree.find(key1)->rank;
            dtotal.slice(loc) += sn * qddo.slice(m);
          } else {
            keys(nddo) = hnod(nod.gams - expn(mp), key1);
            dtotal.slice(nddo) = sn * qddo.slice(m);
            nddo += 1;
          }
        }
      }
    }
  }
}

void blockdeom::remSch(cx_cube &dtotal, const cx_cube &total, const double t) {

  const int nsav = nddo;
  dtotal.slices(0, nddo - 1).zeros();

  cx_cube qddo(nl, nr, nmod);
  cx_cube ddoq(nl, nr, nmod);

  for (int iado = 0; iado < nsav; ++iado) {
    const cx_mat &ado = total.slice(iado);
    if (iado == 0 || is_valid(ado)) {
      const hnod &nod = keys(iado);
      ivec key0(nod.key);
      int tier = tree.find(key0)->tier;

      dtotal.slice(iado) += -deom_ci * (haml * ado - ado * hamr) - nod.gams * ado;
      for (int m = 0; m < nmod; ++m) {
        qddo.slice(m) = qmdl.slice(m) * ado;
        ddoq.slice(m) = ado * qmdr.slice(m);
        if (abs(delt_res(m)) > 1.e-15) {
          dtotal.slice(iado) -= delt_res(m) * (qmdl.slice(m) * (qddo.slice(m) - ddoq.slice(m)) - (qddo.slice(m) - ddoq.slice(m)) * qmdr.slice(m));
        }
      }

      if (tier < lmax) {
        for (int mp = 0; mp < nind; ++mp) {
          ivec key1 = gen_key(key0, mp, 1);
          const int m = modLabel(mp);
          const int n = key1(mp) - 1;
          const cx_double sn = -deom_ci * sqrt((n + 1) / coef_abs(mp));
          const cx_double cl = sn * coef_lft(mp);
          const cx_double cr = sn * coef_rht(mp);
          if (!tree.try_insert(key1, nddo)) {
            int loc = tree.find(key1)->rank;
            dtotal.slice(loc) += cl * qddo.slice(m) - cr * ddoq.slice(m);
          } else {
            keys(nddo) = hnod(nod.gams + expn(mp), key1);
            dtotal.slice(nddo) = cl * qddo.slice(m) - cr * ddoq.slice(m);
            nddo += 1;
          }
        }
      }

      for (int mp = 0; mp < nind; ++mp) {
        ivec key1 = gen_key(key0, mp, -1);
        if (!key1.is_empty()) {
          const int m = modLabel(mp);
          const int n = (key1.n_rows < (unsigned int)(mp + 1)) ? 1 : (key1(mp) + 1);
          const cx_double sn = -deom_ci * sqrt(n * coef_abs(mp));
          if (!tree.try_insert(key1, nddo)) {
            int loc = tree.find(key1)->rank;
            dtotal.slice(loc) += sn * (qddo.slice(m) - ddoq.slice(m));
          } else {
            keys(nddo) = hnod(nod.gams - expn(mp), key1);
            dtotal.slice(nddo) = sn * (qddo.slice(m) - ddoq.slice(m));
            nddo += 1;
          }
        }
      }
    }
  }
}

void blockdeom::remHei(cx_cube &dtotal, const cx_cube &total, const double t) {

  const int nsav = nddo;
  dtotal.slices(0, nddo - 1).zeros();

  cx_cube qddo = zeros<cx_cube>(nl, nr, nmod);
  cx_cube ddoq = zeros<cx_cube>(nl, nr, nmod);

  for (int iado = 0; iado < nsav; ++iado) {
    const cx_mat &ado = total.slice(iado);
    if (iado == 0 || is_valid(ado)) {
      const hnod &nod = keys(iado);
      ivec key0(nod.key);
      int tier = tree.find(key0)->tier;

      dtotal.slice(iado) += -deom_ci * (ado * hamr - haml * ado) - nod.gams * ado;
      for (int m = 0; m < nmod; ++m) {
        qddo.slice(m) = qmdl.slice(m) * ado;
        ddoq.slice(m) = ado * qmdr.slice(m);
        if (abs(delt_res(m)) > 1.e-15) {
          dtotal.slice(iado) -= delt_res(m) * (qmdl.slice(m) * (qddo.slice(m) - ddoq.slice(m)) - (qddo.slice(m) - ddoq.slice(m)) * qmdr.slice(m));
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
          if (!tree.try_insert(key1, nddo)) {
            int loc = tree.find(key1)->rank;
            dtotal.slice(loc) += cl * ddoq.slice(m) - cr * qddo.slice(m);
          } else {
            keys(nddo) = hnod(nod.gams - expn(mp), key1);
            dtotal.slice(nddo) = cl * ddoq.slice(m) - cr * qddo.slice(m);
            nddo += 1;
          }
        }
      }

      if (tier < lmax) {
        for (int mp = 0; mp < nind; ++mp) {
          ivec key1 = gen_key(key0, mp, 1);
          const int m = modLabel(mp);
          const int n = key1(mp) - 1;
          const cx_double sn = -deom_ci * sqrt((n + 1) * coef_abs(mp));
          if (!tree.try_insert(key1, nddo)) {
            int loc = tree.find(key1)->rank;
            dtotal.slice(loc) += sn * (ddoq.slice(m) - qddo.slice(m));
          } else {
            keys(nddo) = hnod(nod.gams + expn(mp), key1);
            dtotal.slice(nddo) = sn * (ddoq.slice(m) - qddo.slice(m));
            nddo += 1;
          }
        }
      }
    }
  }
}
