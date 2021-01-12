/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "deom.hpp"
#include <cmath>


void deom::rem_freeze(cx_cube &dtotal, const cx_cube &total, const double t) {

  dtotal.slices(0, nddo - 1).zeros();

  cx_cube qddo(size(qmd1));
  cx_cube ddoq(size(qmd1));

  for (int iado = 0; iado < nddo; ++iado) {
    const cx_mat &ado = total.slice(iado);
    //if (iado==0 || is_valid (ado)) {
    const hnod &nod = keys(iado);
    ivec key0(nod.key);
    int tier = tree.find(key0)->tier;

    dtotal.slice(iado) += -deom_ci * (ham1 * ado - ado * ham1) - nod.gams * ado;
    for (int m = 0; m < nmod; ++m) {
      qddo.slice(m) = qmd1.slice(m) * ado;
      ddoq.slice(m) = ado * qmd1.slice(m);
      if (abs(delt_res(m)) > 1.e-15) {
        dtotal.slice(iado) -= delt_res(m) * (qmd1.slice(m) * (qddo.slice(m) - ddoq.slice(m)) - (qddo.slice(m) - ddoq.slice(m)) * qmd1.slice(m));
      }
    }

    if (tier < lmax) {
      for (int mp = 0; mp < nind; ++mp) {
        ivec key1 = gen_key(key0, mp, 1);
        int m = modLabel(mp);
        if (flag_rwa != 0) { // RWA
          m = (m % 2 == 0) ? (m + 1) : (m - 1);
        }
        const int n = key1(mp) - 1;
        const cx_double sn = -deom_ci * sqrt((n + 1) / coef_abs(mp));
        const cx_double cl = sn * coef_lft(mp);
        const cx_double cr = sn * coef_rht(mp);
        auto *ptr = tree.find(key1);
        if (ptr && ptr->rank >= 0) {
          int loc = ptr->rank;
          dtotal.slice(loc) += cl * qddo.slice(m) - cr * ddoq.slice(m);
        }
      }
    }

    for (int mp = 0; mp < nind; ++mp) {
      ivec key1 = gen_key(key0, mp, -1);
      if (!key1.is_empty()) {
        const int m = modLabel(mp);
        const int n = (key1.n_rows < (unsigned int)(mp + 1)) ? 1 : (key1(mp) + 1);
        const cx_double sn = -deom_ci * sqrt(n * coef_abs(mp));
        auto *ptr = tree.find(key1);
        if (ptr && ptr->rank >= 0) {
          int loc = ptr->rank;
          dtotal.slice(loc) += sn * (qddo.slice(m) - ddoq.slice(m));
        }
      }
    }
    //}
  }
}

void deom::remSch(cx_cube &dtotal, const cx_cube &total, const double t) {

  const int nsav = nddo;
  dtotal.slices(0, nddo - 1).zeros();

  cx_cube qddo(size(qmd1));
  cx_cube ddoq(size(qmd1));

  for (int iado = 0; iado < nsav; ++iado) {
    const cx_mat &ado = total.slice(iado);
    if (iado == 0 || is_valid(ado)) {
      const hnod &nod = keys(iado);
      ivec key0(nod.key);
      int tier = tree.find(key0)->tier;

      dtotal.slice(iado) += -deom_ci * (ham1 * ado - ado * ham1) - nod.gams * ado;
      for (int m = 0; m < nmod; ++m) {
        qddo.slice(m) = qmd1.slice(m) * ado;
        ddoq.slice(m) = ado * qmd1.slice(m);
        if (abs(delt_res(m)) > 1.e-15) {
          dtotal.slice(iado) -= delt_res(m) * (qmd1.slice(m) * (qddo.slice(m) - ddoq.slice(m)) - (qddo.slice(m) - ddoq.slice(m)) * qmd1.slice(m));
        }
      }

      if (tier < lmax) {
        for (int mp = 0; mp < nind; ++mp) {
          ivec key1 = gen_key(key0, mp, 1);
          int m = modLabel(mp);
          if (flag_rwa != 0) { // RWA
            m = (m % 2 == 0) ? (m + 1) : (m - 1);
          }
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

void deom::remHei(cx_cube &dtotal, const cx_cube &total, const double t) {

  const int nsav = nddo;
  dtotal.slices(0, nddo - 1).zeros();

  cx_cube qddo(size(qmd1));
  cx_cube ddoq(size(qmd1));

  for (int iado = 0; iado < nsav; ++iado) {
    const cx_mat &ado = total.slice(iado);
    if (iado == 0 || is_valid(ado)) {
      const hnod &nod = keys(iado);
      ivec key0(nod.key);
      int tier = tree.find(key0)->tier;

      dtotal.slice(iado) += -deom_ci * (ado * ham1 - ham1 * ado) - nod.gams * ado;
      for (int m = 0; m < nmod; ++m) {
        qddo.slice(m) = qmd1.slice(m) * ado;
        ddoq.slice(m) = ado * qmd1.slice(m);
        if (abs(delt_res(m)) > 1.e-15) {
          dtotal.slice(iado) -= delt_res(m) * (qmd1.slice(m) * (qddo.slice(m) - ddoq.slice(m)) - (qddo.slice(m) - ddoq.slice(m)) * qmd1.slice(m));
        }
      }

      for (int mp = 0; mp < nind; ++mp) {
        ivec key1 = gen_key(key0, mp, -1);
        if (!key1.is_empty()) {
          int m = modLabel(mp);
          if (flag_rwa != 0) { // RWA
            m = (m % 2 == 0) ? (m + 1) : (m - 1);
          }
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

void deom::rem(cx_cube &dtotal, const cx_cube &total, const double t, const int &projection) {

  const int nsav = nddo;
  dtotal.slices(0, nddo - 1).zeros();

  cx_cube qddo(size(qmd1));
  cx_cube ddoq(size(qmd1));

  for (int iado = 0; iado < nsav; ++iado) {
    const cx_mat &ado = total.slice(iado);
    if (iado == 0 || is_valid(ado)) {
      const hnod &nod = keys(iado);
      ivec key0(nod.key);
      int tier = tree.find(key0)->tier;

      dtotal.slice(iado) += -deom_ci * (ham1 * ado - ado * ham1) - nod.gams * ado;
      for (int m = 0; m < nmod; ++m) {
        qddo.slice(m) = qmd1.slice(m) * ado;
        ddoq.slice(m) = ado * qmd1.slice(m);
        if (abs(delt_res(m)) > 1.e-15) {
          dtotal.slice(iado) -= delt_res(m) * (qmd1.slice(m) * (qddo.slice(m) - ddoq.slice(m)) - (qddo.slice(m) - ddoq.slice(m)) * qmd1.slice(m));
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

  if (projection != 0) {
    dtotal.slice(0).zeros();
  }
}

void deom::rem(cx_cube &dtotal, const cx_cube &total, const double t, const ivec &projection) {

  const int nsav = nddo;
  dtotal.slices(0, nddo - 1).zeros();

  cx_cube qddo(size(qmd1));
  cx_cube ddoq(size(qmd1));

  for (int iado = 0; iado < nsav; ++iado) {
    const cx_mat &ado = total.slice(iado);
    if (iado == 0 || is_valid(ado)) {
      const hnod &nod = keys(iado);
      ivec key0(nod.key);
      int tier = tree.find(key0)->tier;

      dtotal.slice(iado) += -deom_ci * (ham1 * ado - ado * ham1) - nod.gams * ado;
      for (int m = 0; m < nmod; ++m) {
        qddo.slice(m) = qmd1.slice(m) * ado;
        ddoq.slice(m) = ado * qmd1.slice(m);
        if (abs(delt_res(m)) > 1.e-15) {
          dtotal.slice(iado) -= delt_res(m) * (qmd1.slice(m) * (qddo.slice(m) - ddoq.slice(m)) - (qddo.slice(m) - ddoq.slice(m)) * qmd1.slice(m));
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

  for (int i = 0; i < nsys; ++i) {
    if (projection(i) != 0) {
      dtotal(i, i, 0) = 0.0;
    }
  }
}

void deom::rem(cx_cube &dtotal, const cx_cube &total, const double t, const mat &sdip, const pulse &p, const mat &damp) {

  cx_mat hamt(ham1);
  if (p.on) {
    const double et = p.et(t);
    hamt -= deom_c1 * sdip * et;
  }
  cx_cube qddo(size(qmd1));
  cx_cube ddoq(size(qmd1));

  const int nsav = nddo;
  dtotal.slices(0, nddo - 1).zeros();

  for (int iado = 0; iado < nsav; ++iado) {
    const cx_mat &ado = total.slice(iado);
    if (iado == 0 || is_valid(ado)) {
      const hnod &nod = keys(iado);
      ivec key0(nod.key);
      int tier = tree.find(key0)->tier;

      dtotal.slice(iado) += -deom_ci * (hamt * ado - ado * hamt) - nod.gams * ado;
      for (int m = 0; m < nmod; ++m) {
        qddo.slice(m) = qmd1.slice(m) * ado;
        ddoq.slice(m) = ado * qmd1.slice(m);
        if (abs(delt_res(m)) > 1.e-15) {
          dtotal.slice(iado) -= delt_res(m) *
                                (qmd1.slice(m) * (qddo.slice(m) - ddoq.slice(m)) - (qddo.slice(m) - ddoq.slice(m)) * qmd1.slice(m));
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

  if (!damp.is_empty()) {
    dtotal.slice(0) -= damp * total.slice(0) + total.slice(0) * damp;
  }
}

void deom::rem(cx_cube &dtotal, const cx_cube &total, const double t, const mat &sdip, const cube &pdip, const vec &bdip, const pulse &p) {

  cx_mat hamt(ham1);
  cx_cube qmdt = zeros<cx_cube>(nsys, nsys, nind);
  if (p.on) {
    const double et = p.et(t);
    hamt -= deom_c1 * sdip * et;
    for (int mp = 0; mp < nind; ++mp) {
      int m = modLabel(mp);
      qmdt.slice(mp) = qmd1.slice(m) - pdip.slice(m) * bdip(mp) * et * deom_c1;
    }
  }
  cx_cube qddo(size(qmdt));
  cx_cube ddoq(size(qmdt));

  const int nsav = nddo;
  dtotal.slices(0, nddo - 1).zeros();

  for (int iado = 0; iado < nsav; ++iado) {
    const cx_mat &ado = total.slice(iado);
    if (iado == 0 || is_valid(ado)) {
      const hnod &nod = keys(iado);
      ivec key0(nod.key);
      int tier = tree.find(key0)->tier;

      dtotal.slice(iado) += -deom_ci * (hamt * ado - ado * hamt) - nod.gams * ado;
      for (int mp = 0; mp < nind; ++mp) {
        int m = modLabel(mp);
        qddo.slice(mp) = qmdt.slice(mp) * ado;
        ddoq.slice(mp) = ado * qmdt.slice(mp);
        if (abs(delt_res(m)) > 1.e-15) {
          dtotal.slice(iado) -= delt_res(m) *
                                (qmdt.slice(mp) * (qddo.slice(mp) - ddoq.slice(mp)) - (qddo.slice(mp) - ddoq.slice(mp)) * qmdt.slice(mp));
        }
      }

      if (tier < lmax) {
        for (int mp = 0; mp < nind; ++mp) {
          ivec key1 = gen_key(key0, mp, 1);
          const int n = key1(mp) - 1;
          const cx_double sn = -deom_ci * sqrt((n + 1) / coef_abs(mp));
          const cx_double cl = sn * coef_lft(mp);
          const cx_double cr = sn * coef_rht(mp);
          if (!tree.try_insert(key1, nddo)) {
            int loc = tree.find(key1)->rank;
            dtotal.slice(loc) += cl * qddo.slice(mp) - cr * ddoq.slice(mp);
          } else {
            keys(nddo) = hnod(nod.gams + expn(mp), key1);
            dtotal.slice(nddo) = cl * qddo.slice(mp) - cr * ddoq.slice(mp);
            nddo += 1;
          }
        }
      }

      for (int mp = 0; mp < nind; ++mp) {
        ivec key1 = gen_key(key0, mp, -1);
        if (!key1.is_empty()) {
          const int n = (key1.n_rows < (unsigned int)(mp + 1)) ? 1 : (key1(mp) + 1);
          const cx_double sn = -deom_ci * sqrt(n * coef_abs(mp));
          if (!tree.try_insert(key1, nddo)) {
            int loc = tree.find(key1)->rank;
            dtotal.slice(loc) += sn * (qddo.slice(mp) - ddoq.slice(mp));
          } else {
            keys(nddo) = hnod(nod.gams - expn(mp), key1);
            dtotal.slice(nddo) = sn * (qddo.slice(mp) - ddoq.slice(mp));
            nddo += 1;
          }
        }
      }
    }
  }
}

void deom::rem(cx_cube &dtotal, const cx_cube &total, const double t, const mat &sdip, const cube &pdip, const vec &bdip, const pulse &p, const ivec &projection) {

  cx_mat hamt(ham1);
  cx_cube qmdt(nsys, nsys, nind);
  if (p.on) {
    const double et = p.et(t);
    hamt -= deom_c1 * sdip * et;
    for (int mp = 0; mp < nind; ++mp) {
      int m = modLabel(mp);
      qmdt.slice(mp) = qmd1.slice(m) - pdip.slice(m) * bdip(mp) * et * deom_c1;
    }
  }
  cx_cube qddo(size(qmdt));
  cx_cube ddoq(size(qmdt));

  const int nsav = nddo;
  dtotal.slices(0, nddo - 1).zeros();

  for (int iado = 0; iado < nsav; ++iado) {
    const cx_mat &ado = total.slice(iado);
    if (iado == 0 || is_valid(ado)) {
      const hnod &nod = keys(iado);
      ivec key0(nod.key);
      int tier = tree.find(key0)->tier;

      dtotal.slice(iado) += -deom_ci * (hamt * ado - ado * hamt) - nod.gams * ado;
      for (int mp = 0; mp < nind; ++mp) {
        int m = modLabel(mp);
        qddo.slice(mp) = qmdt.slice(mp) * ado;
        ddoq.slice(mp) = ado * qmdt.slice(mp);
        if (abs(delt_res(m)) > 1.e-15) {
          dtotal.slice(iado) -= delt_res(m) *
                                (qmdt.slice(mp) * (qddo.slice(mp) - ddoq.slice(mp)) - (qddo.slice(mp) - ddoq.slice(mp)) * qmdt.slice(mp));
        }
      }

      if (tier < lmax) {
        for (int mp = 0; mp < nind; ++mp) {
          ivec key1 = gen_key(key0, mp, 1);
          const int n = key1(mp) - 1;
          const cx_double sn = -deom_ci * sqrt((n + 1) / coef_abs(mp));
          const cx_double cl = sn * coef_lft(mp);
          const cx_double cr = sn * coef_rht(mp);
          if (!tree.try_insert(key1, nddo)) {
            int loc = tree.find(key1)->rank;
            dtotal.slice(loc) += cl * qddo.slice(mp) - cr * ddoq.slice(mp);
          } else {
            keys(nddo) = hnod(nod.gams + expn(mp), key1);
            dtotal.slice(nddo) = cl * qddo.slice(mp) - cr * ddoq.slice(mp);
            nddo += 1;
          }
        }
      }

      for (int mp = 0; mp < nind; ++mp) {
        ivec key1 = gen_key(key0, mp, -1);
        if (!key1.is_empty()) {
          const int n = (key1.n_rows < (unsigned int)(mp + 1)) ? 1 : (key1(mp) + 1);
          const cx_double sn = -deom_ci * sqrt(n * coef_abs(mp));
          if (!tree.try_insert(key1, nddo)) {
            int loc = tree.find(key1)->rank;
            dtotal.slice(loc) += sn * (qddo.slice(mp) - ddoq.slice(mp));
          } else {
            keys(nddo) = hnod(nod.gams - expn(mp), key1);
            dtotal.slice(nddo) = sn * (qddo.slice(mp) - ddoq.slice(mp));
            nddo += 1;
          }
        }
      }
    }
  }

  for (int i = 0; i < nsys; ++i) {
    if (projection(i) != 0) {
      dtotal(i, i, 0) = 0.0;
    }
  }
}

void deom::rem(cx_cube &dtotal, const cx_cube &total, const double t, const double kappa, const double I0) {

  cx_mat opr = zeros<cx_mat>(2, 2);
  opr(0, 0) = 1 + kappa;
  opr(1, 1) = 1 - kappa;

  const int nsav = nddo;
  dtotal.slices(0, nddo - 1).zeros();

  cx_cube qddo(size(qmd1));
  cx_cube ddoq(size(qmd1));

  for (int iado = 0; iado < nsav; ++iado) {
    const cx_mat &ado = total.slice(iado);
    if (iado == 0 || is_valid(ado)) {
      const hnod &nod = keys(iado);
      ivec key0(nod.key);
      int tier = tree.find(key0)->tier;

      dtotal.slice(iado) += -deom_ci * (ham1 * ado - ado * ham1) - nod.gams * ado;
      dtotal.slice(iado) += -0.5 * I0 * (opr * opr * ado + ado * opr * opr - 2.0 * opr * ado * opr);

      for (int m = 0; m < nmod; ++m) {
        qddo.slice(m) = qmd1.slice(m) * ado;
        ddoq.slice(m) = ado * qmd1.slice(m);
        if (abs(delt_res(m)) > 1.0e-15) {
          dtotal.slice(iado) -= delt_res(m) * (qmd1.slice(m) * (qddo.slice(m) - ddoq.slice(m)) - (qddo.slice(m) - ddoq.slice(m)) * qmd1.slice(m));
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

void deom::remHSB(cx_cube &dtotal, const cx_cube &total, const double t) {

  const int nsav = nddo;
  dtotal.slices(0, nddo - 1).zeros();

  cx_cube qddo(size(qmd1));

  for (int iado = 0; iado < nsav; ++iado) {
    const cx_mat &ado = total.slice(iado);
    if (iado == 0 || is_valid(ado)) {
      const hnod &nod = keys(iado);
      ivec key0(nod.key);
      int tier = tree.find(key0)->tier;

      for (int m = 0; m < nmod; ++m) {
        qddo.slice(m) = qmd1.slice(m) * ado;
        if (abs(delt_res(m)) > 1.e-15) {
          dtotal.slice(iado) -= delt_res(m) * (qmd1.slice(m) * qddo.slice(m) - qddo.slice(m) * qmd1.slice(m));
        }
      }

      if (tier < lmax) {
        for (int mp = 0; mp < nind; ++mp) {
          ivec key1 = gen_key(key0, mp, 1);
          int m = modLabel(mp);
          if (flag_rwa != 0) { // RWA
            m = (m % 2 == 0) ? (m + 1) : (m - 1);
          }
          const int n = key1(mp) - 1;
          const cx_double sn = sqrt((n + 1) / coef_abs(mp));
          const cx_double cl = sn * coef_lft(mp);
          if (!tree.try_insert(key1, nddo)) {
            int loc = tree.find(key1)->rank;
            dtotal.slice(loc) += cl * qddo.slice(m);
          }
        }
      }

      for (int mp = 0; mp < nind; ++mp) {
        ivec key1 = gen_key(key0, mp, -1);
        if (!key1.is_empty()) {
          const int m = modLabel(mp);
          const int n = (key1.n_rows < (unsigned int)(mp + 1)) ? 1 : (key1(mp) + 1);
          const cx_double sn = sqrt(n * coef_abs(mp));
          if (!tree.try_insert(key1, nddo)) {
            int loc = tree.find(key1)->rank;
            dtotal.slice(loc) += sn * qddo.slice(m);
          }
        }
      }
    }
  }
}

void deom::rem_diag(cx_cube &dtotal, const cx_cube &total, const double t) {

  cx_cube qddo(size(qmd1));
  cx_cube ddoq(size(qmd1));

  for (int iado = 0; iado < nddo; ++iado) {
    const cx_mat &ado = total.slice(iado);
    const hnod &nod = keys(iado);
    dtotal.slice(iado) = -deom_ci * (ham1 * ado - ado * ham1) - nod.gams * ado;
    for (int m = 0; m < nmod; ++m) {
      qddo.slice(m) = qmd1.slice(m) * ado;
      ddoq.slice(m) = ado * qmd1.slice(m);
      if (abs(delt_res(m)) > 1.e-15) {
        dtotal.slice(iado) -= delt_res(m) * (qmd1.slice(m) * (qddo.slice(m) - ddoq.slice(m)) - (qddo.slice(m) - ddoq.slice(m)) * qmd1.slice(m));
      }
    }
  }
}

void deom::rem_offdiag(cx_cube &dtotal, const cx_cube &total, const double t) {

  const int nsav = nddo;
  dtotal.slices(0, nddo - 1).zeros();

  cx_cube qddo(size(qmd1));
  cx_cube ddoq(size(qmd1));

  for (int iado = 0; iado < nsav; ++iado) {
    const cx_mat &ado = total.slice(iado);
    if (iado == 0 || is_valid(ado)) {
      const hnod &nod = keys(iado);
      ivec key0(nod.key);
      int tier = tree.find(key0)->tier;

      for (int m = 0; m < nmod; ++m) {
        qddo.slice(m) = qmd1.slice(m) * ado;
        ddoq.slice(m) = ado * qmd1.slice(m);
      }

      if (tier < lmax) {
        for (int mp = 0; mp < nind; ++mp) {
          ivec key1 = gen_key(key0, mp, 1);
          int m = modLabel(mp);
          if (flag_rwa != 0) { // RWA
            m = (m % 2 == 0) ? (m + 1) : (m - 1);
          }
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
