/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "blockdeom.hpp"
#include <cmath>

void blockdeom::oprAct(cx_cube &d_ddos, const mat &sdip, const cx_cube &ddos, const char lrc) {
  for (int iado = 0; iado < nddo; ++iado) {
    const cx_mat &ado = ddos.slice(iado);
    if (lrc == 'l') {
      d_ddos.slice(iado) = sdip * ado;
    } else if (lrc == 'r') {
      d_ddos.slice(iado) = ado * sdip;
    } else if (lrc == 'c') {
      d_ddos.slice(iado) = sdip * ado - ado * sdip;
    }
  }
}

void blockdeom::oprAct(cx_cube &d_ddos, const mat &sdip, const cube &pdip, const vec &bdip, const cx_cube &ddos, const char lrc) {

  const int nsav = nddo;
  d_ddos.slices(0, nddo - 1).zeros();

  for (int iado = 0; iado < nsav; ++iado) {

    const cx_mat &ado = ddos.slice(iado);

    if (iado == 0 || is_valid(ado)) {

      const hnod &nod = keys(iado);
      ivec key0(nod.key);
      int tier = tree.find(key0)->tier;

      if (lrc == 'l') {
        d_ddos.slice(iado) += sdip * ado;
      } else if (lrc == 'r') {
        d_ddos.slice(iado) += ado * sdip;
      } else if (lrc == 'c') {
        d_ddos.slice(iado) += sdip * ado - ado * sdip;
      }

      if (tier < lmax) {
        for (int mp = 0; mp < nind; ++mp) {
          ivec key1 = gen_key(key0, mp, 1);
          const int m = modLabel(mp);
          const int n = key1(mp) - 1;
          cx_double cl = sqrt((n + 1) / coef_abs(mp)) * coef_lft(mp) * bdip(mp);
          cx_double cr = sqrt((n + 1) / coef_abs(mp)) * coef_rht(mp) * bdip(mp);
          if (!tree.try_insert(key1, nddo)) {
            int loc = tree.find(key1)->rank;
            if (lrc == 'l') {
              d_ddos.slice(loc) += cl * pdip.slice(m) * ado;
            } else if (lrc == 'r') {
              d_ddos.slice(loc) += cr * ado * pdip.slice(m);
            } else if (lrc == 'c') {
              d_ddos.slice(loc) += cl * pdip.slice(m) * ado - cr * ado * pdip.slice(m);
            }
          } else {
            keys(nddo) = hnod(nod.gams + expn(mp), key1);
            if (lrc == 'l') {
              d_ddos.slice(nddo) = cl * pdip.slice(m) * ado;
            } else if (lrc == 'r') {
              d_ddos.slice(nddo) = cr * ado * pdip.slice(m);
            } else if (lrc == 'c') {
              d_ddos.slice(nddo) = cl * pdip.slice(m) * ado - cr * ado * pdip.slice(m);
            }
            nddo += 1;
          }
        }
      }

      if (iado > 0) {
        for (int mp = 0; mp < nind; ++mp) {
          ivec key1 = gen_key(key0, mp, -1);
          if (!key1.is_empty()) {
            const int m = modLabel(mp);
            const int n = (key1.n_rows < (unsigned int)(mp + 1)) ? 1 : (key1(mp) + 1);
            cx_double sn = sqrt(n * coef_abs(mp)) * bdip(mp);
            if (!tree.try_insert(key1, nddo)) {
              int loc = tree.find(key1)->rank;
              if (lrc == 'l') {
                d_ddos.slice(loc) += sn * pdip.slice(m) * ado;
              } else if (lrc == 'r') {
                d_ddos.slice(loc) += sn * ado * pdip.slice(m);
              } else if (lrc == 'c') {
                d_ddos.slice(loc) += sn * (pdip.slice(m) * ado - ado * pdip.slice(m));
              }
            } else {
              keys(nddo) = hnod(nod.gams - expn(mp), key1);
              if (lrc == 'l') {
                d_ddos.slice(nddo) = sn * pdip.slice(m) * ado;
              } else if (lrc == 'r') {
                d_ddos.slice(nddo) = sn * ado * pdip.slice(m);
              } else if (lrc == 'c') {
                d_ddos.slice(nddo) = sn * (pdip.slice(m) * ado - ado * pdip.slice(m));
              }
              nddo += 1;
            }
          }
        }
      }
    }
  }
}

void blockdeom::iniSch(cx_cube &ddos, const mat &sdip, const cube &pdip, const vec &bdip, const char lr) {

  ddos.slice(0) = deom_c1 * sdip;

  if (nddo != 1) {
    printf("Error! nddo != 1");
  }

  for (int mp = 0; mp < nind; ++mp) {
    int m = modLabel(mp);
    ivec key0 = zeros<ivec>(mp + 1);
    key0(mp) = 1;
    const cx_double eta = (lr == 'r') ? coef_rht(mp) : coef_lft(mp);
    const cx_double sn = bdip(mp) * eta / sqrt(coef_abs(mp));
    if (!tree.try_insert(key0, nddo)) {
      int loc = tree.find(key0)->rank;
      ddos.slice(loc) += sn * pdip.slice(m);
    } else {
      keys(nddo) = hnod(expn(mp), key0);
      ddos.slice(nddo) = sn * pdip.slice(m);
      nddo += 1;
    }
  }
}

void blockdeom::iniHei(cx_cube &oprs, const mat &sdip) {
  oprs.slice(0) = deom_c1 * sdip;
}

void blockdeom::iniHei(cx_cube &oprs, const mat &sdip, const cube &pdip, const vec &bdip) {

  oprs.slice(0) = deom_c1 * sdip;

  if (nddo != 1) {
    printf("Error! nddo != 1");
  }

  for (int mp = 0; mp < nind; ++mp) {
    int m = modLabel(mp);
    ivec key0 = zeros<ivec>(mp + 1);
    key0(mp) = 1;
    const cx_double sn = sqrt(coef_abs(mp)) * bdip(mp);
    if (!tree.try_insert(key0, nddo)) {
      int loc = tree.find(key0)->rank;
      oprs.slice(loc) += sn * pdip.slice(m);
    } else {
      keys(nddo) = hnod(expn(mp), key0);
      oprs.slice(nddo) = sn * pdip.slice(m);
      nddo += 1;
    }
  }
}

static double delta(const int i, const int j) {
  return (i == j) ? 1.0 : 0.0;
}

cx_mat blockdeom::get_superOpr(const cx_mat &Opr) {
  cx_mat superOpr = zeros<cx_mat>(nsys * nsys, nsys * nsys);
  for (int i = 0; i < nsys; ++i) {
    for (int j = 0; j < nsys; ++j) {
      const int ij = i * nsys + j;
      for (int k = 0; k < nsys; ++k) {
        for (int l = 0; l < nsys; ++l) {
          const int kl = k * nsys + l;
          superOpr(ij, kl) = Opr(i, k) * delta(j, l) - Opr(l, j) * delta(i, k);
        }
      }
    }
  }
  return superOpr;
}
