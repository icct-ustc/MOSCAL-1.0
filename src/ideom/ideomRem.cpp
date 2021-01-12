/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "ideom.hpp"
#include <cmath>


void ideom::remSch(cx_cube &dtotal, const cx_cube &total, const double t) {
  printf("irem\n");
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

      dtotal.slice(iado) += -(ham1 * ado - ado * ham1) + deom_ci * nod.gams * ado;
      for (int m = 0; m < nmod; ++m) {
        qddo.slice(m) = qmd1.slice(m) * ado;
        ddoq.slice(m) = ado * qmd1.slice(m);
        if (abs(delt_res(m)) > 1.0e-15) {
          dtotal.slice(iado) += deom_ci * delt_res(m) * (qmd1.slice(m) * (qddo.slice(m) - ddoq.slice(m)) - (qddo.slice(m) - ddoq.slice(m)) * qmd1.slice(m));
        }
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
