/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "blockdeom2.hpp"
#include <cmath>

void blockdeom2::oprAct(cx_cube &d_ddos, const mat &sdip, const cx_cube &ddos, const char lrc) {
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

void blockdeom2::iniSch(cx_cube &ddos, const mat &sdip) {
  ddos.slice(0) = deom_c1 * sdip;
  if (nddo != 1) {
    printf("Error! nddo != 1");
  }
}

void blockdeom2::iniHei(cx_cube &oprs, const mat &sdip) {
  oprs.slice(0) = deom_c1 * sdip;
  if (nddo != 1) {
    printf("Error! nddo != 1");
  }
}
