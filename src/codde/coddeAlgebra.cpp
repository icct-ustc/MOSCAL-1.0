/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "codde.hpp"

void codde::oprAct (cx_cube& d_ddos, const mat& sdip, const cx_cube& ddos, const char lrc) {

    d_ddos.slices(0,nddo-1).zeros();

    for (int iado=0; iado<nddo; ++iado) {
        const cx_mat& ado = ddos.slice(iado);
        if (lrc == 'l') {
            d_ddos.slice(iado) += sdip*ado;
        } else if (lrc == 'r') {
            d_ddos.slice(iado) += ado*sdip;
        } else if (lrc == 'c') {
            d_ddos.slice(iado) += sdip*ado-ado*sdip;
        }
    }

    if (lrc=='l' or lrc=='c') {
        for (int k=0; k<nind; ++k) {
            const iado = k+1;
            d_ddos.slice(iado) += (sdip*qhat.slice(k)-qhat.slice(k)*sdip)*ddos.slice(0);
        }
    }

    if (lrc=='r' or lrc=='c') { 
        for (int k=0; k<nind; ++k) {
            const iado = k+nind+1;
            d_ddos.slice(iado) += ddos.slice(0)*(sdip*qhat.slice(k)-qhat.slice(k)*sdip);
        }
    }
}


void codde::iniHei (cx_cube& oprs, const mat& sdip) {
    oprs.zeros();
    oprs.slice(0) = sdip*deom_c1;
}
