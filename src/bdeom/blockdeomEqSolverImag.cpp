/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "blockdeom.hpp"
#include <cstdio>
#include <cstdlib>

void blockdeom::EqSolverImag(cx_cube &ddos, const double dt, const double err, const int nk) {

  if (nddo != 1) {
    printf("Warning: Hidx is not correctly initialized!\n");
  }

  cx_cube rhot = zeros<cx_cube>(nsys, nsys, nmax);
  rhot.slice(0) = deom_c1 * eye<mat>(nsys, nsys);

  const int nt = 100 * static_cast<int>(1.0 / (temperature * dt));
  const double dtau = 1.0 / (temperature * nt);

  // Imaginary-time propagation
  for (int it = 0; it < nt; ++it) {
    const double tau = it * dtau;
    if (it % nk == 0) {
      printf("Propagation %5.1f%%: nddo=%6d, lddo=%3d\n",
             100 * it / static_cast<double>(nt), nddo, lddo);
    }
    rk4(rhot, tau, dtau, 'i');
  }

  // partition function
  double za = real(trace(rhot.slice(0)));
  printf("Partition function\n");
  printf("za = %16.6e\n", za);

  // Hermitian
  for (int iddo = 0; iddo < nddo; ++iddo) {
    const hnod &nod = keys(iddo);
    ivec key0(nod.key);
    if (abs(imag(nod.gams)) < 1.0e-20) {
      ddos.slice(iddo) = (rhot.slice(iddo) + rhot.slice(iddo).t()) / (2.0 * za);
    } else {
      ivec key1 = gen_nbar(key0);
      auto *ptr = tree.find(key1);
      if (ptr && ptr->rank >= 0) {
        int jddo = ptr->rank;
        ddos.slice(iddo) = (rhot.slice(iddo) + rhot.slice(jddo).t()) / (2.0 * za);
      } else { //TODO: current set the unpaired ddos to zero
        // ddos.slice(iddo) = rhot.slice(iddo)/(2.0*za);
        ddos.slice(iddo).zeros();
      }
    }
  }

  // real-time propagation
  // EqSolverProp (ddos, dt, err, nk);

  // save rho0
  ddos.slice(0).save("rho0.mat", raw_ascii);
}
