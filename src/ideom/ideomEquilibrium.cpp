/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cstdio>
#include <cstdlib>
#include "ideom.hpp"

void ideom::equilibrium (cx_cube& ddos, const double dt, const double err, const int nk) {

    if (nddo != 1) {
        printf("Warning: Hidx is not correctly initialized!\n");
    }

    //ddos.zeros();
    //ddos.slice(0) = deom_c1*eye<mat>(nsys,nsys);
    //cx_cube ddos = zeros<cx_cube>(d.nsys,d.nsys,d.nmax);
    mat eham = expmat_sym(-real(ham1)/temperature);
    ddos.slice(0) = eham/trace(eham)*deom_c1;


    const int nt = floor(1.0/(temperature*dt));
    //const double nt = 1.0/temperature*dt;

    // Imaginary-time propagation
    //double tau = 0;
    FILE *fs = fopen("rhot_img.log","w");
    FILE *fs1 = fopen("zazb.log","w");
    for (int it=0; it<=nt; ++it) {
        const double tau = it*dt;
        printf ("Propagation %5.1f%%: nddo=%6d, lddo=%3d\n",
                100*it/static_cast<double>(nt), nddo, lddo);
        if (it%nk == 0) {
            fprintf(fs,"%26.16e",tau);
	    fprintf (fs, "\n");
	    for (int i=0; i<1; ++i){
	        for (int j=0; j<nsys; ++j)
		    for (int k=0; k<nsys; ++k)
			fprintf (fs, "(%26.16e,%26.16e)\t",ddos(j,k,i).real(),ddos(j,k,i).imag());
		fprintf (fs, "\n\n");
            }
	}
        rk4 (ddos,tau,dt);
        //tau += dt;
        //double za=real(trace(ddos.slice(0)));
        double za = real(trace(ddos.slice(0)));
        double zb = imag(trace(ddos.slice(0)));
        fprintf (fs1, "%26.16e%26.16e%26.16e\n", tau, za, zb);
    }
    fclose(fs);
    fclose(fs1);

    // save rho0
    ddos.slice(0).save("rho0.imag",raw_ascii);
    // save rhon
    fs = fopen("rhonimag.dat","w");
    for (int i=0; i<nddo; ++i){
        for (size_t j=0; j<keys(i).key.n_rows; ++j)
            fprintf(fs,"%d\t",(int)keys(i).key(j));
        fprintf (fs, "\n");
        for (int j=0; j<nsys; ++j)
            for (int k=0; k<nsys; ++k)
                fprintf (fs, "(%26.16e,%26.16e)\t",ddos(j,k,i).real(),ddos(j,k,i).imag());  //%16.6e
        fprintf (fs, "\n\n");
    }
    fclose(fs);

}
