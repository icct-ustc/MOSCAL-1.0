#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "deom.hpp"

static double entropy (const cx_mat& rho) {
    vec eval = eig_sym(rho);
    if (any(eval<0)) {
        return -9520;
    } else {
        return -dot(eval,log(eval));
    }
}

int main () {

    ifstream jsonFile("input.json");
    stringstream strStream;
    strStream << jsonFile.rdbuf();
    string jsonStr = strStream.str();
    string err;

    const Json json = Json::parse(jsonStr,err);
    if (!err.empty()) {
        printf ("Error in parsing input file: %s\n", err.c_str());
        return 0;
    }

    deom d(json["deom"]);

    const int inistate = json["rhot"]["inistate"].int_value();
    const int    nt = json["rhot"]["nt"].int_value();
    const double dt = json["rhot"]["dt"].number_value();
    const int    nk = json["rhot"]["nk"].int_value();

    cx_cube ddos = zeros<cx_cube>(d.nsys,d.nsys,d.nmax);
    ddos(inistate,inistate,0) = 1.0;

    FILE *frho = fopen("prop-rho.dat","w");
    FILE *mdjz = fopen("mdjz.dat","w");
    FILE *fent = fopen("entropy.dat","w");
    FILE *ef = fopen("ef.dat","w");

    cx_cube d_ddos(ddos);

    for (int it=0; it<nt; ++it) {
        const double t = it*dt;
        printf ("Propagation %5.1f%%: nddo=%6d, lddo=%3d\n",
                100*it/static_cast<double>(nt), d.nddo, d.lddo);

        // write rho0_diag
        fprintf (frho, "%16.6e\n", t/deom_fs2unit);
        for (int i=0; i<d.nsys; ++i) {
            fprintf (frho, "%20.10e\n", real(ddos(i,i,0)));
        }
        fprintf (frho, "\n");

        // write rho0
        for (int i=0; i<d.nsys; ++i) {
            for (int j=0; j<d.nsys; ++j) {
                fprintf(mdjz,"(%lf,%lf)\t",real(ddos.slice(0)(i,j)),imag(ddos.slice(0)(i,j)));
            }
            fprintf(mdjz,"\n");
        }
        fprintf(mdjz,"#######################################\n");

        // write entropy
        double S1 = entropy(ddos.slice(0));
        printf("this is entropy: %lf\n", S1);
        if (S1 < 0) {
            fprintf (fent, "%16.6e nan\n", t/deom_fs2unit);
        } else {
            fprintf (fent, "%16.6e%16.6e\n", t/deom_fs2unit,S1);
        }
        // get d\rho/dt
        d.rem_freeze(d_ddos,ddos,t);

        d.rk4 (ddos,t,dt);

        // dS
        // dS1 = -tr[d\rho/dt (\ln\rho)]
        double dS1 = -real(trace(d_ddos.slice(0)*logmat_sympd(ddos.slice(0))));
        // dS2
        double dS2 = real(trace(d.ham1*d_ddos.slice(0)));
        fprintf(ef,"%16.6e%16.6e%16.6e%16.6e\n", t/deom_fs2unit, dS1, dS2, dS1-dS2);

    }
    fclose (fent);
    fclose (frho);
    fclose (ef);
    fclose (mdjz);

    return 0;
}
