/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "deom.hpp"

static const int ntMax = 100000000;

// return (x^n/n!)e^{-x}
static double poisson (const double x, const int n) {
    double result = exp(-x);
    for (int i=1; i<=n; ++i) {
        result *= x/static_cast<double>(i);
    }
    return result;
}

// return (x^n/n!)e^{-y}
static double poisson2 (const double x, const double y, const int n) {
    double result = exp(-y);
    for (int i=1; i<=n; ++i) {
        result *= x/static_cast<double>(i);
    }
    return result;
}

static void collapse (deom& d, cx_cube& ddos, const int n, const double kappa, const double I0, const double tau) {
    double P11 = poisson((1.0+kappa)*(1.0+kappa)*I0*tau,n);
    double P22 = poisson((1.0-kappa)*(1.0-kappa)*I0*tau,n);
    double P12 = poisson2((1.0+kappa)*(1.0-kappa)*I0*tau,(1.0+kappa*kappa)*I0*tau,n);
    double pIv = 1.0/(P11*real(ddos(0,0,0))+P22*real(ddos(1,1,0)));
    // update density operators after weak measure
    for (int iddo=0; iddo<d.nddo; ++iddo) {
        ddos(0,0,iddo) *= pIv*P11;
        ddos(1,0,iddo) *= pIv*P12;
        ddos(0,1,iddo) *= pIv*P12;
        ddos(1,1,iddo) *= pIv*P22;
    }     
}

static double measure (const cx_mat& rho, const int n, const double kappa, const double I0, const double tau) {
    double P11 = poisson((1.0+kappa)*(1.0+kappa)*I0*tau,n);
    double P22 = poisson((1.0-kappa)*(1.0-kappa)*I0*tau,n);
    double pnt = P11*real(rho(0,0))+P22*real(rho(1,1));
    return pnt;
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

    deom d1(json["deom"]);

    const int inistate = json["rhot"]["inistate"].int_value();
    const int    nt = json["rhot"]["nt"].int_value();
    const int    nk = json["rhot"]["nk"].int_value();
    const double dt = json["rhot"]["dt"].number_value();
    const double staticErr = json["rhot"]["staticErr"].number_value();
    const double kappa = json["rhot"]["kappa"].number_value();
    const double tau = json["rhot"]["tau"].number_value();
    const double I0 = json["rhot"]["I0"].number_value();

    cx_cube rho1 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);
    cx_cube rho2 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);
    rho1(inistate,inistate,0) = 1.0;

    d1.equilibrium(rho1, dt, staticErr, 1, "SCI2");

    int nddo_backup(d1.nddo);
    hkey keys_backup(d1.keys);
    cx_cube ddos_backup(rho1);

    int it = 0;
    FILE *frho = fopen("equilibrium.log","w");
    do {
        double t = it*dt;
        printf("Propagate to equilibrium: step %d, nddo=%d, lddo=%d\n", it, d1.nddo, d1.lddo);
        if ((it+1)%nk == 0) {
            fprintf (frho, "%16.6e", t/deom_fs2unit);
            for (int i=0; i<d1.nsys; ++i) {
                fprintf(frho, "%16.6e", real(rho1(i,i,0)));
            }
            fprintf(frho, "\n");

            if (d1.notconverged(nddo_backup,keys_backup,ddos_backup,rho1,staticErr)) {
                nddo_backup = d1.nddo;
                keys_backup.rows(0,d1.nddo-1) = d1.keys.rows(0,d1.nddo-1);
                ddos_backup.head_slices(d1.nddo) = rho1.head_slices(d1.nddo);
            } else {
                break;
            }
        }
        d1.rk4 (rho1,t,dt,kappa,I0);
        it += 1;
    } while (it < ntMax);
    fclose(frho);

    if (it >= ntMax) {
        printf ("Fail to reach equilibrium at error %16.6e\n", staticErr);
    } else {
        printf ("Equilibrium reached with %d steps!\n", it);
    } 

    FILE *fp1 = fopen("pn1.dat","w");
    for (int n1=0; n1<=10; ++n1) {
        // collapse first time
        double pn1 = measure (rho1.slice(0), n1, kappa, I0, tau);
        fprintf (fp1, "%d %16.6e\n", n1, pn1);
        // 
        deom d2(d1);
        rho2.slices(0,d2.nddo-1) = rho1.slices(0,d2.nddo-1);
        collapse (d2, rho2, n1, kappa, I0, tau);
        // save to file
        double t = 0.0;
        stringstream ss;
        ss << "pn1_" << n1 << ".dat";
        FILE *fpn = fopen(ss.str().c_str(),"w");
        for (int it=0; it<nt; ++it) {
            for (int jt=0; jt<nk; ++jt) {
                d2.rk4 (rho2,t,dt,kappa,I0);
                t += dt;
            }
            printf ("Propagation %5.1f%%: n1=%6d, nddo=%6d, lddo=%3d\n",
                     100*it/static_cast<double>(nt), n1, d2.nddo, d2.lddo);
            // collapse first time
            fprintf (fpn, "%16.6e", t);
            for (int n2=0; n2<=10; ++n2) {
                double pn2 = measure (rho2.slice(0), n2, kappa, I0, tau);
                pn2 *= pn1;
                fprintf (fpn, "%16.6e", pn2);
            }
            fprintf (fpn, "\n");
        }
        fclose (fpn);
    }
    fclose(fp1);

    return 0;
}
