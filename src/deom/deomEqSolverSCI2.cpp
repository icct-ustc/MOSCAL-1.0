/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cstdio>
#include <cstdlib>
#include "deom.hpp"

static const int ntMax = 100000000;

static double norm (const int nddo, const cx_cube& x) {
    const int n1 = nddo, n2 = x.n_rows, n3 = x.n_cols;
    double nrm = 0.0;
    for (int k=0; k<n1; ++k) {
    for (int i=0; i<n2; ++i) {
    for (int j=0; j<n3; ++j) {
        double tmp = abs(x(i,j,k));
        nrm+= tmp*tmp;
    } } }
    return sqrt(nrm);
}

void deom::EqSolverSCI2 (cx_cube& ddos, const double dt, const double err, const int nk) {

    if (nddo != 1) {
        printf("Warning: Hidx is not correctly initialized!\n");
    }
    vec eigh = eig_sym(ham1);
    //const double OMG = max(abs(eigh))*100.1;
    const double OMG = max(0.5*max(abs(eigh)),5*sqrt(max(coef_abs)));
    printf ("OMG=%16.6e\n",OMG);

    if (ddos.is_empty()) {
        ddos = zeros<cx_cube>(nsys,nsys,nmax);
        mat eham = expmat_sym(-real(ham1)/temperature);
        ddos.slice(0) = eham/trace(eham)*deom_c1;
    } 
    // else {
    //    ddos.zeros();
    // }
    // mat eham = expmat_sym(-real(ham1)/temperature);
    // ddos.slice(0) = eham/trace(eham)*deom_c1;
    // ddos(0,0,0) = deom_c1;

    // propagation (ddos, dt, 200, 1);

    cx_cube ddos_diff = zeros<cx_cube>(size(ddos));
    cx_mat rho_old(nsys,nsys);
    cx_mat htmp(deom_ci*ham1);
    cx_mat hlft(htmp);
    cx_mat hrht(htmp);
    double max_diff = 0.0;
    FILE *ferr = fopen("sci-error.dat","w");
    int iter = 0;
    do { 
        max_diff = 0.0;
        for (int iddo=0; iddo<nddo; ++iddo) {

            hnod& nod = keys(iddo);
            ivec key0(nod.key);
            const int tier = tree.find(key0)->tier;
            const bool flag = is_valid(ddos.slice(iddo));

            cx_mat rhob = zeros<cx_mat>(nsys,nsys);
            for (int m=0; m<nmod; ++m) {
                rho_old = qmd1.slice(m)*ddos.slice(iddo)-ddos.slice(iddo)*qmd1.slice(m);
                if (abs(delt_res(m))>1.0e-15) {
                    rhob -= delt_res(m)*(qmd1.slice(m)*rho_old-rho_old*qmd1.slice(m));
                }
            }
            if (tier < lmax) {
                for (int mp=0; mp<nind; ++mp) {
                    ivec key1 = gen_key(key0, mp, 1);
                    const int m = modLabel(mp);
                    const int n = key1(mp)-1;
                    const cx_double sn = -deom_ci*sqrt((n+1)*coef_abs(mp));
                    auto *ptr = tree.find(key1);
                    if (ptr && (ptr->rank)>=0) {
                        int loc = ptr->rank;
                        rhob += sn*(qmd1.slice(m)*ddos.slice(loc)-ddos.slice(loc)*qmd1.slice(m));
                    } else if (flag) {
                        tree.try_insert(key1,nddo);
                        keys(nddo) = hnod(nod.gams+expn(mp),key1);
                        ddos.slice(nddo).zeros();
                        nddo += 1;
                    }
                }
            }
            for (int mp=0; mp<nind; ++mp) {
                ivec key1 = gen_key(key0, mp, -1);
                if (!key1.is_empty()) {
                    int m = modLabel(mp);
                    if (flag_rwa != 0) { // RWA
                        m = (modLabel(mp)%2==0)?(modLabel(mp)+1):(modLabel(mp)-1);
                    }
                    const int n = (key1.n_rows<(unsigned int)(mp+1))?1:(key1(mp)+1);
                    const cx_double sn = -deom_ci*sqrt(n/coef_abs(mp));
                    const cx_double cl = sn*coef_lft(mp);
                    const cx_double cr = sn*coef_rht(mp);
                    auto *ptr = tree.find(key1);
                    if (ptr && (ptr->rank)>=0) {
                        int loc = ptr->rank;
                        rhob += cl*qmd1.slice(m)*ddos.slice(loc)-cr*ddos.slice(loc)*qmd1.slice(m); 
                    } else if (flag) {
                        tree.try_insert(key1,nddo);
                        keys(nddo) = hnod(nod.gams-expn(mp),key1);
                        ddos.slice(nddo).zeros();
                        nddo += 1;
                    }
                }
            }

            rhob += OMG*ddos.slice(iddo);

            rho_old = ddos.slice(iddo);
            hlft = htmp;
            hlft.diag() += OMG+nod.gams;
            hrht =-htmp;
            if (iddo != 0) {
                ddos.slice(iddo) = syl(hlft,hrht,-rhob);
            } else {
                cx_mat rhox = syl(hlft,hrht,-rhob);
                ddos.slice(0) = 0.5*(rhox+rhox.t());
            }
            double diff = max(vectorise(abs(rho_old-ddos.slice(iddo))));
            max_diff = max_diff>diff?max_diff:diff;
        }
        
        iter += 1;

        filter(ddos);

        rem_freeze(ddos_diff,ddos,0.0);
        double error = norm(nddo, ddos_diff);

        printf("SCI step %d, nddo=%d, lddo=%d, max_diff=%16.6e\n", iter, nddo, lddo, max_diff);
        fprintf(ferr, "%d %16.6e %16.6e\n", iter, error, max_diff);
        vec pop = real(ddos.slice(0).diag());
        pop.print("population");

        if (max_diff < err) break;

    } while (iter < ntMax);
    fclose(ferr);

    //int nddo_backup(nddo);
    //hkey keys_backup(keys);
    //cx_cube ddos_backup(ddos);
    //int it = 0;
    //FILE *frho = fopen("equilibrium.log","w");
    //do {
    //    double t = it*dt;
    //    printf("Propagate to equilibrium: step %d, nddo=%d, lddo=%d\n", it, nddo, lddo);
    //    if ((it+1)%nk == 0) {
    //        fprintf (frho, "%16.6e", t/deom_fs2unit);
    //        for (int i=0; i<nsys; ++i) {
    //            fprintf(frho, "%16.6e", real(ddos(i,i,0)));
    //        }
    //        fprintf(frho, "\n");

    //        if (notconverged(nddo_backup,keys_backup,ddos_backup,ddos,err)) {
    //            nddo_backup = nddo;
    //            keys_backup.rows(0,nddo-1) = keys.rows(0,nddo-1);
    //            ddos_backup.head_slices(nddo) = ddos.head_slices(nddo);
    //        } else {
    //            break;
    //        }
    //    }
    //    rk4 (ddos,t,dt);
    //    it += 1;
    //} while (it < ntMax);
    //fclose(frho);
    //if (it >= ntMax) {
    //    printf ("Fail to reach equilibrium at error %16.6e\n", err);
    //} else {
    //    printf ("Equilibrium reached with %d steps!\n", it);
    //} 
    ddos.slice(0).save("rho0.real",raw_ascii);
    //// save rhon                                                                        
    FILE *fs = fopen("rhonreal.dat","w");
    for (int i=0; i<nddo; ++i){
        for (size_t j=0; j<keys(i).key.n_rows; ++j)
            fprintf(fs,"%d\t",(int)keys(i).key(j));
        fprintf (fs, "\n");
        for (int j=0; j<nsys; ++j)
            for (int k=0; k<nsys; ++k)
                fprintf (fs, "(%26.16e,%26.16e)\t",ddos(j,k,i).real(),ddos(j,k,i).imag());
        fprintf (fs, "\n\n");
    }
    fclose(fs);

}
