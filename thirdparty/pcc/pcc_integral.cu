/*
 * Doing the 4 dimensional integration in Photon coincidence counting
 */
#include <cstdio>
#include <cstdlib>
#include <thrust/complex.h>
#include "json11.hpp"

using namespace std;
using namespace json11;

const complex<double> _ci = complex<double>(0.0,1.0);

__device__ inline double hstep (const double x) {
    if (x < 0) 
        return 0;
    else if (x > 0) 
        return 1.0;
    return 0.5;
}

__device__ inline complex<double> Dtw (const double tpr, const double tau, const double t, const double w, const double st, const double sw) {
    return 1.0/(2*sw)*hstep(tpr-t)*hstep(tpr+tau-t)*exp(-(_ci*w+st)*tau-2*st*(tpr-t))*(hstep(tau)*exp(-sw*tau)+hstep(-tau)*exp(sw*tau));
}

__device__ inline int at (const int nt, const int it0, const int it1, const int it2, const int it3) {
    return nt*(nt*(nt*it3+it2)+it1)+it0;
} 

__global__ void pccKernel (double* s2, const complex<double>* ft1, const complex<double>* ft2, const int nt, const int nw, const double tf, const double wi, const double wf, const double ws, const double t1, const double t2, const double st1, const double sw1, const double st2, const double sw2) {

    const int bx = blockIdx.x;
    const int tx = threadIdx.x;

    const double dt = tf/nt;
    const double dw = (wf-wi)/nw;
    const double dt4 = dt*dt*dt*dt; 
    const double w1 = wi+dw*bx;
    const double w2 = wi+dw*tx;
    complex<double> signal = 0;
    for (int it0=0; it0<nt; ++it0) { 
    for (int it1=0; it1<nt; ++it1) {
    for (int it2=0; it2<nt; ++it2) {
    for (int it3=0; it3<nt; ++it3) {
        const int loc = at(nt,it0,it1,it2,it3);
        {
            double tau2 = dt*it1;
            double t2pr = dt*it0;
            double t1pr = t2pr+tau2+dt*it2;
            double tau1 = dt*it3;
            singal += Dtw(t1pr,tau1,t1,w1,st1,sw1)*Dtw(t2pr,tau2,t2,w2,st2,sw2)*ft1[loc]*exp(_ci*ws*(tau1+tau2));
        }
        {
            double tau2 =-dt*it1;
            double t2pr = dt*it0-tau2;
            double t1pr = dt*it2+t2pr;
            double tau1 = dt*it3;
            signal += Dtw(t1pr,tau1,t1,w1,st1,sw1)*Dtw(t2pr,tau2,t2,w2,st2,sw2)*ft2[loc]*exp(_ci*ws*(tau1+tau2));
        }
    }}}}
    s2[bx*blockDim.x+tx] = real(signal)*dt4;
}


int main () {

    ifstream jsonFile("input_pcc.json");
    stringstream strStream;
    strStream << jsonFile.rdbuf();
    string jsonStr = strStream.str();
    string err;

    const Json json = Json::parse(jsonStr,err);
    if (!err.empty()) {
        printf ("Error in parsing input file: %s\n", err.c_str());
        return 0;
    } 
    
    const int    nt = json["nt"].int_value();
    const int    nw = json["nw"].int_value();
    const double tf = json["tf"].number_value();
    const double wi = json["wi"].number_value();
    const double wf = json["wf"].number_value();
    const double ws = json["ws"].number_value();
    const double t1 = json["t1"].number_value();
    const double t2 = json["t2"].number_value();
    const double st1 = json["st1"].number_value();
    const double sw1 = json["sw1"].number_value();
    const double st2 = json["st2"].number_value();
    const double sw2 = json["sw2"].number_value();
    
    const double dt = tf/nt;
    const double dw = (wf-wi)/nw;
    const double dt4 = dt*dt*dt*dt; 
    const int nt4 = nt*nt*nt*nt; 
    
    complex<double> *ft1, *ft2;
    double *s2;
    cudaMallocManaged(&ft1, nt4*sizeof(complex<double>));
    cudaMallocManaged(&ft2, nt4*sizeof(complex<double>));
    cudaMallocManaged( &s2, nt4*sizeof(double));

    double re(0), im(0);
    FILE *fs1 = fopen("barePcc_1.dat","r");
    FILE *fs2 = fopen("barePcc_2.dat","r");
    for (int i=0; i<nt4; ++i) {
        fscanf(fs1,"%lf",&re,&im);
        fscanf(fs2,"%lf",&re,&im);
    }

    pccKernel<<<nw,nw>>>(s2, ft1, ft2, nt, nw, tf, wi, wf, ws, t1, t2, st1, sw1, st2, sw2);

    FILE *fs = fopen("s2.mat","w");
    for (int iw1=0; iw1<nw; ++iw1) {
        for (int iw2=0; iw2<nw; ++iw2) {
            fprintf (fs, "%16.6e", s2[iw1*nw+iw2]);
        }
        fprintf (fs, "\n");
    }
    fclose(fs);
    
    cudaFree(s2);
    cudaFree(ft1);
    cudaFree(ft2);

    return 0;
}
