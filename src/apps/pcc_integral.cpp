/*
 * Doing the 4 dimensional integration in Photon coincidence counting
 */
#include <cstdio>
#include <cstdlib>
#include "armadillo"
#include "json11.hpp"

using namespace std;
using namespace arma;
using namespace json11;

const cx_double _ci = cx_double(0.0,1.0);

inline double hstep (const double x) {
    if (x < 0) 
        return 0;
    else if (x > 0) 
        return 1.0;
    return 0.5;
}

inline cx_double Dtw (const double tpr, const double tau, const double t, const double w, const double st, const double sw) {
    return 1.0/(2*sw)*hstep(tpr-t)*hstep(tpr+tau-t)*exp(-(_ci*w+st)*tau-2*st*(tpr-t))*(hstep(tau)*exp(-sw*tau)+hstep(-tau)*exp(sw*tau));
}

static int at (const int nt, const int it0, const int it1, const int it2, const int it3) {
    return nt*(nt*(nt*it3+it2)+it1)+it0;
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
    
    cx_vec ft1, ft2;
    ft1.load("barePcc_1.mat");
    ft2.load("barePcc_2.mat");

    mat s2(nw,nw);
    for (int iw1=0; iw1<nw; ++iw1) {
    for (int iw2=0; iw2<nw; ++iw2) {
        double w1 = wi+dw*iw1;
        double w2 = wi+dw*iw2;
        cx_double signal = 0;
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
                cx_double D1 = Dtw(t1pr,tau1,t1,w1,st1,sw1);
                cx_double D2 = Dtw(t2pr,tau2,t2,w2,st2,sw2);
                signal += D1*D2*ft1(loc)*exp(_ci*ws*(tau1+tau2));
            }
            {
                double tau2 =-dt*it1;
                double t2pr = dt*it0-tau2;
                double t1pr = dt*it2+t2pr;
                double tau1 = dt*it3;
                cx_double D1 = Dtw(t1pr,tau1,t1,w1,st1,sw1);
                cx_double D2 = Dtw(t2pr,tau2,t2,w2,st2,sw2);
                signal += D1*D2*ft2(loc)*exp(_ci*ws*(tau1+tau2));
            }
        }}}}
        printf ("iw1=%d, iw2=%d\n", iw1, iw2);
        s2(iw1,iw2) = real(signal)*dt4;
    }}

    FILE *fs = fopen("s2.mat","w");
    for (int iw1=0; iw1<nw; ++iw1) {
        for (int iw2=0; iw2<nw; ++iw2) {
            fprintf (fs, "%16.6e", s2(iw1,iw2));
        }
        fprintf (fs, "\n");
    }
    fclose(fs);
    
    return 0;
}
