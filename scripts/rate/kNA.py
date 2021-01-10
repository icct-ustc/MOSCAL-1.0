#!/usr/bin/python
# -*- coding: utf-8 -*-
import json
import numpy as np
from cmath import exp as cexp
from math import sqrt, exp
from scipy import integrate
from BoseFermiExpansion import PSD
from unitconverter import kt2au, cm2au

def jwdru (omg,jdru):
    lamd, gamd = jdru['lamd'], jdru['gamd']
    return 2.0*lamd*gamd*omg/(omg**2+gamd**2)

def jwsdr (omg,jsdr):
    lams, omgs, gams = jsdr['lams'], jsdr['omgs'], jsdr['gams']
    return 2.0*lams*omgs**2*gams*omg/((omg**2-omgs**2)**2+(omg*gams)**2)

def fBose (x,pole,resi,rn,tn):
    return 1/x+0.5+rn*x+tn*x**3+sum(2.0*resi[i]*x/(x**2+pole[i]**2) for i in xrange(len(pole)))
    
def init (inidic):
        
    try:
        q_cl = inidic['q_cl']
    except:
        q_cl = 'q'
    npsd = inidic['npsd']
    pade = inidic['pade']
    temp = inidic['temp']*kt2au
    jdru = [{'lamd':x['lamd']*cm2au,'gamd':x['gamd']*cm2au} for x in inidic['jdru']]
    jsdr = [{'lams':x['lams']*cm2au,'omgs':x['omgs']*cm2au,'gams':x['gams']*cm2au} for x in inidic['jsdr']]

    ndru = len(jdru)
    nsdr = len(jsdr)
    nper = ndru+2*nsdr+npsd

    lamT = 0.0
    delr = 0.0
    expn = np.zeros(nper,dtype=complex)
    etal = np.zeros(nper,dtype=complex)

    pole, resi, rn, tn = PSD (npsd,BoseFermi=1,pade=pade)
    
    for idru in xrange(ndru):
        lamd, gamd = jdru[idru]['lamd'], jdru[idru]['gamd']
        iper = idru
        expn[iper] = gamd
        etal[iper] = -2.J*lamd*gamd*fBose(-1.J*gamd/temp,pole,resi,rn,tn)
        if q_cl == 'cl':
            etal[iper] = etal[iper].real
        lamT += lamd
        delr += 2.*lamd*gamd/temp*rn
    
    for isdr in xrange(nsdr):
        lams, omgs, gams = jsdr[isdr]['lams'], jsdr[isdr]['omgs'], jsdr[isdr]['gams']
        iper, jper = ndru+isdr*2, ndru+isdr*2+1
        etaBO = 2.*lams*omgs*omgs*gams
        Delta = omgs*omgs-gams*gams/4.0
        lamT += lams
        delr +=  etaBO*tn/(temp**3)
        if Delta > 0:
            OmgB = sqrt(Delta)
            expn[iper] = 0.5*gams+1.J*OmgB
            expn[jper] = 0.5*gams-1.J*OmgB
        elif Delta < 0:
            OmgB = sqrt(-Delta)
            expn[iper] = 0.5*gams+OmgB
            expn[jper] = 0.5*gams-OmgB           
        else:
            raise ValueError("Not prepared for Delta=0")
        z1, z2 = -1.J*expn[iper], -1.J*expn[jper]
        etal[iper] = -2.J*etaBO*z1/(2.*z1*(z1+z2)*(z1-z2))*fBose(z1/temp,pole,resi,rn,tn)
        etal[jper] = -2.J*etaBO*z2/(2.*z2*(z2+z1)*(z2-z1))*fBose(z2/temp,pole,resi,rn,tn)
        if Delta > 0:
            if q_cl == 'cl':
                etal[iper],etal[jper] = 0.5*(etal[iper]+etal[jper].conj()), \
                        0.5*(etal[jper]+etal[iper].conj())
        elif Delta < 0:
            if q_cl == 'cl':
                etal[iper],etal[jper] = etal[iper].real,etal[jper].real
        else:
            raise ValueError("Not prepared for Delta=0")

    for ipsd in xrange(npsd):
    
        iper = ndru+nsdr*2+ipsd
        zomg = -1.J*pole[ipsd]*temp
        jsum = sum(jwdru(zomg,x) for x in jdru)+sum(jwsdr(zomg,x) for x in jsdr)
        expn[iper] = pole[ipsd]*temp            
        etal[iper] = -2.J*resi[ipsd]*temp*jsum   

    return lamT, delr, expn, etal
                            
def integrant (t,delG,lamT,delr,expn,etal):

    npol = len(expn)
    gt = (1.J*(delG+lamT)+delr)*t  \
           + sum(etal[k]/expn[k] for k in xrange(npol))*t \
           - sum(etal[k]/expn[k]**2*(1-cexp(-expn[k]*t)) for k in xrange(npol))
    result = cexp(-gt).real
    return result

if __name__ == '__main__':

    inidic = {
        "q_cl": 'q',
        "npsd": 1,
        "pade": 2,
        "temp": 1.0,
        "jdru": [{"lamd":3.0,"gamd":1.0}],
        "jsdr": [{"lams":1.0,"omgs":5.0,"gams":1.0}],
        "delG": -5.0,
        "vab" : 1.0,
    }

    with open('rateNA.dat','w') as f:
        for delG in np.linspace(-12.0,2.0,100):

            inidic["delG"] = delG

            delG = inidic["delG"]*cm2au
            vab = inidic["vab"]*cm2au

            lamT, delr, expn, etal = init(inidic)

            exp_gt = lambda t, delG, lamT, delr, expn, etal: integrant(t,delG,lamT,delr,expn,etal)
            y, err = integrate.quad(exp_gt, 0, np.inf, args=(delG,lamT,delr,expn,etal), limit=500)

            k = 2*vab**2*y

            print >> f, 2*'%16.6e'%(delG,k)
