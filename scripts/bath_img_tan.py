#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from math import sqrt, exp, pi
from BoseFermiExpansion import PSD
import armadillo as arma 
from cmath import exp as cexp
import json

def jwdru (omg,jdru):
    lamd, gamd = jdru[0], jdru[1]
    return 2.0*lamd*gamd*omg/(omg**2+gamd**2)

def jwsdr (omg,jsdr):
    lams, omgs, gams = jsdr[0], jsdr[1], jsdr[2]
    return 2.0*lams*omgs**2*gams*omg/((omg**2-omgs**2)**2+(omg*gams)**2)

def fBose (x,pole,resi,rn,tn,sign=0):
    if sign == 0:
        return 1/x+0.5+rn*x+tn*x**3+sum(2.0*resi[i]*x/(x**2+pole[i]**2) for i in xrange(len(pole)))
    else:
        if isinstance(x,complex):
            return 1.0/(1-cexp(-x))
        else:
            return 1.0/(1-exp(-x))
    
def init (inidic):
        
    try:
        exbe = inidic['exbe']
    except:
        exbe = 0

    pade = 0
    npsd = inidic['npsd']
    temp = inidic['temp']
    jomg = inidic['jomg']
    nmod = len(jomg)

    nind = 0
    mode = []
    for m in xrange(nmod):
        try:
            ndru = len(jomg[m]['jdru'])
        except:
            ndru = 0
        nper = ndru+npsd*2
        nind += nper
        mode.extend([m for i in xrange(nper)])
    mode = np.array(mode)
    delr = np.zeros(nmod,dtype=float)
    expn = np.zeros(nind,dtype=complex)
    etal = np.zeros(nind,dtype=complex)
    etar = np.zeros(nind,dtype=complex)
    etaa = np.zeros(nind,dtype=float)

    iind = 0

    for m in xrange(nmod):

        try:
            jdru = jomg[m]['jdru']
        except:
            jdru = []

        ndru = len(jdru)
    
        for idru in xrange(ndru):

            if len(jdru[idru]) != 2:
                raise ValueError('Invalid drude input')

            lamd, gamd = jdru[idru][0], jdru[idru][1]
            expn[iind] = 0.0
            etal[iind] = 2*lamd*temp
            etar[iind] = 0.0
            etaa[iind] = abs(etal[iind])
            iind += 1

        for ipsd in xrange(npsd):
            expn[iind] = 2.0*pi*temp*(ipsd+1)
            expn[iind+1] = -2.0*pi*temp*(ipsd+1)
            etal[iind] = 2*lamd*gamd*temp/(gamd+expn[iind])
            etal[iind+1] = 2*lamd*gamd*temp/(gamd+expn[iind])
            etar[iind] = 0.0
            etar[iind+1] = 0.0
            etaa[iind] = abs(etal[iind])  
            etaa[iind+1] = abs(etal[iind+1])  
            iind += 2

    arma.save (mode,inidic['modeFile'])
    arma.save (etal,inidic['etalFile'])
    arma.save (etar,inidic['etarFile'])
    arma.save (etaa,inidic['etaaFile'])
    arma.save (expn,inidic['expnFile'])
    arma.save (delr,inidic['delrFile'])
            
    return etal, expn
            
if __name__ == '__main__':

    inidic = {
        "q_cl": 'q',
        "exbe": 0,
        "npsd": 2,
        "pade": 2,
        "temp": 1.0,
        "jomg": [{"jdru":[(0.5,0.5)],"jsdr":[(0.1,1.0,0.1)]},{"jdru":[(0.1,0.1)],"jsdr":[(0.1,1.0,0.1)]}],
	"modeFile": "inp_mode.mat",
	"etalFile": "inp_etal.mat",
	"etarFile": "inp_etar.mat",
	"etaaFile": "inp_etaa.mat",
	"expnFile": "inp_expn.mat",
	"delrFile": "inp_delr.mat"
    }
    init(inidic)
    with open('input.json','w') as f:
        json.dump(inidic,f,indent=4) 
    
