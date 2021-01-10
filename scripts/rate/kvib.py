import numpy as np
from math import sqrt, exp, pi
from unitconverter import kt2au, cm2au
import sys
import json

if __name__ == '__main__':

    with open(sys.argv[1]) as f:
        data = json.load(f)

    vij = data['vij']*cm2au
    tem = data['tem']*kt2au
    deG = data['deG']*cm2au
    lam = data['lam']*cm2au
    omg = data['omg']*cm2au
    shr = lam/omg
    fcf = 0.0
    pre = 1.0
    tmp = deG+lam
    for n in xrange(100):
        fcf += pre*exp(-tmp**2/(4.0*lam*tem))
        pre *= shr/(n+1)
        tmp += omg
    kij = vij**2/sqrt(lam*tem/pi)*exp(-shr)*fcf

    print 'rate wt vibration', kij

    kna = vij**2/sqrt(lam*tem/pi)*exp(-(deG+lam)**2/(4.0*lam*tem))
    print 'rate wo vibration', kna

        


