import math

import numpy as np


np.seterr('raise')

def Wihn(XC, YC, XS, YS, delta_i, Q):
    # Number of panels
    numPan = len(XC)  # Number of panels/control points


    Wiqn = np.zeros([numPan, Q])
    for i in range(numPan):
        for q in range(Q):
            rih = math.sqrt((XC[i]-XS[q])**2 + (YC[i]-YS[q])**2)
            Vihx = (YS[q] - YC[i]) / (2 * np.pi * rih**2)
            Vihy = (XC[i] - XS[q]) / (2 * np.pi * rih**2)
            Wihn[i, q] = Vihx * np.cos(delta_i) + Vihy * np.sin(delta_i)


        return Wiqn

