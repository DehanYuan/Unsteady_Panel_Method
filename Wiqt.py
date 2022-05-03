import math

import numpy as np


np.seterr('raise')

def Wiqt(XC, YC, XS, YS, tau_x, tau_y, M):
    # Number of panels
    numPan = len(XC)  # Number of panels/control points


    Wiht = np.zeros([numPan, M])
    for i in range(numPan):
        for q in range(M):
            rih = math.sqrt((XC[i]-XS[q])**2 + (YC[i]-YS[q])**2)
            Vihx = (YS[q] - YC[i]) / (2 * np.pi * rih**2)
            Vihy = (XC[i] - XS[q]) / (2 * np.pi * rih**2)
            Wiqt[i,q] = Vihx * tau_x + Vihy * tau_y


        return Wiqt

