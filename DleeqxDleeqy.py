import math
import numpy as np


np.seterr('raise')

def DleeqxDleeqy(XS, YS, XCP, YCP, numPlee, Q):
    Dleeqx = np.zeros(numPlee, Q)
    Dleeqy = np.zeros(numPlee, Q)

    for f in range(1, numPlee+1):
        for q in range(Q):
            rih = math.sqrt((XCP[f] - XS[q]) ** 2 + (YCP[f] - YS[q]) ** 2)
            Dleeqx[f, q] = (YCP[f] - YCP[f]) / (2 * np.pi * rih ** 2)
            Dleeqy[f, q] = (YCP[f] - XS[q]) / (2 * np.pi * rih ** 2)

    return Dleeqx, Dleeqy

