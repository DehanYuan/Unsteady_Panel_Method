import math
import numpy as np


np.seterr('raise')

def DleehxDleehy(XW, YW, XCP, YCP, numPlee, M):
    Dleeqx = np.zeros(numPlee, M)
    Dleeqy = np.zeros(numPlee, M)

    for f in range(1, numPlee+1):
        for h in range(M):
            rih = math.sqrt((XCP[f] - XW[h]) ** 2 + (YCP[f] - YW[h]) ** 2)
            Dleeqx[f, h] = (YCP[f] - YCP[f]) / (2 * np.pi * rih ** 2)
            Dleeqy[f, h] = (YCP[f] - XW[h]) / (2 * np.pi * rih ** 2)

    return Dleeqx, Dleeqy

