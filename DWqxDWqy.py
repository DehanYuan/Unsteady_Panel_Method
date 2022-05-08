import math
import numpy as np


np.seterr('raise')

def DWqxDWqy(XS, YS, XW, YW, M, Q):
    DWqx = np.zeros(M, Q)
    DWqy = np.zeros(M, Q)

    for h in range(M):
        for q in range(Q):
            rih = math.sqrt((XW[h] - XS[q]) ** 2 + (YW[h] - YS[q]) ** 2)
            DWqx[h, q] = (YS[q] - YW[h]) / (2 * np.pi * rih ** 2)
            DWqy[h, q] = (YW[h] - XS[q]) / (2 * np.pi * rih ** 2)

    return DWqx, DWqy

