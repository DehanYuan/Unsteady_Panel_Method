import math

import numpy as np


np.seterr('raise')


def DShxDShy(XS, YS, XW, YW, M, Q):
    DShx = np.zeros(Q, M)
    DShy = np.zeros(Q, M)
    for q in range(Q):
        for h in range(M):
            rih = math.sqrt((XS[q] - XW[h]) ** 2 + (YS[q] - YW[h]) ** 2)
            DShx[q, h] = (YW[h] - YS[q]) / (2 * np.pi * rih ** 2)
            DShy[q, h] = (YS[q] - XW[h]) / (2 * np.pi * rih ** 2)

    return DShx, DShy

