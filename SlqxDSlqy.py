import math

import numpy as np


np.seterr('raise')

def DSlqxDSlqy(xSl, ySl, XS, YS, M):


    for q in range(M):
        rih = math.sqrt((xSl-XS[q])**2 + (ySl-YS[q])**2)
        DSlqx = (YS[q] - ySl) / (2 * np.pi * rih**2)
        DSlqy = (ySl - XS[q]) / (2 * np.pi * rih**2)



        return DSlqx, DSlqy

