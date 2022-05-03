import math

import numpy as np


np.seterr('raise')

def DSlhxDSlhy(xSl, ySl, XW, YW, M):



    for h in range(M):
        rih = math.sqrt((xSl-XW[h])**2 + (ySl-YW[h])**2)
        DSlhx = (YW[h] - ySl) / (2 * np.pi * rih**2)
        DSlhy = (ySl - XW[h]) / (2 * np.pi * rih**2)



        return DSlhx, DSlhy

