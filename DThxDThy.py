import math

import numpy as np


np.seterr('raise')

def DThxDThy(xT, yT, XW, YW, M):



    for h in range(M):
        rih = math.sqrt((xT-XW[h])**2 + (yT-YW[h])**2)
        DThx = (YW[h] - yT) / (2 * np.pi * rih**2)
        DThy = (yT - XW[h]) / (2 * np.pi * rih**2)



        return DThx, DThy

