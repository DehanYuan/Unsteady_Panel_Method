import math

import numpy as np


np.seterr('raise')

def DTqxDTqy(xT, yT, XS, YS, M):


    for q in range(M):
        rih = math.sqrt((xT-XS[q])**2 + (yT-YS[q])**2)
        DTqx = (YS[q] - yT) / (2 * np.pi * rih**2)
        DTqy = (xT - XS[q]) / (2 * np.pi * rih**2)



        return DTqx, DTqy

