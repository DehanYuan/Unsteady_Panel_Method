import numpy as np


np.seterr('raise')

def VSqxVSqy(lamb, MSjx, MSjy, NSjx, NSjy, gamma, Q, DShx, DShy, DScx, DScy, Gamma_Sq,
               Gamma_wh, AoAR, Vinf, ASy, ASx, gamma_Sl, gamma_T, TSx, TSy):

    VSqx = np.zeros(Q)
    VSqy = np.zeros(Q)
    for q in range(Q):

        term1x = np.dot(MSjx[q], lamb) / (2 * np.pi)

        term2x = -sum(NSjx) * gamma / (2 * np.pi)

        term3x = Vinf * np.cos(AoAR)

        term4x = np.dot(DScx[q], Gamma_Sq)

        term5x = np.dot(DShx[q], Gamma_wh)

        term6x = gamma_Sl * sum(ASx) / (2 * np.pi)

        term7x = gamma_T * sum(TSx) / (2 * np.pi)

    VSqx[q] = term1x + term2x + term3x + term4x + term5x + term6x + term7x

    for q in range(Q):
        term1y = np.dot(MSjy[q], lamb) / (2 * np.pi)

        term2y = -sum(NSjy) * gamma / (2 * np.pi)

        term3y = Vinf * np.sin(AoAR)

        term4y = np.dot(DScy[q], Gamma_Sq)

        term5y = np.dot(DShy[q], Gamma_wh)

        term6y = gamma_Sl * sum(ASy) / (2 * np.pi)

        term7y = gamma_T * sum(TSy) / (2 * np.pi)

        VSqy[q] = term1y + term2y + term3y + term4y + term5y + term6y + term7y

    return VSqx, VSqy

