import numpy as np


np.seterr('raise')

def Vleex(lamb, Mleejx, Nleejx, gamma, numPlee, Dleehx, Gamma_Sq, Dleeqx,
               Gamma_wh, Vinf, Aleex, gamma_Sl, gamma_T, Tleex):

    Vleex = np.zeros(numPlee)
    for f in range(1, numPlee+1):

        term1x = np.dot(Mleejx[f], lamb) / (2 * np.pi)

        term2x = -sum(Nleejx[f]) * gamma / (2 * np.pi)

        term3x = Vinf

        term4x = np.dot(Dleeqx[f], Gamma_Sq)

        term5x = np.dot(Dleehx[f], Gamma_wh)

        term6x = gamma_Sl * sum(Aleex) / (2 * np.pi)

        term7x = gamma_T * sum(Tleex) / (2 * np.pi)

    Vleex[f] = term1x + term2x + term3x + term4x + term5x + term6x + term7x


    return Vleex

