import numpy as np

np.seterr('raise')


def VSqxVSqy(lamb, MWjx, MWjy, NWjx, NWjy, gamma, M, DWqx, DWqy, DWmx, DWmy, Gamma_Sq,
             Gamma_wh, AoAR, Vinf, AWy, AWx, gamma_Sl, gamma_T, TWx, TWy):
    term1x = np.zeros(M)
    term2x = np.zeros(M)
    term3x = np.zeros(M)
    term4x = np.zeros(M)
    term5x = np.zeros(M)
    term6x = np.zeros(M)
    term7x = np.zeros(M)
    term1y = np.zeros(M)
    term2y = np.zeros(M)
    term3y = np.zeros(M)
    term4y = np.zeros(M)
    term5y = np.zeros(M)
    term6y = np.zeros(M)
    term7y = np.zeros(M)
    VSqx = np.zeros(M)
    VSqy = np.zeros(M)
    for h in range(M):
        term1x[h] = np.dot(MWjx, lamb) / (2 * np.pi)

        term2x[h] = -sum(NWjx) * gamma / (2 * np.pi)

        term3x[h] = Vinf * np.cos(AoAR)

        term4x[h] = np.dot(DWmx, Gamma_wh)

        term5x[h] = np.dot(DWqx, Gamma_Sq)

        term6x[h] = gamma_Sl * sum(AWx) / (2 * np.pi)

        term7x[h] = gamma_T * TWx

    VSqx[h] = term1x[h] + term2x[h] + term3x[h] + term4x[h] + term5x[h] + term6x[h] + term7x[h]

    for h in range(M):
        term1y[h] = np.dot(MWjy, lamb) / (2 * np.pi)

        term2y[h] = -sum(NWjy) * gamma / (2 * np.pi)

        term3y[h] = Vinf * np.sin(AoAR)

        term4y[h] = np.dot(DWmy, Gamma_wh)

        term5y[h] = np.dot(DWqy, Gamma_Sq)

        term6y[h] = gamma_Sl * sum(AWy) / (2 * np.pi)

        term7y[h] = gamma_T * TWy

        VSqy[h] = term1y[h] + term2y[h] + term3y[h] + term4y[h] + term5y[h] + term6y[h] + term7y[h]

    return VSqx, VSqy

