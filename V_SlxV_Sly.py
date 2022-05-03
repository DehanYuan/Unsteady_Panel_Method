import numpy as np


np.seterr('raise')

def V_SlxV_Sly(XC,  lamb, MSljx, MSljy, NSljx, NSljy, gamma_T, Q, M, DSlhx, DSlhy, DSlqx, DSlqy, Gamma_Sq,
               Gamma_wh, AoAR, Vinf, TSly, TSlx, gamma):
    # Number of panels

    numPan = len(XC)  # Number of panels/control points


    # Compute term1 of V_Slx
    for j in range(numPan):
        term1x = np.dot(MSljx, lamb) / (2 * np.pi)

    for j in range(numPan):
        term2x = -sum(NSljx) * gamma / (2 * np.pi)

    term3x = Vinf * np.cos(AoAR)

    for h in range(M):
        term4x = np.dot(DSlhx, Gamma_Sq)

    for q in range(Q):
        term5x = np.dot(DSlqx, Gamma_wh)

    term6x = gamma_T * sum(TSlx) / (2 * np.pi)

    V_Slx = term1x + term2x + term3x + term4x + term5x + term6x

    # Compute term1 of V_Sly
    for j in range(numPan):
        term1y = np.dot(MSljy, lamb) / (2 * np.pi)

    for j in range(numPan):
        term2y = -sum(NSljy) * gamma / (2 * np.pi)

    term3y = Vinf * np.sin(AoAR)

    for h in range(M):
        term4y = np.dot(DSlhy, Gamma_Sq)

    for q in range(Q):
        term5y = np.dot(DSlqy, Gamma_wh)

    term6y = gamma_T * sum(TSly) / (2 * np.pi)

    V_Sly = term1y + term2y + term3y + term4y + term5y + term6y

    return V_Slx, V_Sly

