import numpy as np


np.seterr('raise')

def V_TxV_Ty(XC, lamb, MTjx, MTjy, NTjx, NTjy, gamma, Q, M, DThx, DThy, DTqx, DTqy, Gamma_Sq,
               Gamma_wh, AoAR, Vinf, ATy, ATx):
    # Number of panels

    numPan = len(XC)  # Number of panels/control points


    # Compute term1 of V_Slx
    for j in range(numPan):
        term1x = np.dot(MTjx, lamb) / (2 * np.pi)

    for j in range(numPan):
        term2x = -sum(NTjx) * gamma / (2 * np.pi)

    term3x = Vinf * np.cos(AoAR)

    for h in range(M):
        term4x = np.dot(DThx, Gamma_Sq)

    for q in range(Q):
        term5x = np.dot(DTqx, Gamma_wh)

    term6x = gamma * sum(ATx) / (2 * np.pi)

    V_Tx = term1x + term2x + term3x + term4x + term5x + term6x

    # Compute term1 of V_Sly
    for j in range(numPan):
        term1y = np.dot(MTjy, lamb) / (2 * np.pi)

    for j in range(numPan):
        term2y = -sum(NTjy) * gamma / (2 * np.pi)

    term3y = Vinf * np.sin(AoAR)

    for h in range(M):
        term4y = np.dot(DThy, Gamma_Sq)

    for q in range(Q):
        term5y = np.dot(DTqy, Gamma_wh)

    term6y = gamma * sum(ATy) / (2 * np.pi)

    V_Ty = term1y + term2y + term3y + term4y + term5y + term6y

    return V_Tx, V_Ty

