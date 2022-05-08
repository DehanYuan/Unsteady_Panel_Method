import numpy as np

def V_it(XC, Vs_x, tau_x, Vs_y, tau_y, J, lamb, gamma, L, Tit, DCM, gamma_T, Wiht, M, Ait, Gamma_wh, gamma_Sl, Q, Diqt,
         Gamma_Sq):

    numPan = len(XC)
    V_it = np.zeros(numPan)
    for i in range(numPan):
        ACM = 2 * np.pi * (Vs_x[i] * tau_x[i] + Vs_y[i] * tau_y[i])  # Compute A term of tangent velocity to calculate KC
        BCM = np.dot(lamb, J[i]) / (2 * np.pi)   # Compute B term of tangent velocity to calculate KC
        CCM = - sum(L[i] * gamma) / (2 * np.pi)  # Compute C term of tangent velocity to calculate KC
        DCM = np.dot(Wiht, Gamma_wh)  # Compute D term of tangent velocity to calculate KC
        ECM = Tit * gamma_T  # Compute E term of tangent velocity to calculate KC
        FCM = Ait * gamma_Sl  # Compute F term of tangent velocity to calculate KC
        GCM = np.dot(Diqt, Gamma_Sq)  # Compute G term of tangent velocity to calculate KC
        HCM = gamma / 2
        V_it[i] = ACM + BCM + CCM + DCM + ECM + FCM + GCM + HCM  # Compute the tangent velocity of airfoil panel

    return V_it
