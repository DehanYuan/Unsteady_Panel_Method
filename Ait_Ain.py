import numpy as np
import math as math

np.seterr('raise')


def Ait_Ain(XC, YC, x_isl, y_isl, phi, DeltaS_Sl, Theta_Sl):
    # Number of panels
    numPan = len(XC)  # Number of panels

    # Initialize arrays
    Ain = np.zeros(numPan)
    Ait = np.zeros(numPan)

    # Compute integral
    for i in range(numPan):
        A = -(XC[i] - x_isl) * np.cos(Theta_Sl) - (YC[i] - y_isl) * np.sin(Theta_S)  
        B = (XC[i] - x_isl) ** 2 + (YC[i] - y_isl) ** 2  
        Cn = -np.cos(phi[i] - Theta_Sl)  
        Dn = (XC[i] - x_isl) * np.cos(phi[i]) + (YC[i] - y_isl) * np.sin(phi[i])  
        Ct = np.sin(Theta_Sl - phi[i])  
        Dt = (XC[i] - x_isl) * np.sin(phi[i]) - (YC[i] - y_isl) * np.cos(phi[i])  
        E = np.sqrt(B - A ** 2)  
        if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(
                E)):  
            Ain[i] = 0  
            Ait[i] = 0  
        else:

            term1 = 0.5 * Cn * np.log((DeltaS_Sl ** 2 + 2 * A * DeltaS_Sl + B) / B)  
            term2 = ((Dn - A * Cn) / E) * (
                    math.atan2((DeltaS_Sl + A), E) - math.atan2(A, E))  
            Ain[i] = term1 + term2  
            
            
            term1 = 0.5 * Ct * np.log((DeltaS_Sl ** 2 + 2 * A * DeltaS_Sl + B) / B)  
            term2 = ((Dt - A * Ct) / E) * (
                    math.atan2((DeltaS_Sl + A), E) - math.atan2(A, E))  
            Ait[i] = term1 + term2  

            # Zero out any problem values
            if (np.iscomplex(Ain[i]) or np.isnan(Ain[i]) or np.isinf(
                    Ain[i])):  
                Ain[i] = 0  
            if (np.iscomplex(Ait[i]) or np.isnan(Ait[i]) or np.isinf(
                    Ait[i])):  
                Ait[i] = 0  
                
    return Ain, Ait  
