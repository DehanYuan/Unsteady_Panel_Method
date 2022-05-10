
import numpy as np
import math as math

np.seterr('raise')


def AleexAleey(x_isl, y_isl, XCP, YCP, Theta_Sl, DeltaS_Sl, numPlee):
    Aleex = np.zeros(numPlee)
    Aleey = np.zeros(numPlee)
    for f in range(1, numPlee+1):
        A = -(XCP[f] - x_isl) * np.cos(Theta_Sl) - (YCP[f] - y_isl) * np.sin(Theta_Sl)  # A term
        B = (XCP[f] - x_isl) ** 2 + (YCP[f] - y_isl) ** 2  # B term
        Cx = np.sin(Theta_Sl)  # Cx term (X-direction)
        Dx = -(YCP[f] - y_isl)  # Dx term (X-direction)
        Cy = -np.cos(Theta_Sl)  # Cy term (Y-direction)
        Dy = XCP[f] - x_isl  # Dy term (Y-direction)
        E = math.sqrt(B - A ** 2)  # E term

        if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):
            Aleex[f] = 0
            Aleey[f] = 0
        else:

            term1 = 0.5 * Cx * np.log((DeltaS_Sl ** 2 + 2 * A * Theta_Sl + B) / B);
            term2 = ((Dx - A * Cx) / E) * (math.atan2((DeltaS_Sl + A), E) - math.atan2(A, E));
            Aleex[f] = term1 + term2;


            term1 = 0.5 * Cy * np.log((DeltaS_Sl ** 2 + 2 * A * DeltaS_Sl + B) / B);
            term2 = ((Dy - A * Cy) / E) * (
                        math.atan2((DeltaS_Sl + A), E) - math.atan2(A, E));
            Aleey[f] = term1 + term2;

            # Zero out any problem values
        if (np.iscomplex(Aleex[f]) or np.isnan(Aleex[f]) or np.isinf(Aleex[f])):
            Aleex[f] = 0
        if (np.iscomplex(Aleey[f]) or np.isnan(Aleey[f]) or np.isinf(Aleey[f])):
            Aleey[f] = 0

    return Aleex, Aleey
