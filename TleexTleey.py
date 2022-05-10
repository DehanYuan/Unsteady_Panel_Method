
import numpy as np
import math as math

np.seterr('raise')


def TleexTleey(XCP, YCP, xTB, yTB, Theta_T, DeltaS_T, numPlee):
    Tleex = np.zeros(numPlee)
    Tleey = np.zeros(numPlee)
    for f in range(1, numPlee+1):
        # Compute intermediate values
        A = -(XCP[f] - xTB) * np.cos(Theta_T) - (YCP[f] - yTB) * np.sin(Theta_T)  # A term
        B = (XCP[f] - xTB) ** 2 + (YCP[f] - yTB) ** 2  # B term
        Cx = np.sin(Theta_T)  # Cx term (X-direction)
        Dx = -(YCP[f] - yTB)  # Dx term (X-direction)
        Cy = -np.cos(Theta_T)  # Cy term (Y-direction)
        Dy = XCP[f] - xTB  # Dy term (Y-direction)
        E = math.sqrt(B - A ** 2)  # E term

        if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):  # If E term is 0 or complex or a NAN or an INF
            Tleex[f] = 0
            Tleey[f] = 0
        else:

            term1 = 0.5 * Cx * np.log((DeltaS_T ** 2 + 2 * A * Theta_T + B) / B);
            term2 = ((Dx - A * Cx) / E) * (math.atan2((DeltaS_T + A), E) - math.atan2(A, E));
            Tleex[f] = term1 + term2;

            term1 = 0.5 * Cy * np.log((DeltaS_T ** 2 + 2 * A * DeltaS_T + B) / B);
            term2 = ((Dy - A * Cy) / E) * (math.atan2((DeltaS_T + A), E) - math.atan2(A, E));
            Tleey[f] = term1 + term2;

            # Zero out any problem values
        if (np.iscomplex(Tleex[f]) or np.isnan(Tleex[f]) or np.isinf(Tleex[f])):
            Tleex[f] = 0
        if (np.iscomplex(Tleey[f]) or np.isnan(Tleey[f]) or np.isinf(Tleey[f])):
            Tleey[f] = 0


    return Tleex, Tleey
