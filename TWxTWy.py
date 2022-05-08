
import numpy as np
import math as math

np.seterr('raise')


def TWxTWy(XW, YW, xT, yT, Theta_T, DeltaS_T, M):
    TWx = np.zeros(M)
    TWy = np.zeros(M)
    for h in range(M):
        # Compute intermediate values
        A = -(XW[h] - xT) * np.cos(Theta_T) - (YW[h] - yT) * np.sin(Theta_T)  # A term
        B = (XW[h] - xT) ** 2 + (YW[h] - yT) ** 2  # B term
        Cx = np.sin(Theta_T)  # Cx term (X-direction)
        Dx = -(YW[h] - yT)  # Dx term (X-direction)
        Cy = -np.cos(Theta_T)  # Cy term (Y-direction)
        Dy = XW[h] - xT  # Dy term (Y-direction)
        E = math.sqrt(B - A ** 2)  # E term
        if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):  # If E term is 0 or complex or a NAN or an INF
            TWx[h] = 0
            TWy[h] = 0
        else:

            term1 = 0.5 * Cx * np.log((DeltaS_T ** 2 + 2 * A * Theta_T + B) / B);
            term2 = ((Dx - A * Cx) / E) * (math.atan2((DeltaS_T + A), E) - math.atan2(A, E));
            TWx[h] = term1 + term2;

            term1 = 0.5 * Cy * np.log((DeltaS_T ** 2 + 2 * A * DeltaS_T + B) / B);
            term2 = ((Dy - A * Cy) / E) * (math.atan2((DeltaS_T + A), E) - math.atan2(A, E));
            TWy[h] = term1 + term2;

            # Zero out any problem values
        if (np.iscomplex(TWx[h]) or np.isnan(TWx[h]) or np.isinf(TWx[h])):
            TWx[h] = 0
        if (np.iscomplex(TWy[h]) or np.isnan(TWy[h]) or np.isinf(TWy[h])):
            TWy[h] = 0


    return TWx, TWy
