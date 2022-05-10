
import numpy as np
import math as math

np.seterr('raise')


def TSxTSy(XS, YS, xTB, yTB, Theta_T, DeltaS_T, Q):
    TSx = np.zeros(Q)
    TSy = np.zeros(Q)
    for q in range(Q):
        # Compute intermediate values
        A = -(XS[q] - xTB) * np.cos(Theta_T) - (YS[q] - yTB) * np.sin(Theta_T)  # A term
        B = (XS[q] - xTB) ** 2 + (YS[q] - yTB) ** 2  # B term
        Cx = np.sin(Theta_T)  # Cx term (X-direction)
        Dx = -(YS[q] - yTB)  # Dx term (X-direction)
        Cy = -np.cos(Theta_T)  # Cy term (Y-direction)
        Dy = XS[q] - xTB  # Dy term (Y-direction)
        E = math.sqrt(B - A ** 2)  # E term

        if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):  # If E term is 0 or complex or a NAN or an INF
            TSx[q] = 0
            TSy[q] = 0
        else:

            term1 = 0.5 * Cx * np.log((DeltaS_T ** 2 + 2 * A * Theta_T + B) / B);
            term2 = ((Dx - A * Cx) / E) * (math.atan2((DeltaS_T + A), E) - math.atan2(A, E));
            TSx[q] = term1 + term2;

            term1 = 0.5 * Cy * np.log((DeltaS_T ** 2 + 2 * A * DeltaS_T + B) / B);
            term2 = ((Dy - A * Cy) / E) * (math.atan2((DeltaS_T + A), E) - math.atan2(A, E));
            TSy[q] = term1 + term2;

            # Zero out any problem values
        if (np.iscomplex(TSx[q]) or np.isnan(TSx[q]) or np.isinf(TSx[q])):
            TSx[q] = 0
        if (np.iscomplex(TSy[q]) or np.isnan(TSy[q]) or np.isinf(TSy[q])):
            TSy[q] = 0


    return TSx, TSy
