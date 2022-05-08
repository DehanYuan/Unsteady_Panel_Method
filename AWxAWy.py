
import numpy as np
import math as math

np.seterr('raise')


def AWxAWy(xSl, ySl, XW, YW, Theta_Sl, DeltaS_Sl, M):
    AWx = np.zeros(M)
    AWy = np.zeros(M)
    for h in range(M):
        # Compute intermediate values
        A = -(XW[h] - xSl) * np.cos(Theta_Sl) - (YW[h] - ySl) * np.sin(Theta_Sl)  # A term
        B = (XW[h] - xSl) ** 2 + (YW - ySl) ** 2  # B term
        Cx = np.sin(Theta_Sl)  # Cx term (X-direction)
        Dx = -(YW[h] - ySl)  # Dx term (X-direction)
        Cy = -np.cos(Theta_Sl)  # Cy term (Y-direction)
        Dy = XW[h] - xSl  # Dy term (Y-direction)
        E = math.sqrt(B - A ** 2)  # E term

        if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):  # If E term is 0 or complex or a NAN or an INF
            AWx[h] = 0
            AWy[h] = 0
        else:

            term1 = 0.5 * Cx * np.log((DeltaS_Sl ** 2 + 2 * A * Theta_Sl + B) / B);
            term2 = ((Dx - A * Cx) / E) * (math.atan2((DeltaS_Sl + A), E) - math.atan2(A, E));
            AWx[h] = term1 + term2;


            term1 = 0.5 * Cy * np.log((DeltaS_Sl ** 2 + 2 * A * DeltaS_Sl + B) / B);
            term2 = ((Dy - A * Cy) / E) * (
                        math.atan2((DeltaS_Sl + A), E) - math.atan2(A, E));
            AWy[h] = term1 + term2;

            # Zero out any problem values
        if (np.iscomplex(AWx[h]) or np.isnan(AWx[h]) or np.isinf(AWx[h])):
            AWx[h] = 0
        if (np.iscomplex(AWy[h]) or np.isnan(AWy[h]) or np.isinf(AWy[h])):
            AWy[h] = 0


    return AWx, AWy
