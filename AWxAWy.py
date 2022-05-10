
import numpy as np
import math as math

np.seterr('raise')


def AWxAWy(x_isl, y_isl, XW, YW, Theta_Sl, DeltaS_Sl, M):
    AWx = np.zeros(M)
    AWy = np.zeros(M)
    for h in range(M):
        # Compute intermediate values
        A = -(XW[h] - x_isl) * np.cos(Theta_Sl) - (YW[h] - y_isl) * np.sin(Theta_Sl)  # A term
        B = (XW[h] - x_isl) ** 2 + (YW - y_isl) ** 2  # B term
        Cx = np.sin(Theta_Sl)  # Cx term (X-direction)
        Dx = -(YW[h] - y_isl)  # Dx term (X-direction)
        Cy = -np.cos(Theta_Sl)  # Cy term (Y-direction)
        Dy = XW[h] - x_isl  # Dy term (Y-direction)
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
