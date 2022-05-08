
import numpy as np
import math as math

np.seterr('raise')


def ASxASy(xSl, ySl, XS, YS, Theta_Sl, DeltaS_Sl, Q):
    ASx = np.zeros(Q)
    ASy = np.zeros(Q)
    for q in range(Q):
        # Compute intermediate values
        A = -(XS[q] - xSl) * np.cos(Theta_Sl) - (YS[q] - ySl) * np.sin(Theta_Sl)  # A term
        B = (XS[q] - xSl) ** 2 + (YS[q] - ySl) ** 2  # B term
        Cx = np.sin(Theta_Sl)  # Cx term (X-direction)
        Dx = -(YS[q] - ySl)  # Dx term (X-direction)
        Cy = -np.cos(Theta_Sl)  # Cy term (Y-direction)
        Dy = XS[q] - xSl  # Dy term (Y-direction)
        E = math.sqrt(B - A ** 2)  # E term

        if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):
            ASx[q] = 0
            ASy[q] = 0
        else:

            term1 = 0.5 * Cx * np.log((DeltaS_Sl ** 2 + 2 * A * Theta_Sl + B) / B);
            term2 = ((Dx - A * Cx) / E) * (math.atan2((DeltaS_Sl + A), E) - math.atan2(A, E));
            ASx[q] = term1 + term2;


            term1 = 0.5 * Cy * np.log((DeltaS_Sl ** 2 + 2 * A * DeltaS_Sl + B) / B);
            term2 = ((Dy - A * Cy) / E) * (
                        math.atan2((DeltaS_Sl + A), E) - math.atan2(A, E));
            ASy[q] = term1 + term2;

            # Zero out any problem values
        if (np.iscomplex(ASx[q]) or np.isnan(ASx[q]) or np.isinf(ASx[q])):
            ASx[q] = 0
        if (np.iscomplex(ASy[q]) or np.isnan(ASy[q]) or np.isinf(ASy[q])):
            ASy[q] = 0

    return ASx, ASy
