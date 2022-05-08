import numpy as np
import math as math

np.seterr('raise')


def MSjxMSjy(XS, YS, XB, YB, phi, S, Q):
    # Number of panels
    numPan = len(XB) - 1  # Number of panels

    # Initialize arrays
    MSjx = np.zeros(Q, numPan)
    MSjy = np.zeros(Q, numPan)

    # Compute integral
    for q in range(Q):
        for j in range(numPan):  # Loop over all panels
            # Compute intermediate values
            A = -(XS[q] - XB[j]) * np.cos(phi[j]) - (YS[q] - YB[j]) * np.sin(phi[j])  # A term
            B = (XS[q] - XB[j]) ** 2 + (YS - YB[j]) ** 2;  # B term
            Cx = -np.cos(phi[j]);  # C term (X-direction)
            Dx = XS[q] - XB[j];  # D term (X-direction)
            Cy = -np.sin(phi[j]);  # C term (Y-direction)
            Dy = YS[q] - YB[j];  # D term (Y-direction)
            E = math.sqrt(B - A ** 2);  # E term
            if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):
                MSjx[q, j] = 0
                MSjy[q, j] = 0
            else:

                term1 = 0.5 * Cx * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
                term2 = ((Dx - A * Cx) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
                MSjx[j] = term1 + term2;

                term1 = 0.5 * Cy * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
                term2 = ((Dy - A * Cy) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
                MSjy[j] = term1 + term2;

            # Zero out any problem values
            if (np.iscomplex(MSjx[q, j]) or np.isnan(MSjx[q, j]) or np.isinf(MSjx[q, j])):
                MSjx[q, j] = 0
            if (np.iscomplex(MSjy[q, j]) or np.isnan(MSjy[q, j]) or np.isinf(MSjy[q, j])):
                MSjy[q, j] = 0



    return MSjx, MSjy
