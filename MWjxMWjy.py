import numpy as np
import math as math

np.seterr('raise')


def MWjxMWjy(XW, YW, XB, YB, phi, S, M):
    # Number of panels
    numPan = len(XB) - 1  # Number of panels

    # Initialize arrays
    MWjx = np.zeros([M, numPan])
    MWjy = np.zeros([M, numPan])

    # Compute integral
    for h in range(M):
        for j in range(numPan):  # Loop over all panels
            # Compute intermediate values
            A = -(XW[h] - XB[j]) * np.cos(phi[j]) - (YW[h] - YB[j]) * np.sin(phi[j])  # A term
            B = (XW[h] - XB[j]) ** 2 + (YW - YB[j]) ** 2;  # B term
            Cx = -np.cos(phi[j]);  # C term (X-direction)
            Dx = XW[h] - XB[j];  # D term (X-direction)
            Cy = -np.sin(phi[j]);  # C term (Y-direction)
            Dy = YW[h] - YB[j];  # D term (Y-direction)
            E = math.sqrt(B - A ** 2);  # E term
            if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):
                MWjx[h ,j] = 0
                MWjy[h, j] = 0
            else:

                term1 = 0.5 * Cx * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
                term2 = ((Dx - A * Cx) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
                MWjx[h ,j] = term1 + term2;

                term1 = 0.5 * Cy * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
                term2 = ((Dy - A * Cy) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
                MWjy[h, j] = term1 + term2;

            # Zero out any problem values
            if (np.iscomplex(MWjx[h, j]) or np.isnan(MWjx[h, j]) or np.isinf(MWjx[h, j])):
                MWjx[h, j] = 0
            if (np.iscomplex(MWjy[h, j]) or np.isnan(MWjy[h, j]) or np.isinf(MWjy[h, j])):
                MWjy[h, j] = 0

    return MWjx, MWjy
