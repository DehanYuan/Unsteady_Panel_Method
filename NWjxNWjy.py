
import numpy as np
import math as math

np.seterr('raise')


def NWjxNWjy(XW, YW, XB, YB, phi, S, M):
    # Number of panels
    numPan = len(XB) - 1  # Number of panels (control points)

    # Initialize arrays
    NWjx = np.zeros([M, numPan])
    NWjy = np.zeros([M, numPan])

    for h in range(M):
        for j in range(numPan):  # Loop over all panels
            # Compute intermediate values
            A = -(XW[h] - XB[j]) * np.cos(phi[j]) - (YW[h] - YB[j]) * np.sin(phi[j])  # A term
            B = (XW[h] - XB[j]) ** 2 + (YW - YB[j]) ** 2  # B term
            Cx = np.sin(phi[j])  # Cx term (X-direction)
            Dx = -(YW[h] - YB[j])  # Dx term (X-direction)
            Cy = -np.cos(phi[j])  # Cy term (Y-direction)
            Dy = XW[h] - XB[j]  # Dy term (Y-direction)
            E = math.sqrt(B - A ** 2)  # E term
            if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(
                    E)):  # If E term is 0 or complex or a NAN or an INF
                NWjx[h, j] = 0
                NWjy[h, j] = 0
            else:

                term1 = 0.5 * Cx * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
                term2 = ((Dx - A * Cx) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
                NWjx[h, j] = term1 + term2;

                term1 = 0.5 * Cy * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
                term2 = ((Dy - A * Cy) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
                NWjy[h, j] = term1 + term2;

            # Zero out any problem values
            if (np.iscomplex(NWjx[h, j]) or np.isnan(NWjx[h, j]) or np.isinf(NWjx[h, j])):
                NWjx[h, j] = 0
            if (np.iscomplex(NWjy[h, j]) or np.isnan(NWjy[h, j]) or np.isinf(NWjy[h, j])):
                NWjy[h, j] = 0

    return NWjx, NWjy
