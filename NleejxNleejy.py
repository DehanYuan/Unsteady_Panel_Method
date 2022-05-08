import numpy as np
import math as math

np.seterr('raise')


def NleejxNleejy(XCP, YCP, XB, YB, phi, S, numPlee):
    # Number of panels
    numPan = len(XB) - 1  # Number of panels (control points)

    # Initialize arrays
    Nleejx = np.zeros(numPlee, numPan)
    Nleejy = np.zeros(numPlee, numPan)

    for f in range(1, numPlee+1):
        for j in range(numPan):  # Loop over all panels
            # Compute intermediate values
            A = -(XCP[f] - XB[j]) * np.cos(phi[j]) - (YCP[f] - YB[j]) * np.sin(phi[j])  # A term
            B = (XCP[f] - XB[j]) ** 2 + (YCP[f] - YB[j]) ** 2  # B term
            Cx = np.sin(phi[j])  # Cx term (X-direction)
            Dx = -(YCP[f] - YB[j])  # Dx term (X-direction)
            Cy = -np.cos(phi[j])  # Cy term (Y-direction)
            Dy = XCP[f] - XB[j]  # Dy term (Y-direction)
            E = math.sqrt(B - A ** 2)  # E term
            if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(
                    E)):  # If E term is 0 or complex or a NAN or an INF
                Nleejx[f, j] = 0
                Nleejy[f, j] = 0
            else:

                term1 = 0.5 * Cx * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
                term2 = ((Dx - A * Cx) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
                Nleejx[f, j] = term1 + term2;

                term1 = 0.5 * Cy * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
                term2 = ((Dy - A * Cy) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
                Nleejy[f, j] = term1 + term2;

            # Zero out any problem values
            if (np.iscomplex(Nleejx[f, j]) or np.isnan(Nleejx[f, j]) or np.isinf(Nleejx[f, j])):
                Nleejx[f, j] = 0
            if (np.iscomplex(Nleejy[f, j]) or np.isnan(Nleejy[f, j]) or np.isinf(Nleejy[f, j])):
                Nleejy[f, j] = 0

    return Nleejx, Nleejy
