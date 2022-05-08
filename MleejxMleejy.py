import numpy as np
import math as math

np.seterr('raise')


def MleejxMleejy(XCP, YCP, XB, YB, phi, S, numPlee):
    # Number of panels
    numPan = len(XB) - 1  # Number of panels

    # Initialize arrays
    Mleejx = np.zeros(numPlee, numPan)
    Mleejy = np.zeros(numPlee, numPan)

    # Compute integral
    for f in range(1, numPlee+1):
        for j in range(numPan):  # Loop over all panels
            # Compute intermediate values
            A = -(XCP[f] - XB[j]) * np.cos(phi[j]) - (YCP[f] - YB[j]) * np.sin(phi[j])  # A term
            B = (XCP[f] - XB[j]) ** 2 + (YCP[f] - YB[j]) ** 2;  # B term
            Cx = -np.cos(phi[j]);  # C term (X-direction)
            Dx = XCP[f] - XB[j];  # D term (X-direction)
            Cy = -np.sin(phi[j]);  # C term (Y-direction)
            Dy = YCP[f] - YB[j];  # D term (Y-direction)
            E = math.sqrt(B - A ** 2);  # E term
            if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):
                Mleejx[f, j] = 0
                Mleejy[f, j] = 0
            else:

                term1 = 0.5 * Cx * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
                term2 = ((Dx - A * Cx) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
                Mleejx[f, j] = term1 + term2;

                term1 = 0.5 * Cy * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
                term2 = ((Dy - A * Cy) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
                Mleejy[f, j] = term1 + term2;

            # Zero out any problem values
            if (np.iscomplex(Mleejx[f, j]) or np.isnan(Mleejx[f, j]) or np.isinf(Mleejx[f, j])):
                Mleejx[f, j] = 0
            if (np.iscomplex(Mleejy[f, j]) or np.isnan(Mleejy[f, j]) or np.isinf(Mleejy[f, j])):
                Mleejy[f, j] = 0



    return Mleejx, Mleejy
