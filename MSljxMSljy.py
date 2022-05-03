import numpy as np
import math as math

np.seterr('raise')


def MSljxMSljy(xSl, ySl, XB, YB, phi, S):
    # Number of panels
    numPan = len(XB) - 1  # Number of panels

    # Initialize arrays
    MSljx = np.zeros(numPan)
    MSljy = np.zeros(numPan)

    # Compute integral
    for j in range(numPan):  # Loop over all panels
        # Compute intermediate values
        A = -(xSl - XB[j]) * np.cos(phi[j]) - (ySl - YB[j]) * np.sin(phi[j])  # A term
        B = (xSl - XB[j]) ** 2 + (ySl - YB[j]) ** 2;  # B term
        Cx = -np.cos(phi[j]);  # C term (X-direction)
        Dx = xSl - XB[j];  # D term (X-direction)
        Cy = -np.sin(phi[j]);  # C term (Y-direction)
        Dy = ySl - YB[j];  # D term (Y-direction)
        E = math.sqrt(B - A ** 2);  # E term
        if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):
            MSljx[j] = 0
            MSljy[j] = 0
        else:

            term1 = 0.5 * Cx * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
            term2 = ((Dx - A * Cx) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
            MSljx[j] = term1 + term2;


            term1 = 0.5 * Cy * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
            term2 = ((Dy - A * Cy) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
            MSljy[j] = term1 + term2;

        # Zero out any problem values
        if (np.iscomplex(MSljx[j]) or np.isnan(MSljx[j]) or np.isinf(MSljx[j])):
            MSljx[j] = 0
        if (np.iscomplex(MSljy[j]) or np.isnan(MSljy[j]) or np.isinf(MSljy[j])):
            MSljy[j] = 0

    return MSljx, MSljy
