import numpy as np
import math as math

np.seterr('raise')


def MTjxMTjy(xT, yT, XB, YB, phi, S):
    # Number of panels
    numPan = len(XB) - 1  # Number of panels

    # Initialize arrays
    MTjx = np.zeros(numPan)  # Initialize Ix integral array
    MTjy = np.zeros(numPan)  # Initialize Iy integral array

    # Compute integral
    for j in range(numPan):  # Loop over all panels
        # Compute intermediate values
        A = -(xT - XB[j]) * np.cos(phi[j]) - (yT - YB[j]) * np.sin(phi[j])  # A term
        B = (xT - XB[j]) ** 2 + (yT - YB[j]) ** 2;  # B term
        Cx = -np.cos(phi[j]);  # C term (X-direction)
        Dx = xT - XB[j];  # D term (X-direction)
        Cy = -np.sin(phi[j]);  # C term (Y-direction)
        Dy = yT - YB[j];  # D term (Y-direction)
        E = math.sqrt(B - A ** 2);  # E term
        if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):  # If E term is 0 or complex or a NAN or an INF
            MTjx[j] = 0
            MTjy[j] = 0
        else:

            term1 = 0.5 * Cx * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
            term2 = ((Dx - A * Cx) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
            MTjx[j] = term1 + term2;


            term1 = 0.5 * Cy * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
            term2 = ((Dy - A * Cy) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
            MTjy[j] = term1 + term2;

        # Zero out any problem values
        if (np.iscomplex(MTjx[j]) or np.isnan(MTjx[j]) or np.isinf(MTjx[j])):
            MTjx[j] = 0
        if (np.iscomplex(MTjy[j]) or np.isnan(MTjy[j]) or np.isinf(MTjy[j])):
            MTjy[j] = 0

    return MTjx, MTjy
