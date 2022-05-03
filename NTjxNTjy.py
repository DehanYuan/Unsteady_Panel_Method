
import numpy as np
import math as math

np.seterr('raise')


def NTjxNTjy(xT, yT, XB, YB, phi, S):
    # Number of panels
    numPan = len(XB) - 1  # Number of panels (control points)

    # Initialize arrays
    NTjx = np.zeros(numPan)
    NTjy = np.zeros(numPan)


    for j in range(numPan):  # Loop over all panels
        # Compute intermediate values
        A = -(xT - XB[j]) * np.cos(phi[j]) - (xT - YB[j]) * np.sin(phi[j])  # A term
        B = (xT - XB[j]) ** 2 + (yT - YB[j]) ** 2  # B term
        Cx = np.sin(phi[j])  # Cx term (X-direction)
        Dx = -(yT - YB[j])  # Dx term (X-direction)
        Cy = -np.cos(phi[j])  # Cy term (Y-direction)
        Dy = xT - XB[j]  # Dy term (Y-direction)
        E = math.sqrt(B - A ** 2)  # E term
        if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):  # If E term is 0 or complex or a NAN or an INF
            NTjx[j] = 0
            NTjy[j] = 0
        else:

            term1 = 0.5 * Cx * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
            term2 = ((Dx - A * Cx) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
            NTjx[j] = term1 + term2;


            term1 = 0.5 * Cy * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
            term2 = ((Dy - A * Cy) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
            NTjy[j] = term1 + term2;

        # Zero out any problem values
        if (np.iscomplex(NTjx[j]) or np.isnan(NTjx[j]) or np.isinf(NTjx[j])):
            NTjx[j] = 0
        if (np.iscomplex(NTjy[j]) or np.isnan(NTjy[j]) or np.isinf(NTjy[j])):
            NTjy[j] = 0

    return NTjx, NTjy
