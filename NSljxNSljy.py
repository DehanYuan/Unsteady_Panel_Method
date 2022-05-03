
import numpy as np
import math as math

np.seterr('raise')


def NSljxNSljy(xSl, ySl, XB, YB, phi, S):
    # Number of panels
    numPan = len(XB) - 1  # Number of panels (control points)

    # Initialize arrays
    NSljx = np.zeros(numPan)
    NSljy = np.zeros(numPan)


    for j in range(numPan):  # Loop over all panels
        # Compute intermediate values
        A = -(xSl - XB[j]) * np.cos(phi[j]) - (ySl - YB[j]) * np.sin(phi[j])  # A term
        B = (xSl - XB[j]) ** 2 + (ySl - YB[j]) ** 2  # B term
        Cx = np.sin(phi[j])  # Cx term (X-direction)
        Dx = -(ySl - YB[j])  # Dx term (X-direction)
        Cy = -np.cos(phi[j])  # Cy term (Y-direction)
        Dy = xSl - XB[j]  # Dy term (Y-direction)
        E = math.sqrt(B - A ** 2)  # E term
        if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):  # If E term is 0 or complex or a NAN or an INF
            NSljx[j] = 0
            NSljy[j] = 0
        else:

            term1 = 0.5 * Cx * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
            term2 = ((Dx - A * Cx) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
            NSljx[j] = term1 + term2;


            term1 = 0.5 * Cy * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B);
            term2 = ((Dy - A * Cy) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E));
            NSljy[j] = term1 + term2;

        # Zero out any problem values
        if (np.iscomplex(NSljx[j]) or np.isnan(NSljx[j]) or np.isinf(NSljx[j])):
            NSljx[j] = 0
        if (np.iscomplex(NSljy[j]) or np.isnan(NSljy[j]) or np.isinf(NSljy[j])):
            NSljy[j] = 0

    return NSljx, NSljy
