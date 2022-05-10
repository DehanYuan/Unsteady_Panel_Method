
import numpy as np
import math as math

np.seterr('raise')


def ATxATy(x_isl, y_isl, xTB, yTB, Theta_Sl, DeltaS_Sl):

        # Compute intermediate values
    A = -(xTB - x_isl) * np.cos(Theta_Sl) - (yTB - y_isl) * np.sin(Theta_Sl)  # A term
    B = (xTB - x_isl) ** 2 + (yTB - y_isl) ** 2  # B term
    Cx = np.sin(Theta_Sl)  # Cx term (X-direction)
    Dx = -(yTB - y_isl)  # Dx term (X-direction)
    Cy = -np.cos(Theta_Sl)  # Cy term (Y-direction)
    Dy = xTB - x_isl  # Dy term (Y-direction)
    E = math.sqrt(B - A ** 2)  # E term

    if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):  # If E term is 0 or complex or a NAN or an INF
        ATx = 0
        ATy = 0
    else:

        term1 = 0.5 * Cx * np.log((DeltaS_Sl ** 2 + 2 * A * Theta_Sl + B) / B);  # First term in Nx equation
        term2 = ((Dx - A * Cx) / E) * (math.atan2((DeltaS_Sl + A), E) - math.atan2(A, E));
        ATx = term1 + term2;

        # Compute Ny, Ref [1]
        term1 = 0.5 * Cy * np.log((DeltaS_Sl ** 2 + 2 * A * DeltaS_Sl + B) / B);  # First term in Ny equation
        term2 = ((Dy - A * Cy) / E) * (math.atan2((DeltaS_Sl + A), E) - math.atan2(A, E));  # Second term in Ny equation
        ATy = term1 + term2;  # Compute Ny integral

        # Zero out any problem values
    if (np.iscomplex(ATx) or np.isnan(ATx) or np.isinf(ATx)):  # If Nx term is complex or a NAN or an INF
        ATx = 0  # Set Nx value equal to zero
    if (np.iscomplex(ATy) or np.isnan(ATy) or np.isinf(ATy)):  # If Ny term is complex or a NAN or an INF
        ATy = 0

    return ATx, ATy
