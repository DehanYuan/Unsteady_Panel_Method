
import numpy as np
import math as math

np.seterr('raise')


def TSlxTSly(x_isl, y_isl, xTB, yTB, Theta_T, DeltaS_T):
        # Compute intermediate values
    A = -(x_isl - xTB) * np.cos(Theta_T) - (y_isl - yTB) * np.sin(Theta_T)  # A term
    B = (x_isl - xTB) ** 2 + (y_isl - yTB) ** 2  # B term
    Cx = np.sin(Theta_T)  # Cx term (X-direction)
    Dx = -(y_isl - yTB)  # Dx term (X-direction)
    Cy = -np.cos(Theta_T)  # Cy term (Y-direction)
    Dy = x_isl - xTB  # Dy term (Y-direction)
    E = math.sqrt(B - A ** 2)  # E term

    if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):  # If E term is 0 or complex or a NAN or an INF
        TSlx = 0
        TSly = 0
    else:

        term1 = 0.5 * Cx * np.log((DeltaS_T ** 2 + 2 * A * Theta_T + B) / B);
        term2 = ((Dx - A * Cx) / E) * (math.atan2((DeltaS_T + A), E) - math.atan2(A, E));
        TSlx = term1 + term2;


        term1 = 0.5 * Cy * np.log((DeltaS_T ** 2 + 2 * A * DeltaS_T + B) / B);
        term2 = ((Dy - A * Cy) / E) * (math.atan2((DeltaS_T + A), E) - math.atan2(A, E));
        TSly = term1 + term2;

        # Zero out any problem values
    if (np.iscomplex(TSlx) or np.isnan(TSlx) or np.isinf(TSlx)):
        TSlx = 0
    if (np.iscomplex(TSly) or np.isnan(TSly) or np.isinf(TSly)):
        TSly = 0

    return TSlx, TSly
