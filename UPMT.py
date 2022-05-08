import numpy as np
import math as math
import sympy as sym

from XFOIL import XFOIL
from IJ_SPM import IJ_SPM
from KL_VPM import KL_VPM
from Tit_Tin import Tit_Tin
from Ait_Ain import Ait_Ain
from Diqt import Diqt
from Diqn import Diqn
from Wihn import Wihn
from Wiht import Wiht
from MTjxMTjy import MTjxMTjy
from MSljxMSljy import MSljxMSljy
from NSljxNSljy import NSljxNSljy
from DSlhxDSlhy import DSlhxDSlhy
from V_SlxV_Sly import V_SlxV_Sly
from TSlxTSly import TSlxTSly
from DSlqxDSlqy import DSlqxDSlqy
from NTjxNTjy import NTjxNTjy
from DThxDThy import DThxDThy
from DTqxDTqy import DTqxDTqy
from ATxATy import ATxATy
from V_it import V_it
from NSjxNSjy import NSjxNSjy
from ASxASy import ASxASy
from DScxDScy import DScxDScy
from DShxDShy import DShxDShy
from MSjxMSjy import MSjxMSjy
from VSqxVSqy import VSqxVSqy
from TSxTSy import TSxTSy
from MWjxMWjy import MWjxMWjy
from NWjxNWjy import NWjxNWjy
from DWqxDWqy import DWqxDWqy
from DWmxDWmy import DWmxDWmy
from AWxAWy import AWxAWy
from TWxTWy import TWxTWy
from MleejxMleejy import MleejxMleejy
from NleejxNleejy import NleejxNleejy
from TleexTleey import TleexTleey
from AleexAleey import AleexAleey
from DleeqxDleeqy import DleeqxDleeqy
from DleehxDleehy import DleehxDleehy
from Vleex import Vleex
from PotentialC import PotentialC
# Flag to specify creating or loading airfoil
flagAirfoil = [1,  # Create specified NACA airfoil in XFOIL
               0]  # Load Selig-format airfoil from directory

# User-defined knowns
Vinf = 1  # Freestream velocity [] (just leave this at 1)
AoA = 8  # Angle of attack [deg]
NACA = '0012'  # NACA airfoil to load [####]

# Initialize parameters
time_step = 0.02  # Setting time_step
ts = 1  # Setting the time you want to calculate

Kt = ts / time_step  # Compute the number of steps

# Convert angle of attack to radians
AoAR = AoA * (np.pi / 180)  # Angle of attack [rad]

# Plotting flags
flagPlot = [0,  # Airfoil with panel normal vectors
            0,  # Geometry boundary pts, control pts, first panel, second panel
            1,  # Cp vectors at airfoil surface panels
            1,  # Pressure coefficient comparison (XFOIL vs. VPM)
            0,  # Airfoil streamlines
            0]  # Pressure coefficient contour

# %% XFOIL - CREATE/LOAD AIRFOIL

# PPAR menu options
PPAR = ['170',  # "Number of panel nodes"
        '4',  # "Panel bunching paramter"
        '1',  # "TE/LE panel density ratios"
        '1',  # "Refined area/LE panel density ratio"
        '1 1',  # "Top side refined area x/c limits"
        '1 1']  # "Bottom side refined area x/c limits"

# Call XFOIL function to obtain the following:
# - Airfoil coordinates
# - Pressure coefficient along airfoil surface
# - Lift, drag, and moment coefficients
xFoilResults = XFOIL(NACA, PPAR, AoA, flagAirfoil)

# Separate out results from XFOIL function results
afName = xFoilResults[0]  # Airfoil name
xFoilX = xFoilResults[1]  # X-coordinate for Cp result
xFoilY = xFoilResults[2]  # Y-coordinate for Cp result
xFoilCP = xFoilResults[3]  # Pressure coefficient
XB = xFoilResults[4]  # Boundary point X-coordinate
YB = xFoilResults[5]  # Boundary point Y-coordinate
xFoilCL = xFoilResults[6]  # Lift coefficient
xFoilCD = xFoilResults[7]  # Drag coefficient
xFoilCM = xFoilResults[8]  # Moment coefficient

# Number of boundary points and panels
numPts = len(XB)  # Number of boundary points
numPan = numPts - 1  # Number of panels (control points)

# %% CHECK PANEL DIRECTIONS - FLIP IF NECESSARY

# Check for direction of points
edge = np.zeros(numPan)  # Initialize edge value array
for i in range(numPan):  # Loop over all panels
    edge[i] = (XB[i + 1] - XB[i]) * (YB[i + 1] + YB[i])  # Compute edge values

sumEdge = np.sum(edge)  # Sum all edge values

# If panels are CCW, flip them (don't if CW)
if (sumEdge < 0):  # If panels are CCW
    XB = np.flipud(XB)  # Flip the X-data array
    YB = np.flipud(YB)  # Flip the Y-data array

# Initialize variables
XC = np.zeros(numPan)  # Initialize control point X-coordinate array
YC = np.zeros(numPan)  # Initialize control point Y-coordinate array
S = np.zeros(numPan)  # Initialize panel length array
phi = np.zeros(numPan)  # Initialize panel orientation angle array [deg]

# Find geometric quantities of the airfoil
for i in range(numPan):  # Loop over all panels
        XC[i] = 0.5 * (XB[i] + XB[i + 1])  # X-value of control point
        YC[i] = 0.5 * (YB[i] + YB[i + 1])  # Y-value of control point
        dx = XB[i + 1] - XB[i]  # Change in X between boundary points
        dy = YB[i + 1] - YB[i]  # Change in Y between boundary points
        S[i] = (dx ** 2 + dy ** 2) ** 0.5  # Length of the panel
        phi[i] = math.atan2(dy, dx)  # Angle of panel (positive X-axis to inside face)
        if (phi[i] < 0):  # Make all panel angles positive [rad]
            phi[i] = phi[i] + 2 * np.pi

# Compute angle of panel normal w.r.t. horizontal and include AoA
delta = phi + (np.pi / 2)  # Angle from positive X-axis to outward normal vector [rad]
beta = delta - AoAR  # Angle between freestream vector and outward normal vector [rad]
beta[beta > 2 * np.pi] = beta[beta > 2 * np.pi] - 2 * np.pi  # Make all panel angles between 0 and 2pi [rad]

# Separation vortex and Trailing-edge vortex parameters definition
DeltaS_T = 0.1  # Initialize DeltaS_T
Theta_T = 15 / np.pi  # Initialize Theta_T
DeltaS_Sl = 0.1    # Initialize DeltaS_S
Theta_Sl = 15 / np.pi  # Initialize Theta_S

# Parameters initialize
Q = 0  # number of Separate vortex
M = 0  # number of Wake vortex

gammal = 0  # gamma value in previous situation
V_itc = np.zeros(numPan)  # source strength of airfoil (different from panel to panel)
ABM = np.zeros(numPan)
BBM = np.zeros([numPan, numPan])
CBM = np.zeros([numPan, numPan])
DBM = np.zeros(numPan)
EBM = np.zeros(numPan)
FBM = np.zeros(numPan)
GBM = np.zeros(numPan)  # BC parameters

ACM = np.zeros(numPan)
BCM1 = np.zeros([numPan, numPan])
BCM = np.zeros([numPan, numPan])
CCM = np.zeros(numPan)
CCM1 = np.zeros([numPan, numPan])
DCM = np.zeros(numPan)
ECM = np.zeros(numPan)
FCM = np.zeros(numPan)
GCM = np.zeros(numPan)  # KC parameters

l = numPan * S  # compute the perimeter of airfoil
Vs_x = np.zeros(numPan)
Vs_y = np.zeros(numPan)
for i in range(numPan):
    Vs_x[i] = Vinf * np.cos(beta[i])
    Vs_y[i] = Vinf * np.sin(beta[i])

XW = np.zeros(M)
YW = np.zeros(M)
XS = np.zeros(int(Kt))
YS = np.zeros(int(Kt))
Gamma_wh = np.zeros(M)
Gamma_Sq = np.zeros(Q)
xSl = np.zeros(int(Kt))
ySl = np.zeros(int(Kt))  # Wake vortex location, Separation vortex location, LSV and TEV define

tau_x = np.zeros(numPan)
tau_y = np.zeros(numPan)
delta_i = np.zeros(numPan)  # normal vector angle of airfoil panels

Cp = np.zeros([int(Kt), numPan])
Phi = np.zeros(numPan)  # Pressure parameters
gamma = sym.symbols('gamma')  # vortex strength of airfoil
Phil = np.zeros(numPan)
CL = np.zeros(int(Kt))
CM = np.zeros(int(Kt))  # Coefficients definetion
# Determine fixed separation location
nfxsl = 0
fxsl = np.zeros(numPan)
for i in range(numPan):
    fxsl[i] = XB[i] - nfxsl
islindex = np.argmin(fxsl)
x_isl = XB[islindex]
y_isl = YB[islindex]  # Define latest separation vortex associate with the airfoil location

Xle = min(XB)
Yle = min(YB)
YCP = Yle
ile = np.argmin(XB)
numPlee = 100
XBP = np.zeros(numPlee)
XCP = np.zeros(numPlee)
XBP[0] = Xle


# {
#     All the circle code with k write in here:
for k in range(int(Kt)):
    # Compute normal vector angle delta_i and define the tangent vector of foil panel
    for i in range(numPan):
        delta_i[i] = 90 + phi[i]
        tau_x[i] = (XB[i + 1] - XB[i]) / S[i]  # cos(panel incline angle)
        tau_y[i] = (YB[i + 1] - YB[i]) / S[i]  # sin(panel incline angle)

    xT = XB[0] + 0.5 * DeltaS_T * np.cos(Theta_T)
    yT = YB[0] + 0.5 * DeltaS_T * np.sin(Theta_T)  # Compute Trailing edge vortex location

    xSl = x_isl + 0.5 * DeltaS_Sl * np.cos(Theta_Sl)
    ySl = y_isl + 0.5 * DeltaS_Sl * np.sin(Theta_Sl)  # Define latest separation vortex control point location

    # unit TEV strength computation
    gamma_T = (gamma - gammal) * l / DeltaS_T - gamma * DeltaS_Sl

    # Unit LSV strength computation
    gamma_Sl = gamma

    # Geometric integrals for SPM and VPM (normal [I,K] and tangential [J,L])
    I, J = IJ_SPM(XC, YC, XB, YB, phi, S)
    K, L = KL_VPM(XC, YC, XB, YB, phi, S)

    # Compute normal velocity of airfoil panels
    for i in range(numPan):
        ABM = 2 * np.pi * (Vs_x[i] * np.cos(delta_i[i]) + Vs_y[i] * np.sin(delta_i[i]))  # Compute AM term of BC.eq.

    for i in range(numPan):
        for j in range(numPan):
            if (i == j):
                BBM[i, j] = np.pi
            else:
                BBM[i, j] = I[i, j]  # Compute BM term of BC.eq.

    for i in range(numPan):
        for j in range(numPan):
            if (i == j):
                CBM[i, j] = 0
            else:
                CBM[i, j] = K[i, j]  # Compute CM term of BC.eq.


    Wihn = Wihn(XC, YC, XW, YW, delta_i, M)

    DBM = np.dot(Wihn, Gamma_wh)  # Compute DM term of BC

    # Compute EM term of BC
    Tin, Tit = Tit_Tin(XC, YC, xT, yT, phi, DeltaS_T, Theta_T)
    Ain, Ait = Ait_Ain(XC, YC, x_isl, y_isl, phi, DeltaS_Sl, Theta_Sl)
    EBM = Tin * gamma_T / (2 * np.pi)

    # Compute FM term of BC
    FBM = Ain * gamma_Sl / (2 * np.pi)

    # Compute GM term of BC
    Diqn = Diqn(XC, YC, XS, YS, delta_i, Q)

    GBM = np.dot(Diqn, Gamma_Sq)  # Compute DM term of BC

    BBM_inv = np.linalg.inv(BBM)
    lamb = BBM_inv * (
                -CBM * gamma - ABM - DBM - EBM - FBM - GBM)  # Obtain source strength in terms of unit vortex strength
    # V_it = V_it(XC, Vs_x, tau_x, Vs_y, tau_y, J, lamb, gamma, L, Tit, DCM, gamma_T, Wiht, M, Ait, Gamma_wh, gamma_Sl, Q,
    #             Diqt, Gamma_Sq)  # Compute the 1th to the Nth panel tangent velocity
    #
    # Kutta = sym.Eq(V_it[numPan - 1] ** 2 - V_it[0] ** 2, 2 * l * (gamma - gammal) / time_step)
    # solution = sym.solve([V_itc[0], V_itc[numPan - 1], Kutta],
    #                      (V_it[0], V_it[numPan - 1], gamma))  # solve the gamma value
    # gammal = gamma
    # lamb = BBM_inv * (-CBM * gamma - ABM - DBM - EBM - FBM - GBM)  # Recall source strength
    # V_it = V_it(XB, Vs_x, tau_x, Vs_y, tau_y, J, lamb, gamma, L, Tit, DCM, gamma_T, Wiht, M, Ait, Gamma_wh, gamma_Sl, Q,
    #             Diqt,
    #             Gamma_Sq)  # Recall tangent velocity
    #
    # # update latest separation vortex (LSV) and TEV location
    # MSljx, MSljy = MSljxMSljy(xSl, ySl, XB, YB, phi, S)
    # NSljx, NSljy = NSljxNSljy(xSl, ySl, XB, YB, phi, S)
    # DSlhx, DSlhy = DSlhxDSlhy(xSl, ySl, XW, YW, M)
    # DSlqx, DSlqy = DSlqxDSlqy(xSl, ySl, XS, YS, M)
    # TSlx, TSly = TSlxTSly(xSl, ySl, xT, yT, Theta_T, DeltaS_T)
    # V_Slx, V_Sly = V_SlxV_Sly(XC, lamb, MSljx, MSljy, NSljx, NSljy, gamma_T, Q, M, DSlhx, DSlhy, DSlqx, DSlqy, Gamma_Sq,
    #                           Gamma_wh, AoAR, Vinf, TSly, TSlx, gamma)
    # MTjx, MTjy = MTjxMTjy(xT, yT, XB, YB, phi, S)
    # NTjx, NTjy = NTjxNTjy(xT, yT, XB, YB, phi, S)
    # DThx, DThy = DThxDThy(xT, yT, XW, YW, M)
    # DTqx, DTqy = DTqxDTqy(xT, yT, XS, YS, M)
    # ATx, ATy = ATxATy(xSl, ySl, xT, yT, Theta_Sl, DeltaS_Sl)
    # V_Tx, V_Ty = V_SlxV_Sly(XC, lamb, MTjx, MTjy, NTjx, NTjy, gamma, Q, M, DThx, DThy, DTqx, DTqy, Gamma_Sq,
    #                         Gamma_wh, AoAR, Vinf, ATy, ATx)
    # Theta_Sl = math.atan2(V_Sly, V_Slx)
    # DeltaS_Sl = math.sqrt(V_Slx ** 2 + V_Sly ** 2) * time_step
    #
    # # update separate vortex location and Strength
    # XS.append(xSl + V_Slx * time_step)
    # YS.append(ySl + V_Sly * time_step)
    # Q += Q
    # Gamma_Sq.append(gamma * DeltaS_Sl)
    # MSjx, MSjy = MSjxMSjy(XS, YS, XB, YB, phi, S, Q)
    # NSjx, NSjy = NSjxNSjy(XS, YS, XB, YB, phi, S, M)
    # DShx, DShy = DShxDShy(XS, YS, XW, YW, M, Q)
    # DScx, DScy = DScxDScy(XS, YS, Q)
    # ASx, ASy = ASxASy(xSl, ySl, XS, YS, Theta_Sl, DeltaS_Sl, Q)
    # TSx, TSy = TSxTSy(XS, YS, xT, yT, Theta_T, DeltaS_T, Q)
    # VSqx, VSqy = VSqxVSqy(lamb, MSjx, MSjy, NSjx, NSjy, gamma, Q, DShx, DShy, DScx, DScy, Gamma_Sq,
    #                       Gamma_wh, AoAR, Vinf, ASy, ASx, gamma_Sl, gamma_T, TSx, TSy)
    # for q in range(Q):
    #     XS[q] = XS[q] + VSqx[q] * time_step
    #     YS[q] = YS[q] + VSqy[q] * time_step  # Separate vortex location update
    #
    # # update Wake vortex location and Strength
    # XW.append(xT + V_Tx * time_step)
    # YW.append(yT + V_Ty * time_step)
    # M += M
    # Gamma_wh.append(gamma_T * DeltaS_T)
    # MWjx, MWjy = MWjxMWjy(XW, YW, XB, YB, phi, S)
    # NWjx, NWjy = NWjxNWjy(XW, YW, XB, YB, phi, S)
    # DWmx, DWmy = DWmxDWmy(XW, YW, M)
    # DWqx, DWqy = DWqxDWqy(XS, YS, XW, YW, M, Q)
    # AWx, AWy = AWxAWy(xSl, ySl, XW, YW, Theta_Sl, DeltaS_Sl, M)
    # TWx, TWy = TWxTWy(XW, YW, xT, yT, Theta_T, DeltaS_T, M)
    # VWhx, VWhy = VSqxVSqy(lamb, MWjx, MWjy, NWjx, NWjy, gamma, M, DWqx, DWqy, DWmx, DWmy, Gamma_Sq,
    #                       Gamma_wh, AoAR, Vinf, AWy, AWx, gamma_Sl, gamma_T, TWx, TWy)
    # for h in range(M):
    #     XW[h] = XW[h] + VWhx * time_step
    #     YW[h] = YW[h] + VWhy * time_step  # update wake vortex location
    #
    # # Compute all airfoil panels potential
    #
    # for f in range(1, numPlee + 1):
    #     XBP[f] = (XBP[f - 1] - 10) / numPlee  # Obtain all boundary points of extending panels
    #
    # for f in range(1, numPlee + 1):
    #     XCP[f] = (XBP[f] + XBP[f - 1]) / 2  # Obtain all control points of extending panels
    #
    # Mleejx = MleejxMleejy(XCP, YCP, XB, YB, phi, S, numPlee)
    # Nleejx = NleejxNleejy(XCP, YCP, XB, YB, phi, S, numPlee)
    # Dleehx = DleehxDleehy(XW, YW, XCP, YCP, numPlee, M)
    # Dleeqx = DleeqxDleeqy(XS, YS, XCP, YCP, numPlee, Q)
    # Aleex = AleexAleey(xSl, ySl, XCP, YCP, Theta_Sl, DeltaS_Sl, numPlee)
    # Tleex = TleexTleey(XCP, YCP, xT, yT, Theta_T, DeltaS_T, numPlee)
    # Vleex = Vleex(lamb, Mleejx, Nleejx, gamma, numPlee, Dleehx, Gamma_Sq, Dleeqx,
    #               Gamma_wh, Vinf, Aleex, gamma_Sl, gamma_T, Tleex)  # Compute all the LE extending panels velocity
    # lePhi = sum(Vleex * 10)  # Compute LE potential
    # Phi = PotentialC(ile, V_it, S, lePhi, XC)  # Control points Potential
    # # Compute pressure coefficient of airfoil
    # for i in range(numPan):
    #     Cp[k, i] = 1 - (V_it / Vinf) ** 2 - (Phi[i] - Phil[i]) / time_step  # unsteady Bernolli's law
    #     Phil = Phi[i]
    # # Compute lift, drag and moment coefficients
    # CL[k] = sum(-Cp * S * np.sin(beta))  # Lift force coefficient []
    # # CD = sum(-Cp * S * np.cos(beta))  # Drag force coefficient []
    # CM[k] = sum(Cp * (XC - 0.25 * np.cos(AoAR)) * S * np.cos(phi))  # Moment coefficient []
    print(lamb)
# print(CL)
# print(CM)













