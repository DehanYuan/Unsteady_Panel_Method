import numpy as np
import math as math
import matplotlib.pyplot as plt
from matplotlib import path

from XFOIL import XFOIL
from IJ_SPM import IJ_SPM
from KL_VPM import KL_VPM
from Tit_Tin import Tit_Tin
from Ait_Ain import Ait_Ain
from Diqt import Diqt
from Wihn import Wihn
from Wiht import Wiht
# from MnjxMnjy import MxMy
# from NnjxNnjy import NxNy
from Circulation import Circulation

# Flag to specify creating or loading airfoil
flagAirfoil = [1,  # Create specified NACA airfoil in XFOIL
               0]  # Load Selig-format airfoil from directory

# User-defined knowns
Vinf = 1  # Freestream velocity [] (just leave this at 1)
AoA = 8  # Angle of attack [deg]
NACA = '0012'  # NACA airfoil to load [####]
zeta = 0.05   # dimensionless distance between leading edge and pitch axis position


# Initialize parameters
time_step = 0.02 # Setting time_step
ts= 1 # Setting the time you want to calculate

Kt = ts / time_step


# Basic parameters
h_m = 0   # non-dimensional heave amplitude
# frequency
f = 2 #Hz
# heave motion eq. : h(t) = h_m(t) * sin(2 * pi * f * t)
# pitch motion eq. : theta(t) = Theta_m * sin(2 * pi * f * t + varphi)
theta_m = 30 / np.pi  # Pitch Amplitude
varphi = 90 / np.pi  # Phase difference
h = np.zeros(Kt)
theta = np.zeros(Kt)
omega = np.zeros(Kt)


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

for k in range(Kt):
    h[k] = h_m * np.sin(2 * np.pi * f * k * time_step)
    theta[k] = theta_m * np.sin(2 * np.pi * f * k * time_step + varphi)
    omega[k] = 2 * np.pi * f * np.cos(2 * np.pi * k * time_step + varphi) # Compute the angular velocity

# {
#     All the circle code with k write in here:





# Separation vortex and Trailing-edge vortex parameters definition
DeltaS_T = 0.1 # Initialize DeltaS_T
Theta_T = 15 / np.pi  # Initialize Theta_T
DeltaS_Sl = 0.1    # Initialize DeltaS_S
Theta_Sl = 15 / np.pi  # Initialize Theta_S



# compute the perimeter of airfoil
l = numPan * S
# Define the number of wake vortex
numvor = np.zeros(ts)



# Compute relative angle of attack AoAr, update the airfoil panel location
xa = zeta
ya = 0
x_rC = np.zeros(numPan)
y_rC = np.zeros(numPan)
Vs_x = np.zeros(numPan)
Vs_y = np.zeros(numPan)
x_rB = np.zeros(numPan)
y_rB = np.zeros(numPan)
for i in range(numPan):
    x_rC[i] = abs(XC[i] - xa)
    y_rC[i] = abs(YC[i])
    x_rB[i] = abs(XB[i] - xa)
    y_rB[i] = abs(YB[i])

for k in range(Kt):

    for i in range(numPan):  # Foil control points and boundary points update
        XC[i] = XC[i] - y_rC[i] * omega
        YC[i] = YC[i] + h + x_rC * omega
        XB[i] = XC[i] - y_rB[i] * omega
        YB[i] = YB[i] - y_rB[i] * omega

# Wake vortex location and Separation vortex location define
XW = np.zeros([numPan])
YW = np.zeros([numPan])
XS = np.zeros([numPan])
YS = np.zeros([numPan])

# Compute normal vector angle delta_i and define the tangent vector of foil panel
tau_x = np.zeros(numPan)
tau_y = np.zeros(numPan)
delta_i = np.zeros(numPan)
for i in range(numPan):
    delta_i[i] = 90 + phi
    tau_x[i] = (XB[i+1] - XB[i])/S
    tau_y[i] = (YB[i+1] - YB[i])/S

# Define Trailing edge vortex location
xT = XB[0] + 0.5 * DeltaS_T * np.cos(Theta_T)
yT = YB[0] + 0.5 * DeltaS_T * np.sin(Theta_T)

# Determine fixed separation location
nfxsl = 0.05
fxsl = np.zeros(numPan)
for i in range(numPan):
    fxsl[i] = XB[i] - fxsl
islindex = fxsl.index(min(fxsl))
x_isl = XB[islindex]
y_isl = YB[islindex]  # Define latest separation vortex associate with the airfoil location

x_Sl = x_isl + 0.5 * DeltaS_Sl * np.cos(Theta_Sl)
y_Sl = y_isl + 0.5 * DeltaS_Sl * np.sin(Theta_Sl)  # Define latest separation vortex control point location


# Compute relative AoA
for i in range(numPan):
    Vs_x[i] = Vinf * np.cos(AoA) + omega * y_rC[i]
    Vs_y[i] = Vinf * np.sin(AoA) - omega * x_rC[i] - h

cx = np.cos(theta)
cy = np.sin(theta)
AoAri = np.zeros(numPan)
for i in range(numPan):
    AoAri[i] = math.acos((Vs_x * cx + Vs_y * cy) / (math.sqrt(Vs_x**2 + Vs_y**2) + math.sqrt(cx**2 + cy**2)))

sumAoAri = 0
for i in range(numPan):
    sumAoAri += AoAri[i]
AoAr = sumAoAri / numPan




# Determine the separation location
# if -10 <= AoAr <= 10:
#     DeltaS_S = 0
#     Theta_S = 0
# else:
#     DeltaS_S =






# Geometric integrals for SPM and VPM (normal [I,K] and tangential [J,L])
I, J = IJ_SPM(XC, YC, XB, YB, phi, S)
K, L = KL_VPM(XC, YC, XB, YB, phi, S)






# Compute one time-step until DeltaS_T & Theta_T converge

# Compute normal velocity of airfoil panels
ABM = np.zeros(numPan)
for i in range(numPan):
    ABM = 2 * np.pi * (Vs_x[i] * np.cos(delta_i) + Vs_y * np.sin(delta_i))  # Compute AM term of BC.eq.

BBM = np.zeros([numPan, numPan])
BBM1 = np.zeros([numPan, numPan])
for i in range(numPan):
    for j in range(numPan):
        if (i==j):
            BBM1[i, j] = np.pi
        else:BBM1[i, j] = I[i, j]  # Compute BM term of BC.eq.

BBM = np.dot(BBM1, lamb)



CBM = np.zeros(numPan)
for i in range(numPan):
    for j in range(numPan):
        if (i==j):
            CBM[i, j] = 0
        else:CBM[i, j] = K[i, j]  # Compute CM term of BC.eq.


M = np.zeros(Kt)
Gamma_wh = np.zeros(len(M))
Wihn = Wihn(XC, YC, XW, YW, delta_i, M)
DBM = np.zeros(numPan)
for i in range(numPan):
    for h in range(M):
        DBM[i, h] = np.dot(Wihn, Gamma_wh)  # Compute DM term of BC

# Compute EM term of BC
Tin, Tit = Tit_Tin
Ain, Ait = Ait_Ain   # recall coefficients

EBM = np.zeros(numPan)
EBM = Tin

# Compute FM term of BC

FBM = np.zeros(numPan)
FBM = Ain

# Compute GM term of BC

GBM = np.zeros(numPan)
Q = np.zeros(Kt)
Gamma_Sq = np.zeros(len(Q))
Diqt = Diqt(XC, YC, XS, YS, tau_x, tau_y, Q)
GBM = np.zeros([numPan, Q])
for i in range(numPan):
    for q in range(Q):
        GBM[i, h] = np.dot(Wiht, Gamma_Sq)  # Compute DM term of BC

ABM + BBM +CBM + DBM + EBM + FBM +GBM == 0



# Compute the 1th & the Nth panel tangent velocity

ACM = np.zeros(numPan)
for i in range(numPan):
    ACM = 2 * np.pi * (Vs_x[i] * tau_x + Vs_y * tau_y)   # Compute A term of tangent velocity to calculate KC

BCM1 = np.zeros([numPan, numPan])
BCM = np.zeros([numPan, numPan])
for i in range(numPan):
    for j in range(numPan):
        if (i==j):
            BCM1[i, j] = 0
        else:BCM1 = J[i, j]
        BCM = np.dot(BCM1, lamb)  # Compute B term of tangent velocity to calculate KC

CCM = np.zeros(numPan)
CCM1 = np.zeros(numPan, numPan)
for i in range(numPan):
    for j in range(numPan):
        if (i==j):
            CCM[i, j] = 0
        else:CCM1 = L[i, j]  # Compute C term of tangent velocity to calculate KC

CCM = CCM1 * gamma_Sl

DCM = np.zeros(numPan)
for i in range(numPan):
    for h in range(M):
        DCM[i, h] = np.dot(Wiht, Gamma_wh)  # Compute D term of tangent velocity to calculate KC

ECM = np.zeros(numPan)
for i in range(numPan):
    ECM = Tit  # Compute E term of tangent velocity to calculate KC

FCM = np.zeros(numPan)
for i in range(numPan):
    FCM = Ait  # Compute F term of tangent velocity to calculate KC

GCM = np.zeros([numPan, Q])
for i in range(numPan):
    for q in range(Q):
        GCM = np.dot(Diqt[i, q], Gamma_Sq)  # Compute G term of tangent velocity to calculate KC

V_it = np.zeros(numPan)
V_it = ACM + BCM + CCM + DCM + ECM + FCM + GCM







# V_it = np.zeros(numPan)
# Vs_xt = np.zeros(numPan)
# Vs_yt = np.zeros(numPan)
# Vs_xt[i] = Vs_x[i] * np.cos(phi[i])
# Vs_yt[i] = -Vs_y[i] * np.sin(phi[i])
# t1 = np.zeros(numPan)
# t2 = np.zeros(numPan)
# t3 = np.zeros(numPan)
# t4 = np.zeros(numPan)
# t5 = np.zeros(numPan)
# t6 = np.zeros(numPan)
# for i in range(numPan):
#     t1 = Vs_xt[i] + Vs_yt[i]
#     t2 = (1 / (2 * np.pi)) * sum(lamb[k] * J[i, :])
#     t3 = (1 / (2 * np.pi)) * sum(gamma[k] * L[i, :])
#     t4 = gamma[k] / 2
#     t5 = (1 / (2 * np.pi)) * sum((gamma[k - 1] - gamma[k]) * l * Tint[i])
#     for h in range(numvor):
#       t6 = (1 / (2 * np.pi)) * sum((gamma[k - h -1] - gamma[k - h]) * l * Wiht[i, h])
# V_it[i] = t1 + t2 + t3 + t4 + t5 + t6
#
# # Compute KC.eq. results
# gamma[k] = (V_it[0]**2 - V_it[numPan-1]**2) * time_step/ (2 * l) + gamma[k-1]
#
#
# resArr = np.linalg.solve(A, b)


# Compute KC.eq. result









# Call Tinn calculation about Trailing edge vortex
# Tinn = Tinn(XC, YC, Xn, Yn, phi)


# # Judge whether the time of state reaches the critical time t_cr
#
# if 2 * np.pi * f * theta_m * np.cos(2 * np.pi * t_k + varphi) == 0:
#     numvor += numvor
# # compute Circulation at time-step k
# Gamma_k = gamma * numPan
# # Compute V_Txk & V_Tyk
# V_Tx =
# V_Ty =
# # Updated DeltaS_T & Theta_T
# DeltaS_T = time_step * math.sqrt(V_Txk**2+V_Tyk**2)
# # Call Wihn to compute the airfoil panels velocity induced by wake vortex core
# XW[h] = V_Txk * time_step
# Xh = Xh + 1
# Wihn = Wihn(XC, YC, Xh, Yh, phi)


# }
