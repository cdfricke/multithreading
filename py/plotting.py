# PLOTTING SCRIPT FOR rejection_sampling.cpp
# most of this comes from Rejection_Sampling_V6.py by Payton Linton (linton.93@osu.edu)
# Adapted for work with multithreading program by Connor Fricke (fricke.59@osu.edu)
# Latest Rev: 24-Feb-2025

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

# * GLOBAL CONSTS *
k = 0.0
qk = 0.0
GeoConst = 0.0

# * PROGRAM PARAMS *
SHAPE = int(sys.argv[1])
DIST = int(sys.argv[2])
DMIN = float(sys.argv[3])
DMAX = float(sys.argv[4])
DSTEP = float(sys.argv[5])
DEPTH = float(sys.argv[6])
AREA = float(sys.argv[7])

# * INTERPOLATE GLOBAL VARS FROM PARAMS *
shapeStr = ""
distStr = ""
if SHAPE == 0:
    shapeStr = "SPHERE"
    GeoConst = 4.0 / np.pi
elif SHAPE == 1:
    shapeStr = "CUBE"
    GeoConst = 1.0
elif SHAPE == 2:
    shapeStr = "CUBOID"
    GeoConst = 1.0 / (0.8*0.54)
else:
    print("Something went wrong! (Shape Invalid)")
    exit() 

if (DIST == 0): 
    distStr = "CE3"
    k = 0.0125
    qk = 1.743
elif (DIST == 1):
    distStr = "CE4"
    k = 0.0021
    qk = 0.5648 + 0.01258 / k
else:
    print("Something went wrong! (Distribution Invalid)")
    exit()

# * FUNCTIONS WE STILL NEED *
def nDpm2(D):
    return GeoConst*k*qk*np.exp(-qk*D)/D**2
def NgtD4_Di(D):
    from scipy.integrate import quad
    Diam_array = Sampled_D
    N = []
    for i in range(len(Diam_array)):
        N.append(quad(nDpm2, Diam_array[i], np.inf)[0])
    return N
def Fk_Di(D):
    return k*np.exp(-qk*D)

# * GENERATE ARRAYS AND DF *
rockDataFileName = './data/Rock_Data_' + str(int(AREA)) + 'XY_' + str(int(DEPTH)) + 'Z_' + shapeStr + '.csv'
df_RockData = pd.read_csv(rockDataFileName, delimiter=',')
Avg_Sampled_NgtD = np.loadtxt('./data/Avg_Sampled_NgtD.csv', delimiter=',')
Avg_Sampled_FgtD = np.loadtxt('./data/Avg_Sampled_FgtD.csv', delimiter=',')
Sampled_D = np.arange(DMIN, DMAX, DSTEP)

# if C++ generated arrays are larger or smaller than np.arange(DMIN, DMAX, DSTEP), then truncate to match lengths
while (len(Avg_Sampled_NgtD) > len(Sampled_D)):
    Avg_Sampled_NgtD = Avg_Sampled_NgtD[:-1]
while (len(Avg_Sampled_FgtD) > len(Sampled_D)):
    Avg_Sampled_FgtD = Avg_Sampled_FgtD[:-1]
while (len(Avg_Sampled_NgtD) < len(Sampled_D)):
    Sampled_D = Sampled_D[:-1]
while (len(Avg_Sampled_FgtD) < len(Sampled_D)):
    Sampled_D = Sampled_D[:-1]

# * PLOT *
fig, axs = plt.subplots(2,2)
axs[0,0].plot(Sampled_D, NgtD4_Di(DMIN), label = 'Wu et al. (2021)')
axs[0,0].plot(Sampled_D, Avg_Sampled_NgtD, label = "Rejection Sampling using " + distStr)
axs[0,0].legend()
axs[0,0].set_xlabel('D [m]')
axs[0,0].set_ylabel('$Rocks/m^2$')
axs[0,0].set_title("Number of rocks/$m^2$ with diameter > D, Area = " + str(AREA) +', Depth = ' + str(DEPTH) + ': ' + shapeStr)

# BWOKEN
axs[0,1].plot(Sampled_D, Fk_Di(Sampled_D), label = 'Wu et al. (2021)')
axs[0,1].plot(Sampled_D, Avg_Sampled_FgtD, label = "Rejection Sampling using " + distStr)
axs[0,1].legend()
axs[0,1].set_xlabel('D [m]')
axs[0,1].set_ylabel('Fractional Area')
axs[0,1].set_title("Fractional area covered by rocks with diameter > D Area = " + str(AREA) +', Depth = ' + str(DEPTH) + ': ' + shapeStr)

axs[1,0].plot(Sampled_D, NgtD4_Di(DMIN), label = 'Wu et al. (2021)')
axs[1,0].plot(Sampled_D, Avg_Sampled_NgtD, label = "Rejection Sampling using " + distStr)
axs[1,0].set_yscale('log')
axs[1,0].set_xscale('log')
axs[1,0].legend()
axs[1,0].set_xlabel('D [m]')
axs[1,0].set_ylabel('$Rocks/m^2$')
axs[1,0].set_title("Number of rocks/$m^2$ with diameter > D Area = " + str(AREA) +', Depth = ' + str(DEPTH) + ': ' + shapeStr)

# BWOKEN
axs[1,1].plot(Sampled_D, Fk_Di(Sampled_D), label = 'Wu et al. (2021)')
axs[1,1].plot(Sampled_D, Avg_Sampled_FgtD, label = "Rejection Sampling using " + distStr)
axs[1,1].set_yscale('log')
axs[1,1].set_xscale('log')
axs[1,1].legend()
axs[1,1].set_xlabel('D [m]')
axs[1,1].set_ylabel('Fractional Area')
axs[1,1].set_title("Fractional area covered by rocks with diameter > D Area = " + str(AREA) +', Depth = ' + str(DEPTH) + ': ' + shapeStr)
plt.tight_layout()
plt.show()