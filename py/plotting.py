# PLOTTING SCRIPT FOR rejection_sampling.cpp
# most of this comes from Rejection_Sampling_V6.py by Payton Linton (linton.93@osu.edu)
# Adapted for work with multithreading program by Connor Fricke (fricke.59@osu.edu)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

SHAPE = int(sys.argv[1])
DIST = int(sys.argv[2])
DMIN = float(sys.argv[3])
DMAX = float(sys.argv[4])
DSTEP = float(sys.argv[5])
DEPTH = float(sys.argv[6])
AREA = float(sys.argv[7])

shapeStr = ""
if SHAPE == 0:
    shapeStr = "SPHERE"
elif SHAPE == 1:
    shapeStr = "CUBE"
elif SHAPE == 2:
    shapeStr = "CUBOID_A"
elif SHAPE == 3:
    shapeStr = "CUBOID_B"
else:
    print("Something went wrong!")
    exit()

rockDataFileName = './data/Rock_Data_' + str(int(AREA)) + 'XY_' + str(int(DEPTH)) + 'Z_' + shapeStr + '.csv'
df = pd.read_csv(rockDataFileName, delimiter=',')
print(df.head())

# df.to_csv('./Rock_Data_' + str(AREA) + 'XY_' + str(DEPTH) + 'Z_' + shapeStr + '.csv')

# fig, axs = plt.subplots(2,2)
# axs[0,0].plot(Sampled_D, NgtD4_Di(D), label = 'Wu et al. (2021)')
# for i in range(len(Avg_Sampled_NgtD)):
#     axs[0,0].plot(Sampled_D, Avg_Sampled_NgtD[i], label = "Rejection Sampling using " + str(Dist))
# axs[0,0].legend()
# axs[0,0].set_xlabel('D [m]')
# axs[0,0].set_ylabel('$Rocks/m^2$')
# axs[0,0].set_title("Number of rocks/$m^2$ with diameter > D, Area = " + str(Area) +', Depth = ' + str(Depth) + ': ' + Shape)

# axs[0,1].plot(np.arange(), Fk_Di(Sampled_D), label = 'Wu et al. (2021)')
# for i in range(len(Avg_Sampled_FgtD)):
#     axs[0,1].plot(Sampled_D, Avg_Sampled_FgtD[i], label = "Rejection Sampling using " + str(Dist))
# axs[0,1].legend()
# axs[0,1].set_xlabel('D [m]')
# axs[0,1].set_ylabel('Fractional Area')
# axs[0,1].set_title("Fractional area covered by rocks with diameter > D Area = " + str(Area) +', Depth = ' + str(Depth) + ': ' + Shape)



# axs[1,0].plot(Sampled_D, NgtD4_Di(D), label = 'Wu et al. (2021)')
# for i in range(len(Avg_Sampled_NgtD)):
#     axs[1,0].plot(Sampled_D, Avg_Sampled_NgtD[i], label = "Rejection Sampling using " + str(Dist))
# axs[1,0].set_yscale('log')
# axs[1,0].set_xscale('log')
# axs[1,0].legend()
# axs[1,0].set_xlabel('D [m]')
# axs[1,0].set_ylabel('$Rocks/m^2$')
# axs[1,0].set_title("Number of rocks/$m^2$ with diameter > D Area = " + str(Area) +', Depth = ' + str(Depth) + ': ' + Shape)

# axs[1,1].plot(Sampled_D, Fk_Di(Sampled_D), label = 'Wu et al. (2021)')
# for i in range(len(Avg_Sampled_FgtD)):
#     axs[1,1].plot(Sampled_D, Avg_Sampled_FgtD[i], label = "Rejection Sampling using " + str(Dist))
# axs[1,1].set_yscale('log')
# axs[1,1].set_xscale('log')
# axs[1,1].legend()
# axs[1,1].set_xlabel('D [m]')
# axs[1,1].set_ylabel('Fractional Area')
# axs[1,1].set_title("Fractional area covered by rocks with diameter > D Area = " + str(Area) +', Depth = ' + str(Depth) + ': ' + Shape)
# plt.tight_layout()
# plt.show()