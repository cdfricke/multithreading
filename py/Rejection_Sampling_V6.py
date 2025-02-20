###################################################
#
#   Purpose: This code is meant to simulate populating
#            a volume according to the CE4 distribution.
#            The CE4 distribution is an area density and I
#            want to populate a volume SUCH THAT if I take any
#            slice in the z-direction, on average I can get back
#            the fractional area distribution reported by CE4
#            in Wu et al. (2021)
#
#   NOTE: Cuboids are slightly different. Spheres, and Cubes
#         are completely defined by D, but cuboids need a
#         diameter, D, a width, b, and a height, h
#
#   NOTE: I am not assuming anything fancy like giving the shapes
#         some random orientation. So for Cubes and Cuboids, the assumption
#         is that their faces point along the x, y, z directions. Therefore
#         when determining if a particular rock intersects a plane at point z0
#         the only relevant variable is its height, and z position.
#
###################################################

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import quad

def Fk(D):
     ## This function is the cumulative fractional area covered by
     ## rocks with diameter > D per area
     ## Basically says look at the ground, every square area Fk(D) of the
     ## area is covered by rocks that are at least D in diameter
     ## This is the distribution reported in Wu et al. (2021)
     k = 0.0021
     qk = 0.5648 + 0.01258/k
     Fk = k*np.exp(-qk*D)
     return Fk

def Fk_Di(D):
     ## This function is the cumulative fractional area covered by
     ## rocks with diameter > D per area
     ## Basically says look at the ground, every square area Fk(D) of the
     ## area is covered by rocks that are at least D in diameter
     ## This is the distribution reported in Wu et al. (2021)
     k = 0.0125
     qk = 1.743
     Fk = k*np.exp(-qk*D)
     return Fk

def NgtD4_Di(D):
     from scipy.integrate import quad
     Diam_array = np.arange(D, 3, 0.001)
     N = []
     for i in range(len(Diam_array)):
          N.append(quad(nDpm2, Diam_array[i], np.inf)[0])
     return N

def N_SVII(Dmin, Dmax):
     D = np.arange(Dmin, Dmax, 0.001)
     K = 7.9e3/1000**1.8 
     gamma = 1.8
     N = K*D**(-gamma)
     return N

def n_a_SVII(D):
     K = 7.9e3/1000**1.8 
     gamma = 1.8
     n = gamma*K*D**(-(gamma + 1))
     return n

def n_V_SVII(D):
     K = 7.9e3/1000**1.8 
     gamma = 1.8
     n = gamma*K*D**(-(gamma + 2))
     return n

def F_SVII(Dmin, Dmax):
     K = 7.9e3/1000**1.8 
     gamma = 1.8
     coeff = gamma*K/(2 - gamma)
     F = coeff*(Dmax**(2 - gamma) - Dmin**(2 - gamma))
     return F




def expintA(x):
    ## Used in the approximation of integrating Fk(D) to get 
    ## NgtD4
    A = np.log((0.56146/x + 0.65)*(1+x))
    B = (x**4)*np.exp(7.7*x)*(2+x)**3.7
    Ei = ((A)**-7.7+B)**-0.13
    return Ei

def NgtD4(D, Shape = "Sphere"):
    ## This is the cumulative number of rocks with diameter > D
    ## per area. This is derived from Fk(D). Its an integral
    ## but its written in this way to match how this is calculated
    ## in XFdtd, which doesnt have an integrate function so this is
    ## an approximation
    if(Shape == "Sphere" or Shape == "Torus"):
          GeoConst = 4/np.pi
    if(Shape == "Cube" or Shape == "Pyramid"):
          GeoConst = 1
    if(Shape == "Cuboid"):
          GeoConst = 1/0.8
    x = 6.555*D
    N = GeoConst*0.0137608*(np.exp(-x)/D - 6.555*expintA(x))
    return N

def nDpm2(D, Shape = "Sphere", Distribution = "CE3"):
    ## This is the number density of rocks with diameter D per m^2
    ## This is just a derivative of the cumulative number (NgtD4)
    if(Shape == "Sphere" or Shape == "Torus"):
          GeoConst = 4/np.pi
    if(Shape == "Cube" or Shape == "Pyramid"):
          GeoConst = 1
    if(Shape == "Cuboid"):
          GeoConst = 1/0.8
    k = 0.0135
    qk = 1.734
    n = GeoConst*k*qk*np.exp(-qk*D)/D**2
    return n

def DnDpm2(D, Shape = "Cube"):
    ## This is the same as nDpm2 but multiplied by D
    ## used for finding the average rock size
    if(Shape == "Sphere" or Shape == "Torus"):
          GeoConst = 4/np.pi
    if(Shape == "Cube" or Shape == "Pyramid"):
          GeoConst = 1
    if(Shape == "Cuboid"):
          GeoConst = 1/0.8
    k = 0.0021
    qk = 0.5648 + 0.01258/k
    n = GeoConst*k*qk*np.exp(-qk*D)/D
    return n

def nDpm3(D, Shape = "Sphere", Distribution = "CE3"):
    ## This is the number density of rocks with diameter D per m^3. 
    ## This here is a complete guess on how to extend a 2D density into a 3D density
    ## Basically divided by D and replaced the Area constant with a volume constant
    ## This is mainly what I want to check
    if(Shape == "Sphere" or Shape == "Torus"):
          GeoConst = 4/np.pi
    if(Shape == "Cube" or Shape == "Pyramid"):
          GeoConst = 1
    if(Shape == "Cuboid"):
          GeoConst = 1/(0.8*0.54)
    if(Distribution == "CE4"):
        k = 0.0021
        qk = 0.5648 + 0.01258/k
        n = GeoConst*k*qk*np.exp(-qk*D)/D**3
    elif(Distribution == "CE3"):
         k = 0.0125
         qk = 1.743
         n = GeoConst*k*qk*np.exp(-qk*D)/D**3
    return n

def DnDpm3(D, Shape = "Sphere", Distribution = "CE3"):
    ## Same as nDpm3 but multiplied by D
    if(Shape == "Sphere" or Shape == "Torus"):
          GeoConst = 4/np.pi
    if(Shape == "Cube" or Shape == "Pyramid"):
          GeoConst = 1
    if(Shape == "Cuboid"):
          GeoConst = 1/(0.8*0.54)
    if(Distribution == "CE4"):
        k = 0.0021
        qk = 0.5648 + 0.01258/k
        n = GeoConst*k*qk*np.exp(-qk*D)/D**2
    elif(Distribution == "CE3"):
         k = 0.0125
         qk = 1.743
         n = GeoConst*k*qk*np.exp(-qk*D)/D**2
    return n

def CheckNgtD(rocks, Dmin, A, Dmax = 1.25, Shape = "Sphere"):
    ## This function takes in an array of rock sizes and returns the
    ## cumulative number of rocks per A.
    ## This is intended to be passed an array or list of rock sizes
    ## found to intersect some plane at some z0 and return the equivalent
    ## of NgtD4. So this function assumes all rocks in the array 'rocks' are in a single 
    ## plane that has an area A.
    D = np.arange(Dmin, Dmax, 0.001)
    N_rocks = []
    for i in range(len(D)): ## Step through a series of D values
        cuml = 0
        for ii in range(len(rocks)): ## Go through all the rocks in the list given
            if(Shape != "Cuboid"):
                if (rocks[ii] >= D[i]): ## If the rock size is greater than the current
                                        ## D value we are on then count it
                    cuml += 1
            else:
                if (rocks[ii][0] >= D[i]):
                    cuml += 1
        N_rocks.append(cuml/(A))
    return D, N_rocks

def CheckFgtD(rocks, Dmin, A, Dmax = 1.25, Shape = "Sphere"):
    ## This function takes in an array of rock sizes and returns the
    ## cumulative fractional area covered by rocks with size greater than D.
    ## This is intended to be passed an array or list of rock sizes
    ## found to intersect some plane at some z0 and return the equivalent
    ## of NgtD4. So this function assumes all rocks in the array 'rocks' are in a single 
    ## plane that has an area A.

    ## The returns of this function are the diameter list it used to calculate the fractional area
    ## as well as the fractional area greater than each value in the D list
    ## Therefore F_Area[0] will be the cumulative fractional area of all rocks populated into the space
    D = np.arange(Dmin, Dmax, 0.001)
    F_Area = []
    for i in range(len(D)):
        cuml = 0
        for ii in range(len(rocks)):
            if (Shape != "Cuboid"):
                if (rocks[ii] >= D[i]):
                    if (Shape == "Cube"):
                        cuml += rocks[ii]**2
                    if (Shape == "Sphere"):
                        cuml += (np.pi/4)*rocks[ii]**2
            else:
                if (rocks[ii][0] >= D[i]):
                    cuml += rocks[ii][0]*rocks[ii][1]
        F_Area.append(cuml/(A))
    return D, F_Area

def hpD(x):
    ## For cuboids this is a random sampling of the height
    mean = 0.54
    std = 0.03
    expo = 0.5*((x - mean)/std)**2
    return np.exp(-expo)

def bpD(x):
    ## for cuboids this a random sampling of its width
    mean = 0.8
    std = 0.16
    expo = 0.5*((x - mean)/std)**2
    return np.exp(-expo)

#### END FUNCTION DEFINITIONS ####


#### BEGIN MAIN CODE ####

D = 0.05 ## Smallest rock diameter allowed to be sampled
Dmax = 3.0
Distribution = "CE3"
D_ave = quad(DnDpm3, D, Dmax)[0]/quad(nDpm3, D, Dmax)[0] ## Calculated the average rock diameter in the space
                                                             ## according to the assuming volume number density
print("Average Diameter: ", D_ave)


Depth = 10  ## The volume being filled will be a rectangular prism, this is the z-direction
Area = 50*50 ## Cross section of the prism
Shape = "Sphere"  ## Shape that the scatterers will be (This determines constants when calculating volume, area, and number)
OccVol = 0

if (Shape != "Cuboid"): ## Cuboids are special so they need to be treated differently
    #NRocks = NgtD4(D, Shape = Shape)*Area*Depth/D_ave   ## Calculate the number of rocks in the volume from the distribution in literature
    #print("NSphere:", NRocks)
    NRocks = quad(nDpm3, D, Dmax)[0]*Area*Depth
else:   
     NRocks = NgtD4(D, Shape = Shape)*Area*Depth/(0.54*D_ave)
     print("NSphere:", NRocks)
     NRocks = quad(nDpm3, D, np.inf)[0]*Area*Depth
print("NSphere:", NRocks)

functions = [nDpm3]
function_Name = ['CE3']
Avg_Sampled_NgtD = []
Avg_Sampled_FgtD = []


for ii in functions:
    rocks = []
    xyz = []
    Sampled_NgtD = []
    Sampled_FgtD = []
    ## This next loop samples the function and populates the volume with rocks
    for i in range(int(np.floor(NRocks))):
        if i%(1000) == 0:
            print(i)
        fx = -99
        xmax = Dmax
        #fxmax = ii(D, Shape, Distribution)
        fxmax = ii(D)
        W = fxmax
        while W > fx:
            W = fxmax*np.random.random()
            x = (xmax - D)*np.random.random() + D
            #fx = ii(x, Shape, Distribution)
            fx = ii(x)
            xyzr = (np.sqrt(Area)*np.random.random(), np.sqrt(Area)*np.random.random(), Depth*np.random.random())

        if (Shape == "Cuboid"):
                b = 0
                h = 0
                fb = -99
                fh = -99
                bxmax = x
                hxmax = x
                bymax = 1
                hymax = 1
                while(bymax > fb):
                    bymax = np.random.random(); ## Since the functions are normalized gaussian, only need to sample between 0 and 1
                    bx = np.random.random()
                    fb = bpD(bx)

                
                b = bx*x; ## If we break the loop then we have selected a width

                while(hymax > fh):
                    hymax = np.random.random(); ## Since the functions are normalized gaussian, only need to sample between 0 and 1
                    hx = np.random.random()
                    fh = hpD(hx)
                
                h = hx*x ##If we break the loop then we have selected a height
        xyz.append(np.array(xyzr))
        if(Shape != "Cuboid"):
            rocks.append(x)
            if (Shape == 'Cube'):
                 OccVol = OccVol + x**3
            elif (Shape == 'Sphere'):
                 OccVol = OccVol + (4/3)*np.pi*x**3/8
        else:
             rocks.append(np.array((x, b, h)))
             OccVol = OccVol + x*b*h

    ## Now we need to check that if we take any random slice in z, on average the original distributions are preserved
    for iii in range(1000):
         if(iii%10 == 0):
            print(iii)
         z0 = Depth*np.random.random()
         rocks_in_slice = []
         for iii in range(len(rocks)):
            if(Shape != "Cuboid"):
                if(z0 > xyz[iii][2] - rocks[iii]/2 and z0 < xyz[iii][2] + rocks[iii]/2):
                    rocks_in_slice.append(rocks[iii])
            else:
                if(z0 > xyz[iii][2] - rocks[iii][2]/2 and z0 < xyz[iii][2] + rocks[iii][2]/2):
                    rocks_in_slice.append(rocks[iii])
                     
         Sampled_D, SAMPLED_NGTD = CheckNgtD(rocks_in_slice, D, Area, Dmax, Shape = Shape)
         Sampled_D, SAMPLED_FGTD = CheckFgtD(rocks_in_slice, D, Area, Dmax, Shape = Shape)
         Sampled_NgtD.append(SAMPLED_NGTD)
         Sampled_FgtD.append(SAMPLED_FGTD)
    
    Avg_Sampled_NgtD.append(np.mean(Sampled_NgtD, axis = 0))
    Avg_Sampled_FgtD.append(np.mean(Sampled_FgtD, axis = 0))

xyz = np.array(xyz)
print(len(xyz[:,0]))
if (Shape != 'Cuboid'):
    df = pd.DataFrame({'Diameter (m)': rocks, 'x (m)': xyz[:,0], 'y (m)': xyz[:,1], 'z (m)': xyz[:,2]})
else:
    rocks = np.array(rocks)
    df = pd.DataFrame({'Diameter (m)': rocks[:,0], 'b (m)': rocks[:,1], 'h (m)': rocks[:,2], 'x (m)': xyz[:,0], 'y (m)': xyz[:,1], 'z (m)': xyz[:,2]})
df.to_csv('./Rock_Data_' + str(Area) + 'XY_' + str(Depth) + 'Z_' + Shape + '.csv')

fig, axs = plt.subplots(2,2)
axs[0,0].plot(Sampled_D, NgtD4_Di(D), label = 'Wu et al. (2021)')
for i in range(len(Avg_Sampled_NgtD)):
    axs[0,0].plot(Sampled_D, Avg_Sampled_NgtD[i], label = "Rejection Sampling using " + str(function_Name[i]))
axs[0,0].legend()
axs[0,0].set_xlabel('D [m]')
axs[0,0].set_ylabel('$Rocks/m^2$')
axs[0,0].set_title("Number of rocks/$m^2$ with diameter > D, Area = " + str(Area) +', Depth = ' + str(Depth) + ': ' + Shape)

axs[0,1].plot(Sampled_D, Fk_Di(Sampled_D), label = 'Wu et al. (2021)')
for i in range(len(Avg_Sampled_FgtD)):
    axs[0,1].plot(Sampled_D, Avg_Sampled_FgtD[i], label = "Rejection Sampling using " + str(function_Name[i]))
axs[0,1].legend()
axs[0,1].set_xlabel('D [m]')
axs[0,1].set_ylabel('Fractional Area')
axs[0,1].set_title("Fractional area covered by rocks with diameter > D Area = " + str(Area) +', Depth = ' + str(Depth) + ': ' + Shape)



axs[1,0].plot(Sampled_D, NgtD4_Di(D), label = 'Wu et al. (2021)')
for i in range(len(Avg_Sampled_NgtD)):
    axs[1,0].plot(Sampled_D, Avg_Sampled_NgtD[i], label = "Rejection Sampling using " + str(function_Name[i]))
axs[1,0].set_yscale('log')
axs[1,0].set_xscale('log')
axs[1,0].legend()
axs[1,0].set_xlabel('D [m]')
axs[1,0].set_ylabel('$Rocks/m^2$')
axs[1,0].set_title("Number of rocks/$m^2$ with diameter > D Area = " + str(Area) +', Depth = ' + str(Depth) + ': ' + Shape)

axs[1,1].plot(Sampled_D, Fk_Di(Sampled_D), label = 'Wu et al. (2021)')
for i in range(len(Avg_Sampled_FgtD)):
    axs[1,1].plot(Sampled_D, Avg_Sampled_FgtD[i], label = "Rejection Sampling using " + str(function_Name[i]))
axs[1,1].set_yscale('log')
axs[1,1].set_xscale('log')
axs[1,1].legend()
axs[1,1].set_xlabel('D [m]')
axs[1,1].set_ylabel('Fractional Area')
axs[1,1].set_title("Fractional area covered by rocks with diameter > D Area = " + str(Area) +', Depth = ' + str(Depth) + ': ' + Shape)
plt.tight_layout()
plt.show()

print('Occupied volume fraction: ', OccVol/(Area*Depth))
