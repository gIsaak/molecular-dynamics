import numpy as np
import matplotlib.pyplot as plt
import time

##### Parameters #####
sizeOfBox = 100 # size of each periodic unit cell L (in Angstrom)
timestep = 1e-4 # time between iterations of Euler's method

# Lennard-Jones parameters (all in SI, except length (sigma) in Angstrom)
boltzmann = 1.38e-23
epsilon = 119.8*boltzmann
sigma = 7 # playing around for testing, not true value of sigma

# Argon parameters (SI)
mass = 39.948*1.66e-27

numOfParticles = 2 # 2 particles
numOfDimensions = 2 # 2 dimensions
num_t = 100 # store 100 time steps in array

# Create n x d x 2 numpy array of floats "parameterMatrix" to store n particles
# in d dimensions with 2 (position and velocity) parameters, and num_t timesteps stored.

# For example, access particle 2's position in the y-direction at time_step=0 as
# parameterMatrix[1][1][0][0]
# Access particle 1's momentum in the x-direction at time_step=1 as
# parameterMatrix[0][0][1][1]

parameterMatrix = np.zeros((numOfParticles,numOfDimensions,2,num_t), dtype=float)

def getP1Xcoord(ts):
    return parameterMatrix[0][0][0][ts]
def getP2Xcoord(ts):
    return parameterMatrix[1][0][0][ts]
def getP1Ycoord(ts):
    return parameterMatrix[0][1][0][ts]
def getP2Ycoord(ts):
    return parameterMatrix[1][1][0][ts]
def getP1Xvel(ts):
    return parameterMatrix[0][0][1][ts]
def getP2Xvel(ts):
    return parameterMatrix[1][0][1][ts]
def getP1Yvel(ts):
    return parameterMatrix[0][1][1][ts]
def getP2Yvel(ts):
    return parameterMatrix[1][1][1][ts]

def setP1Xcoord(val,ts):
    parameterMatrix[0][0][0][ts] = val
def setP2Xcoord(val,ts):
    parameterMatrix[1][0][0][ts] = val
def setP1Ycoord(val,ts):
    parameterMatrix[0][1][0][ts] = val
def setP2Ycoord(val,ts):
    parameterMatrix[1][1][0][ts] = val
def setP1Xvel(val,ts):
    parameterMatrix[0][0][1][ts] = val
def setP2Xvel(val,ts):
    parameterMatrix[1][0][1][ts] = val
def setP1Yvel(val,ts):
    parameterMatrix[0][1][1][ts] = val
def setP2Yvel(val,ts):
    parameterMatrix[1][1][1][ts] = val

def getdUdr(r):
    # given a distance r, returns the derivative of the Lennard-Jones potential with respect to r evaluated at that point
    dUdr_index = int(round(r/sizeOfBox*n))
    return dUdr[dUdr_index]

def getParticleDistance(ts):
    # function that returns the minimum distance between the two particles at the given timestep "ts"
    # includes periodic boundary conditions
    x_dist = abs(getP1Xcoord(ts)-getP2Xcoord(ts))
    y_dist = abs(getP1Ycoord(ts)-getP2Ycoord(ts))
    if x_dist > sizeOfBox/2:
        # if the x-distance is greater than half the x-length of the box, the shortest x distance wraps over the periodic BCs
        x_dist = sizeOfBox - x_dist
    if y_dist > sizeOfBox/2:
        # same reasoning as for x coordinate
        y_dist = sizeOfBox - y_dist

    # calculate and return distance shortest distance between particles
    return np.sqrt(x_dist**2 + y_dist**2)

def getPotentialEnergy(ts):
    # returns the potential energy (Lennard Jones) at the given timestep
    r = getParticleDistance(ts)
    return 4*epsilon*((sigma/r)**12-(sigma/r)**6)

def getTotalEnergy(ts):
    # calculates sum of potential energy and kinetic energy at the timestep (of particles in a unit cell)
    U = getPotentialEnergy(ts)
    T = 1/2*mass*(getP1Xvel(ts)**2+getP1Yvel(ts)**2+getP2Xvel(ts)**2+getP2Yvel(ts)**2)
    return U + T

def getForce(ts):
    # calculates the force on each particle at a given timestep in the x and y directions (returned as a 4-tuple: f1x,f2x,f1y,f2y)
    rel_x = getP1Xcoord(ts)-getP2Xcoord(ts)
    if abs(rel_x) > sizeOfBox/2:
        if rel_x < 0:
            rel_x = sizeOfBox-abs(rel_x)
        else:
            rel_x = -sizeOfBox+abs(rel_x)
    rel_y = getP1Ycoord(ts)-getP2Ycoord(ts)
    if abs(rel_y) > sizeOfBox/2:
        if rel_y < 0:
            rel_y = sizeOfBox-abs(rel_y)
        else:
            rel_y = -sizeOfBox+abs(rel_y)
    r = np.sqrt(rel_x**2+rel_y**2)
    grad = getdUdr(r)
    force_P1_X = -grad*rel_x/r
    force_P2_X = grad*rel_x/r # change of sign from force_P1_X since forces must be equal and opposite
    force_P1_Y = -grad*rel_y/r
    force_P2_Y = grad*rel_y/r
    return force_P1_X,force_P2_X,force_P1_Y,force_P2_Y

def iterateCoordinates(ts):
    # takes current position and velocity at timestep ts and updates particle coordinates for timestep ts+1
    newP1Xcoord = enforcePeriodicBCs(getP1Xcoord(ts)+getP1Xvel(ts)*timestep)
    newP1Ycoord = enforcePeriodicBCs(getP1Ycoord(ts)+getP1Yvel(ts)*timestep)
    newP2Xcoord = enforcePeriodicBCs(getP2Xcoord(ts)+getP2Xvel(ts)*timestep)
    newP2Ycoord = enforcePeriodicBCs(getP2Ycoord(ts)+getP2Yvel(ts)*timestep)
    
    if ts+1 == num_t:
        next_ts = 0
    else:
        next_ts = ts+1
    setP1Xcoord(newP1Xcoord,next_ts)
    setP1Ycoord(newP1Ycoord,next_ts)
    setP2Xcoord(newP2Xcoord,next_ts)
    setP2Ycoord(newP2Ycoord,next_ts)

def enforcePeriodicBCs(val):
    # accepts a coordinate and makes sure the value is within the interval [0,1)
    return val%sizeOfBox

def iterateVelocities(ts):
    # takes current velocity and force at timestep ts and updates particle velicities for timestep ts+1
    fx1,fx2,fy1,fy2 = getForce(ts)
    newP1Xvel = getP1Xvel(ts) + 1/mass*fx1*timestep 
    newP1Yvel = getP1Yvel(ts) + 1/mass*fy1*timestep
    newP2Xvel = getP2Xvel(ts) + 1/mass*fx2*timestep
    newP2Yvel = getP2Yvel(ts) + 1/mass*fy2*timestep
    
    if ts+1 == num_t:
        next_ts = 0
    else:
        next_ts = ts+1
    setP1Xvel(newP1Xvel,next_ts)
    setP1Yvel(newP1Yvel,next_ts)
    setP2Xvel(newP2Xvel,next_ts)
    setP2Yvel(newP2Yvel,next_ts)
    
#### Precompute Lennard Jones force ####

# since potential is only a function of distance, we can efficiently
# compute the potential for all possible distances (given periodic BCs),
n = 1000 # number of discrete distances with a precomputed gradient
r_vec = np.linspace(sizeOfBox/n,sizeOfBox,n) # 1D array of all possible distances of particles (up to some user-defined precision)

# evaluate the Lennard Jones Potential on the grid
U_lj = 4*epsilon*((sigma/r_vec)**12-(sigma/r_vec)**6)

# evaluate the derivative of the Lennard Jones Potential with respect to r
dUdr = np.gradient(U_lj,sizeOfBox/n) # access dU/dr with dUdr[distance]

################# Begin main program ########################

# set random starting point for particles
setP1Xcoord(sizeOfBox/2,0)
setP1Ycoord(sizeOfBox/3,0)
setP2Xcoord(1,0)
setP2Ycoord(1,0)
setP1Xvel(5000,0)
setP1Yvel(25000,0)
setP2Xvel(10000,0)
setP2Yvel(-5000,0)

energies = np.zeros(100,dtype=float)


fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
p1 = ax.scatter(parameterMatrix[0][0][0][0],parameterMatrix[0][1][0][0], color='r')
p2 = ax.scatter(parameterMatrix[1][0][0][0],parameterMatrix[1][1][0][0], color='b')
ax.set_xlim((0,sizeOfBox))
ax.set_ylim((0,sizeOfBox))
plt.ion()
plt.show()
plt.pause(0.001)

for j in range(100):
    i = j%num_t # don't go over indices of parameterMatrix
    energies[i] = getTotalEnergy(i)
    iterateCoordinates(i)
    iterateVelocities(i)
    p1.remove()
    p2.remove()
    p1 = ax.scatter(getP1Xcoord(i),getP1Ycoord(i),color='r')
    p2 = ax.scatter(getP2Xcoord(i),getP2Ycoord(i),color='b')
    plt.pause(0.001)
    time.sleep(0.05)
    
